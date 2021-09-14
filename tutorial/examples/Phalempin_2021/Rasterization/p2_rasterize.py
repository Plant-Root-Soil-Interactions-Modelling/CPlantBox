import sys
sys.path.append("../../../../..")
import plantbox as pb
import vtk_plot as vp
import numpy as np
import time
import bresenham3D as bres3D
import matplotlib.pyplot as plt
from pyevtk.hl import gridToVTK
import scipy
from scipy import ndimage
import pygorpho as pg
import skimage
from skimage.morphology import ball
import multiprocessing
from multiprocessing import Process, active_children
import psutil
from threading import Thread 
import math


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def soil_cores(x :list, y :list, r :float,  h :list, up: list):
    """
    A lsit of soil core geometries with a fixed location in the field  
 
    @param x     x coordinates of the soil cores (cm)
    @param y     y coordinates of the soil cores (cm)
    @param r     radius of the soil core (cm)
    @param h     height of the soil core (cm)
    """
    assert len(x) == len(y), "coordinate length must be equal"
    core = []
    for i in range(0,len(h)): 
        core.append(pb.SDF_PlantContainer(r, r, h[i], False))
    cores = []
    for i in range(0, len(x)):
        cores.append(pb.SDF_RotateTranslate(core[i], 0., pb.SDF_Axis.xaxis, pb.Vector3d(x[i], y[i], up[i])))  # just translate
    return cores;

def set_all_sd(rs):
    for p in rs.getRootRandomParameter():
        p.las = 0
        p.lbs = 0
        p.rs = 0
        p.lmaxs = 0
        p.thetas = 0
        
def raster_model(number):   
     
    s = segs[number]

    #find start and end points of the segment and use bresenham algorithm to define all voxels in between 
    n1, n2 = nodes[s.x], nodes[s.y]
    (x1, y1, z1) = [np.around(n1.x/res), np.around(n1.y/res), np.around(n1.z/res)]
    (x2, y2, z2) = [np.around(n2.x/res), np.around(n2.y/res), np.around(n2.z/res)]
    ListOfPoints = np.array(bres3D.Bresenham3D(x1, y1, z1, x2, y2, z2))

    allidx_ = [] 
    Czero = np.zeros((nx, ny, nz))
    for j in range(0, len(ListOfPoints)): 
        xidx = np.where(xx==ListOfPoints[j,0])
        yidx = np.where(yy==ListOfPoints[j,1])
        zidx = np.where(zz==ListOfPoints[j,2])
        if (xidx[0].size > 0 and yidx[0].size > 0 and zidx[0].size > 0):
            a = [int(xidx[0]),int(yidx[0]),int(zidx[0])]
            allidx_.append(a)
    allidx = np.array(allidx_)
    if (len(allidx)):
        Czero[allidx[:,0], allidx[:,1], allidx[:,2]] = 1

    #dilate each root segment by the appropriate value
    Czero1 = Czero
    value = int(np.around(radius[number]/res))
    print(value)
    if value > 15:
        value = 15
    if value > 0:
        struct = skimage.morphology.ball(value, np.int8) 
        Czero1 = ndimage.binary_dilation(Czero1, structure = struct).astype(Czero1.dtype)
    Czero1 = np.reshape(Czero1, (nx*ny*nz))
    idx = np.where(Czero1==1)
    indx = idx[0]
    
    return indx;

        
def rasterize(number,q):
    q.put([(raster_model(number))])
#########################################################
mm =1
rs = pb.RootSystem()
path = "../rootsys_params/"
name = "p2_Optimization_diam"
rs.readParameters(path + name + ".xml")

#general info
corenum = 3 #number of soil cores taken 
times = 21 #days after planting
r1 = 1.5  # core radius
bigcyl = pb.SDF_PlantContainer(5, 5, 18, False)

#define the cores (x,y,r,h,up)
A = np.zeros((3,4))
A = [[-2.5,   0,  -3.5,  -6.5],
     [-2.5,    0,  -8.5,  -11.5],
     [-2.5,    0,  -13.5, -16.5]]
A_ = np.array(A)
x = A_[:,0]
y = A_[:,1]
h = np.subtract(A_[:,2],A_[:,3])
up = A_[:,2]
dow = A_[:,3]
cores = soil_cores(x, y, r1, h, up)

#define the grid
width = 3  # cm
depth = 3
res = 0.006  #60microm 0.006
nx = int(width / res);
ny = int(width / res);
nz = int(depth / res);
print("Width", width, "Depth", depth, " at a Resolution", nx, ny, nz)


#initialize and simulate
set_all_sd(rs)
rs.setGeometry(bigcyl)
rs.setSeed(1)
rs.initialize()
rs.simulate(21, False)
rs.write("results/vtp/sims"+str(mm)+"_rootsystem.vtp")


core_RLD = np.zeros((3,))
for z in range(0, len(cores)):

    #mesh
    X = np.linspace(-4, -1, nx+1)
    Y = np.linspace(-1.5, 1.5, ny+1)
    Z = np.linspace(A_[z,3], A_[z,2], nz+1)
    xx = np.around(X[:-1]/res)
    yy = np.around(Y[:-1]/res)
    zz = np.around(Z[:-1]/res)

    ana =  pb.SegmentAnalyser(rs)
    ana.crop(cores[z]);
    ana.pack() #delete all the unneeded nodes and segs
    tl1 = ana.distribution("length", up[z], dow[z], 1, True)
    tl1 = np.array(tl1) / ( r1 * r1 * math.pi * h[z])  
    core_RLD[z] = tl1

    nodes = ana.nodes
    segs = ana.segments
    radius = ana.getParameter("radius")
    
    t = time.time()
    C = np.zeros((nx*ny*nz,))
    
    print(len(segs))
    max = 20 #maximum number of segments that should be processd in parallel 
    tranches = int(np.ceil(len(segs)/max))
    task = np.zeros(tranches)
    

    for j in range(0,tranches): 
        q = multiprocessing.Queue()

        if j== tranches-1: 
            end = len(segs)
        else:
            end = int((j+1)*max)+1
        if j== 0: 
            start = 0
        else:
            start = int(j*max+1)


        for k in range(start,end):
            p = multiprocessing.Process(target=rasterize, args=(k,q,))
            p.start()
        
        for k in range(start,end):
            indx = np.array(q.get())
            C[indx] = 1 


    C = np.reshape(C, (nx,ny,nz))
    elapsed = time.time() - t
    print("Computation took", elapsed, "s")

    univ_name = "sims"+str(mm)+"_core_cropped" + str(z); 
    gridToVTK("./"+'results/vtr/'+univ_name, X, Y, Z, cellData = {"root":C})
    C.astype('int8').tofile('results/raw/'+univ_name+'_'+str(nz)+'x'+str(nx)+'x'+str(ny)+'.raw')
    ana.write('results/vtp/'+univ_name+"_"+str(tl1)+".vtp")
    print('core ', str(z), ' done')

np.savetxt("results/RLD/sims_"+str(mm)+".txt",core_RLD,fmt='%.2f')

