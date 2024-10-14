"""simulates a root system, which is rasterized to a given resolution. To mimic an MRI image, Gaussian noise is additionally added"""
import sys
sys.path.append("../../")
sys.path.append("../../src")
import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np
import functional.bresenham3D as bres3D
import matplotlib.pyplot as plt
from pyevtk.hl import gridToVTK
import scipy
from scipy import ndimage
import subprocess
try:
    import pygorpho as pg
except:
    print("pygorpho library missing. installing pygorpho")    
    subprocess.run(["pip3", "install", "pygorpho"])
    import pygorpho as pg
try:
    import skimage 
except:
    print("skimage library missing. installing scikit-image")    
    subprocess.run(["pip3", "install", "scikit-image"])
    import skimage
from skimage.morphology import ball
try:
    import nibabel as nib 
except:
    print("nibabel library missing. installing nibabel")    
    subprocess.run(["pip3", "install", "nibabel"])
    import nibabel as nib
    

import math

import os

def isConsecutive(A):
    if len(A) <= 1:
        return True
 
    minimum = min(A)
    maximum = max(A)
 
    if maximum - minimum != len(A) - 1:
        return False

    visited = set()
    for i in A:
        if i in visited:
            return False
        visited.add(i)
    return True

def noisy(image, mean, var):
    row,col,ch= image.shape
    sigma = var 
    gauss = np.random.normal(mean,sigma,(row,col,ch))
    gauss = gauss.reshape(row,col,ch)
    noisy = image + gauss
    return noisy
#########################################################
rs = pb.RootSystem()
path = "../../modelparameter/structural/rootsystem/"
name = "Faba_synMRI"
#name = "Zeamays_synMRI" 
SNR = 5 #can be changed to other values

rs.readParameters(path + name + ".xml")

RSage = [5,7,10] 
width = 2.8 #radius of cylindrical container 
depth = 19
resx = 0.05 #choose your resolution 
resy = 0.05
resz = 0.1

#radius histograms for each time point 
fig, axs = plt.subplots(len(RSage))

#define the grid
nx = int(width*2 / resx);
ny = int(width*2 / resy);
nz = int(depth / resz);
print("Width", width, "Depth", depth, " at a Resolution", nx, ny, nz)
cellvol = resx*resy*resz

#general info
bigcyl = pb.SDF_PlantContainer(width-0.05, width-0.05, depth-0.05, False)

#initialize and simulate
rs.setGeometry(bigcyl)
rs.setSeed(1) #use the same seed for all simulations 
rs.initialize()
dummy = 0
for i in range(0,RSage[-1]+1): 
    simtime = 1
    rs.simulate(simtime, False)

    if int(i) in RSage:

        #for visualization in paraview 
        #rs.write("results/"+name+'_day'+str(i)+".vtp")
        
        #mesh
        X = np.linspace(-1*width, -1*width+nx*resx, nx+1)
        Y = np.linspace(-1*width, -1*width+ny*resy, ny+1)
        Z = np.linspace(0, 0-nz*resz, nz+1)
        xx = np.ceil(X[:-1]/resx)
        yy = np.ceil(Y[:-1]/resy)
        zz = np.round(Z[:-1]/resz)

        if not (isConsecutive(xx)):
            xx = np.linspace(np.min(xx), np.max(xx), len(xx))
        if not (isConsecutive(yy)):
            yy = np.linspace(np.min(yy), np.max(yy), len(yy))
        if not (isConsecutive(zz)):
            zz = np.linspace(0, np.min(zz), len(zz))

        ana =  pb.SegmentAnalyser(rs)
        nodes = ana.nodes
        segs = ana.segments
        radius = ana.getParameter("radius")
        seglen = ana.getParameter("length")
        C = np.zeros((nx*ny*nz))

        #plot diameter distribution as histogram
        axs[dummy].hist(radius,15, weights=seglen)
        if dummy == len(RSage)-1: 
            axs[dummy].set_xlabel('Radius (cm)')
            plt.setp(axs, xlim=[0, np.max(radius)])
        axs[dummy].set_ylabel('Length (cm)')
        dummy = dummy+1

        #order segs from small to large radius
        idxrad = np.argsort(radius)

        for k, s in enumerate(segs):
            #find start and end points of the segment and use bresenham algorithm to define all voxels in between
            s1 = segs[idxrad[k]]
            n1, n2 = nodes[s1.x], nodes[s1.y]
            (x1, y1, z1) = [np.around(n1.x/resx), np.around(n1.y/resy), np.around(n1.z/resz)]
            (x2, y2, z2) = [np.around(n2.x/resx), np.around(n2.y/resy), np.around(n2.z/resz)]
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

            if (np.round(radius[idxrad[k]]*2/resx)>1):
                if (len(allidx)):
                    Czero[allidx[:,0], allidx[:,1], allidx[:,2]] = 1
                Czero1 = Czero
                if radius[idxrad[k]]> 0.1:
                    radius[idxrad[k]]= 0.1
                value = int(np.around(radius[idxrad[k]]/resx))
                struct = skimage.morphology.ball(value, np.int8)
                Czero1 = ndimage.binary_dilation(Czero1, structure = struct).astype(Czero1.dtype)
                Czero1 = np.reshape(Czero1, (nx*ny*nz))
                idx = np.where(Czero1==1)
                C[idx[0]] = int(255)
            else:
                estlen = resx #estimated segment length within voxel: very rough estimation
                rootvol = radius[idxrad[k]]**2*math.pi*estlen
                frac = rootvol/cellvol
                if frac>1:
                    frac = 1
                #print(rootvol, cellvol, rootvol/cellvol, radius[idxrad[k]]**2*math.pi)
                if (len(allidx)):
                    Czero[allidx[:,0], allidx[:,1], allidx[:,2]] = 1

                #set root voxels to the appropriate value
                Czero = np.reshape(Czero, (nx*ny*nz))
                idx = np.where(Czero==1)
                C[idx[0]] = int(np.floor(frac*255))


        C = np.reshape(C, (nx,ny,nz))


        #add Gaussian noise for given snr: SNR = 0.66 * mean(signal) / std(air), SNR ~ +/-80 - (Firbank et al, 1999)
        meanC = np.mean(C[C>0])
        n_var = 0.66/SNR*meanC
        C_ = noisy(C,meanC,n_var)
        C_ = C_-255
        C_[C_<0]= 0
        C_[C_>255] = 255

        namep = name.split('_')

        #for visualization in paraview 
        #gridToVTK("./"+'results/'+namep[0]+'_add_n_'+str(SNR)+'_day'+str(i), X, Y, Z, cellData = {"root":C_})

        C_ = np.swapaxes(C_, 0,2)
        C_ = C_[::-1]
        C_.astype('int8').tofile('results/'+namep[0]+'_SNR_'+str(SNR)+'_day'+str(i)+'_'+str(nx)+'x'+str(ny)+'x'+str(nz)+'.raw')
        print('done')

plt.show()

