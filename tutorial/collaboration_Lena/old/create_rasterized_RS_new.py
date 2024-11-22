"""simulates a root system, which is rasterized to a given resolution. """


import sys
sys.path.append("../../")
sys.path.append("../../src")
import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np
import functional.bresenham3D as bres3D
import matplotlib.pyplot as plt
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
try:
    from pyevtk.hl import gridToVTK
except:
    print("pyevtk library missing. installing pyevtk")
    subprocess.run(["pip3", "install", "pyevtk"])
    from pyevtk.hl import gridToVTK
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
    row, col, ch = image.shape
    sigma = var
    gauss = np.random.normal(mean, sigma, (row, col, ch))
    gauss = gauss.reshape(row, col, ch)
    noisy = image + gauss
    return noisy


def voxel_lookup(x :list, y :list, w :float,  h :list, up: list):
    """
    A lsit of soil core geometries with a fixed location in the field  
 
    @param x     x coordinates of the voxel (cm)
    @param y     y coordinates of the voxel(cm)
    @param w     width of the voxel(cm)
    @param h     height of the soil core (cm)
    """
    assert len(x) == len(y), "coordinate length must be equal"
    voxel = pb.SDF_PlantContainer(w/2, w/2, h, True)
    voxels = []
    for i in range(0, len(x)):
        voxels.append(pb.SDF_RotateTranslate(voxel, 0., pb.SDF_Axis.xaxis, pb.Vector3d(x[i], y[i], up[i])))  # just translate
    return voxels;

#########################################################
rs = pb.RootSystem()
path = "../../modelparameter/structural/rootsystem/"
name = "Faba_synMRI"
# name = "Zeamays_synMRI"
wcroot = 0.85

rs.readParameters(path + name + ".xml")

RSage = [5,7]
width = 2.8  # radius of cylindrical container
depth = 19
resx = 0.25  # choose your resolution
resy = 0.25
resz = 0.25

# radius histograms for each time point
fig, axs = plt.subplots(len(RSage))

# define the grid
nx = int(width * 2 / resx);
ny = int(width * 2 / resy);
nz = int(depth / resz);
print("Width", width, "Depth", depth, " at a Resolution", nx, ny, nz)
cellvol = resx * resy * resz

# general info
bigcyl = pb.SDF_PlantContainer(width - resx/2, width - resy/2, depth - resz/2, False) #some space in order to not touch the edges

# initialize and simulate
rs.setGeometry(bigcyl)
rs.setSeed(1)  # use the same seed for all simulations
rs.initialize()
dummy = 0
for i in range(0, RSage[-1] + 1):
    simtime = 1
    rs.simulate(simtime, False)

    if int(i) in RSage:

        # for visualization in paraview
        # rs.write("results/"+name+'_day'+str(i)+".vtp")

        # mesh
        X = np.linspace(-1 * width, -1 * width + nx * resx, nx + 1)
        Y = np.linspace(-1 * width, -1 * width + ny * resy, ny + 1)
        Z = np.linspace(0, 0 - nz * resz, nz + 1)
        xx = np.ceil(X[:-1] / resx)
        yy = np.ceil(Y[:-1] / resy)
        zz = np.round(Z[:-1] / resz)

        if not (isConsecutive(xx)):
            xx = np.linspace(np.min(xx), np.max(xx), len(xx))
        if not (isConsecutive(yy)):
            yy = np.linspace(np.min(yy), np.max(yy), len(yy))
        if not (isConsecutive(zz)):
            zz = np.linspace(0, np.min(zz), len(zz))

        ana = pb.SegmentAnalyser(rs)
        nodes = ana.nodes
        segs = ana.segments
        radius = ana.getParameter("radius")
        seglen = ana.getParameter("length")
        C = np.zeros((nx * ny * nz))

        # plot diameter distribution as histogram
        axs[dummy].hist(radius, 15, weights = seglen)
        if dummy == len(RSage) - 1:
            axs[dummy].set_xlabel('Radius (cm)')
            plt.setp(axs, xlim = [0, np.max(radius)])
        axs[dummy].set_ylabel('Length (cm)')
        dummy = dummy + 1

        # order segs from small to large radius
        idxrad = np.argsort(radius)

        for k, s in enumerate(segs):

            #cropping modifies segs info
            ana = pb.SegmentAnalyser(rs)
            nodes = ana.nodes
            segs = ana.segments
            radius = ana.getParameter("radius")
            seglen = ana.getParameter("length")

            
            # find start and end points of the segment and use bresenham algorithm to define all voxels in between
            s1 = segs[idxrad[k]]
            n1, n2 = nodes[s1.x], nodes[s1.y]

            allidx_ = []
            for z in range(0,2):
                if z == 0:
                    (x1, y1, z1) = [np.floor((n1.x ) / resx), np.floor((n1.y ) / resy), np.floor((n1.z ) / resz)]
                    (x2, y2, z2) = [np.floor((n2.x ) / resx), np.floor((n2.y ) / resy), np.floor((n2.z ) / resz)]
                else: 
                    (x1, y1, z1) = [np.ceil((n1.x ) / resx), np.ceil((n1.y ) / resy), np.ceil((n1.z ) / resz)]
                    (x2, y2, z2) = [np.ceil((n2.x ) / resx), np.ceil((n2.y ) / resy), np.ceil((n2.z ) / resz)]

                ListOfPoints = np.array(bres3D.Bresenham3D(x1, y1, z1, x2, y2, z2))

                Czero = np.zeros((nx, ny, nz))
                for j in range(0, len(ListOfPoints)):
                    xidx = np.where(xx == ListOfPoints[j, 0])
                    yidx = np.where(yy == ListOfPoints[j, 1])
                    zidx = np.where(zz == ListOfPoints[j, 2])
                    if (xidx[0].size > 0 and yidx[0].size > 0 and zidx[0].size > 0):
                        a = [int((xidx[0]).item()), int((yidx[0]).item()), int((zidx[0]).item())]
                        allidx_.append(a)
            allidx = np.array(allidx_)

            if (np.round(radius[idxrad[k]] * 2 / resx) > 1):
                if (len(allidx)):
                    Czero[allidx[:, 0], allidx[:, 1], allidx[:, 2]] = 1
                Czero1 = Czero
                if radius[idxrad[k]] > 0.1:
                    radius[idxrad[k]] = 0.1
                value = int(np.around(radius[idxrad[k]] / resx))
                struct = skimage.morphology.ball(value, np.int8)
                Czero1 = ndimage.binary_dilation(Czero1, structure = struct).astype(Czero1.dtype)
                Czero1 = np.reshape(Czero1, (nx * ny * nz))
                idx = np.where(Czero1 == 1)
                C[idx[0]] = wcroot
            else:
                if (len(allidx)):
                    x_vox = X[allidx[:, 0]]
                    y_vox = Y[allidx[:, 1]]
                    up = Z[allidx[:, 2]]
                    dow = Z[allidx[:, 2]+1]
                    #print(x_vox, y_vox)
                    #print(up, dow)
                    voxels = voxel_lookup(x_vox, y_vox, resx, resz, up)

                    for p in range(0,len(allidx)):
                        ana1 = ana
                        ana1.crop(voxels[p]);
                        ana1.pack()
                        tl = ana1.distribution("length", up[p], dow[p], 1, True)[0]

                        print(tl) 
                        estlen = tl #resx  # estimated segment length within voxel: very rough estimation
                        rootvol = radius[idxrad[k]] ** 2 * math.pi * estlen
                        frac = rootvol / cellvol
                        if frac > 1:
                            frac = 1
                        Czero[allidx[p, 0], allidx[p, 1], allidx[p, 2]] = frac * wcroot

                # set root voxels to the appropriate value
                Czero = np.reshape(Czero, (nx * ny * nz))
                idx = np.where(Czero>0)
                C[idx[0]] = Czero[idx[0]]

        C = np.reshape(C, (nx, ny, nz))

        namep = name.split('_')

        #save original RS
        rs.write('results/'+namep[0]+'_day'+str(i)+'.vtp')

        # for visualization in paraview
        gridToVTK("./"+'results/'+namep[0]+'_res_'+str(resx)+'_day'+str(i)+'_final', X, Y, Z, cellData = {"root":C})

        C = np.swapaxes(C, 0, 2)
        C = C[::-1]
        #C.astype('int8').tofile('results/' + namep[0] + '_day' + str(i) + '_' + str(nx) + 'x' + str(ny) + 'x' + str(nz) + '.raw')
        print('done')

plt.show()

