""" map root segments to a soil grid """
import sys;
sys.path.append("../..");
sys.path.append("../../src/")
import plantbox as pb
import visualisation.vtk_plot as vp
import math
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors
from scipy import ndimage

try:
    from pyevtk.hl import gridToVTK
except:
    print("pyevtk library missing. installing pyevtk")
    subprocess.run(["pip3", "install", "pyevtk"])
    from pyevtk.hl import gridToVTK

""" root system """
rs = pb.MappedPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Faba_synMRI" 
rs.readParameters(path + name + ".xml")

#general info
RSage = 10
width = 5  # container width
depth = 19 #container depth 
res = 0.5 #try different resolutions 

fig = plt.figure(figsize=plt.figaspect(0.5))


""" root system """
rs = pb.MappedPlant()
rs.readParameters(path + name + ".xml")


resx = res 
resy = res
resz = res
cyl = pb.SDF_PlantContainer(width - resx/2, width - resy/2, depth - resz/2, False) #some space in order to not touch the edges

# define the grid
nx = int(width / resx);
ny = int(width / resy);
nz = int(depth / resz);
print("Width", width, "Depth", depth, " at a Resolution", nx, ny, nz)
cellvol = resx * resy * resz
RLD = np.zeros((nx * ny * nz))
RVD = np.zeros((nx * ny * nz))

# initialize and simulate
rs.setGeometry(cyl)
rs.setSeed(1)  # use the same seed for all simulations
rs.initialize()
rs.simulate(RSage, False)

""" soil """
min_ = np.array([-width, -width, -depth])
max_ = np.array([width, width, 0])
res_ = np.array([nx, ny, nz])

rs.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments

""" add segment indices """
ana = pb.SegmentAnalyser(rs.mappedSegments())
segs = ana.segments
nodes = ana.nodes
radius = ana.getParameter("radius")
seglen = ana.getParameter("length")

x = np.zeros(len(segs))
for i, s in enumerate(segs):
    try:
        x[i] = rs.seg2cell[i]
    except:  # in case the segment is not within the domain
        x[i] = -1

sc = np.unique(x).astype(int) #all voxels with a segment 
for i in range(0,len(sc)):
    sg = rs.cell2seg[sc[i]] #numbers of all segments within that voxel

    rootvol = 0
    rootlen = 0
    for j in range(0,len(sg)):
        rootlen = rootlen + seglen[sg[j]]
        rootvol = rootvol + radius[sg[j]] ** 2 * math.pi * seglen[sg[j]]
    RLD[sc[i]] = rootlen/cellvol
    RVD[sc[i]] = rootvol/cellvol


RLD = np.reshape(RLD, (nz, ny, nx))
RLD = np.swapaxes(RLD, 0, 2)
RLD_ = RLD[:,int(ny/2),:]

RVD = np.reshape(RVD, (nz, ny, nx))
RVD = np.swapaxes(RVD, 0, 2)
RVD_ = RVD[:,int(ny/2),:]

#PLOT
ax2 = fig.add_subplot(1, 3, 2)
colors = plt.cm.seismic_r(RLD_)
#print(np.min(RLD_),np.max(RLD_))
norm = matplotlib.colors.Normalize(vmin=np.min(RLD_), vmax=np.max(RLD_))
RLDrot90 = ndimage.rotate(RLD_, 90)
ax2.imshow(RLDrot90, cmap = 'seismic_r')
m = cm.ScalarMappable(cmap=plt.cm.seismic_r, norm=norm)
m.set_array([])
plt.colorbar(m, ax = ax2, label = 'Root length density (cm / cm³)')
ax2.set_aspect('equal')
ax2.axis('off')

ax3 = fig.add_subplot(1, 3, 3)
colors = plt.cm.seismic_r(RVD_)
norm = matplotlib.colors.Normalize(vmin=np.min(RVD_), vmax=np.max(RVD_))
RVDrot90 = ndimage.rotate(RVD_, 90)
ax3.imshow(RVDrot90, cmap = 'seismic_r')
m = cm.ScalarMappable(cmap=plt.cm.seismic_r, norm=norm)
m.set_array([])
plt.colorbar(m, ax = ax3, label = 'Root volume density (cm³ / cm³)')
ax3.set_aspect('equal')
ax3.axis('off')

#save as .raw
#namep = name.split('_')
#C = np.swapaxes(C, 0, 2)
#C = C*int(255) #int(8)
#C = C[::-1] #up-down
#C.astype('int8').tofile('results/' + namep[0] + '_day' + str(RSage) + '_' + str(nx) + 'x' + str(ny) + 'x' + str(nz) + '.raw')
#print('done')

ax1 = fig.add_subplot(1, 3, 1)
cmap = cm.seismic
norm = matplotlib.colors.Normalize(vmin=np.min(radius), vmax=np.max(radius))
fc = cmap(norm(radius))
for k, s in enumerate(segs):
    s1 = segs[k]
    n1, n2 = nodes[s1.x], nodes[s1.y]
    ax1.plot([n1.x, n2.x], [n1.z, n2.z], color = fc[k,:])
m = cm.ScalarMappable(cmap=plt.cm.seismic)
m.set_array([])
plt.colorbar(m, ax = ax1, label = 'Radius (cm)')
ax1.set_aspect('equal')
ax1.axis('off')
plt.show()

