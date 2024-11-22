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
res = [0.25, 0.5, 1] #try different resolutions 

fig = plt.figure(figsize=plt.figaspect(0.5))
for m in range(0,3):

    """ root system """
    rs = pb.MappedPlant()
    rs.readParameters(path + name + ".xml")

    
    resx = res[m]  
    resy = res[m]
    resz = res[m]
    cyl = pb.SDF_PlantContainer(width - resx/2, width - resy/2, depth - resz/2, False) #some space in order to not touch the edges
    wcroot = 0.85

    # define the grid
    nx = int(width / resx);
    ny = int(width / resy);
    nz = int(depth / resz);
    print("Width", width, "Depth", depth, " at a Resolution", nx, ny, nz)
    cellvol = resx * resy * resz
    C = np.zeros((nx * ny * nz))

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
        for j in range(0,len(sg)):
            rootvol = rootvol + radius[sg[j]] ** 2 * math.pi * seglen[sg[j]]
        frac = rootvol / cellvol
        C[sc[i]] = frac * wcroot

    C[C>1] = 1
    C = np.reshape(C, (nz, ny, nx))
    C = np.swapaxes(C, 0, 2)
    
    #PLOT
    ax2 = fig.add_subplot(2, 2, m+2, projection='3d')
    colors = plt.cm.seismic_r(C)
    norm = matplotlib.colors.Normalize(vmin=np.min(C), vmax=np.max(C))
    ax2.voxels(C, facecolors=colors, edgecolor='none')
    m = cm.ScalarMappable(cmap=plt.cm.seismic_r, norm=norm)
    m.set_array([])
    plt.colorbar(m, ax = ax2, label = 'Water content (-)')
    ax2.set_aspect('equal')
    ax2.axis('off')

    #save as .raw
    namep = name.split('_')
    C = np.swapaxes(C, 0, 2)
    C = C*int(255) #int(8)
    C = C[::-1] #up-down
    C.astype('int8').tofile('results/' + namep[0] + '_day' + str(RSage) + '_' + str(nx) + 'x' + str(ny) + 'x' + str(nz) + '.raw')
    print('done')


ax1 = fig.add_subplot(2, 2, 1, projection='3d')
cmap = cm.seismic
norm = matplotlib.colors.Normalize(vmin=np.min(radius), vmax=np.max(radius))
fc = cmap(norm(radius))
for k, s in enumerate(segs):
    s1 = segs[k]
    n1, n2 = nodes[s1.x], nodes[s1.y]
    ax1.plot3D([n1.x, n2.x], [n1.y, n2.y], [n1.z, n2.z], color = fc[k,:])
m = cm.ScalarMappable(cmap=plt.cm.seismic)
m.set_array([])
plt.colorbar(m, ax = ax1, label = 'Radius (cm)')
ax1.set_aspect('equal')
ax1.axis('off')
#ax1.axis('equal')
plt.show()


sys.exit()
namep = name.split('_')

#save original RS
rs.write('results/'+namep[0]+'_day'+str(RSage)+'.vtp')

# for visualization in paraview
X = np.linspace(-1 * width/2, -1 * width/2 + nx * resx, nx + 1)
Y = np.linspace(-1 * width/2, -1 * width/2 + ny * resy, ny + 1)
Z = np.linspace(0, 0 - nz * resz, nz + 1)
gridToVTK("./"+'results/'+namep[0]+'_res_'+str(resx)+'_day'+str(RSage), X, Y, Z, cellData = {"root":C})

C = np.swapaxes(C, 0, 2)
C = C[::-1]
#C.astype('int8').tofile('results/' + namep[0] + '_day' + str(i) + '_' + str(nx) + 'x' + str(ny) + 'x' + str(nz) + '.raw')
print('done')


#print(x)
#sys.exit()


""" vizualise roots """
# ana = pb.SegmentAnalyser(rs)  # <---- wrong!
ana = pb.SegmentAnalyser(rs.mappedSegments())
ana.addData("linear_index", x)
pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "linear_index"])
rootActor, rootCBar = vp.plot_roots(pd, "linear_index", "segment index plot", False)

"""  vizualise soil  """
grid = vp.uniform_grid(min_, max_, res_)  # for visualization
meshActor, meshCBar = vp.plot_mesh(grid, "", "", False)
vp.render_window([meshActor[0], rootActor], "Test mapping", rootCBar, grid.GetBounds()).Start()

