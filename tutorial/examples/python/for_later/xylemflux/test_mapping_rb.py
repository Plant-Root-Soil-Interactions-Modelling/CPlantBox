""" 
Tests if MappedSegments works for PlantBox.RootSystem under a periodic setting  
"""
import sys
sys.path.append("../../../..")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
import vtk_plot as vp

from math import *
import numpy as np


def vector_3d(a):
    return pb.Vector3d(a[0], a[1], a[2])


# grid to test cutting
min_ = np.array([-2, -2, -20])
max_ = np.array([2, 2, -3])

""" pick one of the resolutions """
# res_ = np.array([1, 1, 10])
res_ = np.array([1, 5, 1])
# res_ = np.array([10, 1, 1])


def periodic_soil_index(x, y, z):
    """ not in C++, because this should be handeld by mesh.pick(x,y,z) implemented for dumux in solverbase.cc """
    pos = np.array([x, y, z])
    for i in range(0, 2):  # periodic mapping for x and y, not z
        minx = min_[i]
        xx = max_[i] - minx
        if not isinf(xx):  # periodic in x
            pos[i] -= minx;  # start at 0
            if pos[i] >= 0:
                pos[i] = pos[i] - int(pos[i] / xx) * xx;
            else:
                pos[i] = pos[i] + int((xx - pos[i]) / xx) * xx;
            pos[i] += minx;
    # map to linear index
    w = max_ - min_
    p0 = pos - min_
    i = np.array([p0[0] / w[0] * res_[0], p0[1] / w[1] * res_[1], p0[2] / w[2] * res_[2] ])
    for k in range(0, 3):
        if ((i[k] < 0) or (i[k] >= res_[k])):
            return -1;

    return floor(i[0]) * res_[1] * res_[2] + floor(i[1]) * res_[2] + floor(i[2])  # a linear index not periodic


rs = pb.MappedRootSystem()
path = "../../../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")
for p in rs.getRootRandomParameter():  # Modify axial resolution
    p.dx = 5  # [cm] adjust resolution
rs.initialize()
rs.simulate(20)  # days

rs.setSoilGrid(periodic_soil_index)
rs.setRectangularGrid(vector_3d(min_), vector_3d(max_), vector_3d(res_))  # cut and map segments
rs.sort()

""" Assign soil cell index to root segments """
segs = rs.segments

for i, s in enumerate(segs):
    if s.y - 1 != i :
        print("ohno!!!!")
        input()

nodes = rs.nodes
x = np.zeros(len(segs))
for c, s in enumerate(segs):
#     i, j = s.x, s.y
#     n1 = np.array(nodes[i])
#     n2 = np.array(nodes[j])
#     mid = 0.5 * (n1 + n2)
#     x[c] = rs.soil_index(mid[0], mid[1], mid[2])
    try:
        x[c] = rs.seg2cell[s.y - 1]
    except:
        x[c] = -1

y = np.zeros(x.shape[0])
for i in range(0, x.shape[0]):
   y[i] = i % 2 if x[i] >= 0 else i % 2 - 2  #

for i in range(0, x.shape[0]):
   x[i] = x[i] if x[i] >= 0 else -1

print("RS number of segments (mapped)", len(rs.segments), "(rs)", rs.getNumberOfSegments())  # shoot segments are mapped
ana = pb.SegmentAnalyser(rs.mappedSegments())  # rs is MappedSegments and RootSystem. So passing rs it is not unique which constructor is called.
print("Number of segments", len(ana.segments), "x", len(x), "and", len(y))
ana.addData("linear_index", x)  # node data are converted to segment data
ana.addData("zebra", y)
ana.mapPeriodic(max_[0] - min_[0], max_[1] - min_[1])  # data are also split

pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "linear_index", "zebra"])
rootActor, rootCBar = vp.plot_roots(pd, "linear_index", False)

""" Mesh """
width_ = max_ - min_
ind_ = res_ / width_
print("min    ", min_)
print("max    ", max_)
print("width  ", width_)
print("cuboids", 1 / ind_)
grid = vp.uniform_grid(min_, max_, res_)
meshActor, meshCBar = vp.plot_mesh(grid, "", "", False)

# vp.render_window([meshActor, rootActor], "Test mapping", rootCBar).Start()
grid = vp.uniform_grid(min_, max_, res_)
actors = vp.plot_mesh_cuts(pd, "linear_index")
print(actors)
print(len(actors))

# actors.extend([rootActor])
vp.render_window(actors, "Test mapping", rootCBar).Start()

