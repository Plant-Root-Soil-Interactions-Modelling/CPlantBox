""" map root segments to a soil grid """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np

""" root system """
rs = pb.MappedRootSystem()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml")
rs.initialize()
rs.simulate(10., False)

""" soil """
min_ = np.array([-5, -5, -15])
max_ = np.array([9, 4, -5])
res_ = np.array([3, 1, 5])
rs.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments

""" add segment indices """
segs = rs.segments
x = np.zeros(len(segs))
for i, s in enumerate(segs):
    try:
        x[i] = rs.seg2cell[i]
    except:  # in case the segment is not within the domain
        x[i] = -1

""" infos on a specific cell"""
ci = rs.soil_index(0, 0, -7)
print("Cell at [0,0,-7] has index", ci)
try:
    print(len(rs.cell2seg[ci]), "segments in this cell:")
    print(rs.cell2seg[ci])
except:
    print("There are no segments in this cell")

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

