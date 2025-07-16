""" map root segments to a soil grid """
import sys; sys.path.append("../.."); sys.path.append("../../src/");
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/"); sys.path.append("../../../dumux-rosi/python/modules/");

import plantbox as pb
import visualisation.vtk_plot as vp
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import numpy as np

""" Root system """  # |\label{l61m:root_system_start}|
plant = pb.MappedPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
plant.readParameters(path + name + ".xml")
plant.initialize()  # |\label{l61m:root_system_end}|

""" Macroscopic soil grid """  # |\label{l61m:grid_start}|
min_b = np.array([-2, -2, -15])  # [cm]
max_b = np.array([3, 4, -5])  # [cm]
cell_number = np.array([3, 4, 11])  # [1]
periodic = False
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)
s.setVGParameters([[0.08, 0.43, 0.04, 1.6, 50]])
s.setHomogeneousIC(-300, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.initializeProblem()  # |\label{l61m:grid_start}|

""" Coupling """
picker = lambda x, y, z: s.pick([x, y, z])
plant.setSoilGrid(picker)

""" Simulate """
plant.simulate(10., False)

""" find grid cell index for segment """
segs = plant.segments  # |\label{l61m:segments}|
x = np.array([plant.seg2cell[i] for i in range(0, len(segs))])

""" find segment indices in a grid cell"""
ci = plant.soil_index(0, -3, -7)  # grid cell index
print("Cell at [0,0,-7] has index", ci)
try:
    print(len(plant.cell2seg[ci]), "segments in this cell:")
    print(plant.cell2seg[ci])
except:
    print("There are no segments in this cell")

""" Visualize"""
ana = pb.SegmentAnalyser(plant.mappedSegments())  # ana = pb.SegmentAnalyser(rs)  # <---- wrong!
ana.addData("linear_index", x)
pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "linear_index"])
rootActor, rootCBar = vp.plot_roots(pd, "linear_index", "segment index plot", False)
grid = vp.uniform_grid(min_b, max_b, cell_number)  # for visualization
meshActor, meshCBar = vp.plot_mesh(grid, "", "", False)
vp.render_window([meshActor[0], rootActor], "Test mapping", rootCBar, grid.GetBounds()).Start()

