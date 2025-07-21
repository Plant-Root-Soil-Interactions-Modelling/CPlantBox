""" map root segments to a soil grid """
import sys; sys.path.append("../.."); sys.path.append("../../src/");
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/"); sys.path.append("../../../dumux-rosi/python/modules/");

import plantbox as pb
import visualisation.vtk_plot as vp
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import numpy as np

""" Root system """  # |\label{l71m:root_system_start}|
plant = pb.MappedPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010, Anagallis_femina_Leitner_2010
plant.readParameters(path + name + ".xml")
plant.setSeed(4)
plant.initialize()  # |\label{l71m:root_system_end}|

""" Macroscopic soil grid """  # |\label{l71m:grid_start}|
min_b = np.array([-2, -2, -15])  # [cm]
max_b = np.array([2, 2, -5])  # [cm]
cell_number = np.array([2, 3, 6])  # [1]
s = RichardsWrapper(RichardsSP())
s.initialize()
periodic = True
s.createGrid(min_b, max_b, cell_number, periodic)
s.setVGParameters([[0.08, 0.43, 0.04, 1.6, 50]])
s.setHomogeneousIC(-300, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.initializeProblem()  # |\label{l71m:grid_end}|

""" Coupling """
picker = lambda x, y, z: s.pick([x, y, z])  # |\label{l71m:picker}|
plant.setSoilGrid(picker)  # |\label{l71m:picker_end}|

""" Simulate """
plant.simulate(10., False)  # |\label{l71m:simulate}|

""" Find segment indices in a grid cell"""
ci = plant.soil_index(0, 0, -7)  # grid cell index |\label{l71m:soil_index}|
print("Cell at [0,0,-7] has index", ci)
try:
    print(len(plant.cell2seg[ci]), "segments in this cell:")  # |\label{l71m:cell2seg}|
    print(plant.cell2seg[ci])
except:
    print("There are no segments in this cell")

""" Find grid cell index for segment """
segs = plant.segments  # |\label{l71m:segments}|
x = np.array([plant.seg2cell[i] for i in range(0, len(segs))])  # |\label{l71m:seg2cell}|

""" Visualize"""  # |\label{l71m:visualize}|
ana = pb.SegmentAnalyser(plant.mappedSegments())  # |\label{l71m:mappedSegments}|
ana.addData("linear_index", x)
pd = vp.segs_to_polydata(ana, 1., ["radius", "linear_index"])
rootActor, rootCBar = vp.plot_roots(pd, "linear_index", "Segment index", False)
grid = vp.uniform_grid(min_b, max_b, cell_number)
meshActor, meshCBar = vp.plot_mesh(grid, "", "", False)
vp.render_window([meshActor[0], rootActor], "Test mapping", rootCBar, grid.GetBounds()).Start()

