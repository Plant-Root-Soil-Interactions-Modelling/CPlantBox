"""map root segments to a soil grid"""

import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from rosi.richards import RichardsWrapper  # Python part
from rosi.rosi_richards import RichardsSP  # C++ part (Dumux binding)

# Root system |\label{l71m:root_system_start}|
plant = pb.MappedPlant()
path = "../../modelparameter/structural/rootsystem/"
filename = "Anagallis_femina_Leitner_2010"
plant.readParameters(path + filename + ".xml")
plant.setSeed(4)  # |\label{l71m:random}|
plant.initialize()  # |\label{l71m:root_system_end}|

# Macroscopic soil grid |\label{l71m:grid_start}|
box_min = np.array([-2, -2, -15])  # cm
box_max = np.array([2, 2, -5])  # cm
cell_number = np.array([2, 3, 6])  # 1
s = RichardsWrapper(RichardsSP())
s.initialize()
periodic = True  # |\label{l71m:periodic}|
s.createGrid(box_min, box_max, cell_number, periodic)
s.setVGParameters([[0.08, 0.43, 0.04, 1.6, 50]])
s.setHomogeneousIC(-300, True)  # cm pressure head
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.initializeProblem()  # |\label{l71m:grid_end}|


# Coupling
def picker(x, y, z):
    """soil grid cell index for position (x, y, z)"""
    return s.pick([x, y, z])  # |\label{l71m:picker}|


plant.setSoilGrid(picker)  # |\label{l71m:picker_end}|

# Simulate
plant.simulate(10.0, False)  # |\label{l71m:simulate}|

# Find segment indices in a grid cell
ci = plant.soil_index(0, 0, -7)  # grid cell index |\label{l71m:soil_index}|
print("Cell at [0,0,-7] has index", ci)
try:
    print(len(plant.cell2seg[ci]), "segments in this cell:")  # |\label{l71m:cell2seg}|
    print(plant.cell2seg[ci])
except Exception:
    print("There are no segments in this cell")

# Find grid cell index for segment
segs = plant.segments  # |\label{l71m:segments}|
x = np.array([plant.seg2cell[i] for i in range(0, len(segs))])  # |\label{l71m:seg2cell}|

# Visualize |\label{l71m:visualize}|
ana = pb.SegmentAnalyser(plant.mappedSegments())  # |\label{l71m:mappedSegments}|
ana.addData("linear_index", x)
pd = vp.segs_to_polydata(ana, 1.0, ["radius", "linear_index"])
rootActor, rootCBar = vp.plot_roots(pd, "linear_index", "Segment index", False)
grid = vp.uniform_grid(box_min, box_max, cell_number)
meshActor, meshCBar = vp.plot_mesh(grid, "", "", False)
vp.render_window([meshActor[0], rootActor], "Test mapping", rootCBar, grid.GetBounds()).Start()
