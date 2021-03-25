import sys; sys.path.append("../../python/modules/"); sys.path.append("../../../CPlantBox/")
sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../build-cmake/cpp/python_binding/")

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import vtk_plot as vp
import vtk_tools as vt

import numpy as np

"""
1. plot DuMux .vtu output

2. pull out data from vtu

3. make own grid, put in data, plot again

TODO currently, the DuMux grid of the dumux-rosi binding is used, a CPlantBox grid implementation is missing (pick for SoilLookUp)
"""

name = "soybean_Honly-00001"

""" 1. """
# Open .vtu
pd = vp.read_vtu(name + ".vtu")
print(pd.GetBounds())  # xmin, xmax, ymin, ymax, zmin, zmax
print("Number of cells", vt.np_cells(pd).shape[0])

# Convert m -> cm
np_points = vt.np_points(pd)
points = vt.vtk_points(np_points * 100)  # m -> cm
pd.SetPoints(points)

# Plot
# vp.plot_mesh(pd, "pressure head") # useless plot
# vp.plot_mesh_cuts(pd, "pressure head", 5)

""" 2. """
data, _ = vt.np_data(pd, 9, True)  # grid, data_index, cell data
print("Data range from {:g} to {:g}".format(np.min(data), np.max(data)))

""" 3. """
min_ = np.array([-18.5, -3, -150])
max_ = np.array([18.5, 3, 0.])
res_ = np.array([18, 3, 75])
periodic = True

s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_, max_, res_, periodic)  # [cm]
loam = [0.08, 0.43, 0.04, 1.6, 50]  # we do not plan to calculate a thing, but we need parameters for initialisation
s.setVGParameters([loam])
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.initializeProblem()
s.setInitialCondition(data)  # put data to the grid

pd = vp.solver_to_polydata(s, min_, max_, res_)
vp.plot_mesh_cuts(pd, "pressure head", 5)  # should look like the first plot

print("fin")
