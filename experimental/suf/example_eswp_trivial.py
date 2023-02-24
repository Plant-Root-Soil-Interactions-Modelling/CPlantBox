import sys; sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src/python_modules")
sys.path.append("../../python/modules/")
sys.path.append("../../build-cmake/cpp/python_binding/")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
import vtk_plot as vp

from math import *
import numpy as np
import matplotlib.pyplot as plt

""" 
ESWP can be more negative than the root collar potential (?)
"""

""" 1. SUF """

""" Parameters """
kz = 4.32e-2  # axial conductivity [cm^3/day]
kr = 1.728e-4  # radial conductivity [1/day]
simtime = 154  # [day] 

""" root system """
rs = pb.MappedRootSystem()
p_s = np.linspace(-500, -200, 3001)  # 3 meter down, from -200 to -500, resolution in mm
soil_index = lambda x, y, z : int(-10 * z)  # maps to p_s 
rs.setSoilGrid(soil_index)

soilbox = pb.SDF_PlantBox(1e6, 1e6, 300)
rs.setGeometry(soilbox)

path = "../../../CPlantBox//modelparameter/rootsystem/"
name = "Glycine_max"  
rs.setSeed(1)
rs.readParameters(path + name + ".xml")
rs.getRootSystemParameter().seedPos.z = -0.1
rs.initialize()
rs.simulate(simtime, False)

""" set up xylem parameters """
r = XylemFluxPython(rs)
r.setKr([kr])  # or use setKrTables, see XylemFlux.h
r.setKx([kz])

""" numerical solution of transpiration -10000 cm3/day - very high to avoid numerical errors"""
rx = r.solve_neumann(simtime, -10000, p_s, True)  # True: matric potential given per cell (not per segment)
print("solved")

fluxes = r.segFluxes(simtime, rx, p_s, False, True)  # cm3/day (double simTime,  rx,  sx,  approx, cells
print("Transpiration", r.collar_flux(simtime, rx, p_s), np.sum(fluxes), "cm3/day")
suf = np.array(fluxes) / -10000.  # [1]

# # """ vtk plot """
# ana = pb.SegmentAnalyser(r.rs)
# ana.addData("SUF", suf) # cut off for vizualisation np.minimum(suf, 1.e-4)
# vp.plot_roots(ana, "SUF", "Soil uptake fraction (cm3 day)")  # "fluxes"

print("Sum of suf", np.sum(suf), "from", np.min(suf), "to", np.max(suf))

""" 2. Equivalent soil water potential """
#nodes = r.get_nodes()
p_s = np.linspace(-14000, -10000, 3001)
p_s[0:10]=-300;  # 3 meter down, resolution in mm, dry with moist top
soil_index = lambda x, y, z : int(-10 * z)
r.rs.setSoilGrid(soil_index)

""" Numerical solution"""
rx = r.solve_neumann(simtime, -0., p_s, True)  # True: matric potential given per cell (not per segment)

print("Transpiration", r.collar_flux(0., rx, [p_s]), "cm3/day")
eswp = 0.
n = len(r.rs.segments)
seg2cell = r.rs.seg2cell
for i in range(0, n):
    eswp += suf[i] * p_s[seg2cell[i]]

print()
print("Equivalent soil water potential", eswp)

z = r.rs.nodes[1].z
print("Root collar potential          ", rx[1], "node z", z, "psi", p_s[soil_index(0,0,z)])
print()

""" Additional vtk plot """
ana = pb.SegmentAnalyser(r.rs)
ana.addData("rx", rx)  # node data are converted to segment data
pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "length", "rx","suf"])
vp.plot_roots(pd, "rx")

