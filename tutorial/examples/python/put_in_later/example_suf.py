""" 
soil uptake fraction of a root system (soil is in hydrostatic equilibrium) 
"""
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """
kz = 4.32e-2  # axial conductivity [cm^3/day]
kr = 1.728e-4  # radial conductivity [1/day]
simtime = 14  # [day] for task b

""" root system """
rs = pb.MappedRootSystem()
p_s = np.linspace(-200, -400, 2001)  # 2 meter down, from -200 to -400, resolution in mm
soil_index = lambda x, y, z : int(-10*z) # maps to p_s (hydrostatic equilibirum)
rs.setSoilGrid(soil_index)

path = "../../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.setSeed(1)
rs.readParameters(path + name + ".xml")
rs.initialize()
rs.simulate(simtime, False)
p_s = np.linspace(-200, -400, 2001)  # 2 meter down, resolution in mm

""" set up xylem parameters """
r = XylemFluxPython(rs)
r.setKr([kr]) # or use setKrTables, see XylemFlux.h
r.setKx([kz])

""" numerical solution of transpiration -1 cm3/day"""
rx = r.solve_neumann(simtime, -1, p_s, True) # True: matric potential given per cell (not per segment)
print("solved")

fluxes = r.segFluxes(simtime, rx, p_s, False, True)  # cm3/day (double simTime,  rx,  sx,  approx, cells
print("Transpiration", r.collar_flux(simtime, rx, p_s), np.sum(fluxes), "cm3/day")
suf = np.array(fluxes) / -1. # [1]

""" Additional vtk plot """
ana = pb.SegmentAnalyser(r.rs)
ana.addData("SUF", np.minimum(suf, 1.e-2))  # cut off for vizualisation
vp.plot_roots(ana, "SUF", "Soil uptake fraction (cm3 day)")  # "fluxes"
