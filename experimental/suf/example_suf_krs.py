""" 
soil uptake fraction of a root system (soil is in hydrostatic equilibrium) 
"""
import sys; sys.path.append("../../python/modules/"); sys.path.append("../../../CPlantBox/"); 
sys.path.append("../../../CPlantBox/src/python_modules")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """
kz = 4.32e-2  # axial conductivity [cm^3/day]
kr = 1.728e-4  # radial conductivity [1/day]
simtime = 20  # [day] for task b

""" root system """
rs = pb.MappedRootSystem()

path = "../../../CPlantBox//modelparameter/rootsystem/"
name = "Glycine_max"  
rs.setSeed(1)
rs.readParameters(path + name + ".xml")
rs.getRootSystemParameter().seedPos.z = -0.1
rs.initialize()
rs.simulate(simtime, False)

print()
print("Shoot segments: ", [str(s) for s in rs.getShootSegments()])
print()

""" set up xylem parameters """
r = XylemFluxPython(rs)
r.setKr([kr])  # or use setKrTables, see XylemFlux.h
r.setKx([kz])
r.test()

""" numerical solution of transpiration -1 cm3/day"""
suf = r.get_suf(0.)
print("Sum of SUF", np.sum(suf), "from", np.min(suf), "to", np.max(suf), "summed positive", np.sum(suf[suf >= 0]))
krs = r.get_krs(0.)
print("Krs: ", krs)

""" Additional vtk plot """
ana = pb.SegmentAnalyser(r.rs)
ana.addData("SUF", suf)  # cut off for vizualisation
vp.plot_roots(ana, "SUF", "Soil uptake fraction (cm3 day)")  # "fluxes"
