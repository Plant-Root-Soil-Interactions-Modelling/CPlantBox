"""root system length over time"""
import sys
sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
import plantbox as pb

import numpy as np
import matplotlib.pyplot as plt
import rsml_reader as rsml
import math
from xylem_flux import XylemFluxPython  # Python hybrid solver
import vtk_plot as vp

""" Read RSML file and test whether it is OK """
r = XylemFluxPython("RootSystem.rsml")
r.test()  # here you could add the addition of artificial shoot and creation time (using the viewer as template)

"""Retrieve information from the RSML file"""
polylines, props, functions, metadata = rsml.read_rsml("RootSystem.rsml")
print(len(polylines), "roots")
# print(props["parent-poly"])
# print(props["parent-node"])
print("Tap root number of nodes", len(polylines[0]))

""" compute Krs and SUF """
kz = 4.32e-2  # axial conductivity [cm^3/day]
kr = 1.728e-4  # radial conductivity [1/day]
r.setKr([kr])  # or use setKrTables according to any function, see XylemFlux.h
r.setKx([kz])

suf = r.get_suf(0.)
print("Sum of SUF", np.sum(suf), "from", np.min(suf), "to", np.max(suf), "summed positive", np.sum(suf[suf >= 0]))
krs = r.get_krs(0.)
print("Krs: ", krs)

""" eswp for specific soil scenario """
p_s = np.linspace(-500, -100, 300)   # 3 m down, resolution in cm
soil_index = lambda x, y, z : int(-1* z)
r.rs.setSoilGrid(soil_index)
eswp = r.get_eswp(0.,p_s)
print("eswp: ",eswp)

""" vtk plot """
ana = pb.SegmentAnalyser(r.rs)
ana.addData("suf", suf) 
#vp.plot_roots(ana, "suf", "Soil uptake fraction (cm3 day)")  # "fluxes"


