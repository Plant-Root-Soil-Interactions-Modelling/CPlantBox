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

""" root problem """
r = XylemFluxPython("RootSystem.rsml")
r.test()

polylines, props, funcs = rsml.read_rsml("RootSystem.rsml")
print(len(polylines), "roots")
# print(props["parent-poly"])
# print(props["parent-node"])
print("Tap root number of nodes", len(polylines[0]))
for i in range(1, len(polylines)):
    pp = int(props["parent-poly"][i])
    pn = int(props["parent-node"][i])
    try:
        n0 = polylines[pp][pn]
        n1 = polylines[i][0]
        print("root", i, np.linalg.norm(np.array(n1) - np.array(n0)))
    except:
        print("root", i, "index out of range", "parent-node", pn, "length", len(polylines[pp]))

ana = pb.SegmentAnalyser(r.rs)
print("surface area", math.sqrt(ana.getSummed("volume")/ana.getSummed("length")/math.pi)*2*math.pi*ana.getSummed("length"))

