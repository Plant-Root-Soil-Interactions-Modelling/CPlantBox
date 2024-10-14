"""root system length over time"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import visualisation.vtk_plot as vp
import rsml.rsml_reader as rsml

import numpy as np
import matplotlib.pyplot as plt

""" root problem """
r = XylemFluxPython("../../../dumux-rosi/grids/RootSystem.rsml")
r.test()

polylines, props, funcs, metadata = rsml.read_rsml("../../../dumux-rosi/grids/RootSystem.rsml")
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
print("surface area", np.sqrt(ana.getSummed("volume") / ana.getSummed("length") / np.pi) * 2 * np.pi * ana.getSummed("length"))

