import sys; sys.path.append("../../..");  sys.path.append("../../../src")

import plantbox as pb
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import rsml.rsml_reader as rsml
import visualisation.vtk_tools as vt
import visualisation.vtk_plot as vp

import matplotlib.pyplot as plt
import numpy as np

""" 
opens and vtk plots the rsml
"""

file_name = "RootSys_verysimple.rsml"

""" root problem """
r = XylemFluxPython(file_name)
r.test()

polylines, props, funcs, metadata = rsml.read_rsml(file_name)
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

""" Additional vtk plot """
ana = pb.SegmentAnalyser(r.rs)
pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "length"])
vp.plot_roots(pd, "subType")

