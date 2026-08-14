import matplotlib.pyplot as plt
import numpy as np
import plantbox as pb
import plantbox.rsml.rsml_reader as rsml
import plantbox.visualisation.vtk_plot as vp
import plantbox.visualisation.vtk_tools as vt
from plantbox.functional.PlantHydraulicModel import (
    HydraulicModel_Doussan,
    HydraulicModel_Meunier,
)

""" 
opens and vtk plots the rsml
"""

file_name = "../../grids/RootSys_verysimple.rsml"

r = HydraulicModel_Doussan(file_name, None, cached=False)  # or HydraulicModel_Doussan
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
ana = pb.SegmentAnalyser(r.ms)
pd = vp.segs_to_polydata(ana, 1.0, ["radius", "subType", "creationTime", "length"])
vp.plot_roots(pd, "subType")
