"""root system length over time"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import visualisation.vtk_plot as vp
import rsml.rsml_reader as rsml
import numpy as np
import matplotlib.pyplot as plt
import os
failed = []
for year in [1,2,48]: #,
    path = "./rsml/RSML_year"+str(year)+"/"
    files = os.listdir(path)
    for name in files:
        print(name)
        name = name[:-5]
        try:
            r = XylemFluxPython(path + name + ".rsml")
            ana = pb.SegmentAnalyser(r.rs)
            os.makedirs("./results/RSML_year"+str(year), exist_ok=True)
            vp.plot_roots(ana, "subType", win_title = "RSML_year"+str(year)+"/" +name, render = False)
        except:
            failed.append(name)

print('failed',failed)