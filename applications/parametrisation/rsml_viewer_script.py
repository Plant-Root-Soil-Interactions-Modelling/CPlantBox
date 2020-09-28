import sys;  sys.path.append("../..")
""" determines the (more easy) parameters la, lb, ln, a, theta by order and prints mean and std"""

import numpy as np
import matplotlib.pyplot as plt

import rsml_reader as rsml
import estimate_root_params as es
from xylem_flux import XylemFluxPython
import vtk_plot as vp
import plantbox as pb

time = range(1, 3)  # measurement times (not in the rsml)
name = ["RSML/m1/dicot/lupin/lupin_d{:g}.rsml".format(a) for a in range(1, 10)]

# time = range(1, 3)  # measurement times (not in the rsml)
# name = ["RSML/m1/monocot/maize/PL0{:g}_DAS0{:g}.rsml".format(1, a) for a in range(1, 8)]

# time = [75]  # measurement times (not in the rsml)
# name = ["RSML/Maize_Kutschera.rsml"]

rs = XylemFluxPython(name[0])  # parses rsml, XylemFluxPython.rs is of type MappedRootSegments
ana = pb.SegmentAnalyser(rs.rs)  # convert MappedRootSegments to a SegmentAnalyser

# radii = ana.data["radius"]  # DOES NOT WORK (why?)
# for i in range(0, len(radii)):
#     radii[i] /= 116.93
# ana.data["radius"] = radii

pd = vp.segs_to_polydata(ana, 1. / 116.93)  # makes a vtkPolydata (to save as .vtp, or visualize with vtk)
vp.plot_roots(pd, "radius")  # plots vtkPolydata into an interactive window
