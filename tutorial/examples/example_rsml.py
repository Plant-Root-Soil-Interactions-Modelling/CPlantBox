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

ana = pb.SegmentAnalyser(r.rs)
