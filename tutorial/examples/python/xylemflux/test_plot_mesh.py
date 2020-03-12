import sys
sys.path.append("../../../..")

# from xylem_flux import XylemFluxPython  # Python hybrid solver
# import plantbox as pb
# import rsml_reader as rsml
import vtk_plot as vp

# from math import *
# import numpy as np
# import matplotlib.pyplot as plt

#
# def vector_3d(a):
#     return pb.Vector3d(a[0], a[1], a[2])

""" 
Tests if we manage to plot a mesh
"""

ug = vp.read_vtu("benchmark3d_2-00001.vtu")

p_name = "water content"  # "pressure head, "S_liq" "water content"
meshActor, scalarBar = vp.plot_mesh_wireframe(ug, p_name, False)

# vp.plot_mesh_cuts(ug, p_name)

renWin = vp.render_window([meshActor], "Wireframe mesh", scalarBar)

vp.write_png(renWin, "test_plot_mesh")
