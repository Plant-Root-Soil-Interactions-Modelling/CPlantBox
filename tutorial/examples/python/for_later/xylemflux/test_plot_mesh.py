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
vp.plot_mesh(ug, "water content", "Mesh")  # "pressure head, "S_liq" "water content"
# vp.plot_mesh_cuts(ug, p_name)

# p         write png (added it)
# s         as solid (built in default)
# w         as wireframe (built in default)
