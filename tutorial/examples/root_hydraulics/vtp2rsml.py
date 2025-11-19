import sys; sys.path.append("../../..");  sys.path.append("../../../src")

import plantbox as pb
from plantbox.functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox.rsml.rsml_writer as rsml
import plantbox.visualisation.vtk_tools as vt
import plantbox.visualisation.vtk_plot as vp

import matplotlib.pyplot as plt
import numpy as np

""" 
Converts a DuMux output vtp to a RSML
"""

file_in = "test.vtp"  # ../../grids/RootSystem8.vtp"
file_out = "test.rsml"  # "../../grids/RootSystem8.rsml"

""" read vtp """
pd = vt.read_vtp(file_in)

meta = rsml.Metadata()
meta.unit = "m"
meta.add_property(rsml.Property("radius [m]", "float", "m", None))

order_id = 4

vt.write_rsml(file_out, pd, order_id, meta)  # meta is optional now

print("fin")
