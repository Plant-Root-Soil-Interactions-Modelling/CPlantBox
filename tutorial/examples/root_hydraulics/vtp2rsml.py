import matplotlib.pyplot as plt
import numpy as np
import plantbox as pb
import plantbox.rsml.rsml_writer as rsml
import plantbox.visualisation.vtk_plot as vp
import plantbox.visualisation.vtk_tools as vt
from plantbox.functional.PlantHydraulicModel import (
    HydraulicModel_Doussan,
    HydraulicModel_Meunier,
)

""" 
Converts a DuMux output vtp to a RSML
"""

file_in = "../../grids/RootSystem8.vtp"
file_out = "../../grids/RootSystem8_vtp.rsml"

""" read vtp """
pd = vt.read_vtp(file_in)

meta = rsml.Metadata()
meta.unit = "m"
meta.add_property(rsml.Property("radius [m]", "float", "m", None))

order_id = 4

vt.write_rsml(file_out, pd, order_id, meta)  # meta is optional now

print("fin")
