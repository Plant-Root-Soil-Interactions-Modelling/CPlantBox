import py_rootbox as rb
import math
import numpy as np
from rb_tools import *

rs = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010"
rs.openFile(name)

scale_elongation = rb.EquidistantGrid1D(0, -50, 100)  # for root elongation from 0 cm to -50 cm, 100 layers

soil_strength = np.ones((99,)) * 0.5  # some data

scales = np.exp(-0.4 * soil_strength)  # scales from some equation (TODO)

scale_elongation.data = a2v(scales)  # set proportionality factors

print("value at -3 cm", scale_elongation.getValue(rb.Vector3d(0, 0, -3)))

# Manually set scale elongation function
for p in rs.getRootTypeParameter():
    p.se = scale_elongation

# Simulation
rs.initialize()
simtime = 120.
dt = 1.
N = 120 / dt
for i in range(0, round(N)):

    # update soil model (e.g. soil_strength)

    # update scales (e.g. from water content, soil_strength)
    scales = np.exp(-0.4 * soil_strength)  # scales from some equation (TODO)

    # copy scales into scaling funciton
    scale_elongation.data = a2v(scales)

    rs.simulate(dt, True)

rs.write("../results/example_4c_2.vtp")

print("done.")
