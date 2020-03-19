"""Scale root elongation based on SoilLookup"""  # todo improve example
import sys
sys.path.append("../../..")
import plantbox as pb
import math
import numpy as np

rs = pb.RootSystem()
path = "../../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")

scale_elongation = pb.EquidistantGrid1D(0, -50, 100)  # for root elongation from 0 cm to -50 cm, 100 layers
soil_strength = np.ones((99,)) * 0.5  # some data
scales = np.exp(-0.4 * soil_strength)  # scales from some equation (TODO)
scale_elongation.data = scales  # set proportionality factors

print("value at -3 cm", scale_elongation.getValue(pb.Vector3d(0, 0, -3)))

# Manually set scale elongation function
for p in rs.getRootRandomParameter():
    p.f_se = scale_elongation

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
    scale_elongation.data = scales

    rs.simulate(dt, True)

rs.write("../results/example_6b.vtp")

