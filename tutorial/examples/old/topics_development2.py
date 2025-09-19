"""how to find root tips and bases"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import numpy as np
import matplotlib.pyplot as plt

path = "../../modelparameter/structural/rootsystem/"
name = "Brassica_napus_a_Leitner_2010"

plant = pb.RootSystem()
plant.readParameters(path + name + ".xml")
plant.initialize()

simtime = 30
dt = 0.5
N = round(simtime / dt)

bases, tips = [], []

for i in range(0, N):

    plant.simulate(dt)

    polylines = plant.getPolylines()  # polyline representation
    for i, r in enumerate(polylines):
        tips.append([r[-1].x, r[-1].y, r[-1].z])  # last index is the tip
        bases.append([r[0].x, r[0].y, r[0].z])  # first index is the base

bases = np.array(bases)  # convert to numpy arrays
tips = np.array(tips)

# Plot results
plt.title("Top view")
plt.xlabel("cm")
plt.ylabel("cm")
plt.scatter(tips[:, 0], tips[:, 1], c = "g", label = "organ tips", alpha = 0.2)
plt.scatter(bases[:, 0], bases[:, 1], c = "r", label = "organ bases", alpha = 0.2)
plt.legend()
plt.savefig("results/topics_development2.png")
plt.show()
