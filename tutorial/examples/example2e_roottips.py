"""find root tips and bases (two approaches)"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb

import numpy as np
import matplotlib.pyplot as plt

path = "../../modelparameter/structural/rootsystem/"
name = "Brassica_napus_a_Leitner_2010"

rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")
rs.initialize()
rs.simulate(21, True)

print(rs.getNumberOfNodes(), "nodes")
print(rs.getNumberOfSegments(), "segments")

# Use polyline representation of the roots
polylines = rs.getPolylines()
bases = np.zeros((len(polylines), 3))
tips = np.zeros((len(polylines), 3))
for i, r in enumerate(polylines):
    bases[i,:] = [r[0].x, r[0].y, r[0].z]  # first index is the base
    tips[i,:] = [r[-1].x, r[-1].y, r[-1].z]  # last index is the tip

# Or, use node indices to find tip or base nodes
nodes = np.array((list(map(np.array, rs.getNodes()))))
tipI = rs.getRootTips()
baseI = rs.getRootBases()

# Plot results (1st approach)
plt.title("Top view")
plt.xlabel("cm")
plt.ylabel("cm")
plt.scatter(nodes[baseI, 0], nodes[baseI, 1], c = "g", label = "root bases")
plt.scatter(nodes[tipI, 0], nodes[tipI, 1], c = "r", label = "root tips")
plt.legend()
plt.savefig("results/example_2e.png")
plt.show()

 # check if the two approaches yield the same result
uneq = np.sum(nodes[baseI,:] != bases) + np.sum(nodes[tipI,:] != tips)
print("Unequal tips and basals:", uneq)
