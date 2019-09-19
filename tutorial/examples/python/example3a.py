"""root system length over time, root tip distriubtion"""
import py_rootbox as rb
from rb_tools import *
import numpy as np
import matplotlib.pyplot as plt

rs = rb.RootSystem()
name = "Brassica_oleracea_Vansteenkiste_2014"
rs.readParameters("modelparameter/" + name + ".xml")
rs.initialize()

simtime = 60.  # days
dt = 1.
N = round(simtime / dt)  # steps

# Plot some scalar value over time
stype = "length"
stype_str = "length (cm)"
v_ = np.zeros(N)
v1_ = np.zeros(N)
v2_ = np.zeros(N)
v3_ = np.zeros(N)
for i in range(0, N):
    rs.simulate(dt)
    t = v2a(rs.getParameter("type"))
    v = v2a(rs.getParameter(stype))
    v_[i] = np.sum(v)
    v1_[i] = np.sum(v[t == 1])
    v2_[i] = np.sum(v[t == 2])
    v3_[i] = np.sum(v[t == 3])

t_ = np.linspace(dt, N * dt, N)
plt.plot(t_, v_)
plt.plot(t_, v1_)
plt.plot(t_, v2_)
plt.plot(t_, v3_)
plt.xlabel("time (days)")
plt.ylabel(stype_str)
plt.legend(["total", "tap root", "lateral", "2. order lateral"])
plt.savefig("../results/example_3a.png")
plt.show()

# Find root tips and bases (two approaches)
rs.initialize()  # <------------
rs.simulate(7)  # 7 days young....

print(rs.getNumberOfNodes(), "nodes")
print(rs.getNumberOfSegments(), "segments")

# Use polyline representation of the roots
polylines = rs.getPolylines()
bases = np.zeros((len(polylines), 3))
tips = np.zeros((len(polylines), 3))
for i, r in enumerate(polylines):
    bases[i, :] = [r[0].x, r[0].y, r[0].z]
    tips[i, :] = [r[-1].x, r[-1].y, r[-1].z]

# Or, use node indices to find tip or base nodes
nodes = vv2a(rs.getNodes())
tipI = rs.getRootTips()
baseI = rs.getRootBases()

# Plot results (1st approach)
plt.title("Top view")
plt.xlabel("cm")
plt.ylabel("cm")
plt.scatter(nodes[baseI, 0], nodes[baseI, 1], c = "g", label = "root bases")
plt.scatter(nodes[tipI, 0], nodes[tipI, 1], c = "r", label = "root tips")
plt.savefig("../results/example_3a2.png")
plt.show()

 # check if the two approaches yield the same result
uneq = np.sum(nodes[baseI, :] != bases) + np.sum(nodes[tipI, :] != tips)
print("Unequal tips and basals:", uneq)

print("done.")
