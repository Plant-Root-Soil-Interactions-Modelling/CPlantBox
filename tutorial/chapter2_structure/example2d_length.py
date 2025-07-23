"""root system length over time"""
<<<<<<< HEAD
import sys; sys.path.append("../.."); sys.path.append("../../src/")
=======
import sys; sys.path.append("../.."); sys.path.append("../../src/") # |\label{l2_2d:importStart}|
>>>>>>> master

import plantbox as pb

import numpy as np
<<<<<<< HEAD
import matplotlib.pyplot as plt

path = "../../modelparameter/structural/rootsystem/"
=======
import matplotlib.pyplot as plt # |\label{l2_2d:importEnd}|

path = "../../modelparameter/structural/rootsystem/"  # |\label{l2_2d:defineStart}|
>>>>>>> master
name = "Brassica_napus_a_Leitner_2010"

rs = pb.Plant()
rs.readParameters(path + name + ".xml")
rs.initialize()

simtime = 60.  # days
dt = 1.
<<<<<<< HEAD
N = round(simtime / dt)  # steps

# Plot some scalar value over time
stype = "length"
v_, v1_, v2_, v3_ = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
for i in range(0, N):
    rs.simulate(dt)
    t = np.array(rs.getParameter("type"))
    v = np.array(rs.getParameter(stype))
=======
N = round(simtime / dt)  # steps # |\label{l2_2d:defineEnd}|

# Plot some scalar value over time
stype = "length" # |\label{l2_2d:plotStart}|
v_, v1_, v2_, v3_ = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
for i in range(0, N):
    rs.simulate(dt)
    t = np.array(rs.getParameter("subType")) # |\label{l2_2d:getParam1}|
    v = np.array(rs.getParameter(stype))  # |\label{l2_2d:getParam2}|
>>>>>>> master
    v_[i] = np.sum(v)
    v1_[i] = np.sum(v[t == 1])
    v2_[i] = np.sum(v[t == 2])
    v3_[i] = np.sum(v[t == 3])

t_ = np.linspace(dt, N * dt, N)
plt.plot(t_, v_, t_, v1_, t_, v2_, t_, v3_)
plt.xlabel("time (days)")
plt.ylabel(stype + " (cm)")
plt.legend(["total", "tap root", "lateral", "2. order lateral"])
plt.savefig("results/example_2d.png")
<<<<<<< HEAD
plt.show()
=======
plt.show() # |\label{l2_2d:plotEnd}|
>>>>>>> master
