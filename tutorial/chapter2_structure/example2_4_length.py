"""root system length over time"""

import matplotlib.pyplot as plt  # |\label{l2_2d:importEnd}|
import numpy as np
import plantbox as pb  # |\label{l2_2d:importStart}|

path = "../../modelparameter/structural/rootsystem/"  # |\label{l2_2d:defineStart}|
name = "Brassica_napus_a_Leitner_2010"

rs = pb.Plant()
rs.readParameters(path + name + ".xml")
rs.initialize()

simtime = 60.0  # days
dt = 1.0
N = round(simtime / dt)  # steps # |\label{l2_2d:defineEnd}|

# Plot some scalar value over time
stype = "length"  # |\label{l2_2d:plotStart}|
v_, v1_, v2_, v3_ = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
for i in range(0, N):
    rs.simulate(dt)
    t = np.array(rs.getParameter("subType"))  # |\label{l2_2d:getParam1}|
    v = np.array(rs.getParameter(stype))  # |\label{l2_2d:getParam2}|
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
plt.show()  # |\label{l2_2d:plotEnd}|
