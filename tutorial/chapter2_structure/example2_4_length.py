"""root system length over time"""

import matplotlib.pyplot as plt  # |\label{l2_2d:importStart}|
import numpy as np

import plantbox as pb  # |\label{l2_2d:importEnd}|
from plantbox.visualisation import figure_style

path = "../../modelparameter/structural/rootsystem/"  # |\label{l2_2d:defineStart}|
filename = "Brassica_napus_a_Leitner_2010"

rs = pb.Plant()
rs.readParameters(path + filename + ".xml")
rs.initialize()

sim_time = 60.0  # days
dt = 1.0  # days
n_steps = round(sim_time / dt)  # number of iterations  |\label{l2_2d:defineEnd}|

# Plot some scalar value over time
stype = "length"  # |\label{l2_2d:plotStart}|
v_, v1_, v2_, v3_ = np.zeros(n_steps), np.zeros(n_steps), np.zeros(n_steps), np.zeros(n_steps)
for i in range(0, n_steps):
    rs.simulate(dt)
    t = np.array(rs.getParameter("subType"))  # |\label{l2_2d:getParam1}|
    v = np.array(rs.getParameter(stype))  # |\label{l2_2d:getParam2}|
    v_[i] = np.sum(v)
    v1_[i] = np.sum(v[t == 1])
    v2_[i] = np.sum(v[t == 2])
    v3_[i] = np.sum(v[t == 3])
t_ = np.linspace(dt, n_steps * dt, n_steps)

fig, ax = figure_style.subplots12(1, 1)
ax.plot(t_, v_, t_, v1_, t_, v2_, t_, v3_)
ax.set_xlabel("time (days)")
ax.set_ylabel(stype + " (cm)")
ax.legend(["total", "tap root", "lateral", "2. order lateral"])
plt.savefig("results/example_2_4_length.png")
plt.show()  # |\label{l2_2d:plotEnd}|
