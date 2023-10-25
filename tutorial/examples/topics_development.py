"""plant organ lengths over time"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import numpy as np
import matplotlib.pyplot as plt

path = "../../modelparameter/structural/plant/"
name = "hello_world"

plant = pb.Plant()
plant.readParameters(path + name + ".xml")
plant.initialize()

simtime = 60.  # days
dt = 1.
N = round(simtime / dt)  # steps

# calculate organ lengths over time
total, roots, stems, leafs = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
roots_ = [[], [], [], []]

# Simualtion loop
for i in range(0, N):

    plant.simulate(dt)

    t = np.array(plant.getParameter("organType"))
    st = np.array(plant.getParameter("subType"))
    v = np.array(plant.getParameter("length"))  # surface, volume.
    total[i] = np.sum(v)
    roots[i] = np.sum(v[t == 2])  # root
    stems[i] = np.sum(v[t == 3])  # stem
    leafs[i] = np.sum(v[t == 4])  # leaf
    for j in range(0, 4):
        roots_[j].append(np.sum(v[np.logical_and(t == 2, st == j + 1)]))  # roots per sub type

# plot results
t_ = np.linspace(dt, N * dt, N)
plt.plot(t_, total, t_, roots, t_, stems, t_, leafs)
for r in roots_:
    plt.plot(t_, r, '--')
plt.xlabel("time (days)")
plt.ylabel("Organ length (cm)")
plt.legend(["total", "roots", "stem", "leafs", "tap root", "first order", "second order", "basal"])
plt.savefig("results/topics_development.png")
plt.show()
