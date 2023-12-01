"""shows the influence of tropism paramters"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import matplotlib.pyplot as plt
import numpy as np

N_ = [0, 0.5, 1, 2]  # strength [1]
sigma_ = [0.1, 0.2, 0.4, 0.8]  # flexibility [1/m]
dx = 0.5  # axial resolution
theta = 70 / 180 * np.pi  # insertion angle [1]
simtime = 65

plant = pb.Plant()
srp = pb.SeedRandomParameter(plant)
srp.firstB, srp.delayB, srp.maxB = 3, 3 , 100
plant.setOrganRandomParameter(srp)

p0 = pb.RootRandomParameter(plant)
p0.name, p0.subType, p0.lmax, p0.r, p0.dx, p0.theta = "taproot", 1, 100, 1, dx, theta
p0.tropismT = 1  # gravitropism
plant.setOrganRandomParameter(p0)

fig, axes = plt.subplots(4, 4, figsize = (15, 10))
for i, n in enumerate(N_):
    for j, sigma in enumerate(sigma_):
        a = axes[i][j]
        plant.reset()  # does not delete parameters

        p0.tropismN = n
        p0.tropismS = sigma

        plant.initializeLB(True)
        plant.simulate(65, True)

        nodes = plant.getNodes()
        segs = plant.getSegments()
        for s in segs:
            n1, n2 = nodes[s.x], nodes[s.y]
            a.plot([n1.x, n2.x], [n1.z, n2.z], "r")
        a.set_title("$N$ = {}, $\sigma$ = {}".format(n, sigma))  #
        a.axis('equal')
        a.set_xlim([-30, 30.])
        a.set_ylim([-40., 0.])

fig.tight_layout()
plt.savefig("results/topics_tropism.png")
plt.show()
