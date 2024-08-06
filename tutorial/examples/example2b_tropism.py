"""shows the influence of tropism paramters"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import matplotlib.pyplot as plt
import math

fig, axes = plt.subplots(4, 4, figsize = (15, 10))

N_ = [0, 1, 2, 4]  # strength [1]
sigma_ = [0, 0.1, 0.3, 0.6]  # flexibility [1/m]
dx = 1  # axial resolution
theta = 70 / 180 * math.pi  # insertion angle [1]

rs = pb.RootSystem()
srp = pb.SeedRandomParameter(rs)
srp.firstB, srp.delayB, srp.maxB = 3, 3 , 100
rs.setRootSystemParameter(srp)

p0 = pb.RootRandomParameter(rs)
p0.name, p0.subType, p0.lmax, p0.r, p0.dx, p0.theta = "taproot", 1, 100, 1, dx, theta
p0.tropismT = 1  # gravitropism
rs.setOrganRandomParameter(p0)

for i, n in enumerate(N_):
    for j, sigma in enumerate(sigma_):
        a = axes[i][j]
        rs.reset()  # does not delete parameters

        p0.tropismN = n
        p0.tropismS = sigma

        print("*")
        rs.initializeLB(1, 1)
        print("*")
        rs.simulate(50, False)

        nodes = rs.getNodes()
        segs = rs.getSegments()
        for s in segs:
            n1, n2 = nodes[s.x], nodes[s.y]
            a.plot([n1.x, n2.x], [n1.z, n2.z], "r")
        a.set_title("$N$ = {}, $\sigma$ = {}".format(n, sigma))  #
        a.axis('equal')
        a.set_xlim([-30, 30.])
        a.set_ylim([-40., 0.])

fig.tight_layout()
fig.canvas.manager.set_window_title("Gravitropism parameters")
plt.show()
