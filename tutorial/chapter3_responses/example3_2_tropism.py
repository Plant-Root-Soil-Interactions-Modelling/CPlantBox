"""shows the influence of tropism paramters"""

import matplotlib.pyplot as plt
import numpy as np

import plantbox as pb

fig, axes = plt.subplots(4, 4, figsize=(15, 10))

N_ = [0, 1, 2, 4]  # strength [1] # |\label{l3_2_tropism:valuesstart}|
sigma_ = [0, 0.1, 0.3, 0.6]  # flexibility [1/m]  # |\label{l3_2_tropism:valuesend}|
dx = 1  # axial resolution # |\label{l3_2_tropism:paramsstart}|
theta = 70 / 180 * np.pi  # insertion angle [1]

rs = pb.Plant()
srp = pb.SeedRandomParameter(rs)
srp.firstB, srp.delayB, srp.maxB = 3, 3, 100
rs.setOrganRandomParameter(srp)

p0 = pb.RootRandomParameter(rs)
p0.name, p0.subType, p0.lmax, p0.r, p0.dx, p0.theta = "taproot", 1, 100, 1, dx, theta
p0.tropismT = 1  # gravitropism
rs.setOrganRandomParameter(p0)  # |\label{l3_2_tropism:paramsend}|

for i, n in enumerate(N_):  # |\label{l3_2_tropism:loopstart}|
    for j, sigma in enumerate(sigma_):
        a = axes[i][j]
        rs.reset()  # does not delete parameters

        p0.tropismN = n
        p0.tropismS = sigma

        rs.initializeLB()
        rs.simulate(50, False)  # |\label{l3_2_tropism:loopend}|

        nodes = rs.getNodes()  # |\label{l3_2_tropism:readstart}|
        segs = rs.getSegments()
        for s in segs:
            n1, n2 = nodes[s.x], nodes[s.y]  # |\label{l3_2_tropism:readend}|
            a.plot([n1.x, n2.x], [n1.z, n2.z], "r")  # |\label{l3_2_tropism:plotstart}|
        a.set_title(rf"$N$ = {n}, $\sigma$ = {sigma}")
        a.axis("equal")
        a.set_xlim([-30, 30.0])
        a.set_ylim([-40.0, 0.0])

fig.tight_layout()
fig.canvas.manager.set_window_title("Gravitropism parameters")
plt.show()  # |\label{l3_2_tropism:plotend}|
