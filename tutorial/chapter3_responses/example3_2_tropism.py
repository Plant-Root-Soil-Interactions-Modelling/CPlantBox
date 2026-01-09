"""shows the influence of tropism paramters"""

import matplotlib.pyplot as plt
import numpy as np

import plantbox as pb
from plantbox.visualisation import figure_style

N_ = [0, 1, 2, 4]  # strength  |\label{l3_2_tropism:valuesstart}|
sigma_ = [0, 0.1, 0.3, 0.6]  # flexibility (1/cm)   |\label{l3_2_tropism:valuesend}|
dx = 1  # axial resolution (cm) |\label{l3_2_tropism:paramsstart}|
theta = 70 / 180 * np.pi  # insertion angle

plant = pb.Plant()
srp = pb.SeedRandomParameter(plant)
srp.firstB, srp.delayB, srp.maxB = 3, 3, 100
plant.setOrganRandomParameter(srp)

p0 = pb.RootRandomParameter(plant)
p0.name, p0.subType, p0.lmax, p0.r, p0.dx, p0.theta = "taproot", 1, 100, 1, dx, theta
p0.tropismT = 1  # gravitropism
plant.setOrganRandomParameter(p0)  # |\label{l3_2_tropism:paramsend}|

fig, ax_ = figure_style.subplots11large(4, 4)

for i, n in enumerate(N_):  # |\label{l3_2_tropism:loopstart}|
    for j, sigma in enumerate(sigma_):
        ax = ax_[i][j]
        plant.reset()  # does not delete parameters

        p0.tropismN = n
        p0.tropismS = sigma

        plant.initializeLB()
        plant.simulate(50, False)  # |\label{l3_2_tropism:loopend}|

        nodes = plant.getNodes()  # |\label{l3_2_tropism:readstart}|
        segs = plant.getSegments()
        for s in segs:
            n1, n2 = nodes[s.x], nodes[s.y]  # |\label{l3_2_tropism:readend}|
            ax.plot([n1.x, n2.x], [n1.z, n2.z], "r")  # |\label{l3_2_tropism:plotstart}|
        ax.set_title(rf"$N$ = {n}, $\sigma$ = {sigma}")
        ax.axis("equal")
        ax.set_xlim([-30, 30.0])
        ax.set_ylim([-40.0, 0.0])

fig.tight_layout()
fig.canvas.manager.set_window_title("Gravitropism parameters")
plt.show()  # |\label{l3_2_tropism:plotend}|
