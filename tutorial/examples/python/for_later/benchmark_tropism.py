"""
Shows the influence of tropism paramter on basal roots
with 4x4 subfigures, different axial resolutions, should lead
to similar results
"""
import sys
sys.path.append("..")
import plantbox as pb
import matplotlib.pyplot as plt
from rb_tools import *
from rsml import *
from cmath import pi

dx = 1
theta = 70 / 180 * pi
N = [0, 1, 2, 4]
sigma = [0, 0.2, 0.4, 0.6]

fig, axes = plt.subplots(4, 4)

c = 0
for i, n in enumerate(N):
    for j, s in enumerate(sigma):

        c += 1
        rs = pb.RootSystem()
        rs.setSeed(c)  # random seed

        srp = pb.SeedRandomParameter(rs)
        srp.firstB, srp.delayB, srp.maxB = 10, 3 , 100
        rs.setRootSystemParameter(srp)

        p0 = pb.RootRandomParameter(rs)
        p0.name, p0.subType, p0.la, p0.nob, p0.ln, p0.r, p0.dx, p0.theta = "taproot", 1, 10, 20, 89. / 19., 1, dx, theta
        p0.tropismT = 1  # GRAVI
        p0.tropismN = n
        p0.tropismS = s

        rs.setOrganRandomParameter(p0)

        rs.initialize()
        rs.simulate(50)

        axes[i][j].set_title("n: " + str(n) + ", sigma: " + str(s))  #
        axes[i][j].axis('equal')

        nodes = vv2a(rs.getNodes())
        segs = seg2a(rs.getSegments())
        for s in segs:
            axes[i][j].plot([nodes[s[0], 1], nodes[s[1], 1]], [nodes[s[0], 2], nodes[s[1], 2]], "r")

        # todo fix axis -30+ 30, -40 + 10

print("Figure for dx = ", dx)

plt.show()
