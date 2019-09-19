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

dx = 0.5
theta = 70 / 180 * pi
N = [0, 1, 2, 4]
sigma = [0, 0.2, 0.4, 0.6]

fig, axes = plt.subplots(4, 4)

c = 0
for i, n in enumerate(N):
    for j, s in enumerate(sigma):

        c += 1
        rs = pb.RootSystem()
        rs.setSeed(c)
        maxB, firstB, delayB = 100, 10., 3
        rsp = pb.RootSystemParameter()
        rsp.set(-3., firstB, delayB, maxB, 0, 1.e9, 1.e9, 1.e9, 0., 0.)
        rs.setRootSystemParameter(rsp)
        p0 = pb.RootTypeParameter(rs)
        p0.name, p0.type, p0.la, p0.nob, p0.ln, p0.r, p0.dx = "taproot", 1, 10, 20, 89. / 19., 1, 0.5

        p0.tropismT = 1  # GRAVI
        p0.tropismN = n
        p0.tropismS = s
        p0.dx = dx
        p0.theta = theta

        rs.setOrganTypeParameter(p0)

        rs.initialize()
        rs.simulate(50)

        axes[i][j].set_title("n: " + str(n) + ", sigma: " + str(s))  #
        axes[i][j].axis('equal')

        nodes = vv2a(rs.getNodes())
        segs = seg2a(rs.getSegments())
        for s in segs:
            axes[i][j].plot([nodes[s[0], 1], nodes[s[1], 1]], [nodes[s[0], 2], nodes[s[1], 2]], "r")

plt.show()
