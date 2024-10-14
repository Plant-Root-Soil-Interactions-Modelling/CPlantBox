""" shows inter lateral spacing (ln) and how a linear slope (lnk) can modify them """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import matplotlib.pyplot as plt
import numpy as np

fig, axes = plt.subplots(2, 3, figsize = (15, 7))

rs = pb.RootSystem()
srp = pb.SeedRandomParameter(rs)
srp.firstB, srp.delayB, srp.maxB = 1., 1., 0
rs.setRootSystemParameter(srp)

dx0, dx1 = 100, 1  # very large resolution for taproot
theta = 70 / 180 * np.pi

p0 = pb.RootRandomParameter(rs)
p0.name, p0.subType, p0.lmax, p0.r, p0.dx, p0.theta = "taproot", 1, 50., 2., dx0, 0.  # parameters as before
p0.tropismT, p0.tropismN, p0.tropismS = pb.TropismType.gravi, 2, 0.1
p0.successor, p0.successorP = [[2]], [[1.]]  # set up successors
p0.lns = 0.  # test with other values...
p0.lb, p0.la = 0., 0.
rs.setOrganRandomParameter(p0)

p1 = pb.RootRandomParameter(rs)
p1.name, p1.subType, p1.lmax, p1.r, p1.dx, p1.theta = "lateral", 2, 30., 1., dx1, theta
p1.tropismT, p1.tropismN, p1.tropismS = pb.TropismType.exo, 2, 0.2
rs.setOrganRandomParameter(p1)

ln_ = [4., 2.]  # inter lateral distance
lnk_ = [-2. / 45., 0, 2. / 45]  # slope

for i, ln in enumerate(ln_):
    for j, lnk in enumerate(lnk_):
        rs.reset()  # does not delete parameters

        p0.ln = ln
        p0.lnk = ln * lnk  # set up linearly altered spaces

        rs.initializeLB(1, 1)
        rs.simulate(100, False)

        a = axes[i][j]
        nodes = rs.getNodes()
        segs = rs.getSegments()
        for k, s in enumerate(segs):
            n1, n2 = nodes[s.x], nodes[s.y]
            a.plot([n1.x, n2.x], [n1.z, n2.z], 'g')

        pl = rs.getPolylines()
        for n in pl[0]:
             axes[i][j].plot([n.x], [n.z], "r*")

        axes[i][j].set_title("$ln$ = {:.2f}, $lnk$ = {:.2f}".format(p0.ln, p0.lnk))
        axes[i][j].axis('equal')

fig.tight_layout()
fig.canvas.manager.set_window_title("Inter lateral distances")
plt.show()
