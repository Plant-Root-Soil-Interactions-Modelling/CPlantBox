"""
Shows the influence of tropism paramter on basal roots
with 4x4 subfigures, different axial resolutions, should lead
to similar results
"""
import sys
sys.path.append("../../..")
import plantbox as pb
import matplotlib.pyplot as plt
import math

fig, axes = plt.subplots(4, 4, figsize=(20,10))

dx = 0.3 # A change in axial resolution will not qualitatively change the resulting image
theta = 70 / 180 * math.pi
N = [0, 1, 2, 4]
sigma = [0, 0.2, 0.4, 0.6]

for i, n in enumerate(N):
    for j, s in enumerate(sigma):        
        rs = pb.RootSystem()

        srp = pb.SeedRandomParameter(rs)
        srp.firstB, srp.delayB, srp.maxB = 10, 3 , 100
        rs.setRootSystemParameter(srp)

        p0 = pb.RootRandomParameter(rs)
        p0.name, p0.subType, p0.la, p0.lmax, p0.ln, p0.r, p0.dx, p0.theta = "taproot", 1, 10, 20, 100, 1, dx, theta
        p0.tropismT = 1  # gravitropism
        p0.tropismN = n
        p0.tropismS = s
        rs.setOrganRandomParameter(p0)

        rs.initializeLB(1,1)
        rs.simulate(50, False)

        axes[i][j].set_title("$N$ = " + str(n) + ", $\sigma$ = " + str(s))  #
        axes[i][j].axis('equal')

        nodes = rs.getNodes()
        segs = rs.getSegments()
        for s in segs:
            n1, n2 = nodes[s.x], nodes[s.y] 
            axes[i][j].plot([n1.x, n2.x], [n1.z, n2.z], "r")


fig.tight_layout()         
fig.canvas.set_window_title("Gravitropism parameters")
plt.show()
