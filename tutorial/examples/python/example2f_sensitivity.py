"""sensitivity analysis: insertion anlge on root tip distribution"""
import sys
sys.path.append("../../..")
import plantbox as pb

import math
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt


# sets all standard deviation to value*s
def set_all_sd(rs, s):
    for p in rs.getRootRandomParameter():
        p.lbs = p.lb * s
        p.las = p.la * s
        p.lns = p.ln * s
        # p.nobs = p.nob * s TODO
        p.rs = p.r * s
        p.a_s = p.a * s


# Parameters
path = "../../../modelparameter/rootsystem/"
name = "Triticum_aestivum_a_Bingham_2011"
simtime = 20
N = 25  # resolution of paramter
runs = 20  # iterations
theta0_ = np.linspace(0, math.pi / 2, N)


# One simulation
def simulate(i):
    rs = pb.RootSystem()
    rs.readParameters(path + name + ".xml")
    set_all_sd(rs, 0.)  # set all sd to zero
    # vary parameter
    p1 = rs.getRootRandomParameter(1)  # tap and basal root type
    p1.theta = theta0_[i]
    # simulation
    rs.initialize(1, 1, False)
    rs.simulate(simtime, False)
    # calculate target
    depth = 0.  # mean depth
    rad_dist = 0.  # mean raidal distance
    roots = rs.getPolylines()
    for r in roots:
        depth += r[-1].z
        rad_dist += math.hypot(r[-1].x, r[-1].y)
    depth /= len(roots)
    rad_dist /= len(roots)
    return depth, rad_dist


depth_ = np.zeros(N)
rad_dist_ = np.zeros(N)

for r in range(0, runs):

    print("run", r + 1)

    # Parallel execution
    param = []  # param is a list of tuples
    for i in range(0, N):
        param.append((i,))
    pool = Pool()
    output = pool.starmap(simulate, param)
    pool.close()

    # Copy results
    for i, o in enumerate(output):
        depth_[i] += (o[0] / runs)
        rad_dist_[i] += (o[1] / runs)

# Figure
fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (10, 8))
axes[0].set_xlabel('Insertion angle theta (-)')
axes[1].set_xlabel('Insertion angle theta (-)')
axes[0].set_ylabel('Mean tip depth (cm)')
axes[1].set_ylabel('Mean tip radial distance (cm)')
axes[0].plot(theta0_, depth_)
axes[1].plot(theta0_, rad_dist_)
fig.subplots_adjust()
plt.savefig("results/example_4b.png")
plt.show()
