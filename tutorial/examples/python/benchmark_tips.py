"""
Analysis of root tip distribution

works in parallel using multiprocessing

TODO this is not working...
"""

#
# plots root tip radial distance
# plots root tip depth
# x-y tip density map
#
# computes in parallel to enable a lot of runs
#
import sys
sys.path.append("../../..")
import plantbox as pb

from multiprocessing import Pool
from itertools import product

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from matplotlib.dates import mx2num


def simOnce(name, simtime, lbins, lrange, zbins, zrange, dx, dt):

    abins = 100
    arange = (-math.pi, math.pi)
    hl,hz,ha = None, None, None # histograms
    
    # Simulation
    rs = pb.RootSystem() 
    path = "../../../modelparameter/rootsystem/"
    rs.readParameters(path + name + ".xml")
    for p in rs.getRootRandomParameter():
        p.dx = dx
    rs.getRootRandomParameter(4).theta = 80. / 180.*math.pi  # fix insertion angle of the basal roots
    rs.initialize()

    N = round(simtime / dt)
    for i in range(0, N):
        rs.simulate(dt, True);

    # Analysis
    img = np.zeros((2 * lbins, 2 * lbins))
    nodes = rs.getNodes();
    tips = rs.getRootTips()
    z_ = np.zeros(len(tips))
    l_ = np.zeros(len(tips))
    alpha_ = np.zeros(len(tips))    
    c = 0 # counter
    mx, my = 0, 0 # mean tip position
    for t in tips:
        x, y, z = nodes[t].x, nodes[t].y, nodes[t].z

        # tip top view
        i = math.floor(((x + lrange[1]) / (2.*lrange[1])) * 2.*float(lbins))  # np.around((x/lrange[1])*lbins+lbins)
        j = math.floor(((y + lrange[1]) / (2.*lrange[1])) * 2.*float(lbins))  # np.around((y/lrange[1])*lbins+lbins)
        i = min(i, 2 * lbins - 1)
        j = min(j, 2 * lbins - 1)
        i = int(max(i, 0))
        j = int(max(j, 0))
        img[j, i] += 1.

        # depth and length distribution
        z_[c] = z
        l_[c] = math.sqrt(x * x + y * y)
        alpha_[c] = math.atan2(x, y)
        c += 1
        hl, bins = np.histogram(l_, bins = lbins, range = lrange)
        hz, bins = np.histogram(z_, bins = zbins, range = zrange)
        ha, bins = np.histogram(alpha_, bins = abins, range = arange)
        mx += nodes[t].x
        my += nodes[t].y
    
    if (len(tips)==0):
        print("No tips?")
        return hl, hz, ha, img, mx, my 
    else:
        return hl, hz, ha, img, mx / len(tips), my / len(tips)


# Params
dx = 0.5
dt = 5
simtime = 60

runs = 100
 #name = "Lupinus_albus_Leitner_2014"
# name = "Zea_mays_1_Leitner_2010"
name = "Anagallis_femina_Leitner_2010"
# name = "Triticum_aestivum_a_Bingham_2011"

# Histogram params
lbins = 40
lrange = (0., 20.) # cm
zbins = 120
zrange = (-120., 0.)

pool = Pool()  # defaults to number of available CPU's
chunksize = 100  # this may take some guessing ... take a look at the docs to decide
# simOnce(name, simtime, lbins, lrange, zbins, zrange, dx, dt)
output = pool.starmap(simOnce, [(name, simtime, lbins, lrange, zbins, zrange, dx, dt)] * runs)

mmx = 0
mmy = 0
allL, allZ, allA, tiptop, mx, my = output[0];
for i in range(1, runs):
    L, Z, A, img, mx, my = output[i];
    allL = np.vstack((allL, L))
    allZ = np.vstack((allZ, Z))
    allA = np.vstack((allA, A))
    tiptop += img
    mmx += mx
    mmy += my

meanZ = np.mean(allZ, 0)
semZ = stats.sem(allZ, 0)
meanL = np.mean(allL, 0)
semL = stats.sem(allL, 0)
meanA = np.mean(allA, 0)
semA = stats.sem(allA, 0)

plt.figure(1)
plt.errorbar(np.linspace(lrange[0], lrange[1], lbins), meanL, semL, linestyle = 'None', marker = '^')
plt.title("Root tip radial distance")
plt.show(False)

# plt.figure(2)
# plt.errorbar(np.linspace(zrange[0],zrange[1],zbins), meanZ, semZ, linestyle='None', marker='^')
# plt.title("Root tip depth")

plt.figure(3)
abins = 100
arange = (-math.pi, math.pi)
plt.errorbar(np.linspace(arange[0], arange[1], abins), meanA, semA, linestyle = 'None', marker = '^')
plt.title("Angular distribution")

print("\nnumber of tips")
print("total \t", np.sum(tiptop[:, :]))
print("top \t", np.sum(tiptop[0:lbins, :]))
print("bot \t", np.sum(tiptop[lbins:2 * lbins, :]))
print("left \t", np.sum(tiptop[:, 0:lbins]))
print("right \t", np.sum(tiptop[:, lbins:2 * lbins]))

# tiptop[0:lbins,:] = 0
# tiptop[lbins:2*lbins,:] = 0
# tiptop[0:lbins,:] += 1
# tiptop[lbins:2*lbins,:] += 2
plt.figure(4)
plt.title("Tip top view")
plt.imshow(tiptop)

print("\nmysterious drift (cm)")
print(mmx / runs)
print(mmy / runs)

plt.show()

