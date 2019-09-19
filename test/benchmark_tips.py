''' More like a case study, Should be moved to a jupyter notebook '''

#
# Compares numerical approximations with altering resolutions dt, dx (L86 - )
#
# plots root tip radial distance
# plots root tip depth
# x-y tip density map
#
# computes in parallel to enable a lot of runs
#
import ../rootbox as rb

from multiprocessing import Pool
from itertools import product

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from matplotlib.dates import mx2num


def v2a(vd):  # rb.std_vector_double_ to numpy array
    l = np.zeros(len(vd))
    for i in range(0, len(vd)):
        l[i] = vd[i]
    return l


def a2v(a):  #  numpy array to rb.std_vector_double
    l = rb.std_vector_double_()
    for d in a:
        l.append(d)
    return l


def a2i(a):  #  numpy array to rb.std_int_double
    l = rb.std_vector_int_()
    for i in a:
        l.append(i)
    return l


def vv2a(vd):  # rb.std_vector_Vector3_ to numpy array
    N = len(vd)
    l = np.zeros((N, 3))
    for i in range(0, N):
        l[i, :] = [vd[i].x, vd[i].y, vd[i].z]
    return l


def simOnce(name, simtime, lbins, lrange, zbins, zrange, dx, dt):
    abins = 100
    arange = (-math.pi, math.pi)
    # simulation
    rs = rb.RootSystem()
    rs.openFile(name, "modelparameter/")
    for p in rs.getRootTypeParameter():
        p.dx = dx

    rs.initialize()
    rs.getRootTypeParameter(4).theta = 80. / 180.*math.pi  # fix insertion angle of the basal roots
    rs.initialize()

    N = round(simtime / dt)
    for i in range(0, N):
        rs.simulate(dt, True);
    # analysis
    img = np.zeros((2 * lbins, 2 * lbins))
    nodes = vv2a(rs.getNodes());
    tips = rs.getRootTips()
    notips = len(tips)
    z_ = np.zeros(notips)
    l_ = np.zeros(notips)
    alpha_ = np.zeros(notips)
    c = 0;
    mx = 0
    my = 0
    for t in tips:
        x = nodes[t, 0]  # -1.e-9
        y = nodes[t, 1]  # -1.e-9
        z = nodes[t, 2]

        # tip top view
        i = math.floor(((x + lrange[1]) / (2.*lrange[1])) * 2.*float(lbins))  # np.around((x/lrange[1])*lbins+lbins)
        j = math.floor(((y + lrange[1]) / (2.*lrange[1])) * 2.*float(lbins))  # np.around((y/lrange[1])*lbins+lbins)

        i = min(i, 2 * lbins - 1)
        j = min(j, 2 * lbins - 1)
        i = int(max(i, 0))
        j = int(max(j, 0))
        # print(x,y,i,j)

        img[j, i] += 1.

        # depth, and length distribution
        z_[c] = z
        l_[c] = math.sqrt(x * x + y * y)
        alpha_[c] = math.atan2(x, y)
        c += 1
        hl, bins = np.histogram(l_, bins = lbins, range = lrange)
        hz, bins = np.histogram(z_, bins = zbins, range = zrange)
        ha, bins = np.histogram(alpha_, bins = abins, range = arange)
        mx += nodes[t, 0]
        my += nodes[t, 1]

    return hl, hz, ha, img, mx / len(tips), my / len(tips)


# Params
dx = 0.5
dt = 5

runs = 100
# name = "Lupinus_albus_Leitner_2014"
name = "Zea_mays_1_Leitner_2010"
# name = "Anagallis_femina_Leitner_2010"
# name = "Triticum_aestivum_a_Bingham_2011"
simtime = 60

# Histogram params
lbins = 40
lrange = (0., 20.)
zbins = 120
zrange = (-120., 0.)

pool = Pool()  # defaults to number of available CPU's
chunksize = 20  # this may take some guessing ... take a look at the docs to decide
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

