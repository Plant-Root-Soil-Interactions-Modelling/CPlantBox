import sys
sys.path.append("../..")
sys.path.append("../../src/python_modules/")
import time
import numpy as np
import plantbox as rb
import matplotlib
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import math


def uniqueIndexes(l):
    seen = set()
    res = []
    for i, n in enumerate(l):
        if n not in seen:
            res.append(i)
            seen.add(n)
    return res

#############################################################


path = 'dgf_files'

# #
# read in dgf
# #
name = 'RootSys_verysimple.dgf'
with open(path + '/' + name) as f:
    content = f.readlines()

i = 3
coordlist = []
while ('#' in content[i]) == False:
    line = content[i]
    nums = line.split()
    nums = map(float, nums)
    nums = list(nums)
    coordlist.extend((nums))
    i += 1
A = np.asarray(coordlist)
nodes = np.reshape(A, (int(len(A) / 3), 3))

i = i + 3
paramlist = []
while ('#' in content[i]) == False:
    line = content[i]
    nums = line.split()
    nums = map(float, nums)
    nums = list(nums)
    paramlist.extend((nums))
    i += 1
B = np.asarray(paramlist)
params = np.reshape(B, (int(len(B) / 10), 10))

# #
# parameters...
# #
x = nodes[2:-1, 0] * 100;
y = nodes[2:-1, 1] * 100;
z = nodes[2:-1, 2] * 100;
leng = params[:, 5];
rad = params[:, 6];
etime = params[:, 9];
branch = params[:, 3];
prev = params[:, 0];
order = params[:, 2];

# #
# compute max root age
# #
etime = np.ceil(etime)
maxage = int(np.max(etime))

# #
# check dgf
# #

branches = np.unique(branch)
branches = branches[1:-1]
idxbr = uniqueIndexes(branch)
idxbr = idxbr[1:-1]

# #
# are the connections correct??
# #
plt.figure()
ax = plt.axes(projection = '3d')
for ii in range(0, len(branches)):
    idx1 = int(prev[idxbr[ii]])
    idx2 = np.where(branch == branches[ii])
    ax.plot3D(np.append(x[idx1], x[idx2[0] - 1]), np.append(y[idx1], y[idx2[0] - 1]), np.append(z[idx1], z[idx2[0] - 1]))
    ax.set_title("Connections", fontsize = 20, pad = 250)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_zlabel("z (cm)")

# #
# are the orders correct??
# #
n = np.max(order) + 1
colors = plt.cm.jet(np.linspace(0, 1, n))
plt.figure()
ax = plt.axes(projection = '3d')
for ii in range(0, len(branches)):
    idx1 = int(prev[idxbr[ii]])
    idx2 = np.where(branch == branches[ii])
    idxcol = int(order[idxbr[ii]])
    print(idx1)
    print(idx2)
    print(idxcol)
    ax.plot3D(np.append(x[idx1], x[idx2[0] - 1]), np.append(y[idx1], y[idx2[0] - 1]), np.append(z[idx1], z[idx2[0] - 1]), color = colors[idxcol])
    ax.set_title("Orders", fontsize = 20, pad = 250)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_zlabel("z (cm)")

# #
# is the emergence time correct??
# #
n = maxage + 1
colors = plt.cm.viridis(np.linspace(0, 1, n))
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
for ii in range(0, len(branches)):
    idx1 = int(prev[idxbr[ii]])
    idx2 = np.where(branch == branches[ii])
    colidx = etime[idx2[0] - 1]
    colidx = colidx.astype(int)
    cs = ax.scatter3D(x[idx2[0] - 1], y[idx2[0] - 1], z[idx2[0] - 1], c = colors[colidx])
    ax.set_title("Emergence times (days)", fontsize = 20, pad = 250)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_zlabel("z (cm)")

start = np.min(etime)
mid = np.min(etime) + (np.max(etime) - np.min(etime)) / 2
end = np.max(etime)
cax = fig.add_axes([.918, 0.175, 0.025, 0.4])
cb = fig.colorbar(cs, cax = cax, ticks = [0, .5, 1])
cb.ax.set_yticklabels([str(start), str(mid), str(end)])

# #
# is the radius correct??
# #
fg = np.ceil((rad - np.min(rad)) * 10 ** 2);
colors = plt.cm.viridis(np.linspace(0, 1, 100))

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
for ii in range(0, len(branches)):
    idx1 = int(prev[idxbr[ii]])
    idx2 = np.where(branch == branches[ii])
    colidx = fg[idx2[0] - 1]
    colidx = colidx.astype(int)
    cs = ax.scatter3D(x[idx2[0] - 1], y[idx2[0] - 1], z[idx2[0] - 1], c = colors[colidx])
    ax.set_title("Radii", fontsize = 20, pad = 250)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_zlabel("z (cm)")

start = np.min(rad)
mid = np.min(rad) + (np.max(rad) - np.min(rad)) / 2
end = np.max(rad)
cax = fig.add_axes([.918, 0.175, 0.025, 0.4])
cb = fig.colorbar(cs, cax = cax, ticks = [0, .5, 1])
cb.ax.set_yticklabels([str(start), str(mid), str(end)])

plt.show()

