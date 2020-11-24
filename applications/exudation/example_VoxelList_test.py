import sys; sys.path.append("../..")

# Computation took 27.435727834701538 s
# max 29.96163225689554 min 0.0 sum 1018.807994059238
# volume of concentration above threshold:  189.46530339293474
# this is 16.86677722498618 % of the overall volume

# Computation took 28.135575532913208 s
# max 406.35834594892395 min 0.0 sum 6237.840381419113
# volume of concentration above threshold:  34.79822890056393
# this is 3.097844112769486 % of the overall volume

import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
import plantbox as rb
from pyevtk.hl import gridToVTK

#
# Root system
#

# ExudationModel              old model
# ExudationModel2             with voxelList and maps

rs = rb.RootSystem()

path = "../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # "Anagallis_straight_simple"
rs.readParameters(path + name + ".xml")

for p in rs.getRootRandomParameter():
    p.gf = 2  # linear growth function

rs.setSeed(0)
rs.initialize()
simtime = 2
rs.simulate(simtime, True);
rs.write(name + ".vtp")

#
# Grid parameter
#
nodes = np.array([np.array(n) for n in rs.getNodes()])
boxmin = nodes.min(axis = 0); boxmax = nodes.max(axis = 0);  # cm
width = abs(max(boxmax[0], boxmax[1]) - min(boxmin[0], boxmin[1])) + 6;  # cm
depth = abs(boxmin[2]) + 3
xres = 0.2;
yres = 0.2;
zres = 0.2;
nx = int(width / xres);
ny = int(width / yres);
nz = int(depth / zres);
print("Width", width, "Depth", depth, " at a Resolution", nx, ny, nz)

#
# Model parameter
#
model = rb.ExudationModel(width, width, depth, nx, ny, nz, rs)
model2 = rb.ExudationModel2(width, width, depth, nx, ny, nz, rs)

for m in [model, model2]:
    m.Q = 4  # µg/d/tip
    m.Dl = 2.43e-6 * 3600 * 24  # cm2/d
    m.theta = 0.3
    m.R = 16.7  # 16.7  # -
    m.k = 2.60e-6 * 3600 * 24  # d-1
    m.l = 4  # cm (for line source only)
    
    #
    # Numerical parameter
    #
    m.type = rb.IntegrationType.mls;  # mps, mps_straight, mls
    m.n0 = 10  # integration points per cm
    m.thresh13 = 1.e-15;  # threshold to neglect diffusing g (eqn 13)
    m.calc13 = True;  # turns Eqn 13  on (True) and off (False)
    m.observationRadius = 2;  # limits computational domain around roots [cm]

roots = rs.getPolylines()
print("Number of roots", len(roots))


# t = time.time()
# C = model.calculate(simtime)
# elapsed = time.time() - t
# print("Computation took", elapsed, "s")

t = time.time()
print("make voxel lists")
# C = np.zeros((nx*ny*nz,))
model2.makeVoxelLists() # equals model.makeVoxelLists(0, len(roots)) 
C = np.array(model2.calculate(simtime)) # equals model.calculate(simtime, 0, len(roots))  
print(C.shape)
# C = model2.addResults(C)
elapsed = time.time() - t
print("Computation took", elapsed, "s")

#
# post processing...
#
C_ = np.zeros((nx,ny,nz))
C = np.reshape(C, (nz, ny, nx))  
for i in range(0,np.shape(C)[0]):
    for j in range(0,np.shape(C)[1]):
        for k in range(0,np.shape(C)[2]):
            C_[k,j,i] = C[i,j,k]
            
            
del C
C = C_

X = np.linspace(-width / 2, width / 2, nx)
Y = np.linspace(-width / 2, width / 2, ny)
Z = np.linspace(-depth, 0, nz)
X_,Y_, Z_ = np.meshgrid(X, Y, Z, indexing = "ij")

print("max" , np.max(C.flat), "min", np.min(C.flat), "sum", np.sum(C.flat))
num_th = (C > 0).sum()  # number of points for which concentration is larger than threshold
print("volume of concentration above threshold: ", num_th * width / nx * width / ny * depth / nz)
print("this is", num_th / (nx * ny * nz) * 100, "% of the overall volume")

gridToVTK("./Exudates", X, Y, Z, pointData = {"Exudates":C})

fig1 = plt.figure()
ax = plt.axes()
C_ = C[:, int(nx/2), :]
# print(np.max(C_))
# levels = np.logspace(np.log10(np.max(C_+1)) - 5, np.log10(np.max(C_)), 10)
# print(levels)
levels = np.linspace(np.min(C_[:]), np.max(C_[:]))
print(levels)
cs = ax.contourf(X_[:, int(nx / 2), :], Z_[:, int(nx / 2), :], C_, levels = levels, cmap = 'jet')  # , locator = ticker.LogLocator(),
ax.set_xlabel('x')
ax.set_ylabel('z')
plt.axis('equal')
cbar = fig1.colorbar(cs)

fig2 = plt.figure()
node_tip = nodes[-1]
idy = (np.abs(Y - node_tip[1])).argmin()  # y-index of soil element closest to the point on the root axis
idx = (np.abs(X - node_tip[0])).argmin()  # x-index of soil element closest to the point on the root axis
C_ = C[idx, idy, :]  # 1d array
plt.plot(Z, C_, 'k-')
plt.xlabel('z')
plt.ylabel('c')

plt.show()
