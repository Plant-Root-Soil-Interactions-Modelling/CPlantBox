import sys;
sys.path.append("../..")
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
from pyevtk.hl import gridToVTK
import plantbox as rb                                                   

#
# Root system
#
rs = rb.RootSystem()

path = "../../modelparameter/rootsystem/"
name = "Faba_exudation" 
rs.readParameters(path + name + ".xml")

for p in rs.getRootRandomParameter():
    p.gf = 2  # linear growth function
    
srp = rb.SeedRandomParameter(rs)  # with default values
srp.seedPos = rb.Vector3d(0., 0., -3.)  # [cm] seed position

#set geometry 
width = 20  # cm
depth = 45    
soilcore = rb.SDF_PlantContainer(width, width, depth, True)
rs.setGeometry(soilcore)  

rs.setSeed(0)
rs.initialize()

for ii in range(0,21): 
    simtime = 1
    rs.simulate(simtime, True);
#
# Grid parameter
#
nodes = np.array([np.array(n) for n in rs.getNodes()])

#np.save("../nodes/day"+str(ii+1), nodes)
xres = 0.3;
yres = 0.3;
zres = 0.3;
nx = int(width / xres);
ny = int(width / yres);
nz = int(depth / zres);
print("Width", width, "Depth", depth, " at a Resolution", nx, ny, nz)


X = np.linspace(-width / 2, width / 2, nx)
Y = np.linspace(-width / 2, width / 2, ny)
Z = np.linspace(-depth, 0, nz)
X_,Y_, Z_ = np.meshgrid(X, Y, Z, indexing = "ij")
nodes = (np.vstack([X_.flatten(), Y_.flatten(), Z_.flatten()])).T

#
# Model parameter
#
model = rb.ExudationModel2(width, width, depth, nx, ny, nz, rs)
model.Q = 33.38 # µg/d/tip
model.Dl = 1.04e-3  # cm2/d - with impedance factor
model.theta = 0.3 #-
model.R = 1  # -
model.k = 0.22  # d-1

#
# Numerical parameter
#
model.type = rb.IntegrationType.mps;  # mps, mps_straight, mls
model.n0 = 20  # integration points per cm
model.thresh13 = 1.e-15;  # threshold to neglect diffusing g (eqn 13)
model.calc13 = True;  # turns Eqn 13  on (True) and off (False)
model.observationRadius = 0.5;  # limits computational domain around roots [cm]


#
#set threshold 
#
thresh = int(1300)
rn = 344 #problematic root

t = time.time()
print("make voxel lists")
model.makeVoxelLists(rn,rn+1)


print('start check')
C = np.zeros((nx*ny*nz,))
roots = rs.getPolylines()
for i in [rn]:
    C1 = model.calculate(ii+1,i,i+1)
    print(np.sum(C1), i)
    C = np.add(C1,C)
    
print('end check')
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

gridToVTK("./Exudates_day21_problemroot", X, Y, Z, pointData = {"Exudates":C})




