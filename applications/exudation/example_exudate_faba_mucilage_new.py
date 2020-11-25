import sys;
sys.path.append("../..")
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
from pyevtk.hl import gridToVTK
import plantbox as rb                                              
    
############################
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

for ii in range(0,2): 
    simtime = 1
    rs.simulate(simtime, True);
#rs.write("../vtp/day"+("{:02d}".format(ii+1)) + ".vtp")

#
# Grid parameter
#
nodes = np.array([np.array(n) for n in rs.getNodes()])
#np.save("../nodes/day"+str(ii+1), nodes)
xres = 0.1;
yres = 0.1;
zres = 0.1;
nx = int(width / xres);
ny = int(width / yres);
nz = int(depth / zres);
print("Width", width, "Depth", depth, " at a Resolution", nx, ny, nz)

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
model.n0 = 10  # integration points per cm
model.thresh13 = 1.e-15;  # threshold to neglect diffusing g (eqn 13)
model.calc13 = True;  # turns Eqn 13  on (True) and off (False)
model.observationRadius = 0.5;  # limits computational domain around roots [cm]

t = time.time()
roots = rs.getPolylines()
order = np.array(rs.getParameter("type"))
thresh = int(1300)

print("make voxel lists")
model.makeVoxelLists()

C = np.zeros((nx*ny*nz,))
C2lat = np.zeros((nx*ny*nz,))
Cthresh = np.zeros((nx*ny*nz,))

for i in range(0,len(roots)):
    model.calculate(ii+1,i,i+1)
    C = model.addResults(C)
    if order[i]<3:
        C2lat = model.addResults(C2lat)
    Cthresh = np.array(model.addResults(Cthresh))
    Cthresh[Cthresh<thresh] = 0

elapsed = time.time() - t
print("Computation took", elapsed, "s")

#
# post processing...
#
C_ = np.zeros((nx,ny,nz))
C2lat_ = np.zeros((nx,ny,nz))
Cthresh_ = np.zeros((nx,ny,nz))

C = np.reshape(C, (nz, ny, nx))  
C2lat = np.reshape(C2lat, (nz, ny, nx))  
Cthresh = np.reshape(Cthresh, (nz, ny, nx))  
for i in range(0,np.shape(C)[0]):
    for j in range(0,np.shape(C)[1]):
        for k in range(0,np.shape(C)[2]):
            C_[k,j,i] = C[i,j,k]
            C2lat_[k,j,i] = C2lat[i,j,k]
            Cthresh_[k,j,i] = Cthresh[i,j,k]
            
del C, Cthresh, C2lat
C = C_
C2lat = C2lat_
Cthresh = Cthresh_

#
# post processing...
#
#np.savez_compressed("../concentration/day"+("{:02d}".format(ii+1)), C, C2lat, Cthresh)

X = np.linspace(-width / 2, width / 2, nx)
Y = np.linspace(-width / 2, width / 2, ny)
Z = np.linspace(-depth, 0, nz)

gridToVTK("./Exudates_day"+("{:02d}".format(ii+1)), X, Y, Z, pointData = {"Exudates":C})


