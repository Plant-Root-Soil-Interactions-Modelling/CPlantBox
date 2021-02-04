import sys; 
sys.path.append("../../../../../..")
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
from pyevtk.hl import gridToVTK
import plantbox as rb
import multiprocessing
from multiprocessing import Process, active_children
import psutil
from threading import Thread                                                       

def calc_model(number,q):
    q.put([(model.calculate(ii+1,number,number+1))])
    
############################
#
# Root system
#
rs = rb.RootSystem()

path = "../../"
name = "Zea_mays_exudation"  
rs.readParameters(path + name + ".xml")

srp = rb.SeedRandomParameter(rs)  # with default values
srp.seedPos = rb.Vector3d(0., 0., -3.)  # [cm] seed position

#set geometry 
width = 40  # cm
depth = 35    
soilcore = rb.SDF_PlantContainer(width, width, depth, True)
rs.setGeometry(soilcore)  

rs.setSeed(0)
rs.initialize()

for ii in range(0,20): 
    simtime = 1
    rs.simulate(simtime, True);

#
# Grid parameter
#
nodes = np.array([np.array(n) for n in rs.getNodes()])
np.save("../nodes/day"+str(ii+1), nodes)
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
model.Q = 3.7  # µg/d/tip
model.Dl = 0.171  # cm2/d
model.theta = 0.3 #-
model.R = 16.7  # -
model.k = 1.42  # d-1
model.l = 5  # cm (for line source only)

#
# Numerical parameter
#
model.type = rb.IntegrationType.mls;  # mps, mps_straight, mls
model.n0 = 10  # integration points per cm
model.thresh13 = 1.e-15;  # threshold to neglect diffusing g (eqn 13)
model.calc13 = True;  # turns Eqn 13  on (True) and off (False)
model.observationRadius = 0.8;  # limits computational domain around roots [cm]

#
#set threshold 
#
thresh = int(55)    
t = time.time()
print("make voxel lists")
model.makeVoxelLists()

order = np.array(rs.getParameter("type"))                       
order[order==4]=1
roots = rs.getPolylines()

#make lists of the different orders
numo = np.unique(order)
numo = numo.astype(int)
ords= []
for i in range(0,len(numo)):
    ords.append(np.where(order==numo[i])[0])
# preallocate concentration matrizes
C = np.zeros((nx*ny*nz,))
C1lat = np.zeros((nx*ny*nz,))
C2lat = np.zeros((nx*ny*nz,))
Cthresh = np.zeros((nx*ny*nz,))
for r in range(0,len(numo)):
    segnum = np.zeros(len(ords[r]))
    for i in range(0,len(ords[r])):
       segnum[i] = len(roots[ords[r][i]])
    max = 50 #maximum number of segments that should be processd in parallel 
    tranches = int(np.ceil(np.sum(segnum)/max))
    task = np.zeros(tranches)

    nn = 0
    i=0
    while i<len(ords[r]):
        segtest = 0
        while (segtest<=max):
            segtest = segtest+segnum[i]
            i = i+1
            if i==len(ords[r]):
                break
        task[nn] = int(i)-np.sum(task)
        nn=nn+1

    for j in range(0,tranches): 
        q = multiprocessing.Queue()
        if j== tranches-1: 
            end = len(ords[r])
        else:
            end = int(np.sum(task[0:(j+1)]))
        
        for i in range(int(np.sum(task[0:j])),end):
            p = multiprocessing.Process(target=calc_model, args=(ords[r][i],q,))
            p.start()
            
        for i in range(int(np.sum(task[0:j])),end):
            add = np.array(q.get())
            C = np.add(add, C) 
            if numo[r]<2:
                C1lat = np.add(add, C1lat) 
            if numo[r]<3:
                C2lat = np.add(add, C2lat) 
            Cthresh = np.add(add, Cthresh) 
            Cthresh[Cthresh<thresh] = 0    

elapsed = time.time() - t
print("Computation took", elapsed, "s")

#
# post processing...
#
C_ = np.zeros((nx,ny,nz))
C1lat_ = np.zeros((nx,ny,nz))
C2lat_ = np.zeros((nx,ny,nz))
Cthresh_ = np.zeros((nx,ny,nz))

C = np.reshape(C, (nz, ny, nx))
C1lat = np.reshape(C1lat, (nz, ny, nx))  
C2lat = np.reshape(C2lat, (nz, ny, nx))  
Cthresh = np.reshape(Cthresh, (nz, ny, nx))  
for i in range(0,np.shape(C)[0]):
    for j in range(0,np.shape(C)[1]):
        for k in range(0,np.shape(C)[2]):
            C_[k,j,i] = C[i,j,k]
            C1lat_[k,j,i] = C1lat[i,j,k]
            C2lat_[k,j,i] = C2lat[i,j,k]
            Cthresh_[k,j,i] = Cthresh[i,j,k]
            
del C, Cthresh, C2lat, C1lat
C = C_
C1lat = C1lat_
C2lat = C2lat_
Cthresh = Cthresh_

np.savez_compressed("../concentration/day"+("{:02d}".format(ii+1)), C, C2lat, Cthresh, C1lat)

X = np.linspace(-width / 2, width / 2, nx)
Y = np.linspace(-width / 2, width / 2, ny)
Z = np.linspace(-depth, 0, nz)


gridToVTK("../exud/./Exudates_day"+("{:02d}".format(ii+1)), X, Y, Z, pointData = {"Exudates":C})
gridToVTK("../exud/./Exudates_C1lat_day"+("{:02d}".format(ii+1)), X, Y, Z, pointData = {"Exudates":C1lat})
gridToVTK("../exud/./Exudates_C2lat_day"+("{:02d}".format(ii+1)), X, Y, Z, pointData = {"Exudates":C2lat})
gridToVTK("../exud/./Exudates_thresh_day"+("{:02d}".format(ii+1)), X, Y, Z, pointData = {"Exudates":Cthresh})

rs.write("../vtp/day"+("{:02d}".format(ii+1)) + ".vtp")