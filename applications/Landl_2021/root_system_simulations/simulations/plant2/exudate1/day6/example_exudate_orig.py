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

#set geometry 
width = 40  # cm
depth = 70    
soilcore = rb.SDF_PlantContainer(width, width, depth, True)
rs.setGeometry(soilcore)  

rs.setSeed(0)
rs.initialize()

for ii in range(0,6): 
    simtime = 1
    rs.simulate(simtime, True);
rs.write("../vtp/day"+("{:02d}".format(ii+1)) + ".vtp")

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
model = rb.ExudationModel(width, width, depth, nx, ny, nz, rs)
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

t = time.time()
root = rs.getPolylines()
maxt = 40
tasks = int(np.ceil(len(root)/maxt)) #should not launch more than maxt tasks at once 

C_ = np.zeros(nx*ny*nz)
for j in range(0,tasks): 
    m = multiprocessing.Manager()
    q = m.Queue()
    if j== tasks-1: 
        end = len(root)
    else:
        end = (j+1)*maxt
    
    for i in range(j*maxt,end):
        p = multiprocessing.Process(target=calc_model, args=(i,q,))
        p.start()
        
    for i in range(j*maxt,end):
        C_ = np.add(C_,np.array(q.get()))

elapsed = time.time() - t
print("Computation took", elapsed, "s")

#
# post processing...
#
C = np.reshape(C_, (nx, ny, nz))
np.save("../concentration/day"+("{:02d}".format(ii+1)), C)

X = np.linspace(-width / 2, width / 2, nx)
Y = np.linspace(-width / 2, width / 2, ny)
Z = np.linspace(-depth, 0, nz)

#print("max" , np.max(C[:]), "min", np.min(C[:]))
#num_th = (C > 0).sum()  # number of points for which concentration is larger than threshold
#print("volume of concentration above threshold: ", num_th * width / nx * width / ny * depth / nz)
#print("this is", num_th / (nx * ny * nz) * 100, "% of the overall volume")

gridToVTK("../exud/./Exudates_day"+("{:02d}".format(ii+1)), X, Y, Z, pointData = {"Exudates":C})
