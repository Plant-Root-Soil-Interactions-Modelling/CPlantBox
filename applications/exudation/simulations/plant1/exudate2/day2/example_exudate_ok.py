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
import psutil

def calc_model(number,q):
    q.put([(model.calculate(ii+1,number,number+1))])
    
############################
#
# Root system
#
rs = rb.RootSystem()

path = "../../"
name = "Faba_exudation" 
rs.readParameters(path + name + ".xml")

for p in rs.getRootRandomParameter():
    p.gf = 2  # linear growth function

#set geometry 
width = 20  # cm
depth = 90    
soilcore = rb.SDF_PlantContainer(width, width, depth, True)
rs.setGeometry(soilcore)  

rs.setSeed(0)
rs.initialize()

for ii in range(0,2): 
    simtime = 1
    rs.simulate(simtime, True);
rs.write("../vtp/day"+str(ii+1) + ".vtp")

#
# Grid parameter
#
nodes = np.array([np.array(n) for n in rs.getNodes()])
np.save("../nodes/day"+str(ii+1), nodes)
xres = 0.3;
yres = 0.3;
zres = 0.3;
nx = int(width / xres);
ny = int(width / yres);
nz = int(depth / zres);
print("Width", width, "Depth", depth, " at a Resolution", nx, ny, nz)

#
# Model parameter
#
model = rb.ExudationModel(width, width, depth, nx, ny, nz, rs)
model.Q = 33.38 # µg/d/tip
model.Dl = 3.64e-3  # cm2/d
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
model.observationRadius = 5;  # limits computational domain around roots [cm]

t = time.time()
root = rs.getPolylines()
print(psutil.virtual_memory())
sys.exit()

q = multiprocessing.Queue()
for i in range(len(root)):
    p = multiprocessing.Process(target=calc_model, args=(i,q,))
    p.start()
    
C_ =  np.array(q.get())
for i in range(1,len(root)):
    C_ = np.add(C_,np.array(q.get()))

elapsed = time.time() - t
print("Computation took", elapsed, "s")

#
# post processing...
#
C = np.reshape(C_, (nx, ny, nz))
np.save("../concentration/day"+str(ii+1), C)

X = np.linspace(-width / 2, width / 2, nx)
Y = np.linspace(-width / 2, width / 2, ny)
Z = np.linspace(-depth, 0, nz)

#print("max" , np.max(C[:]), "min", np.min(C[:]))
#num_th = (C > 0).sum()  # number of points for which concentration is larger than threshold
#print("volume of concentration above threshold: ", num_th * width / nx * width / ny * depth / nz)
#print("this is", num_th / (nx * ny * nz) * 100, "% of the overall volume")

gridToVTK("../exud/./Exudates_day"+str(ii+1), X, Y, Z, pointData = {"Exudates":C})


