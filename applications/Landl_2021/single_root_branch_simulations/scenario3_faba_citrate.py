import sys;
sys.path.append("../../..")
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

def calc_model(number,q):
    q.put([(model.calculate(ii+1,number,number+1))])
    
############################

name_ = ['branch1', 'branch2']
for jj in range(0,len(name_)):

    #
    # Root system
    #
    rs = rb.RootSystem()

    path = "../../../modelparameter/rootsystem/" 
    rs.readParameters(path + name_[jj] + ".xml")

    #set geometry 
    width = 4  # cm
    depth = 15    
    soilcore = rb.SDF_PlantContainer(width, width, depth, True)
    rs.setGeometry(soilcore)  
    rs.setSeed(0)

    rs.initialize()

    for ii in range(0,10): 
        simtime = 1
        rs.simulate(simtime, True);
    rs.write("vtp/faba_citrate/"+name_[jj] +"_day"+("{:02d}".format(ii+1))+".vtp")

    #
    # Grid parameter
    #
    nodes = np.array([np.array(n) for n in rs.getNodes()])
    np.save("nodes/faba_citrate/"+name_[jj] + "_day"+str(ii+1), nodes)
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
    model.Q = 18.4  # µg/d/tip
    model.Dl = 0.171  # cm2/d
    model.theta = 0.3	 #-
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
    model.observationRadius = 5;  # limits computational domain around roots [cm]

    t = time.time()
    model.makeVoxelLists()                               
    roots = rs.getPolylines()
    segnum = np.zeros(len(roots))
    for i in range(0,len(roots)):
        segnum[i] = len(roots[i])
    max = 500 #maximum number of segments that should be processd in parallel 
    tranches = int(np.ceil(np.sum(segnum)/max))
    task = np.zeros(tranches)

    nn = 0
    i=0
    while i<len(roots):
        segtest = 0
        while (segtest<=max):
            segtest = segtest+segnum[i]
            i = i+1
            if i==len(roots):
                break
        task[nn] = int(i)-np.sum(task)
        nn=nn+1

    # preallocate concentration matrizes
    C = np.zeros((nx*ny*nz,))  
    for j in range(0,tranches): 
        q = multiprocessing.Queue()
        if j== tranches-1: 
            end = len(roots)
        else:
            end = int(np.sum(task[0:(j+1)]))
        
        for i in range(int(np.sum(task[0:j])),end):
            p = multiprocessing.Process(target=calc_model, args=(i,q,))
            p.start()
            
        for i in range(int(np.sum(task[0:j])),end):
            C = np.add(C,q.get()) 
        
    elapsed = time.time() - t
    print("Computation took", elapsed, "s")

    C_ = np.zeros((nx,ny,nz))
    C = np.reshape(C, (nz, ny, nx))
    for i in range(0,np.shape(C)[0]):
        for j in range(0,np.shape(C)[1]):
            for k in range(0,np.shape(C)[2]):
                C_[k,j,i] = C[i,j,k]
                
    del C
    C = C_   
    
    C = np.reshape(C, (nx, ny, nz))
    np.save("concentration/faba_citrate/"+name_[jj] +"_day"+("{:02d}".format(ii+1)), C)

    X = np.linspace(-width / 2, width / 2, nx)
    Y = np.linspace(-width / 2, width / 2, ny)
    Z = np.linspace(-depth, 0, nz)

    gridToVTK("exud/faba_citrate/./Exudates_"+name_[jj] + "_day"+("{:02d}".format(ii+1)), X, Y, Z, pointData = {"Exudates":C})
