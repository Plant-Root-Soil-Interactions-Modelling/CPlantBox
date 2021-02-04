import sys;
sys.path.append("../../..")
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
from pyevtk.hl import gridToVTK
import plantbox as rb

simtime = [13,20]
for st in simtime:
    #
    # Root system
    #
    rs = rb.RootSystem()

    path = "../../../modelparameter/rootsystem/"
    name = "straight_root"  
    rs.readParameters(path + name + ".xml")

    rg = 2
    for p in rs.getRootRandomParameter():
        p.r = rg  # growth rate

    #set geometry 
    width = 3  # cm
    depth = 30    
    soilcore = rb.SDF_PlantContainer(width, width, depth, True)
    rs.setGeometry(soilcore)  
    rs.setSeed(0)

    rs.initialize()

    for ii in range(0,st): 
        simtime = 1
        rs.simulate(simtime, True);
    rs.write("vtp/faba_citrate/rg"+str(rg) +"_day"+("{:02d}".format(ii+1))+".vtp")
    #
    # Grid parameter
    #
    nodes = np.array([np.array(n) for n in rs.getNodes()])
    np.save("nodes/faba_citrate/day"+str(ii+1), nodes)
    xres = 0.05;
    yres = 0.05;
    zres = 0.05;
    nx = int(width / xres);
    ny = int(width / yres);
    nz = int(depth / zres);
    print("Width", width, "Depth", depth, " at a Resolution", nx, ny, nz)

    #
    # Model parameter
    #
    model = rb.ExudationModel(width, width, depth, nx, ny, nz, rs)
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
    model.observationRadius = 2;  # limits computational domain around roots [cm]

    t = time.time()
    C = model.calculate(ii+1)
    elapsed = time.time() - t
    print("Computation took", elapsed, "s")

    C = np.reshape(C, (nx, ny, nz))
    np.save("concentration/faba_citrate/day"+("{:02d}".format(ii+1)), C)

    X = np.linspace(-width / 2, width / 2, nx)
    Y = np.linspace(-width / 2, width / 2, ny)
    Z = np.linspace(-depth, 0, nz)

    gridToVTK("exud/faba_citrate/./Exudates_day"+("{:02d}".format(ii+1)), X, Y, Z, pointData = {"Exudates":C})
