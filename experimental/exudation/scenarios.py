import sys;
sys.path.append("../..")
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
from pyevtk.hl import gridToVTK
import plantbox as rb


simtime = np.linspace(5, 15, num=11)
for st in range(0, len(simtime)):
    
    print("\nIteration", st)
    
    #
    # Root system
    #
    rs = rb.RootSystem()

    path = "../../modelparameter/rootsystem/"
    name = "test_root"  
    rs.readParameters(path + name + ".xml")

    #set geometry 
    width = 4  # cm
    depth = 15    
    soilcore = rb.SDF_PlantContainer(width, width, depth, True)
    rs.setGeometry(soilcore)  
    rs.setSeed(0)

    rs.initialize()

    rs.simulate(simtime[st], True);
    rs.write("vtp/day_"+str(simtime[st]) +".vtp")
    
    #
    # Grid parameter
    #
    nodes = np.array([np.array(n) for n in rs.getNodes()])
    np.save("nodes/day"+str(simtime[st]), nodes)
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
    model.thresh13 = 0.# 1.e-15;  # threshold to neglect diffusing g (eqn 13)
    model.calc13 = True;  # turns Eqn 13  on (True) and off (False)
    model.observationRadius = 5;  # limits computational domain around roots [cm]

    t = time.time()
    C = model.calculate(simtime[st])
    elapsed = time.time() - t
    print("Computation took", elapsed, "s")

    C = np.reshape(C, (nx, ny, nz))
    np.save("concentration/day"+str(simtime[st]), C)
    print("concentration maximum", np.max(C.flat))

    X = np.linspace(-width / 2, width / 2, nx)
    Y = np.linspace(-width / 2, width / 2, ny)
    Z = np.linspace(-depth, 0, nz)

    gridToVTK("exud/./Exudates_day"+str(simtime[st]), X, Y, Z, pointData = {"Exudates":C})
