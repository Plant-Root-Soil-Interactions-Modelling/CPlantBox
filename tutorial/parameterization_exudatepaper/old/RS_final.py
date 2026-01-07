"""small example"""
import sys
sys.path.append("../../../../CPlantBox")
sys.path.append("../../../../CPlantBox/src")
import plantbox as pb
import visualisation.vtk_plot as vp
import vtk 
import math
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import pandas as pd
import timeit

rs = pb.RootSystem()

# Open plant and root parameter from a file
path = "rootsystem_params/"
name = "Zea_mays_5_Leitner_Streuber_2014_mod2"
rs.readParameters(path + name + ".xml")
min_b = [-10., -22, -74.] 
max_b = [10., 22, 0.]
rs.setGeometry(pb.SDF_PlantBox(max_b[0]*2, max_b[1]*2, np.abs(min_b[2])))

# Initialize
rs.initializeLB(5,4)
rs.setSeed(1)
times = [42, 63, 98]
simtime = 1
dummy = 0
for j in range(0,times[-1]+1): 
    rs.simulate(simtime, True);
    rs.write('vtps/day'+str(j)+'.vtp')
