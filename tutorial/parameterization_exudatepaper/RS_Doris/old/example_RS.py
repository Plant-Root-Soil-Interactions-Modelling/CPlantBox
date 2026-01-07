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
path = "rootsystem/"
name = "RS_optimized_L_WT"
rs.readParameters(path + name + ".xml")

# Initialize
rs.initializeLB(5,4)

times = [42, 63, 98]
comp_length = np.zeros((len(times)))
comp_diam = np.zeros((len(times)))
simtime = 1
dummy = 0
radiishare = np.zeros((len(times), 5))
diffradii = [0.01,  0.013, 0.075,0.1,  0.2]

for j in range(0,times[-1]+1): 
    rs.simulate(simtime, True);

    if j in times: 
        length = np.array(rs.getParameter("length"))
        radius = np.array(rs.getParameter("radius"))
        meanrad = np.sum(length*radius)/np.sum(length)
        #plt.hist(radius, bins='auto')
        #plt.show()

        #diffradii = np.sort(np.unique(radius))
        diffradii = np.sort(np.unique(radius))
        for i in range(0,len(diffradii)):
            radiishare[dummy,i] = np.round(np.sum(length[radius == diffradii[i]])/np.sum(length)*100,2)
            
        
        comp_length[dummy] =  np.sum(length)
        comp_diam[dummy] =  meanrad*20
        rs.write("test_ex/day_"+str(j)+".vtp")
        dummy = dummy+1
    

print('computed length', comp_length) 
print('real length', np.array([1259, 12946, 55529, 51365]))

sys.exit()

print('computed diameter', comp_diam) 
print('real diameter', np.array([0.41, 0.32, 0.32, 0.30]))


print('these are the radii: ', diffradii)
print(radiishare) 
