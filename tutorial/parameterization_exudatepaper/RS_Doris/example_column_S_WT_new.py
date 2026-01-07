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


fig, ax = plt.subplots(2,1, figsize = (18, 10))

rs = pb.RootSystem()

# Open plant and root parameter from a file
path = "rootsystem/"
name = "RS_optimized_column_S_WT_new"
rs.readParameters(path + name + ".xml")

#data
times_ = [42, 63, 98, 154]
real_length = np.array([1272.200082, 15346.47835, 55587.97799, 54147.01866])
real_length_SE = np.array([272.2051323, 2879.332343, 4303.265028, 4974.350854])
real_diam = np.array([0.380962036, 0.262333815, 0.283614617, 0.265687526])
real_diam_SE = np.array([0.027869171, 0.004690655, 0.010294733, 0.002938094])

# Initialize
rs.initializeLB(5,4)

#times = [42, 63, 98, 154]
times = np.linspace(0,154,31)
times = times.astype(int)
#print(times)   
comp_length = np.zeros((len(times)))
comp_length_test = np.zeros((len(times)))
comp_diam = np.zeros((len(times)))
simtime = 1
dummy = 0
radiishare = np.zeros((len(times), 5))
diffradii = [0.01,  0.013, 0.075,0.1,  0.2]

for j in range(0,times[-1]+1):

    if j<=98: 
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
        rs.write("test_ex/S_WT_day_"+str(j)+".vtp")
        dummy = dummy+1
    

print('computed length', comp_length) 
print('real length', real_length)

ax[0].plot(times, comp_length, linestyle = '--')
ax[0].errorbar(times_, real_length, real_length_SE, color ='r', fmt='o')

ax[1].plot(times, comp_diam, linestyle = '--')
ax[1].errorbar(times_, real_diam, real_diam_SE, color ='r', fmt='o')
ax[1].set_ylim([0,0.5])

plt.show()
sys.exit()

print('computed diameter', comp_diam) 
print('real diameter', real_diam)


print('these are the radii: ', diffradii)
print(radiishare) 
