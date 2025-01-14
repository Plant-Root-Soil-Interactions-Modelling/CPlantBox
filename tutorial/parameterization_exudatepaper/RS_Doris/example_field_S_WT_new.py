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
name = "RS_optimized_field_S_WT_new"
rs.readParameters(path + name + ".xml")

#data
times_ = [42, 63, 98, 154]
real_length = np.array([3036, 40550, 106061, 137579])
real_length_SE = np.array([1237, 7572, 14495, 10183])
real_diam = np.array([0.31, 0.27, 0.25, 0.24])
real_diam_SE = np.array([0.019,0.010, 0.007, 0.009])

# Initialize
rs.initializeLB(5,4)

#times = [42, 63, 98, 154]
times = np.linspace(0,154,31)
times = times.astype(int)
print(times)
comp_length = np.zeros((len(times)))
comp_diam = np.zeros((len(times)))
simtime = 1
dummy = 0
radiishare = np.zeros((len(times), 5))
diffradii = [0.01,  0.013, 0.075,0.1,  0.2]
depth = 60

fig, ax = plt.subplots(2,1, figsize = (18, 10))

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
            
        ana = pb.SegmentAnalyser(rs)
        comp_length[dummy] = np.sum(ana.distribution("length", 0., -depth, 1, True))
        #comp_length[dummy] =  np.sum(length)
        comp_diam[dummy] =  meanrad*20
        rs.write("test_ex/S_WT_day_"+str(j)+".vtp")
        dummy = dummy+1
    

print('computed length', comp_length) 
print('real length', real_length)

ax[0].plot(times, comp_length, linestyle = '--')
#ax[0].plot(times_, real_length, 'o')
ax[0].errorbar(times_, real_length, real_length_SE, color ='r', fmt='o')
ax[1].plot(times, comp_diam, linestyle = '--')
ax[1].errorbar(times_, real_diam, real_diam_SE, color ='r', fmt='o')
plt.show()
sys.exit()

print('computed diameter', comp_diam) 
print('real diameter', np.array([0.41, 0.32, 0.32, 0.30]))


print('these are the radii: ', diffradii)
print(radiishare) 
