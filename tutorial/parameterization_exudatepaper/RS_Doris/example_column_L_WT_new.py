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


fig, ax = plt.subplots(3,1, figsize = (18, 10))

rs = pb.RootSystem()

# Open plant and root parameter from a file
path = "rootsystem/"
name = "RS_optimized_column_L_WT_new"
rs.readParameters(path + name + ".xml")

#data
times_ = [42, 63, 98, 154]

df = pd.read_csv("data/column_experiment_mean.csv")
data = df[(df['substrate']=='L') & (df['genotype']=="WT")]
DAS = data["DAS"].loc[:].values
real_length = data['RL_mean']
real_length_SE = data['RL_SE']
real_diam = data['RD_mean']
real_diam_SE = data['RD_SE']
real_RSA =data['RSA_mean']
real_RSA_SE = data['RSA_SE']


# Initialize
rs.initializeLB(5,4)

#times = [42, 63, 98, 154]
times = np.linspace(0,154,31)
times = times.astype(int)
#print(times)   
comp_length = np.zeros((len(times)))
comp_length_test = np.zeros((len(times)))
comp_diam = np.zeros((len(times)))
comp_RSA = np.zeros((len(times)))
comp_RSA2 = np.zeros((len(times)))
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
        surf = np.sum(2*radius*math.pi*length)
        #plt.hist(radius, bins='auto')
        #plt.show()

        #diffradii = np.sort(np.unique(radius))
        diffradii = np.sort(np.unique(radius))
        for i in range(0,len(diffradii)):
            radiishare[dummy,i] = np.round(np.sum(length[radius == diffradii[i]])/np.sum(length)*100,2)
            
        comp_length[dummy] =  np.sum(length)
        comp_diam[dummy] =  meanrad*20
        comp_RSA[dummy] = surf
        comp_RSA2[dummy] = 2*meanrad*math.pi*np.sum(length)
        rs.write("test_ex/L_WT_day_"+str(j)+".vtp")
        dummy = dummy+1
    

print('computed length', comp_length) 
print('real length', real_length)
print('computed RSA', comp_RSA, comp_RSA2)
print('real RSA', real_RSA)

ax[0].plot(times, comp_length, linestyle = '--')
ax[0].errorbar(times_, real_length, real_length_SE, color ='r', fmt='o')

ax[1].plot(times, comp_diam, linestyle = '--')
ax[1].errorbar(times_, real_diam, real_diam_SE, color ='r', fmt='o')
ax[1].set_ylim([0,0.5])

ax[2].plot(times, comp_RSA, linestyle = '--')
ax[2].plot(times, comp_RSA2, linestyle = ':')
ax[2].errorbar(times_, real_RSA, real_RSA_SE, color ='r', fmt='o')

plt.show()
sys.exit()

print('computed diameter', comp_diam) 
print('real diameter', real_diam)


print('these are the radii: ', diffradii)
print(radiishare) 
