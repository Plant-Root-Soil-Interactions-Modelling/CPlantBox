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


font = {'size'   : 18}
plt.rc('font', **font)

#import exudate data
df = pd.read_csv("../data_magda/mean_measures_Eva_Michael.csv")
#print(df.columns.values.tolist()) 
DAS = df["DAS"].loc[:].values
a = df['Phenolics'].loc[:].values
b = df['Sugars'].loc[:].values
c = df['AminoAcids'].loc[:].values
exudates = np.vstack((a,b,c))
ex_types = ['Phenolics', 'Sugars', 'Amino Acids']

lat2 = np.zeros((len(exudates),len(DAS)-1))
lat1 = np.zeros((len(exudates),len(DAS)-1))
basal = np.zeros((len(exudates),len(DAS)-1))
crown = np.zeros((len(exudates),len(DAS)-1))
prim = np.zeros((len(exudates),len(DAS)-1))
numtips = np.zeros((len(DAS)-1))
cols = ['b', 'g', 'y', 'c', 'm']

#create root system 
rs = pb.RootSystem()
path = "rootsystem_params/"
name = "Zea_mays_5_Leitner_Streuber_2014_mod2"
rs.readParameters(path + name + ".xml")

# Initialize
rs.initializeLB(5,4)

times = [42, 63, 98]
comp_length = np.zeros((len(times)))
comp_diam = np.zeros((len(times)))
simtime = 1
dummy = 0
radiishare = np.zeros((len(times), 5))
for j in range(0,times[-1]+1): 
    rs.simulate(simtime, True);

    if j in times:
        
        radius = np.array(rs.getParameter("radius"))
        radii = np.sort(np.unique(radius))
        print(radii)
        numtips[dummy] = len(radius) 
        weight = 0
        for i in range(0,len(radii)):
            weight = weight + radii[i]*len(np.where(radius == radii[i])[0])

        for i in range(0,len(exudates)): 
            lat2[i,dummy] = exudates[i,dummy]/weight*radii[0]
            lat1[i,dummy] = exudates[i,dummy]/weight*radii[1]
            basal[i,dummy] = exudates[i,dummy]/weight*radii[2]
            crown[i,dummy] = exudates[i,dummy]/weight*radii[2]
            prim[i,dummy] = exudates[i,dummy]/weight*radii[3]
            
        dummy = dummy+1


lat2_ = np.zeros((len(exudates),len(DAS)))
lat1_ = np.zeros((len(exudates),len(DAS)))
basal_ = np.zeros((len(exudates),len(DAS)))
crown_ = np.zeros((len(exudates),len(DAS)))
prim_ = np.zeros((len(exudates),len(DAS)))
for i in range(0,len(exudates)):
    lat2_[i,:] = np.hstack((lat2[i,:],lat2[i,-1]))
    lat1_[i,:] = np.hstack((lat1[i,:],lat1[i,-1]))
    basal_[i,:] = np.hstack((basal[i,:],basal[i,-1]))
    crown_[i,:] = np.hstack((crown[i,:],crown[i,-1]))
    prim_[i,:] = np.hstack((prim[i,:],prim[i,-1]))

DAS_ = df['DAS']
#DAS_ = np.hstack((0,DAS[:]))
numtips_ = np.hstack((numtips, numtips[-1]))   

fig, ax = plt.subplots(2,2, figsize = (18, 10))

for i in range(0,len(exudates)):
    ax[int(np.floor(i/2)), int(i%2)].plot(DAS_, lat2_[i,:]*10**3, color = cols[0], label = '2nd order lateral')
    ax[int(np.floor(i/2)), int(i%2)].plot(DAS_, lat1_[i,:]*10**3, color = cols[1], label = '1st order lateral')
    ax[int(np.floor(i/2)), int(i%2)].plot(DAS_, basal_[i,:]*10**3, color = cols[2], label = 'seminal root')
    ax[int(np.floor(i/2)), int(i%2)].plot(DAS_, crown_[i,:]*10**3, color = cols[3], label = 'crown root')
    ax[int(np.floor(i/2)), int(i%2)].plot(DAS_, prim_[i,:]*10**3, color = cols[4], label = 'primary root')
    if i == len(exudates)-1: 
        ax[int(np.floor(i/2)), int(i%2)].set_xlabel('Time (d)')
    ax[int(np.floor(i/2)), int(i%2)].set_ylabel(ex_types[i] +'\n ($\mu$mol/ h/ root tip)')
    #ax[int(np.floor(i/2)), int(i%2)].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    
    ax2 = ax[int(np.floor(i/2)), int(i%2)].twinx()
    ax2.plot(DAS_, numtips_,label = 'number of root tips', linestyle = '--', color = 'r')
    ax2.set_ylabel("Number of root tips (-)")
    ax2.spines['right'].set_color('red')
    ax2.tick_params(axis='y', colors='red')
    ax2.yaxis.label.set_color('red')
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(4,4))
    

handles, labels = ax[0,0].get_legend_handles_labels()
order = [4,3,2,1,0]
ax[0,0].legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.tight_layout()

ax[1,1].axis('off')
plt.show()
