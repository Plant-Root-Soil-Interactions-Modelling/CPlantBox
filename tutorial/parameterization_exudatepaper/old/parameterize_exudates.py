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
exudate = 'Phenolics'
df = pd.read_csv("../data_magda/mean_measures_Eva_Michael.csv")
DAS = df["DAS"].loc[:].values
mean = df[exudate]*12*24/10**3
std = df[exudate+' std']*12*24/10**3

#create root system 
#rs = pb.RootSystem()
rs = pb.MappedRootSystem()
path = "rootsystem_params/"
name = "Zea_mays_5_Leitner_Streuber_2014_mod2"
rs.readParameters(path + name + ".xml")

# Initialize
rs.initializeLB(5,4)

times = [42, 63, 98, 154]
plant_exud = np.zeros((len(times)))
exu_prop1 = 0.00045
exu_prop2 = 0.00008
exu_prop3 = 0.00041
exu_prop4 = 0.0003
dummy = 0

#exudation rates per root type
#[age,age] [value,  value] for each time 
kex = np.array([[[0., 2.], [exu_prop1, 0.]],[[0., 2.], [exu_prop2, 0.]],[[0., 2.], [exu_prop3, 0.]], [[0., 2.], [exu_prop4, 0.]]])



for j in range(0,times[-1]+1):

    if j<=times[2]: 
        rs.simulate(1, True);

    if j in times:

        segs = rs.segments
        ages =  np.asarray(j - np.array([rs.nodeCTs[int(seg.y)] for seg in segs]), int)
        types = np.asarray(rs.subTypes, int)
        a = rs.radii
        l = rs.segLength()
        sf = np.zeros(len(segs),)

        for i, s in enumerate(segs):
            if types[i] == 0:
                kexu_ = 0
            else: 
                kex_ = kex[dummy]
                age = ages[i]
                if (age<0) and (age>kex_[0][1]):
                    kexu_ = 0
                else:
                    kexu_ = (kex_[1][1]-kex_[1][0])/(kex_[0][1]-kex_[0][0])*age+kex_[1][0]
                kexu = max(0,kexu_)

                if kexu!=0: 
                    sf[i] = 2 * np.pi * a[i] * l[i] * 1.e-4 * kexu * 1.e3 # g/day/root

        
        plant_exud[dummy] = np.sum(sf) # g/day/plant
        dummy = dummy+1

#plant_exud[-1] = plant_exud[-2]
print('computed exudation', plant_exud)
print('measured exudation', mean) 

#PLOT
plt.plot(df['DAS'], mean, color = 'b', label = 'Measurement')
plt.fill_between(df['DAS'], mean-std, mean+std, color = 'b', alpha = 0.2)
plt.plot(df['DAS'], plant_exud, color = 'b',linestyle = '--', label = 'Simulation')
plt.ylabel(exudate + '\n (g / plant / d)')
plt.legend()
#plt.text(50, 33, '(c)', fontsize=16,  va='top', ha='right')
plt.tight_layout()
plt.show()

sys.exit()
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
    ax2.ylabel("Number of root tips (-)")
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
