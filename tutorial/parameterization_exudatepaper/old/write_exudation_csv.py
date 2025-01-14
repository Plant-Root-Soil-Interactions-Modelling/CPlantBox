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
import csv


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
cols = ['b', 'g', 'r', 'c', 'm']

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

DAS = np.hstack((0, df['DAS'].values))

for i in range(0,len(exudates)):
    prims = np.hstack((prim[i,0]*10**3,prim[i,:]*10**3,prim[i,-1]*10**3))
    lat1s = np.hstack((lat1[i,0]*10**3,lat1[i,:]*10**3,lat1[i,-1]*10**3))
    lat2s = np.hstack((lat2[i,0]*10**3,lat2[i,:]*10**3,lat2[i,-1]*10**3))
    crowns = np.hstack((crown[i,0]*10**3,crown[i,:]*10**3,crown[i,-1]*10**3))
    basals = np.hstack((basal[i,0]*10**3,basal[i,:]*10**3,basal[i,-1]*10**3))

    df = pd.DataFrame(zip(DAS, prims, lat1s, lat2s, crowns, basals), columns=['DAS','prim','lat1', 'lat2', 'crown', 'basal'])
    df.to_csv('../data_magda/'+ex_types[i]+'_2019.csv') 


sys.exit()
    
dicts = [
    {
        "DAS": DAS,
     "prim": prims,
        "lat1": lat1s,
        "lat2": lat2s,
        "crown": crowns,
        "basal": basals,
    }
    for DAS, prims, lat1s, lat2s, crowns, basals in zip(DAS, prims, lat1s, lat2s, crowns, basals)
]
print(dicts[:]) 
with open('../data_magda/'+ex_types[i]+'_2019.csv', 'w') as f:
    writer = csv.DictWriter(f, list[dicts[:]])
    writer.writeheader()
    writer.writerows(dicts)

        

sys.exit()
fig, ax = plt.subplots(2,2, figsize = (18, 10))

for i in range(0,len(exudates)):
    ax[int(np.floor(i/2)), int(i%2)].plot(df['DAS'], lat2[i,:]*10**3, color = cols[0], label = '2nd order lateral')
    ax[int(np.floor(i/2)), int(i%2)].plot(df['DAS'], lat1[i,:]*10**3, color = cols[1], label = '1st order lateral')
    ax[int(np.floor(i/2)), int(i%2)].plot(df['DAS'], basal[i,:]*10**3, color = cols[2], label = 'basal root')
    ax[int(np.floor(i/2)), int(i%2)].plot(df['DAS'], crown[i,:]*10**3, color = cols[3], label = 'crown root')
    ax[int(np.floor(i/2)), int(i%2)].plot(df['DAS'], prim[i,:]*10**3, color = cols[4], label = 'primary root')
    if i == len(exudates)-1: 
        ax[int(np.floor(i/2)), int(i%2)].set_xlabel('Time (d)')
    ax[int(np.floor(i/2)), int(i%2)].set_ylabel(ex_types[i] +'\n ($\mu$mol/ h/ root tip)')


handles, labels = ax[0,0].get_legend_handles_labels()
order = [4,3,2,1,0]
ax[0,0].legend([handles[idx] for idx in order],[labels[idx] for idx in order])
#plt.tight_layout()

ax[1,1].axis('off')
plt.show()
