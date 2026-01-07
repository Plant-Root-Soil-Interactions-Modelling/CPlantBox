"""small example"""
import sys
sys.path.append("../../../../../CPlantBox")
sys.path.append("../../../../../CPlantBox/src")
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

#colors & labels 
a = [0, 0, 1]
b = [106/255, 166/255, 1]
c = [1, 0, 0]
d = [1, 128/255, 0]
cols = [a,b,c,d]
label_root = ["L_WT", "S_WT", "L_rth3","S_rth3"]
soils = ['L', 'S']
gt = ['WT', 'rth3']
linst = ['-', '--','-.',':']
height = [20, 20, 35]
DAS = [42, 63, 98, 154]
area = 20*45 #cm 


#import RLD data
df = pd.read_csv("data/RLD_Doris.csv")
DAS = df["DAS"].loc[:].values
soils = df["substrate"].loc[:].values
gt = df["genotype"].loc[:].values
RL_mean = df["RL_mean"].loc[:].values
RL_SE = df["RL_SE"].loc[:].values

soils_ = np.unique(soils)
gt_ = np.unique(gt)
#BBCH_ = np.unique(BBCH)

fig, ax = plt.subplots(2,1, figsize = (18, 10))

for i in range(0, len(soils_)):
    for j in range(0, len(gt_)):

        data = df[(df['substrate']==soils_[i]) & (df['genotype']==gt_[j])]
        x = data["DAS"].loc[:].values
        y = data["RL_mean"].loc[:].values
        y1 = data["RL_mean"].loc[:].values-data["RL_SE"].loc[:].values
        y2 = data["RL_mean"].loc[:].values+data["RL_SE"].loc[:].values
            
        ax[0].plot(data["DAS"],data["RL_mean"]/100,color = cols[i+2*j], label = label_root[i+2*j])
        ax[0].fill_between(data['DAS'], (data['RL_mean']-data['RL_SE'])/100, (data['RL_mean']+data['RL_SE'])/100, color = cols[i+2*j], alpha = 0.2)
        ax[0].set_ylabel("Total root length (m)")
        ax[0].set_ylim([0,2000])
        ax[0].legend()

        if j == 0: 
            ax[1].plot(data["DAS"],data["RL_mean"]/100,color = cols[i+2*j])
            ax[1].fill_between(data['DAS'], data['RL_mean']/100*0.6, data['RL_mean']/100*1.4, color = cols[i+2*j], alpha = 0.2, label = 'rRMSE = 40%')
            #ax[1].fill_between(np.nan, np.nan, np.nan, color = cols[i+2*j], alpha = 0.2, label = 'rRMSE = 40%')

        ax[1].legend()
plt.ylabel("Total root length (m)")
plt.xlabel("DAS (d) ")
plt.show()
