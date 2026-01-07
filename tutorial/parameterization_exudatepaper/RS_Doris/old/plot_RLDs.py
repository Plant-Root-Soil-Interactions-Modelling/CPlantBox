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
label_root = ["L_WT", "L_rth3", "S_WT", "S_rth3"]
soils = ['L', 'S']
gt = ['WT', 'rth3']
linst = ['-', '--','-.',':']
depth = [0, -20, -20, -40, -40, -75]
DAS = [42, 63, 98, 154]


#import RLD data
df = pd.read_csv("data/RLD.csv")
BBCH = df["BBCH"].loc[:].values
soils = df["substrate"].loc[:].values
gt = df["genotype"].loc[:].values
#depth = df["depth"].loc[:].values
RLD_mean = df["RLD_mean"].loc[:].values
RLD_SE = df["RLD_SE"].loc[:].values

soils_ = np.unique(soils)
gt_ = np.unique(gt)
BBCH_ = np.unique(BBCH)

fig, axs = plt.subplots(2,2)

dummy = -1
for i in range(0, len(soils_)):
    for j in range(0, len(gt_)):
        dummy = dummy+1
        for k in range(0, len(BBCH_)):

            RLD = df[(df['substrate']==soils_[i]) & (df['genotype']==gt_[j]) & (df['BBCH']==BBCH_[k])]
            x_ = RLD["RLD_mean"].loc[:].values
            x1_ = RLD["RLD_mean"].loc[:].values-RLD["RLD_SE"].loc[:].values
            x2_ = RLD["RLD_mean"].loc[:].values+RLD["RLD_SE"].loc[:].values
            x = np.repeat(x_, 2)
            x1 = np.repeat(x1_, 2)
            x2 = np.repeat(x2_, 2)
            y = depth[:len(x)]

            axs[i,j].plot(x,y,color = cols[dummy], linestyle = linst[k])
            axs[i,j].fill_betweenx(y, x1, x2, facecolor=cols[dummy], alpha=.2)
            if i == 1 & j == 1:
                axs[i,j].plot(np.nan, np.nan, color = 'k', linestyle = linst[k], label = DAS[k])
                axs[i,j].legend()
            axs[i,j].set_xlim([0,4])
            axs[i,j].set_ylim([-75,0])
            axs[i,j].set_title(label_root[dummy], color = cols[dummy])

        if i == 0:
            axs[i,j].set_xticks([])
        else: 
            axs[i,j].set_xlabel("RLD, soil cores (cm / cm^3)")
        if j == 0: 
            axs[i,j].set_ylabel("soil depth (cm)")
        else:
            axs[i,j].set_yticks([])
plt.show()
