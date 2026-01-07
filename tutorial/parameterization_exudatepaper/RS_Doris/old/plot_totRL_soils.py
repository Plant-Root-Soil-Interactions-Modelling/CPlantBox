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
cols = [a,d]
label_root = ["loam", "sand"]
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
DAS_ = np.unique(DAS)


for i in range(0, len(soils_)):

    x = DAS_
    y = np.zeros(len(DAS_))
    SE = np.zeros(len(DAS_))
                                         
    for j in range(0, len(DAS_)): 
    
        data = df[(df['substrate']==soils_[i]) & (df['DAS']==DAS_[j])]
        y[j] = np.mean(data["RL_mean"].loc[:].values)
        SE[j] = np.mean(data["RL_SE"].loc[:].values)

    plt.plot(x,y/100,color = cols[i], label = label_root[i])
    plt.fill_between(x, (y-SE)/100, (y+SE)/100, color = cols[i], alpha = 0.2)


plt.legend()
plt.ylabel("Total root length (m)")
plt.xlabel("DAS (d) ")
plt.show()
