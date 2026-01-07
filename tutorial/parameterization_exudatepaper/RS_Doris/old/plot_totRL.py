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
height = [20, 20, 20]
DAS = [42, 63, 98, 154]
#area = 20*45 #cm
area = 1052.5


#import RLD data
df = pd.read_csv("data/RLD.csv")
BBCH = df["BBCH"].loc[:].values
soils = df["substrate"].loc[:].values
gt = df["genotype"].loc[:].values
#depth = df["depth"].loc[:].values
RLD_mean = df["RLD_mean"].loc[:].values
RLD_SE = df["RLD_SE"].loc[:].values

soils_ = np.unique(soils)
gt_ = np.unique(gt)[::-1]
BBCH_ = np.unique(BBCH)


dummy = -1
for i in range(0, len(soils_)):
    for j in range(0, len(gt_)):
        dummy = dummy+1
        RL = np.zeros(len(BBCH_))
        for k in range(0, len(BBCH_)):

            RLD = df[(df['substrate']==soils_[i]) & (df['genotype']==gt_[j]) & (df['BBCH']==BBCH_[k])]
            x = RLD["RLD_mean"].loc[:].values
            RL_ = x*height[:len(x)]
            RL[k] = np.sum(RL_)*area
            
        plt.plot(DAS,RL/100,color = cols[i+2*j], linestyle = '--', label = soils_[i]+'_'+gt_[j])


df = pd.read_csv("../../data_magda/diam_length_exudates.csv")
print(df.columns.values.tolist()) 
DAS = df["DAS"].loc[:].values

soils = ['L', 'S']
gt = ['WT', 'rth3']


for i in range(0,len(soils)):
    for j in range(0, len(gt)):

        y = df[(df['Substrate']== soils[i]) & (df['Genotype']==gt[j])]

        plt.plot(y['DAS'], y['Length']/100, color = cols[i+2*j], label = soils[i]+'_'+gt[j])
        plt.fill_between(y['DAS'], (y['Length']-y['Length std']/2)/100, (y['Length']+y['Length std']/2)/100, color = cols[i+2*j], alpha = 0.2)

plt.plot(np.nan, np.nan, color = 'k', linestyle = '-', label = 'RLD, Eva Michael')
plt.plot(np.nan, np.nan, color = 'k', linestyle = '--', label = 'RLD, Doris') 

plt.legend()
plt.ylabel("Total root length (m)")
plt.xlabel("DAS (d) ")
plt.show()
