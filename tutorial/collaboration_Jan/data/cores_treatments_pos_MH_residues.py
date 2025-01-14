import sys
sys.path.append("../../../..")
sys.path.append("../../../../src")

import plantbox as pb
import visualisation.vtk_plot as vp
import vtk 
import math
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import pandas as pd
import timeit


font = {'size'   : 16}
plt.rc('font', **font)

#import RLD data Maxime 
df = pd.read_csv("df_rld_Maxime.csv")
depth = df['depth'].values
treatment = df["treatment"].values
year = df["year"].values
BBCH = df["growth.stage"].values
RLD = df["rld"].values
res_v = df["root_res_volume"].values


#conversion BBCH to DAS
#conversion BBCH->days: BBCH14, BBCH19, BBCH59 and BBCH83 = days 42, 63, 98 and 154 
DAS = np.zeros((len(RLD)))
DAS[BBCH == 'BBCH14'] = 42
DAS[BBCH == 'BBCH19'] = 63
DAS[BBCH == 'BBCH59'] = 98
DAS[BBCH == 'BBCH83'] = 154

#import data Helena
df2 = pd.read_csv("RLD_2019_Helena.csv")
depth2 = df2['depth'].values
BBCH2 = df2["BBCH"].values
treatment2 = df["treatment"].values
RLD2 = df2["RLD"].values


#colors & labels 
a = [0, 0, 1]
b = [106/255, 166/255, 1]
c = [1, 0, 0]
d = [1, 128/255, 0]
cols = [a,c,b,d]
linst = ['-', '--','-.',':']
diam_ = [0.025, 0.03]

depth_ = np.unique(depth)
treatment_ = np.unique(treatment)
BBCH_ = np.unique(BBCH)
DAS_ = np.unique(DAS)
year_ = np.unique(year)
BBCH2_ = np.unique(BBCH2)
depth2_ = np.unique(depth2)

data = np.zeros((len(treatment_), len(DAS_),len(depth_)))
data2 = np.zeros((len(treatment_), len(DAS_),len(depth2_)))
data_SE = np.zeros((len(treatment_), len(DAS_),len(depth_)))
data2_SE = np.zeros((len(treatment_), len(DAS_),len(depth2_)))
for i in range(0, len(treatment_)):
    if i <2:
        diam = diam_[0]
    else:
        diam = diam_[1]
            
    for j in range(0, len(DAS_)):
        for k in range(0, len(depth_)): 
            data_ = df[(df['treatment']==treatment_[i])&(df['growth.stage']==BBCH_[j])&(df['depth']==depth_[k])&(df['year']==2019)&(df['position']=='IR')]

            if (data_["rld"].loc[:].values).size > 0:
                    data[i,j,k] = np.nanmean(data_["rld"].loc[:].values)+np.nanmean(data_["root_res_volume"].loc[:].values)/(math.pi*(diam/2)**2)
                    data_SE[i,j,k] = np.nanstd(data_["rld"].loc[:].values)/len(data_["rld"].loc[:].values)**0.5

        for k in range(0, len(depth2_)):
            data2_ = df2[(df2['treatment']==treatment_[i])&(df2['BBCH']==BBCH2_[j])&(df2['depth']==depth2_[k])]
            #print(data2_["RLD"].loc[:].values)
            if (data2_["RLD"].loc[:].values).size > 0:
                data2[i,j,k] = data2_["RLD"].loc[:].values
                data2_SE[i,j,k] = data2_["RLD_SE"].loc[:].values
                    

#PLOT
fig, ax = plt.subplots(2,2, figsize = (18, 10))     
for i in range(0, len(treatment_)):
    for j in range(0, len(DAS_)):
        
        xind = int(j/2)
        yind = int(j%2!=0)
        
        for k in range(0, len(depth_)):
            x = data[i,j,k] #Maxime
            ax[xind,yind].plot([x,x],[depth_[k]-2.5,depth_[k]+2.5], color = cols[i],linewidth = 2)

        for k in range(0, len(depth2_)): 
            x2 = data2[i,j,k] #Helena
            ax[xind,yind].plot([x2,x2],[depth2_[k]-20, depth2_[k]], color = cols[i],linewidth = 0.8)

            
        ax[xind,yind].yaxis.set_inverted(True)
        ax[xind,yind].set_xlim([0,2.5])
        if xind == 1: 
            ax[xind,yind].set_xlabel('RLD [cm/cm³]')
        if yind == 0: 
            ax[xind,yind].set_ylabel('Depth [cm]')
        ax[xind,yind].set_title(BBCH_[j], color = 'k')
        ax[xind,yind].set_xlim([-0.1,2.5])


for j in range(0, len(treatment_)):
    ax[0,0].plot(np.nan,np.nan, color = cols[j], label = treatment_[j])
ax[0,0].plot(np.nan,np.nan, color = 'k', linewidth = 2, label = 'Maxime')
ax[0,0].plot(np.nan,np.nan, color = 'k', linewidth = .8, label = 'Helena')
ax[0,0].legend(loc='lower right')


plt.axis('tight')
plt.show()



