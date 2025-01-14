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
BBCH = df["growth.stage"].values
RLD = df["rld"].values

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
linst = ['-', '--']


depth_ = np.unique(depth)
treatment_ = np.unique(treatment)
BBCH_ = np.unique(BBCH)
DAS_ = np.unique(DAS)
BBCH2_ = np.unique(BBCH2)
depth2_ = np.unique(depth2)

data = np.zeros((len(treatment_), len(DAS_),len(depth_)))
data2 = np.zeros((len(treatment_), len(DAS_),len(depth2_)))
data2_SE = np.zeros((len(treatment_), len(DAS_),len(depth2_)))
for i in range(0, len(treatment_)):
    for j in range(0, len(DAS_)):
        for k in range(0, len(depth_)): 
            data_ = df[(df['treatment']==treatment_[i])&(df['growth.stage']==BBCH_[j])&(df['depth']==depth_[k])]

            if (data_["rld"].loc[:].values).size > 0:
                    data[i,j,k] = np.nanmean(data_["rld"].loc[:].values)

        for k in range(0, len(depth2_)):
            data2_ = df2[(df2['treatment']==treatment_[i])&(df2['BBCH']==BBCH2_[j])&(df2['depth']==depth2_[k])]
            if (data2_["RLD"].loc[:].values).size > 0:
                data2[i,j,k] = data2_["RLD"].loc[:].values
                data2_SE[i,j,k] = data2_["RLD_SE"].loc[:].values
                    

#PLOT
fig, ax = plt.subplots(2,2, figsize = (18, 10))     
for i in range(0, len(treatment_)):
    for j in range(0, len(DAS_)):
        x = data[i,j,:] #Maxime
        x2 = data2[i,j,:] #Helena
        error2 = data2_SE[i,j,:] #Helena
        xind = int(j/2)
        yind = int(j%2!=0)
        #ax[xind,yind].plot(x,depth_, color = cols[i], linestyle = linst[j],linewidth = 2)
        ax[xind,yind].plot(x2,depth2_-10, color = cols[i])
        ax[xind,yind].fill_betweenx(depth2_-10, x2-error2, x2+error2, color = cols[i],alpha = 0.3)
        ax[xind,yind].yaxis.set_inverted(True)
        ax[xind,yind].set_xlim([0,2.5])
        if xind == 1: 
            ax[xind,yind].set_xlabel('RLD [cm/cm³]')
        if yind == 0: 
            ax[xind,yind].set_ylabel('Depth [cm]')
        ax[xind,yind].set_title(BBCH_[j], color = 'k')
        ax[xind,yind].set_xlim([-0.1,2.5])

for j in range(0, len(treatment_)):
    plt.plot(np.nan,np.nan, color = cols[j], label = treatment_[j])
#plt.plot(np.nan,np.nan, color = 'k', linewidth = 2, label = 'Maxime')
#plt.plot(np.nan,np.nan, color = 'k', linewidth = .8, label = 'Helena')


plt.legend(loc='lower right')
plt.axis('tight')
plt.show()



