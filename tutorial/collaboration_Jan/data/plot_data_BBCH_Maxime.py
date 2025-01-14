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

#import RLD data
df = pd.read_csv("df_rld_Maxime.csv")
depth = df['depth'].values
treatment = df["treatment"].values
year = df["year"].values
BBCH = df["growth.stage"].values
RLD = df["rld"].values

#conversion BBCH to DAS
#conversion BBCH->days: BBCH14, BBCH19, BBCH59 and BBCH83 = days 42, 63, 98 and 154 
DAS = np.zeros((len(RLD)))
DAS[BBCH == 'BBCH14'] = 42
DAS[BBCH == 'BBCH19'] = 63
DAS[BBCH == 'BBCH59'] = 98
DAS[BBCH == 'BBCH83'] = 154


#colors & labels 
a = [0, 0, 1]
b = [106/255, 166/255, 1]
c = [1, 0, 0]
d = [1, 128/255, 0]
cols = [a,c,b,d]
linst = ['-', '--','-.',':']


depth_ = np.unique(depth)
treatment_ = np.unique(treatment)
BBCH_ = np.unique(BBCH)
DAS_ = np.unique(DAS)
year_ = np.unique(year)

data = np.zeros((len(treatment_), len(DAS_),len(depth_)))
data_SE = np.zeros((len(treatment_), len(DAS_),len(depth_)))
for i in range(0, len(treatment_)):
    for j in range(0, len(DAS_)):
        for k in range(0, len(depth_)): 
            #print(treatment_[i], BBCH_[j], depth_[k])
            
            data_ = df[(df['treatment']==treatment_[i])&(df['growth.stage']==BBCH_[j])&(df['depth']==depth_[k])&(df['year']==2019)]
            if treatment_[i]=='L_WT':
                if BBCH_[j] == 'BBCH83':
                    print(depth_[k])
                    #print((data_["rld"].loc[:].values))
                    #print(np.nanmean(data_["rld"].loc[:].values))
            if (data_["rld"].loc[:].values).size > 0:
                    data[i,j,k] = np.nanmean(data_["rld"].loc[:].values)
                    data_SE[i,j,k] = np.nanstd(data_["rld"].loc[:].values)/(len(data_["rld"].loc[:].values)**.5)



#PLOT
fig, ax = plt.subplots(2,2, figsize = (18, 10))     
for i in range(0, len(treatment_)):
    for j in range(0, len(DAS_)):
        x = data[i,j,:]
        error = data_SE[i,j,:]
        xind = int(j/2)
        yind = int(j%2!=0)
        ax[xind,yind].plot(x,depth_, color = cols[i])
        ax[xind,yind].fill_betweenx(depth_, x-error, x+error, color = cols[i],alpha = 0.3)
        ax[xind,yind].yaxis.set_inverted(True)
        ax[xind,yind].set_xlim([0,2.5])
        if xind == 1: 
            ax[xind,yind].set_xlabel('RLD [cm/cm³]')
        if yind == 0: 
            ax[xind,yind].set_ylabel('Depth [cm]')
        ax[xind,yind].set_title(BBCH_[j], color = 'k')
        ax[xind,yind].set_xlim([-0.1,2.5])

for j in range(0, len(treatment_)):
    plt.plot(np.nan,np.nan, color = cols[j],  label = treatment_[j])


plt.legend(loc='upper left')
plt.axis('tight')
plt.show()



