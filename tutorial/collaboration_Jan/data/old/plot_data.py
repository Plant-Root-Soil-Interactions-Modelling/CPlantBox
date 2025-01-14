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
BBCH = df["growth.stage"].values
RLD = df["rld"].values

#conversion BBCH to DAS
#conversion BBCH->days: BBCH14, BBCH19, BBCH59 and BBCH83 = days 42, 63, 98 and 154 
DAS = np.zeros((len(RLD)))
DAS[BBCH == 'BBCH14'] = 42
DAS[BBCH == 'BBCH19'] = 63
DAS[BBCH == 'BBCH59'] = 98
DAS[BBCH == 'BBCH83'] = 154


#plot

#colors & labels 
a = [0, 0, 1]
b = [106/255, 166/255, 1]
c = [1, 0, 0]
d = [1, 128/255, 0]
cols = [a,b,c,d]
labels = ["L_WT", "S_WT", "L_rth3","S_rth3"]
linst = ['-', '--','-.',':']


depth_ = np.unique(depth)
treatment_ = np.unique(treatment)
BBCH_ = np.unique(BBCH)
DAS_ = np.unique(DAS)

data = np.zeros((len(treatment_), len(DAS_),len(depth_)))
for i in range(0, len(treatment_)):
    for j in range(0, len(DAS_)):
        for k in range(0, len(depth_)): 

                data_ = df[(df['treatment']==treatment_[i])&(df['growth.stage']==BBCH_[j])&(df['depth']==depth_[k])]

                if (data_["rld"].loc[:].values).size > 0:
                    data[i,j,k] = np.nanmean(data_["rld"].loc[:].values)

                
for i in range(0, len(treatment_)):
    for j in range(0, len(DAS_)):
        x = data[i,j,:]
        plt.plot(x,depth_, color = cols[i], linestyle = linst[j])

for i in range(0, len(treatment_)):
    plt.plot(np.nan,np.nan, color = cols[i], label = labels[i])
for j in range(0, len(DAS_)):
    plt.plot(np.nan,np.nan, color = 'k', linestyle = linst[j], label = BBCH_[j])
        
plt.gca().invert_yaxis()
plt.xlabel('RLD [cm/cm³]')
plt.ylabel('Depth [cm]')
plt.legend()
plt.axis('tight')
plt.show()
        
