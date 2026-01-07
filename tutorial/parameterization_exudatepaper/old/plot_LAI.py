import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import *
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import os

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title


#colors
a = [0, 0, 1]
b = [106/255, 166/255, 1]
c = [1, 0, 0]
d = [1, 128/255, 0]
cols = [a,b,c,d]
label_root = ["L_WT", "L_rth3", "S_WT", "S_rth3"]
soils = ['L', 'S']
gt = ['WT', 'rth3']

"""load LAI """
sim_time = 154
year = 2019
genotype_ = ['WT', 'RTH3']
soil_type_ = ['loam', 'sand']
df3 = pd.read_csv("../data_magda/LAI_"+str(year)+".csv")  # LAI
t_ = np.linspace(0, sim_time, (sim_time+1) * 24)  # relative time in hours


for i in range(0,len(soil_type_)):
    for j in range(0,len(genotype_)):

        soil_type = soil_type_[i]
        genotype = genotype_[j]

        if soil_type == 'loam':
            st = 'L'
        else:
            st = 'S'

            
        LAI =  df3[st+"_"+genotype].loc[0:sim_time].values
        f = interpolate.interp1d(np.linspace(0,sim_time, sim_time+1), LAI)

        plt.plot(t_, [f(t) for t in t_], color = cols[i+2*j], label = soils[i]+'_'+gt[j])

plt.xlabel('Time (DAS)')
plt.ylabel('LAI (-)')
plt.legend()
plt.show()
