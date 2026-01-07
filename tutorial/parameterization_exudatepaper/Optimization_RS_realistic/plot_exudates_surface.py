"""small example"""
import sys
sys.path.append("../../../../CPlantBox")
sys.path.append("../../../../CPlantBox/src")
import plantbox as pb
from functional.xylem_flux import XylemFluxPython
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


#import exudate data
exudate = 'Exudation surface'
df = pd.read_csv("data/column_experiment_mean.csv")
DAS = np.unique(df["DAS"].loc[:].values)
soils = df["substrate"].loc[:].values
soils_ = np.unique(soils)
ssoils = ['Loam', 'Sand']

cols = ['b', 'r']

#create root system 
rs = pb.MappedRootSystem()
path = "rootsystem/"
names = ["RS_optimized_column_L_WT", "RS_optimized_column_S_WT"]


for p in range(0, len(names)):

    name = names[p]
    times = [42, 63, 98, 154]
    cols = ['b', 'r', 'g']

    #measured data
    data = df[(df['substrate']==soils_[p]) & (df['genotype']=="WT")]
    mean = data[exudate].values*12*24/10**3 #mmol/cm²/h --> g/cm²/day
    SE = data[exudate+' SE'].values*12*24/10**3
    
    #PLOT
    plt.plot(DAS, mean, color = cols[p], label = ssoils[p])
    plt.fill_between(DAS, mean-SE, mean+SE, color = cols[p], alpha = 0.2)

plt.ylabel('Exudation rate per' + '\n root surface area  (g / $cm^2$ / d)')
plt.xlabel('Time (d)') 
plt.legend()
#plt.text(50, 33, '(c)', fontsize=16,  va='top', ha='right')
plt.tight_layout()
plt.show()

