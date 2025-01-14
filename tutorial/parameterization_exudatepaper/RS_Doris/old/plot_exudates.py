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
from scipy import interpolate

font = {'size'   : 18}
plt.rc('font', **font)

#import exudate data
names = ['column experiment', 'unconfined growth experiment: loam', 'unconfined growth experiment: sand']
times = [0, 42, 63, 98, 154]
exu_rates = np.array([0.002, 0.002,0.00065,0.00045,0.00045]) #[kg/(m2 day)]
data = np.load("plant_exud_rates.npy")

fig, ax = plt.subplots(1,2, figsize = (18, 8))

#PLOT 1
ax[0].plot(times, exu_rates*10**3/10**4,  label = 'tip exudation rate')
ax[0].plot(times, exu_rates/2*10**3/10**4,  label = 'basal exudation rate')
ax[0].set_ylabel('Exudation rate' + '\n $(g / cm^2 / d)$')
ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax[0].set_xlabel('Time (d)')

for i in range(0,len(names)): 

    #PLOT 2
    ax[1].plot(np.linspace(0,154,155), data[i])
    ax[1].set_ylabel('Exudation rate' + '\n $(g / plant / d)$')
    ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax[1].set_xlabel('Time (d)')
    
plt.legend()
plt.tight_layout()
plt.show()
