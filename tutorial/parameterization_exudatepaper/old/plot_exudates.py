"""small example"""
import sys
sys.path.append("../../../../CPlantBox")
sys.path.append("../../../../CPlantBox/src")
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
exudates = ['Phenolics', 'Sugars', 'AminoAcids']
times = [0, 42, 63, 98, 154]
exu_rates = np.array([[0.00012, 0.00012, 0.00003, 0.000019, 0.000015],[0.0006, 0.0006, 0.00015, 0.00017, 0.00016],[0.0001, 0.0001, 0.000025, 0.0000168, 0.0000135]]) #[kg/(m2 day)]
data = np.load("plant_exud_rates.npy")

fig, ax = plt.subplots(1,2, figsize = (18, 8))
for i in range(0,len(exudates)): 

    #PLOT 1
    ax[0].plot(times, exu_rates[i]*10**3/10**4,  label = exudates[i])
    ax[0].ylabel('Exudation rate' + '\n $(g / cm^2 / d)$')
    ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax[0].xlabel('Time (d)')

    #PLOT 2
    ax[1].plot(np.linspace(0,154,155), data[i])
    ax[1].ylabel('Exudation rate' + '\n $(g / plant / d)$')
    ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax[1].xlabel('Time (d)')
    
plt.legend()
plt.tight_layout()
plt.show()
