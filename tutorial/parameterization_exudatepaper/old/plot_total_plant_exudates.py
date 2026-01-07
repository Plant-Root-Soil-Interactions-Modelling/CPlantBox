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

def get_axis_limits(ax, scale=.9):
    return ax.get_xlim()[1]*scale, ax.get_ylim()[1]*scale

#import exudate data
exudates = ['Phenolics', 'Sugars', 'AminoAcids']
times = [0, 42, 63, 98, 153]
exu_rates = np.array([[0.00012, 0.00012, 0.00003, 0.000019, 0.000015],[0.0006, 0.0006, 0.00015, 0.00017, 0.00016],[0.0001, 0.0001, 0.000025, 0.0000168, 0.0000135]]) #[kg/(m2 day)]
#data = np.load("plant_exud_rates.npy")
data = np.load("plant_exud_rates.npz", allow_pickle=True)
exud = data['arr_0']
surf = data['arr_1']
numtips = data['arr_2']

#import measured exudate data
df = pd.read_csv("../data_magda/mean_measures_Eva_Michael.csv")
DAS = df["DAS"].loc[:].values

fig, ax = plt.subplots(2,1, figsize = (18, 8))
for i in range(0,len(exudates)): 

    mean = df[exudates[i]].values*12*24/10**3 #mmol/plant/h--> g/plant/d
    std = df[exudates[i]+' std'].values*12*24/10**3

    #plt.plot(times[1:], mean, 'ro')
    #plt.plot(np.linspace(0,153,154), exud[:,i], label = exudates[i])
    #if i == len(exudates)-1:
    #    plt.plot(times[1:], mean, 'ro', label = 'Measurement')
    #plt.ylabel('Exudation rate' + '\n (g / plant / d)')
    #plt.ylim([0,1e-1])
    #plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    #plt.xlabel('Time (d)')

    #PLOT
    ax[0].plot(times, exu_rates[i]*10**6/10, label = exudates[i])
    ax[0].set_ylabel('Maximum exudation rate' + '\n $(\mu g / cm^2 / d)$')
    #ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax[0].set_xlabel('Time (d)')
    
    #PLOT
    p = ax[1].plot(np.linspace(0,153,154), exud[:,i]*10**3, label = exudates[i])
    ax[1].errorbar(times[1:], mean*10**3, yerr = std/2*10**3, fmt = 'o',color = p[0].get_color())
    ax[1].set_ylabel('Total root system exudation' + '\n $(mg / plant / d)$')
    #ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax[1].set_xlabel('Time (d)')

    if i == len(exudates)-1:
        ax[1].plot(np.nan, np.nan, 'ko', label = 'Mean of measurements +/- SE')
        ax[1].plot(np.nan, np.nan, 'k-', label = 'Simulation')
        ax[0].annotate('(a)', xy=get_axis_limits(ax[0]))
        ax[1].annotate('(b)', xy=get_axis_limits(ax[1]))

    ex_min = np.min(exud[:,i]*10**6/numtips[:,i]/1) #microg / cm root / day
    ex_max = np.max(exud[:,i]*10**6/numtips[:,i]/1) #microg / cm root / day
    print('Exudation of '+exudates[i]+' is between '+str(ex_min)+' and '+str(ex_max)+'$(\mu g / cm root / day)$' )
        
    #PLOT
    #ax[1].plot(np.linspace(0,153,154), surf[:,i], label = exudates[i])
    #ax[1].set_ylabel('Exuding root surface' + '\n $(cm^2)$')
    #ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    #ax[1].set_xlabel('Time (d)')

    
plt.legend()
plt.tight_layout()
plt.show()
