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

def get_axis_limits(ax, scale=.9):
    return ax.get_xlim()[1]*scale, ax.get_ylim()[1]*scale


#colors & labels
a = [1,0,0]
b = [0, 0, 1]
cols = [a,b,a,b]
linst = ['-', '-', '--', '--']
names = ['column experiment: loam', 'column experiment: sand', 'unconfined growth experiment: loam', 'unconfined growth experiment: sand']

times = [0, 42, 63, 98, 153]
exu_rates = [np.array([0.0018,0.0018,0.0007,0.0004,0.00043]),
              np.array([0.003,0.003,0.00046,0.00031,0.00033])]#[kg/(m2 day)]
data = np.load("plant_exud_rates.npz", allow_pickle=True)
exud = data['arr_0']
surf = data['arr_1']
numtips = data['arr_2']

#import measured exudate data
df = pd.read_csv("data/column_experiment_mean.csv")
DAS = np.unique(df["DAS"].loc[:].values)
soils = df["substrate"].loc[:].values
soils_ = np.unique(soils)

fig, ax = plt.subplots(2,1, figsize = (18, 8))

#PLOT1
ax[0].plot(times, exu_rates[0]*10**6/10, color = cols[0], label = 'Tip exudation rate, loam')
ax[0].plot(times, exu_rates[0]/2*10**6/10, color = cols[1], label = 'Base exudation rate, loam')
ax[0].plot(times, exu_rates[1]*10**6/10, linestyle = '--', color = cols[0], label = 'Tip exudation rate, sand')
ax[0].plot(times, exu_rates[1]/2*10**6/10, linestyle = '--', color = cols[1], label = 'Base exudation rate, sand')
ax[0].set_ylabel('Exudation rate' + '\n $(\mu g / cm^2 / d)$')
ax[0].set_xlabel('Time (d)')
ax[0].legend(loc = 'lower left')

#print(np.max(exu_rates*10**6/10))
#print(np.min(exu_rates/2*10**6/10))

for i in range(0,len(names)): 

    #PLOT2
    p = ax[1].plot(np.linspace(0,153,154), exud[:,i]*10**3/surf[:,i], color = cols[i], linestyle = linst[i])
    if i <=1:
        data = df[(df['substrate']==soils_[i]) & (df['genotype']=="WT")]
        mean = data["Exudation surface"].values*12*24/10**6 #µmol cm-2 h-1--> g/cm^2/d
        SE = data['Exudation surface SE'].values*12*24/10**6 
        ax[1].errorbar(times[1:], mean*10**3, yerr = SE*10**3, fmt = 'o',color = p[0].get_color())
    ax[1].set_ylabel('Total root system exudation per root surface area' + '\n $(mg / cm^2 / d)$')
    #ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax[1].set_xlabel('Time (d)')

    if i == len(names)-1:
        ax[1].plot(np.nan, np.nan, color = cols[0], label = 'loam')
        ax[1].plot(np.nan, np.nan, color = cols[1], label = 'sand')
        ax[1].plot(np.nan, np.nan, color = 'k', linestyle = '-', label = 'column experiment')
        ax[1].plot(np.nan, np.nan, color = 'k', linestyle = '--', label = 'unconfined growth experiment')
        ax[1].errorbar(np.nan,np.nan,np.nan, color = 'k', fmt = 'o', label = 'Mean of measurements +/- SE')
        ax[0].annotate('(a)', xy=get_axis_limits(ax[0]))
        ax[1].annotate('(b)', xy=get_axis_limits(ax[1]))

    #ex_min = np.min(exud[:,i]*10**6/numtips[:,i]/1) #microg / cm root / day
    #ex_max = np.max(exud[:,i]*10**6/numtips[:,i]/1) #microg / cm root / day
    #print('Exudation is between '+str(ex_min)+' and '+str(ex_max)+'$(\mu g / cm root / day)$' )
        
#ax[1].set_ylim([0,300])
ax[1].legend()
plt.tight_layout()
plt.show()
