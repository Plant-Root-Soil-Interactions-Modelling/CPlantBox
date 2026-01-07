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


#colors & labels
b = [1,0,0]
a = [0, 0, 1]
cols = [a,b,a,b]
cols_ = ['k', 'g']
linst = ['-', '-', '--', '--']
names = ["Field column experiment, L_WT", "Field column experiment, S_WT", "Unconfined growth experiment, L_WT", "Unconfined growth experiment, S_WT"]
times = np.linspace(0,153,154)
data = np.load("plant_exud_rates_layers.npz", allow_pickle=True)
exud = data['arr_0']
x = [0,0,1,1]
y = [0,1,0,1]
depths = ['0-20cm depth', '20-40cm depth', '40-60cm depth']

fig, ax = plt.subplots(2,2, figsize = (18, 8))
for i in range(0,len(names)): 

    ax[0,0].plot(times, np.cumsum(exud[:,i,0]), color = cols[i], linestyle = linst[i], label = names[i]) #(g C / plant)
    ax[0,1].plot(times, np.cumsum(exud[:,i,1]), color = cols[i], linestyle = linst[i]) #(g C / m^2)
    ax[1,0].plot(times, np.cumsum(exud[:,i,2]), color = cols[i], linestyle = linst[i]) #(g C / cm^3)



ax[0,0].set_ylabel('Cumulative total root system' + '\n exudation $(g C / plant)$')
ax[0,1].set_ylabel('Cumulative total root system' + '\n exudation $(g C / plant)$')
ax[1,0].set_ylabel('Cumulative total root system' + '\n exudation $(g C / plant)$')


        
# ax[0].set_ylim([0,400])
for i in range(0,3): 
    q = get_axis_limits(ax[x[i],y[i]])
    ax[x[i],y[i]].annotate(depths[i], xy=(0, q[1]*0.9))

legend = ax[0,0].legend(loc = 'lower left')

# ax[1].set_ylim([0,25])
# ax[1].annotate('(b)', xy=get_axis_limits(ax[1]))

for i in range(0,3):
    ax[x[i],y[i]].set_xlabel('Time (d)')
    
ax[1,1].axis("off")
plt.tight_layout()
plt.show()
