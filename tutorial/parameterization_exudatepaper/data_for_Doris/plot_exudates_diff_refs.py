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
plt.rcParams['mathtext.default'] = 'regular'

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
data = np.load("plant_exud_rates_new.npz", allow_pickle=True)
exud = data['arr_0']
x = [0,0,1,1]
y = [0,1,0,1]
rho = [1.4,1.5,1.4,1.5]

fig, ax = plt.subplots(2,2, figsize = (18, 8))
for i in range(0,len(names)): 

    ax[0,0].plot(times, np.cumsum(exud[:,i]), color = cols[i], linestyle = linst[i], label = names[i]) #(g C / plant)
    ax[0,1].plot(times, np.cumsum(exud[:,i])/(0.2*0.45), color = cols[i], linestyle = linst[i]) #(g C / m^2)
    ax[1,0].plot(times, np.cumsum(exud[:,i])/(20*45*75), color = cols[i], linestyle = linst[i]) #(g C / cm^3)
    ax[1,1].plot(times, np.cumsum(exud[:,i])/(20*45*75)*rho[i], color = cols[i], linestyle = linst[i]) #(g C / g soil)


ax[0,0].set_ylabel('Cumulative total root system' + '\n exudation (g C / plant)')
ax[0,1].set_ylabel('Cumulative total root system' + '\n exudation (g C / $m^2$ soil surface)')
ax[1,0].set_ylabel('Cumulative total root system' + '\n exudation (g C / $cm^3$ soil)')
ax[1,1].set_ylabel('Cumulative total root system' + '\n exudation (g C / g soil)')

ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        
# ax[0].set_ylim([0,400])
# ax[0].annotate('(a)', xy=get_axis_limits(ax[0]))
ax[0,0].legend()

# ax[1].set_ylim([0,25])
# ax[1].annotate('(b)', xy=get_axis_limits(ax[1]))

for i in range(2,4):
    ax[x[i],y[i]].set_xlabel('Time (d)')

plt.tight_layout()
plt.show()
