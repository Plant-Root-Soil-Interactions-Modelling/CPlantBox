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
names = ['column experiment: loam', 'column experiment: sand', 'unconfined growth experiment: loam', 'unconfined growth experiment: sand']

times = [0, 42, 63, 98, 153]
exu_rates = [np.array([0.0015, 0.0015,0.0009,0.00045,0.00044]),
              np.array([0.002, 0.002,0.0006,0.00034,0.00033])]#[kg/(m2 day)]

data = np.load("plant_exud_rates_new.npz", allow_pickle=True)
exud = data['arr_0']
surf = data['arr_1']
numtips = data['arr_2']

#field trap experiment (mean, min, max)
L_tip_f = [441, 146, 1155]
S_tip_f = [299, 143, 702]
L_seg_f = [213, 121, 292]
S_seg_f = [142, 76, 267]

L_tip_f_err = np.atleast_2d([L_tip_f[0]-L_tip_f[1], L_tip_f[2]-L_tip_f[0]]).T
S_tip_f_err = np.atleast_2d([S_tip_f[0]-S_tip_f[1], S_tip_f[2]-S_tip_f[0]]).T
L_seg_f_err = np.atleast_2d([L_seg_f[0]-L_seg_f[1], L_seg_f[2]-L_seg_f[0]]).T
S_seg_f_err = np.atleast_2d([S_seg_f[0]-S_seg_f[1], S_seg_f[2]-S_seg_f[0]]).T

#lab trap experiments (mean, min, max)
L_tip_l = [490, 128, 941]
L_seg_l = [194, 87, 342]

L_tip_l_err = np.atleast_2d([L_tip_l[0]-L_tip_l[1], L_tip_l[2]-L_tip_l[0]]).T
L_seg_l_err = np.atleast_2d([L_seg_l[0]-L_seg_l[1], L_seg_l[2]-L_seg_l[0]]).T


#import measured exudate data
df = pd.read_csv("data/column_experiment_mean.csv")
DAS = np.unique(df["DAS"].loc[:].values)
soils = df["substrate"].loc[:].values
soils_ = np.unique(soils)

fig, ax = plt.subplots(2,1, figsize = (18, 8))


for i in range(0,len(names)): 

    #Plot1
    p = ax[0].plot(np.linspace(0,153,154), exud[:,i]*10**3, color = cols[i], linestyle = linst[i])
    if i <=1:
        data = df[(df['substrate']==soils_[i]) & (df['genotype']=="WT")]
        mean = data["Total exudation"].values*12*24/10**3 #mmol/plant/h--> mg/plant/d
        SE = data['Total exudation SE'].values*12*24/10**3 
        ax[0].errorbar(times[1:], mean*10**3, yerr = SE*10**3, fmt = 'o',color = p[0].get_color())
    ax[0].set_ylabel('Total root system exudation' + '\n $(mg / plant / d)$')
    #ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax[0].set_xlabel('Time (d)')

    if i == len(names)-1:
        ax[0].plot(np.nan, np.nan, color = cols[0], label = 'Loam, WT')
        ax[0].plot(np.nan, np.nan, color = cols[1], label = 'Sand, WT')
        ax[0].plot(np.nan, np.nan, color = 'k', linestyle = '-', label = 'Field column experiment')
        ax[0].plot(np.nan, np.nan, color = 'k', linestyle = '--', label = 'Unconfined growth experiment')
        ax[0].errorbar(np.nan,np.nan,np.nan, color = 'k', fmt = 'o', label = 'Mean of measurements +/- SE')
    
    #Plot2
    q = ax[1].plot(np.linspace(0,153,154), np.cumsum(exud[:,i]), color = cols[i], linestyle = linst[i])
    print(np.cumsum(exud[-1,i]))
    ax[1].annotate(np.around(np.sum(exud[:,i]),2), (153, np.sum(exud[:,i])),fontsize = 12)
    ax[1].set_ylabel('Cumulative total root system' + '\n exudation $(g C / plant)$')
    ax[1].set_xlabel('Time (d)')
        
ax[0].set_ylim([0,400])
ax[0].annotate('(a)', xy=get_axis_limits(ax[0]))
ax[0].legend()

ax[1].set_ylim([0,25])
ax[1].annotate('(b)', xy=get_axis_limits(ax[1]))


plt.tight_layout()
plt.show()
