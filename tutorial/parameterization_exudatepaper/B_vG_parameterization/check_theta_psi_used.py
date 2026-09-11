import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy import interpolate
import pandas as pd
from os import walk
import re
import csv
from pathlib import Path
import os
from matplotlib.colors import ListedColormap
from scipy import interpolate

def get_axis_limits(ax, scalex=-0.01, scaley=0.9):
    return ax.get_xlim()[1]*scalex, ax.get_ylim()[1]*scaley    

fontsize = 24
mpl.rcParams.update({
    'font.size': fontsize,
    'axes.titlesize': fontsize,
    'axes.labelsize': fontsize,
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,
    'legend.fontsize': fontsize,
    'figure.titlesize': fontsize,
    'mathtext.default': 'regular'
})


left  = 0.1  # the left side of the subplots of the figure
right = 0.85    # the right side of the subplots of the figure
bottom = 0.12   # the bottom of the subplots of the figure
top = 0.85     # the top of the subplots of the figure
wspace = 0.18  # the amount of width reserved for blank space between subplots
hspace = 0.3   # the amount of height reserved for white space between subplots


def vG(vG, pHead): 
    
    theta_r = vG[0]
    theta_s = vG[1]
    alpha = vG[2]
    n = vG[3]
    m = 1-1/n
    wc = []
    for i in range(0, len(pHead)): 
        wc.append(theta_r + (theta_s - theta_r) / (1 + (alpha * abs(pHead[i]))** n)** m) 
    return np.asarray(wc)

path2exp = "data/"
cols = ['b', 'r', 'g']
soiltype = ['loam', 'sand']
depth_points = np.array([-10, -20, -40, -60])
exp_labels = [10,20,40,60]
cols = ['r', 'b', 'g', 'c', 'm']
vG_params = [[0.062, 0.337, 0.0182, 2.733, 1174],[0.1, 0.411, 0.05, 1.267, 441]]
soil_types = {
    "loam": {"j": [0], "vG":[0.1, 0.411, 0.05, 1.267, 441]},
    "sand": {"j": [1], "vG":[0.062, 0.337, 0.0182, 2.733, 1174]},
}

fig, ax = plt.subplots(1,2,figsize=(16, 8))
for soil, settings in soil_types.items():
    j = settings['j'][0]

    #experiment 
    df_theta = pd.read_csv(path2exp+"WC_"+soiltype[j]+".csv")
    time_theta = df_theta['time'].values
    theta = []
    for k in range(0, len(depth_points)): 
        theta.append(df_theta['cm'+str(exp_labels[k])].values)
    theta = np.asarray(theta).T

    df_psi = pd.read_csv(path2exp+"WP_"+soiltype[j]+".csv")
    time_psi = df_psi['time'].values
    psi = []
    for k in range(0, len(depth_points)): 
        psi.append(df_psi['cm'+str(exp_labels[k])].values)
    psi = np.asarray(psi).T

    for k in range(0, np.shape(theta)[1]): 
        ax[j].scatter(np.log10(abs(psi[:,k])),theta[:,k],color = cols[k],label= depth_points[k])
    
    #used retention curves 
    pHead = np.arange(-15000,0,1)
    wc = vG(np.asarray(settings['vG']), np.asarray(pHead))
    ax[j].plot(np.log10(abs(pHead)), wc, color = 'k') 
    
    
    #annotations
    # ax[i].set_xlim([0,50])
    # ax[i].set_ylim([0,0.4])
    # ax[j].set_ylabel('Soil water content ($cm^3$ $cm^{-3}$)', x = 0.47, y = 0.03, fontsize=fontsize)
    # ax[j].set_xlabel('abs Soil water potential  ($cm$)',  fontsize=fontsize)    
    ax[j].legend()

# #add legend
# handles1, labels1 = ax[0].get_legend_handles_labels()
# fig.legend(handles1,labels1,loc='upper left',bbox_to_anchor=(0.05, 1.02),ncol=2,frameon=False,fontsize=fontsize-2)
        

fig.supylabel('Soil water content ($cm^3$ $cm^{-3}$)', fontsize=fontsize)
fig.supxlabel('pf  ($-$)',  fontsize=fontsize)    
plt.show()

