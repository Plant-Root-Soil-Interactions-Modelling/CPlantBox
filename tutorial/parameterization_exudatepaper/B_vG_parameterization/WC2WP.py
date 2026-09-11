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


def vG(vG, theta): 
    
    theta_r = vG[0]
    theta_s = vG[1]
    alpha = vG[2]
    n = vG[3]
    m = 1-1/n
    pHead = []
    for i in range(0, len(theta)): 
        if theta[i] >= theta_s:
            pHead_ = 0.0
        elif theta[i] <= theta_r:
            pHead_ = np.nan
        else:
            pHead_ = -1 / alpha * (((theta_s - theta_r) / (theta[i] - theta_r)) ** (1 / m) - 1) ** (1 / n)
        pHead.append(pHead_)
    return np.asarray(pHead)

path2exp = "data/"
soiltype = ['loam', 'sand']
depth_points = np.array([0,-10, -20, -40, -60,-75])
exp_labels = [0,10,20,40,60,75]
headers = ['time', 'cm0', 'cm10', 'cm20', 'cm40', 'cm60', 'cm75']
# soil_types = {
    # "loam": {"j": [0], "vG":[0.1, 0.411, 0.05, 1.267, 441]},
    # "sand": {"j": [1], "vG":[0.062, 0.337, 0.0182, 2.733, 1174]},
# }
# soil_types = {
    # "loam": {"j": [0], "vG":[0.1, 0.411, 0.05, 1.3315, 441]},
    # "sand": {"j": [1], "vG":[0.038, 0.337, 0.1357, 1.442, 1174]},
# }
soil_types = {
    "loam": {"j": [0], "vG":[0.1, 0.411, 0.05, 1.389, 441]},
    "sand": {"j": [1], "vG":[0.038, 0.337, 0.1, 1.459, 1174]},
}

fig, ax = plt.subplots(1,2,figsize=(16, 8))
for soil, settings in soil_types.items():
    j = settings['j'][0]

    #experiment 
    df_theta = pd.read_csv(path2exp+"WC_"+soiltype[j]+".csv")
    time = df_theta['time'].values
    pHead = []
    for k in range(0, len(depth_points)): 
        theta = (df_theta['cm'+str(exp_labels[k])].values)
        pHead.append(vG(np.asarray(settings['vG']), np.asarray(theta)))
    pHead = np.asarray(pHead).T
    data = np.column_stack((time, pHead))

    # Write to CSV
    np.savetxt(
        "data/WPfromWC_"+soil+".csv",
        data,
        delimiter=",",
        header=",".join(headers),
        comments=""
    )

