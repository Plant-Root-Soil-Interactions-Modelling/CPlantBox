import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy import interpolate
import pandas as pd
from os import walk
import re
import csv
from pathlib import Path

fontsize = 20
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


left  = 0.15  # the left side of the subplots of the figure
right = 0.85    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.18  # the amount of width reserved for blank space between subplots
hspace = 0.5   # the amount of height reserved for white space between subplots


#data
# kr0_loam = np.array([[-1e4, 0.], [-0.1, 0.], [0., 7e-4], [12,1e-5],[300, 1e-5]])
# kr1_loam = np.array([[-1e4, 0.], [-0.1, 0.], [0., 1e-4], [12,1e-5],[300, 1e-5]])

# kr0_sand = np.array([[-1e4, 0.], [-0.1, 0.], [0., 1e-3], [8., 1e-5], [300, 1e-5]])
# kr1_sand = np.array([[-1e4, 0.], [-0.1, 0.], [0., 7e-4], [10., 1e-5],[300, 1e-5]])

# kx0 = np.array([[0., 0.000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [300., 0.432]])
# kx1 = np.array([[0., 0.000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])

kr0_loam = np.array([[0., 7e-4], [10., 7e-4], [12,1e-5],[30, 1e-5]])
kr1_loam = np.array([[0., 1e-4], [10., 1e-4], [12,1e-5],[30, 1e-5]])

kr0_sand = np.array([[0., 1e-3],[7., 1e-3], [9., 1e-5], [30, 1e-5]])
kr1_sand = np.array([[0., 7e-5],[8., 7e-5],[10., 1e-5],[30, 1e-5]])

kx0 = np.array([[0., 0.000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [30., 0.432]])
kx1 = np.array([[0., 0.0000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [30., 0.00173]])

k_loam = dict(kr0 = kr0_loam, kr1 = kr1_loam, kx0 = kx0, kx1= kx1)
k_sand = dict(kr0 = kr0_sand, kr1 = kr1_sand, kx0 = kx0, kx1= kx1)

soiltype = ['Loam', 'Sand']
roottype = ['Axial roots', 'Lateral roots']
conds = ['Radial conductivity\n(cm $jPa^{-1}$ $d^{-1}$)','Axial conductance\n($cm^4$ $jPa^{-1}$ $d^{-1}$)']
var = np.array([["kr0","kr1"],["kx0","kx1"]])
fig, ax = plt.subplots(len(roottype),len(conds))
for i in range(0, len(roottype)): 
    for j in range(0, len(conds)): 
        print((k_loam["kr0"][:,0]))
        ax[j,i].plot(k_loam[var[i,j]][:,0], k_loam[var[i,j]][:,1],'b',label = 'Loam')
        ax[j,i].plot(k_sand[var[i,j]][:,0], k_sand[var[i,j]][:,1],'r',label = 'Sand')
ax[0,0].legend(loc = 'lower right')        
fig.supxlabel('Root segment age (d)', y = 0.01, fontsize=fontsize)
subfigs1 = fig.subfigures(1, 2)
fig.text(
    0.04,          # x-position
    0.5,           # y-position
    conds[0],
    va='center',
    rotation='vertical',
    fontsize=fontsize
)
fig.text(
    0.52,          # x-position
    0.5,           # y-position
    conds[1],
    va='center',
    rotation='vertical',
    fontsize=fontsize
)
# subfigs1.supylabel(conds[0], x = 0.03, fontsize=fontsize)
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show()

