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



#data
# kr0_loam = np.array([[-1e4, 0.], [-0.1, 0.], [0., 7e-4], [12,1e-5],[300, 1e-5]])
# kr1_loam = np.array([[-1e4, 0.], [-0.1, 0.], [0., 1e-4], [12,1e-5],[300, 1e-5]])

# kr0_sand = np.array([[-1e4, 0.], [-0.1, 0.], [0., 1e-3], [8., 1e-5], [300, 1e-5]])
# kr1_sand = np.array([[-1e4, 0.], [-0.1, 0.], [0., 7e-4], [10., 1e-5],[300, 1e-5]])

# kx0 = np.array([[0., 0.000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [300., 0.432]])
# kx1 = np.array([[0., 0.0000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])

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

fig, ax = plt.subplots(
    len(roottype),
    len(conds),
    figsize=(14, 9)
)

for i in range(0, len(roottype)):
    for j in range(0, len(conds)):

        # ------------------------------------------------
        # colors
        # second column black
        # ------------------------------------------------

        if i == 0:
            color_loam = 'b'
            color_sand = 'r'
        else:
            color_loam = 'black'
            color_sand = 'black'

        # ------------------------------------------------
        # plotting
        # ------------------------------------------------

        ax[j,i].plot(
            k_loam[var[i,j]][:,0],
            k_loam[var[i,j]][:,1],
            linestyle='-',
            color=color_loam,
            linewidth=2,
            label='Loam'
        )

        ax[j,i].plot(
            k_sand[var[i,j]][:,0],
            k_sand[var[i,j]][:,1],
            linestyle='-',
            color=color_sand,
            linewidth=2,
            label='Sand'
        )

        # ------------------------------------------------
        # scientific notation
        # ------------------------------------------------

        ax[j,i].ticklabel_format(
            axis='x',
            style='sci',
            scilimits=(0,0)
        )

        ax[j,i].ticklabel_format(
            axis='y',
            style='sci',
            scilimits=(0,0)
        )

        # nicer scientific notation placement
        ax[j,i].xaxis.get_offset_text().set_fontsize(fontsize-2)
        ax[j,i].yaxis.get_offset_text().set_fontsize(fontsize-2)
        
        #text 
        ylim = ax[j,i].get_ylim()
        ax[j,i].text(15, ylim[1]*0.7, roottype[j],bbox=dict(facecolor='white',edgecolor='none'))
# --------------------------------------------------------
# legend
# --------------------------------------------------------

ax[0,0].legend(loc='lower right')

# --------------------------------------------------------
# spacing
# --------------------------------------------------------

plt.subplots_adjust(
    left=0.12,
    right=0.95,
    bottom=0.15,
    top=0.93,
    wspace=0.38,
    hspace=0.30
)

# --------------------------------------------------------
# global x label
# --------------------------------------------------------

fig.supxlabel(
    'Root segment age (d)',
    fontsize=fontsize,
    y=0.04
)

# --------------------------------------------------------
# y labels for BOTH columns
# automatically aligned with subplot positions
# --------------------------------------------------------

# left column label
x_left = ax[0,0].get_position().x0 - 0.07

fig.text(
    x_left,
    0.5,
    conds[0],
    va='center',
    ha='center',
    rotation='vertical',
    fontsize=fontsize
)

# right column label
x_right = ax[0,1].get_position().x0 - 0.07

fig.text(
    x_right,
    0.5,
    conds[1],
    va='center',
    ha='center',
    rotation='vertical',
    fontsize=fontsize
)

plt.show()