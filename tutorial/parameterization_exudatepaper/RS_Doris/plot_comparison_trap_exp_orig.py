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

#exudation rates in (µg C/cm2*d)

#colors & labels
b = [1,0,0]
a = [0,0,1]
cols = [a,b]
hatches = [' ', '//']

#read in data
df = pd.read_csv("data/comp_trap_simulation.csv")

#plot
fig, (ax1, ax2) = plt.subplots(1, 2)

#subplot lab
data = df[(df['Day']==22)]
DAS = data["Day"]
Origin = data["Origin"]
Substrate = data["Soil"]
Location = data["Location"]
mean = data["mean exu rate"]
SE = data["SE exu rate"]

Origin_ = np.unique(Origin)[::-1] 
Substrate_ = np.unique(Substrate)
Location_ = np.unique(Location)[::-1]

x = np.arange(len(np.unique(Origin)))  # the label locations
width = 1  # the width of the bars
multiplier = 0

for h in range(0, len(Origin_)): 
    for i in range(0, len(Substrate_)):
        for j in range(0, len(Location_)):
            data1 = data[(Origin==Origin_[h]) & (Substrate==Substrate_[i]) & (Location==Location_[j])]
            print(data1)
            if data1.empty:
                print('dataframe is empty')
            else: 
                offset = width * multiplier
                ax1.bar(offset+h, data1['mean exu rate'], color = cols[i], hatch = hatches[j])
                if h == 0: 
                    ax1.errorbar(offset+h, data1['mean exu rate'], yerr = data1['SE exu rate'], fmt = '.', color = 'k')
                print(multiplier, 'multiplier') 
                multiplier += 1


ax1.set_xticks([0.5,4.5])
ax1.set_xticklabels(Origin_)
ax1.bar(np.nan, np.nan, color = cols[0], label = 'Loam')
ax1.bar(np.nan, np.nan, color = cols[1], label = 'Sand')
ax1.bar(np.nan, np.nan, color = 'w', edgecolor = 'k',  hatch = hatches[0], label = 'Tip exudation rate')
ax1.bar(np.nan, np.nan, color = 'w', edgecolor = 'k', hatch = hatches[1], label = 'Base exudation rate')
ax1.legend()
ax1.set_ylabel(u'Exudation rate (\u03bcg C $cm^{-2}d^{-1}$)')
ax1.set_xlabel('DAS 22')


#subplot field
data = df[(df['Day']==63)]
DAS = data["Day"]
Origin = data["Origin"]
Substrate = data["Soil"]
Location = data["Location"]
mean = data["mean exu rate"]
SE = data["SE exu rate"]

Origin_ = np.unique(Origin)[::-1] 
Substrate_ = np.unique(Substrate)
Location_ = np.unique(Location)[::-1]

x = np.arange(len(np.unique(Origin)))  # the label locations
width = 1  # the width of the bars
multiplier = 0

for h in range(0, len(Origin_)): 
    for i in range(0, len(Substrate_)):
        for j in range(0, len(Location_)):
            data1 = data[(Origin==Origin_[h]) & (Substrate==Substrate_[i]) & (Location==Location_[j])]
            print(data1)
            if data1.empty:
                print('dataframe is empty')
            else: 
                offset = width * multiplier
                ax2.bar(offset+h, data1['mean exu rate'], color = cols[i], hatch = hatches[j])
                if h == 0: 
                    ax2.errorbar(offset+h, data1['mean exu rate'], yerr = data1['SE exu rate'], fmt = '.', color = 'k')
                print(multiplier, 'multiplier') 
                multiplier += 1


ax2.set_xticks([1.5,6.5])
ax2.set_xticklabels(Origin_)
ax2.set_xlabel('DAS 63')


plt.show()
