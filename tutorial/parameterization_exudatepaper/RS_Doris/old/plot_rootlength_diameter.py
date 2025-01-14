import sys
sys.path.append("../../../../../CPlantBox")
sys.path.append("../../../../../CPlantBox/src")
import plantbox as pb
import visualisation.vtk_plot as vp

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import pandas as pd
import numpy as np

font = {'size'   : 18}
matplotlib.rc('font', **font)

def get_axis_limits(ax, scale=.9):
    return ax.get_xlim()[1]*scale, ax.get_ylim()[1]*scale

#colors & labels 
a = [0, 0, 1]
b = [106/255, 166/255, 1]
c = [1,0,0]
cols = [a,b,c]
label_root = ["unconfined growth experiment: loam", "unconfined growth experiment: sand", "column experiment"]


df = pd.read_csv("data/RLD_Doris.csv")
df1 = pd.read_csv("data/mean_exudation.csv")
print(df.columns.values.tolist()) 
DAS = df["DAS"].loc[:].values
soils = df["substrate"].loc[:].values
soils_ = np.unique(soils)
depth = 60

#create root system 
rs = pb.MappedRootSystem()
path = "rootsystem/"
names = ["RS_optimized_L_WT", "RS_optimized_S_WT", "RS_optimized_exudates"]

fig, ax = plt.subplots(2,1, figsize = (18, 10))
for i in range(0, len(names)):
    
    rs.readParameters(path + names[i] + ".xml")

    # Initialize
    rs.initializeLB(5,4)

    RL = np.zeros((154))
    test = np.zeros((154))
    rad = np.zeros((154))
    times = np.linspace(1,154,154)
    for j in range(0,154):

        if i==1 or j<=98: 
            rs.simulate(1, True);

        length = np.array(rs.getParameter("length"))
        radius = np.array(rs.getParameter("radius"))
        ana = pb.SegmentAnalyser(rs)
        RL[j] = np.sum(ana.distribution("length", 0., -depth, 1, True))
        test[j] = np.sum(length)/RL[j]
        rad[j] = np.sum(length*radius)/np.sum(length)
    #np.savez('RS_params', RL, rad) 

    ax[0].plot(times, RL/100, color = cols[i], linestyle = '--')
    ax[1].plot(times, rad*20, color = cols[i], linestyle = '--')

    if i<2: 
        data = df[(df['substrate']==soils_[i]) & (df['genotype']=="WT")]
        ax[0].errorbar(data['DAS'], data['RL_mean']/100, yerr=data['RL_SE']/100, color = cols[i], fmt='o')
        #ax[1].errorbar(data['DAS'], data['RootDiameter'], yerr=data['RootDiameter std'], color = cols[i], fmt='o')
    else: 
        ax[0].errorbar(df1['DAS'], df1['Length']/100, yerr=df1['Length std']/100/2, color =cols[i], fmt='o')
        ax[1].errorbar(df1['DAS'], df1['RootDiameter'], yerr=df1['RootDiameter std']/2, color = cols[i], fmt='o')
        
ax[0].plot(np.nan, np.nan, color = cols[0], label = label_root[0])
ax[0].plot(np.nan, np.nan, color = cols[1], label = label_root[1])
ax[0].plot(np.nan, np.nan, color = cols[2], label = label_root[2])
ax[0].plot(np.nan, np.nan, 'k--', label = 'Simulation')
ax[0].errorbar(np.nan, np.nan, np.nan, color = 'k', fmt = 'o', label = 'Mean of measurements +/- SE')
ax[0].set_ylabel('Total root length\n (m / plant)')
ax[1].set_ylabel('Root diameter\n (mm)')
ax[0].annotate('(a)', xy=get_axis_limits(ax[0]))
ax[1].annotate('(b)', xy=get_axis_limits(ax[1]))
ax[0].legend(loc = 'upper left')

for i in range(0,1):
    ax[i].set_xlabel('Time (d)')
plt.tight_layout()
plt.show()

print(test) 
