import sys
sys.path.append("../../../../CPlantBox")
sys.path.append("../../../../CPlantBox/src")
import plantbox as pb
import visualisation.vtk_plot as vp
import math
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
b = [1, 0, 0]
#b = [106/255, 166/255, 1]
cols = [a,b,a,b]

df = pd.read_csv("data/field_experiment_mean.csv")
df1 = pd.read_csv("data/column_experiment_mean.csv")
print(df1.columns.values.tolist())
DAS = df["DAS"].loc[:].values
soils = df["substrate"].loc[:].values
soils_ = np.unique(soils)
depth = 60
linst = ['-', '-', '--', '--']
times_ = [42, 63, 98, 153]
rtd = 0.25 #g/cm^3
#create root system 
rs = pb.MappedRootSystem()
path = "rootsystem/"
names = ["RS_optimized_column_L_WT_new", "RS_optimized_column_S_WT_new", "RS_optimized_field_L_WT_new", "RS_optimized_field_S_WT_new"]

fig, ax = plt.subplots(2,1, figsize = (10, 10))
#plt.subplots_adjust(left=None, bottom=None, right=0.8)
for i in range(0, len(names)):
    
    rs.readParameters(path + names[i] + ".xml")

    # Initialize
    rs.initializeLB(5,4)

    RL = np.zeros((154))
    RL1 = np.zeros((154))
    rad = np.zeros((154))
    RV = np.zeros((154))
    times = np.linspace(1,154,154)
    for j in range(0,154):

        if i==3 or j<=98: 
            rs.simulate(1, True);

        length = np.array(rs.getParameter("length"))
        radius = np.array(rs.getParameter("radius"))

        if i>1: #for field experiment RL only until 60cm soil depth 
            ana = pb.SegmentAnalyser(rs)
            RL[j] = np.sum(ana.distribution("length", 0., -depth, 1, True))
            RL1[j] = np.sum(length)
        else:
            RL[j] = np.sum(length)
        rad[j] = np.sum(length*radius)/np.sum(length)
        RV[j] = rad[j]**2*math.pi*RL[j]
        #RV[j] = np.sum(length*radius**2*math.pi)

    ax[0].plot(times, RV, color = cols[i], linestyle = linst[i]) #cm^3
    ax[1].plot(times, RV*rtd, color = cols[i], linestyle = linst[i]) #root biomass in g

    if i>1: 
        data = df[(df['substrate']==soils_[i-2]) & (df['genotype']=="WT")]
        ax[0].scatter(data['DAS'], data['RL_mean']*(data['RD_mean']/20)**2*math.pi, color = cols[i], marker='x')

    else:
        data = df1[(df1['substrate']==soils_[i]) & (df1['genotype']=="WT")]
        ax[0].scatter(data['DAS'], data['RL_mean']*(data['RD_mean']/20)**2*math.pi, color = cols[i], marker='o')

        
ax[0].plot(np.nan, np.nan, color = cols[0], label = 'Loam, WT')
ax[0].plot(np.nan, np.nan, color = cols[1], label = 'Sand, WT')
ax[0].plot(np.nan, np.nan, 'k-', label = 'Simulation, field column experiment')
ax[0].plot(np.nan, np.nan, 'k--', label = 'Simulation, unconfined growth experiment')
ax[0].plot(np.nan, np.nan, 'ko', label = 'Measurement, field column experiment')
ax[0].plot(np.nan, np.nan, 'kx', label = 'Measurement, unconfined growth experiment')
#ax[0].errorbar(np.nan, np.nan, np.nan, color = 'k', fmt = 'o', label = 'Mean of measurements +/- SE')
ax[0].set_ylabel('Root volume\n ($cm^3$ / plant)')
ax[1].set_ylabel('Root biomass (g)')

#ax[0].set_ylim([0,70])
ax[1].set_ylim([0,20])
ax0 = get_axis_limits(ax[0])
ax1 = get_axis_limits(ax[1])


ax[0].annotate('(a)', xy=(ax0[0], ax0[1]*0.9))
ax[1].annotate('(b)', xy=(ax1[0], ax1[1]*0.9))

legend = ax[0].legend(loc = 'upper left', bbox_to_anchor=(0, 1.2),facecolor = 'white', framealpha = 1)


for i in range(1,2):
    ax[i].set_xlabel('Time (d)')
#plt.tight_layout()
plt.show()

