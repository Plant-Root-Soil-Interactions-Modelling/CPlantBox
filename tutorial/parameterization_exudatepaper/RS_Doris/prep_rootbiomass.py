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
RV = np.zeros((154,len(names)))
for i in range(0, len(names)):
    
    rs.readParameters(path + names[i] + ".xml")

    # Initialize
    rs.initializeLB(5,4)

    RL = np.zeros((154))
    RL1 = np.zeros((154))
    rad = np.zeros((154))
    
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
        RV[j,i] = rad[j]**2*math.pi*RL[j]
        #RV[j] = np.sum(length*radius**2*math.pi)
        
np.savez('root_volume', RV)
