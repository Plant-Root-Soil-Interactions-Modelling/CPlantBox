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


times = np.linspace(0,153,154)
rtd = 0.25 #root tissue density (g/cm^3), Birouste et al. 2013
rootC = 0.41 #41% of root biomass can be considered as C , Hirte et al. 2018
names1 = ["RS_optimized_column_L_WT_new", "RS_optimized_column_S_WT_new", "RS_optimized_field_L_WT_new", "RS_optimized_field_S_WT_new"]
names = ["Field column experiment, L_WT", "Field column experiment, S_WT", "Unconfined growth experiment, L_WT", "Unconfined growth experiment, S_WT"]

data1 = np.load("root_volume.npz", allow_pickle=True)
RV = data1['arr_0']
data = np.load("plant_exud_rates_new.npz", allow_pickle=True)
exud = data['arr_0']

x = [0,0,1,1]
y = [0,1,0,1]

fig, ax = plt.subplots(2,2, figsize = (10, 10))

for i in range(0, len(names)): 

    ax[x[i],y[i]].bar(times, RV[:,i]*rtd*rootC, color = 'r', label = 'Root biomass C') #root biomass in g
    ax[x[i],y[i]].bar(times, np.cumsum(exud[:,i]), bottom = RV[:,i]*rtd*rootC, color = 'b', label = 'Exudate C')
    ax[x[i],y[i]].set_ylim([0,160])
    ax[x[i],y[i]].set_ylim([0,30])
    q = get_axis_limits(ax[x[i],y[i]])
    ax[x[i],y[i]].annotate(names[i], xy=(0, q[1]*0.9))

ax[0,0].annotate('Root tissue density = 0.25 $g/cm^3$ \nRoot biomass C content = 41% \n ', xy=(0, q[1]*1.1))
legend = ax[0,0].legend(loc = 'lower left')


for i in range(2,4):
    ax[x[i],y[i]].set_xlabel('Time (d)')
    
for i in [0,2]:
    ax[x[i],y[i]].set_ylabel('Root biomass C (g) \n Cumulative exudate C (g)')
    
    
#plt.tight_layout()
plt.show()

