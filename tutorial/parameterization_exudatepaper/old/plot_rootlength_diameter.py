import sys
sys.path.append("../../../../CPlantBox")
sys.path.append("../../../../CPlantBox/src")
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


df = pd.read_csv("../data_magda/mean_measures_Eva_Michael.csv")
print(df.columns.values.tolist()) 
DAS = df["DAS"].loc[:].values

#create root system 
rs = pb.MappedRootSystem()
path = "rootsystem_params/"
name = "Zea_mays_5_Leitner_Streuber_2014_mod2"
rs.readParameters(path + name + ".xml")

# Initialize
rs.initializeLB(5,4)

RL = np.zeros((154))
rad = np.zeros((154))
times = np.linspace(1,154,154)
for j in range(0,154):

    if j<=98: 
        rs.simulate(1, True);

    length = np.array(rs.getParameter("length"))
    radius = np.array(rs.getParameter("radius"))
    RL[j] = np.sum(length)
    rad[j] = np.sum(length*radius)/np.sum(length)
np.savez('RS_params', RL, rad) 


fig, ax = plt.subplots(2,1, figsize = (18, 10))
optimized2 = np.array([1223.83740793, 13125.69800546, 55855.62879929,55855.62879929 ])
optimized_diam2 = np.array([0.41008175, 0.3496782,  0.29133455, 0.29133455])


#ax[0].plot(df['DAS'], df['Length']/100, color = 'b', label = 'mean of measured values')
#ax[0].fill_between(df['DAS'], (df['Length']-df['Length std'])/100, (df['Length']+df['Length std'])/100, color = 'b', alpha = 0.2, label = 'standard deviation of measured values')
ax[0].errorbar(df['DAS'], df['Length']/100, yerr=df['Length std']/100, color ='r', fmt='.', label = 'Measurement')
ax[0].plot(times, RL/100, color = 'b', linestyle = '--', label = 'Simulation')
ax[0].set_ylabel('Total root length\n (m / plant)')
#ax[0].text(50, 600, '(a)', fontsize=16, va='top', ha='right')

#ax[1].plot(df['DAS'], df['RootDiameter'], color = 'b')
#ax[1].fill_between(df['DAS'], df['RootDiameter']-df['RootDiameter std'], df['RootDiameter']+df['RootDiameter std'], color = 'b', alpha = 0.2)
ax[1].errorbar(df['DAS'], df['RootDiameter'], yerr=df['RootDiameter std'], color = 'r', fmt='.')
ax[1].plot(times, rad*20, color = 'b', linestyle = '--')
ax[1].set_ylabel('Root diameter\n (mm)')
#ax[1].text(50, 0.45, '(b)', fontsize=16, va='top', ha='right')

        
for i in range(0,1):
    ax[i].set_xlabel('Time (d)')

ax[0].annotate('(a)', xy=get_axis_limits(ax[0]))
ax[1].annotate('(b)', xy=get_axis_limits(ax[1]))

ax[0].legend(loc = 'lower right')
plt.tight_layout()
plt.show()



