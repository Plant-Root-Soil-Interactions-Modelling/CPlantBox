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

#create root system 
rs = pb.MappedRootSystem()
path = "rootsystem/"
names = ["RS_optimized_field_L_WT", "RS_optimized_field_S_WT"]
cols = ['b', 'r']
linst = ['-', ':', '-.', '--']
titles = ['Loam', 'Sand']

#data
times= [42, 63, 98, 154]

plt.figure(figsize=(9, 6))
ax1 = plt.subplot(2,2,1)
ax2 = plt.subplot(2,2,2)
ax3 = plt.subplot(2,1,2)
ax = [ax1, ax2, ax3]

depth = 75
layers = 10
area = 20 * 45
soilvolume = (depth / layers) * area
z = np.linspace(0, -depth, layers)  # z - axis


for i in range(0, len(names)):
    dummy = 0
    rs.readParameters(path + names[i] + ".xml")
    rs.initializeLB(5,4)
    
    numtips = []
    nodes = []
    x = np.linspace(0, times[-1], times[-1]+1)
    for j in range(0,times[-1]+1):
        if i==1 or (i==0 and j<=98): 
            rs.simulate(1, True);
            
        if j in times: 
            if i==1 or (i==0 and j<=98): 
                ana = pb.SegmentAnalyser(rs)
                rld_ = np.array(ana.distribution("length", 0., -depth, layers, True))
                rld = rld_/ soilvolume
                ax[i].plot(rld, z, color = cols[i], linestyle = linst[dummy], label = 'Day '+str(int(j)))
                dummy = dummy+1
        nodes.append(len(rs.getNodes()))
        numtips.append(len(rs.getPolylines()))

    ax[2].plot(x, nodes, color = cols[i], label = titles[i])
    ax[2].legend()
    ax[2].set_xlim([0,40])
    ax[2].set_ylim([0,10000])

    ax[i].set_xlabel("RLD (cm / cm^3)")
    ax[i].set_ylabel("Soil depth (cm)")
    ax[i].set_title(titles[i])
    ax[i].legend()
plt.show()
