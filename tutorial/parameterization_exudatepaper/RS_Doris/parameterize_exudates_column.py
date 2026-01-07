"""small example"""
import sys
sys.path.append("../../../../CPlantBox")
sys.path.append("../../../../CPlantBox/src")
import plantbox as pb
from functional.xylem_flux import XylemFluxPython
import visualisation.vtk_plot as vp
import vtk 
import math
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import pandas as pd
import timeit


font = {'size'   : 18}
plt.rc('font', **font)

#tip exudation rates per time step - to be optimized
kexx = [np.array([0.0015,0.0009,0.00045,0.00044]),
       np.array([0.002,0.0006,0.00034,0.00033])]#[kg/(m2 day)]

#import exudate data
exudate = 'Total exudation'
df = pd.read_csv("data/column_experiment_mean.csv")
DAS = np.unique(df["DAS"].loc[:].values)
soils = df["substrate"].loc[:].values
soils_ = np.unique(soils)

cols = ['b', 'r']

#create root system 
rs = pb.MappedRootSystem()
path = "rootsystem/"
names = ["RS_optimized_column_L_WT_new", "RS_optimized_column_S_WT_new"]


for p in range(0, len(names)):

    name = names[p]

    rs.readParameters(path + name + ".xml")

    #get lmax of the different root types
    lmax = []
    for pp in rs.getRootRandomParameter():
        lmax.append(pp.lmax)


    # Initialize
    rs.initializeLB(5,4)

    times = [42, 63, 98, 154]
    plant_exud = np.zeros((len(times)))
    dummy = 0
    cols = ['b', 'r', 'g']


    for j in range(0,times[-1]+1):

        if j<=times[2]: 
            rs.simulate(1, True);

        if j in times:
            polylengths = rs.getParameter("length")
            radii = rs.getParameter("radius")
            types = rs.getParameter("type")
            polylines = rs.getPolylines()

            #rs.write('test_exudate_RS/day'+str(j)+'.vtp')
            #ax = plt.axes(projection='3d')

            kex = kexx[p]
            kex_ = kex[dummy]
            sf = []
            for i in range(0, len(polylines)):
                a = radii[i]
                roottype = int(types[i])
                l_ = 0
                for k in range(0, len(polylines[i])-1):

                    m = polylines[i][-1-k]
                    n = polylines[i][-2-k]
                    p0 = np.array([m.x, m.y, m.z])
                    p1 = np.array([n.x, n.y, n.z])
                    l  = np.linalg.norm(p0 - p1)
                    l_ = l_+l
                    #tip exudation rate 
                    if l_<3.5:
                        kexu = kex_
                        c = 2
                    #base exudation rate 
                    else:
                        kexu = kex_/2
                        c = 1
                    #if growth has already stopped (95% of total length reached) 
                    if lmax[roottype]*0.99>= polylengths[i] or j>times[2]:
                        #print('REACHED')
                        kexu = kex_/2
                        c = 1
                    #if artificial shoot 
                    if roottype == 0:
                        kexu = 0
                        c = 0

                    sf.append(2 * np.pi * a * l * 1.e-4 * kexu * 1.e3) # g/day/root

                    #test plot 
                    #x, y, z = [p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]]
                    #ax.plot3D(x, y, z, cols[c])

            plant_exud[dummy] = np.sum(sf) # g/day/plant
            dummy = dummy+1
            #plt.show()

    #measured data
    data = df[(df['substrate']==soils_[p]) & (df['genotype']=="WT")]
    mean = data[exudate].values*12*24/10**3 #mmol/plant/h --> g/plant/day
    print(mean)
    SE = data[exudate+' SE'].values*12*24/10**3

    print('computed exudation', plant_exud)
    print('measured exudation', mean)
    
    #PLOT
    plt.plot(DAS, mean, color = cols[p])
    plt.fill_between(DAS, mean-SE, mean+SE, color = cols[p], alpha = 0.2)
    plt.plot(DAS, plant_exud, color = cols[p],linestyle = '--')

plt.plot(np.nan,np.nan,color = 'b', label = 'Loam')
plt.plot(np.nan,np.nan,color = 'r', label = 'Sand')
plt.plot(np.nan,np.nan,color = 'k',linestyle = '-', label = 'Measurement')
plt.plot(np.nan,np.nan,color = 'k',linestyle = '--', label = 'Simulation')
plt.ylabel(exudate + '\n (g / plant / d)')
plt.legend()
#plt.text(50, 33, '(c)', fontsize=16,  va='top', ha='right')
plt.tight_layout()
plt.show()

