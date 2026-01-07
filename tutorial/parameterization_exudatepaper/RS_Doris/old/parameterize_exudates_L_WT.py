"""small example"""
import sys
sys.path.append("../../../../../CPlantBox")
sys.path.append("../../../../../CPlantBox/src")
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
kex = np.array([0.0019,0.0009,0.00055,0.00055]) #[kg/(m2 day)]

#import exudate data
exudate = 'Total exudation'
df = pd.read_csv("data/Total_exudation.csv")
data = df[(df['Substrate']=='L') & (df['Genotype']=="WT")]
DAS = data["DAS"].loc[:].values
mean = data[exudate].values*12*24/10**3 #mmol/plant/h --> g/plant/day
SE = data[exudate+' SE '].values*12*24/10**3

#create root system 
rs = pb.MappedRootSystem()
path = "rootsystem/"
name = "RS_optimized_L_WT"
rs.readParameters(path + name + ".xml")

#get lmax of the different root types
lmax = np.zeros((6))
for i in range(0,6):
    p = rs.getRootRandomParameter(i)
    lmax[i] = p.lmax

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
                if lmax[roottype]*0.99>= polylengths[i] or j == times[-1]:
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

print('computed exudation', plant_exud)
print('measured exudation', mean) 

#PLOT
plt.plot(DAS, mean, color = 'b', label = 'Measurement')
plt.fill_between(DAS, mean-SE, mean+SE, color = 'b', alpha = 0.2)
plt.plot(DAS, plant_exud, color = 'b',linestyle = '--', label = 'Simulation')
plt.ylabel(exudate + '\n (g / plant / d)')
plt.legend()
#plt.text(50, 33, '(c)', fontsize=16,  va='top', ha='right')
plt.tight_layout()
plt.show()

