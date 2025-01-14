"""small example"""
import sys
sys.path.append("../../../../CPlantBox")
sys.path.append("../../../../CPlantBox/src")
import plantbox as pb
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

#import exudate data
exudate = 'Phenolics'
df = pd.read_csv("../data_magda/mean_measures_Eva_Michael.csv")
DAS = df["DAS"].loc[:].values
mean = df[exudate].values*12*24/10**3 #g/day/plant
std = df[exudate+' std'].values*12*24/10**3 #g/day/plant

#create root system 
#rs = pb.RootSystem()
rs = pb.MappedRootSystem()
path = "rootsystem_params/"
name = "Zea_mays_5_Leitner_Streuber_2014_mod2"
rs.readParameters(path + name + ".xml")

# Initialize
rs.initializeLB(5,4)

times = [42, 63, 98, 154]
plant_exud = np.zeros((len(times)))
exu_prop1 = 0.00012 #[kg/(m2 day)]
exu_prop2 = 0.00003 #[kg/(m2 day)]
exu_prop3 = 0.000019 #[kg/(m2 day)]
exu_prop4 = 0.000015 #[kg/(m2 day)]
dummy = 0

#exudation rates per root type
#[age,age] [value,  value] for each time 
kex = np.array([[[0., 5.], [exu_prop1, 0.]],[[0., 5.], [exu_prop2, 0.]],[[0., 5.], [exu_prop3, 0.]], [[0., 5.], [exu_prop4, 0.]]])



for j in range(0,times[-1]+1):

    if j<=times[2]: 
        rs.simulate(1, True);

    if j in times:

        segs = rs.segments
        tipI = np.array(rs.getRootTips())-1
        polylengths = rs.getParameter("length")
        types = rs.getParameter("type")
        
        
        #types = np.asarray(rs.subTypes, int)
        a = rs.radii
        l = rs.segLength()
        print('number of segs', len(l))
        sf = []

        for i, s in enumerate(tipI):
            if types[i] == 0:
                kexu_ = 0
            else: 
                kex_ = kex[dummy]
                l_ = 0
                s_ = s
                while l_<= kex_[0][1]:
                    #print(s_, l_,'huhu') 
                    l_ = l_+l[s_]
                    kexu_ = (kex_[1][1]-kex_[1][0])/(kex_[0][1]-kex_[0][0])*l_+kex_[1][0]
                    kexu = max(0,kexu_)

                    sf.append(2 * np.pi * a[s_] * l[s_] * 1.e-4 * kexu * 1.e3) # g/day/root
                    s_ = s_-1
                    if l_>polylengths[i]:
                        #print('here its done')
                        break

        
        plant_exud[dummy] = np.sum(sf) # g/day/plant
        dummy = dummy+1


print('computed exudation', plant_exud)
print('measured exudation', mean) 

#PLOT
plt.plot(df['DAS'], mean, color = 'b', label = 'Measurement')
plt.fill_between(df['DAS'], mean-std, mean+std, color = 'b', alpha = 0.2)
plt.plot(df['DAS'], plant_exud, color = 'b',linestyle = '--', label = 'Simulation')
plt.ylabel(exudate + '\n (g / plant / d)')
plt.legend()
#plt.text(50, 33, '(c)', fontsize=16,  va='top', ha='right')
plt.tight_layout()
plt.show()

