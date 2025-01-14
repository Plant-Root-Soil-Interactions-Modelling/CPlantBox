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
from scipy import interpolate

font = {'size'   : 18}
plt.rc('font', **font)

#import exudate data
exudates = ['Phenolics', 'Sugars', 'AminoAcids']
times = [0, 42, 63, 98, 154]
exu_rates = np.array([[0.00012, 0.00012, 0.00003, 0.000019, 0.000015],[0.0006, 0.0006, 0.00015, 0.00017, 0.00016],[0.0001, 0.0001, 0.000025, 0.0000168, 0.0000135]]) #[kg/(m2 day)]
plant_exud = np.zeros((154,len(exudates)))
surf = np.zeros((154,len(exudates)))
numtips = np.zeros((154,len(exudates)))


for i in range(0,len(exudates)): 

    #interpolate exudation rates
    f = interpolate.interp1d(times, exu_rates[i])
    
    #create root system 
    rs = pb.MappedRootSystem()
    path = "rootsystem_params/"
    name = "Zea_mays_5_Leitner_Streuber_2014_mod2"
    rs.readParameters(path + name + ".xml")

    # Initialize
    rs.initializeLB(5,4)

    for j in range(0,154):

        kex = np.array([[0., 5.], [f(j), 0.]])

        if j<=98: 
            rs.simulate(1, True);
        segs = rs.segments
        tipI = np.array(rs.getRootTips())-1
        numtips[j,i] = len(tipI)
        polylengths = rs.getParameter("length")
        types = rs.getParameter("type")
        a = rs.radii
        l = rs.segLength()
        print('number of segs', len(l))
        sf = []
        surf_ = []
        for k, s in enumerate(tipI):
            if types[k] != 0:
                l_ = 0
                s_ = s
                while l_<= kex[0][1]:
                    l_ = l_+l[s_]
                    kexu_ = (kex[1][1]-kex[1][0])/(kex[0][1]-kex[0][0])*l_+kex[1][0]
                    kexu = max(0,kexu_)
                    sf.append(2 * np.pi * a[s_] * l[s_] * 1.e-4 * kexu * 1.e3) # g/day/root
                    surf_.append(2 * np.pi * a[s_] * l[s_])
                    s_ = s_-1
                    if l_>polylengths[k]:
                        break

        surf[j,i] = np.sum(surf_) 
        plant_exud[j,i] = np.sum(sf) # g/day/plant
    del rs


np.savez('plant_exud_rates', plant_exud, surf,numtips)
