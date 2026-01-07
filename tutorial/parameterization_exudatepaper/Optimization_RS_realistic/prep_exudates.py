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

path = "rootsystem/"
names = ["RS_optimized_column_L_WT", "RS_optimized_column_S_WT","RS_optimized_field_L_WT", "RS_optimized_field_S_WT"]

#import exudate data
times = [0, 42, 63, 98, 154]
exu_rates_ = [np.array([0.001, 0.001,0.00055,0.00039,0.00045]),
              np.array([0.0011, 0.0011,0.0005,0.0003,0.00033]),
              np.array([0.001, 0.001,0.00055,0.00039,0.00045]),
              np.array([0.0011, 0.0011,0.0005,0.0003,0.00033])]#[kg/(m2 day)]
              

plant_exud = np.zeros((154,len(names)))
surf = np.zeros((154,len(names)))
numtips = np.zeros((154,len(names)))

for z in range(0,len(names)): 

    exu_rates = exu_rates_[z]
    f = interpolate.interp1d(times, exu_rates)
    
    #create root system 
    rs = pb.MappedRootSystem()
    rs.readParameters(path + names[z] + ".xml")

    #get lmax of the different root types
    lmax = np.zeros((6))
    for i in range(0,6):
        p = rs.getRootRandomParameter(i)
        lmax[i] = p.lmax

    # Initialize
    rs.initializeLB(5,4)

    for j in range(0,154):

        kex = np.array([[0., 5.], [f(j), 0.]])

        if z <=2: 
            if j<=98: 
                rs.simulate(1, True);
        elif z>2: 
            rs.simulate(1, True);

        polylengths = rs.getParameter("length")
        radii = rs.getParameter("radius")
        types = rs.getParameter("type")
        polylines = rs.getPolylines()
        
        kex_ = f(j)
        sf = []
        surf_ = []
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
                #if growth has already stopped (99% of total length reached) 
                if lmax[roottype]*0.99>= polylengths[i]:
                    #print('REACHED')
                    kexu = kex_/2
                    c = 1
                if z <=2 and j>98:
                    kexu = kex_/2
                    c = 1
                #if artificial shoot 
                if roottype == 0:
                    kexu = 0
                    c = 0

                sf.append(2 * np.pi * a * l * 1.e-4 * kexu * 1.e3) # g/day/root
                surf_.append(2 * np.pi * a * l)

        surf[j,z] = np.sum(surf_) 
        plant_exud[j,z] = np.sum(sf) # g/day/plant

    del rs

np.savez('plant_exud_rates', plant_exud, surf,numtips)
