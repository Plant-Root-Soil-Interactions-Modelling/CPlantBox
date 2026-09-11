"""optimize root length / diameter"""
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
import matplotlib as mpl 
from scipy import interpolate

WRITE_ROOTSYS = False

fontsize = 24
mpl.rcParams.update({
    'font.size': fontsize,
    'axes.titlesize': fontsize,
    'axes.labelsize': fontsize,
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,
    'legend.fontsize': fontsize,
    'figure.titlesize': fontsize,
    'mathtext.default': 'regular'
})

soiltypes = ['L', 'S']
titles = ['Loam', 'Sand']
cols = ['r', 'b']
fig, ax = plt.subplots(3,1, figsize = (15, 10))
for i in range(0, len(soiltypes)): 
    
    #Experimental data
    df = pd.read_csv("data/column_experiment_mean.csv")
    data = df[(df['substrate']==soiltypes[i]) & (df['genotype']=="WT")]
    DAS = data["DAS"].loc[:].values
    real_length_ = data['RL_mean'].loc[:].values
    real_length_SE_ = data['RL_SE'].loc[:].values
    real_diam_ = data['RD_mean'].loc[:].values
    real_diam_SE_ = data['RD_SE'].loc[:].values
    real_exu_ = data["Total exudation"].values*24/10**3 #mmol/plant/h --> mol/plant/day
    real_exu_SE_ = data["Total exudation SE"].values*24/10**3 #mmol/plant/h --> mol/plant/day

    #target times 
    times_ = DAS[:2]
    real_length = real_length_[:len(times_)]
    real_length_SE = real_length_SE_[:len(times_)]
    real_diam = real_diam_[:len(times_)] 
    real_diam_SE = real_diam_SE_[:len(times_)] 
    real_exu = real_exu_[:len(times_)]
    real_exu_SE = real_exu_SE_[:len(times_)]


    # plot simulated and measured total root length 
    path = "rootsystem/"
    name = "RS_optimized_"+soiltypes[i]+"_WT"
    rs = pb.RootSystem()
    rs.readParameters(path + name + ".xml")
    rs.setSeed(0)

    #simulations without tube 
    rs.initializeLB(5,4)
    times = np.linspace(0,times_[-1],times_[-1]+1)
    times = times.astype(int)
    print('time', times) 
    comp_length = np.zeros((len(times)))
    comp_diam = np.zeros((len(times)))

    for day in range(1,times[-1]+1):
        rs.simulate(1, True)
        if WRITE_ROOTSYS: 
            rs.write("RS_visualisation/notube_"+soiltypes[i]+"_WT_day_"+str(day)+".vtp")


    #simulations with tube 
    tube_ = pb.SDF_PlantContainer(10, 10, 80, False) #tube with diameter 20, length 60 cm 
    tube = pb.SDF_RotateTranslate(tube_, 0, pb.SDF_Axis.zaxis, pb.Vector3d(0, 0, 10))
    rs.setGeometry(tube)
    rs.initializeLB(5,4)
    times = np.linspace(0,times_[-1],times_[-1]+1)
    times = times.astype(int)
    print('time', times) 
    comp_length = np.zeros((len(times)))
    comp_diam = np.zeros((len(times)))
    comp_exu = np.zeros((len(times)))
    
    #get lmax of the different root types
    lmax = []
    for pp in rs.getRootRandomParameter():
        lmax.append(pp.lmax)
        
    df = pd.read_csv("data/exudation_rates.csv")
    DAS = df["DAS"].values
    DAS = np.insert(DAS, 0, 0)
    exu_rates = df[titles[i]].values
    exu_rates= np.insert(exu_rates, 0,exu_rates[0])
    f = interpolate.interp1d(DAS, exu_rates)
    tip = 3.5 #cm

    for day in range(1,times[-1]+1):
        rs.simulate(1, True)
        
        #Length and diameter
        length = np.array(rs.getParameter("length"))
        radius = np.array(rs.getParameter("radius"))
        meanrad = np.sum(length*radius)/np.sum(length)

        comp_length[day] =  np.sum(length) #cm
        comp_diam[day] =  meanrad*20 #mm
        
        #write to vtp 
        if WRITE_ROOTSYS: 
            rs.write("RS_visualisation/withtube_"+soiltypes[i]+"_WT_day_"+str(day)+".vtp")
        
        #compute exudation 
        kex_tip = f(day)
        kex_base = kex_tip/2
        
        polylengths = np.asarray(rs.getParameter("length"))
        radii = np.asarray(rs.getParameter("radius"))
        types = np.asarray(rs.getParameter("type"))
        polylines = rs.getPolylines()
        

        sf = []
        for j in range(0, len(polylines)):
            a = radii[j]
            roottype = int(types[j])
            l_ = 0
            for k in range(0, len(polylines[j])-1):

                m = polylines[j][-1-k]
                n = polylines[j][-2-k]
                p0 = np.array([m.x, m.y, m.z])
                p1 = np.array([n.x, n.y, n.z])
                l  = np.linalg.norm(p0 - p1)
                l_ = l_+l
                #tip exudation rate 
                if (polylengths[j]-l_)<tip:
                    kexu = kex_tip
                    c = 2
                #base exudation rate 
                else:
                    kexu = kex_base
                    c = 1
                #if growth has already stopped (95% of total length reached) 
                if polylengths[j]>=lmax[roottype]*0.99:
                    #print('REACHED')
                    kexu = kex_base
                    c = 1
                #if artificial shoot 
                if roottype == 0:
                    kexu = 0
                    c = 0

                sf.append(2 * np.pi * a * l * kexu) # mol/root segment / day

        comp_exu[day] = np.sum(sf) # mol/plant/day
        

    ax[0].plot(times, comp_length/100, linestyle = '--', color = cols[i])
    ax[0].errorbar(times_, real_length/100, real_length_SE/100, color = cols[i], fmt='o')
    if i == 1: 
        ax[0].plot(np.nan, np.nan, linestyle = '-', color = cols[0], label = titles[0])
        ax[0].plot(np.nan, np.nan, linestyle = '-', color = cols[1], label = titles[1])
        ax[0].errorbar(np.nan, np.nan, np.nan, marker = 'o', color = 'k', label = 'Experiment')
        ax[0].plot(np.nan, np.nan, color = 'k', linestyle = '--', label = 'Simulation')
    ax[0].legend(loc = 'upper left')
    ax[0].set_ylabel('Total root \nlength (m)') 

    ax[1].plot(times, comp_diam, color = cols[i], linestyle = '--')
    ax[1].errorbar(times_, real_diam, real_diam_SE, color =cols[i], fmt='o')
    ax[1].set_ylim([0,1.5])
    ax[1].set_ylabel('Mean root \ndimeter (mm)') 
    
    ax[2].plot(times, comp_exu, color = cols[i], linestyle = '--')
    ax[2].errorbar(times_, real_exu, real_exu_SE, color =cols[i], fmt='o')
    ax[2].set_xlabel('Time') 
    ax[2].set_ylabel('Total plant exudation \n(mol/day/plant)')     
    

plt.show()

