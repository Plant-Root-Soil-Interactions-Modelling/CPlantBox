#!/usr/bin/python3

# import libs
import sys
import os
import numpy as np
import pandas as pd
import datetime
from pathlib import Path

# get home dir from ~
home_dir = os.path.expanduser('~')

# append CPlantBox module paths
sys.path.append(home_dir+'/workspace/CPlantBox')
sys.path.append(home_dir+'/workspace/CPlantBox/src/') 
sys.path.append(home_dir+'/workspace/simplace_run/simulation/milena/PyPlantBox/misc')

# load CPlantBox modules and SIMPLACE
import plantbox as pb
import simplace # https://simplace.net/doc/python_wrapper.htm

# custom functions
import visualisation.vtk_plot as vp
from SimplacePlantbox.simplace import lintulslim_interact as sp_int
from SimplacePlantbox.util import getRootOutputs
from SimplacePlantbox.util import checkRuns

# Simulation configuration
gram_per_cm = .000035       # specific root length density used to calculate the maximum root increment in a timestep [g/cm]
area = 6*12.5               # plant area [cm * cm]
simtime = 600               # maximum simulation time-steps [days]
dt = 1                      # simulation timestep [days]
plot = False                # plot root system at the end ?

# Feddes parameters for root elongation restriction due to soil water potential
# Eq. 10 of https://doi.org/10.3389/fpls.2022.865188
h1 = 0; h2 = -1; h3 = -5; h4 = -150

# SIMPLACE initialization |\label{l7_5_simplace:InitStart}|

# ModelSolution
solutiondir = 'milena/PyPlantBox/solution/Lintul5SlimPB.sol.xml'
rainscale = 1.0 # scale rain from 0 to 1
elongationrestriction = int(1) # elongation restriction by soil/water - 0:no, 1:yes

# Paths
jd  = home_dir+'/workspace/'
wd  = jd + 'simplace_run/simulation/'
od  = jd + 'simplace_run/output/'
sol = wd + solutiondir

# read soil information from the ModelSolution
sol_lines = sp_int.ReadFile(sol)
vSoilFile = sp_int.get_SolutionVariable(sol_lines,'variables','vSoilFile')
vSoilName = sp_int.get_SolutionVariable(sol_lines,'variables','vSoilName') 
vSoilFile_path = sp_int.get_SolutionVariable(sol_lines,'interface','vSoilFile')
layerthickness = float(sp_int.get_SolutionVariable(sol_lines,'transform','layerthickness')) * 100 # cm

# read soil depth and layers for vertical grid from the soil file used in the simplace solution
vSoilFile_path = str(vSoilFile_path).replace('${_WORKDIR_}/', wd).replace('${vSoilFile}',vSoilFile)
SoilFile = pd.read_csv(vSoilFile_path, sep=';')
soildepth = max(SoilFile.loc[SoilFile['soilname'] == vSoilName,'depth']) * 100 # cm
layers = soildepth / layerthickness # number of soil layers [#]

soildepth = int(soildepth)
layers = int(layers)

# Instantiate simplace for the ModelSolution
sim = simplace.SimplaceInstance(jd, wd, od, wd, wd)
sim.setLogLevel('ERROR')
sim.openProject(sol)
par = {'startdate':"01.01.2022","vRainScale":rainscale,"projectid":str(rainscale)}
sim.createSimulation(par)
ids = sim.getSimulationIDs()

# CPlantBox initialization

# read root Parameters
rootparname = wd + 'milena/PyPlantBox/data/modelparameter/sbarley_cPlantBox_Sabine_modified.xml'

# Set scale elongaton according to tutorial:
# tutorial/chapter3_responses/example3_1_carbon.py
scale_elongation = pb.EquidistantGrid1D(0, -soildepth, layers)
scale_elongation.data = np.ones((layers))
se = pb.ProportionalElongation()
se.setBaseLookUp(scale_elongation)

rs = pb.Plant()
rs.setSeed(0)
rs.readParameters(rootparname)
for p in rs.getOrganRandomParameter(pb.root):
    p.f_se = se  # set scale elongation function

# intialize root system
rs.initialize()

# initialize length and rld states
ol = 0
vRLD=np.zeros(layers)

# switcher to initialize PB [only True when first root growth happens]
init_pb_switch = False # |\label{l7_5_simplace:InitEnd}|

# logs
warns = [str(datetime.datetime.now())]

# Simulation loop
for s in range(0,simtime): # |\label{l7_5_simplace:LoopStart}|
    
    # simulate simplace step to get maxinc and re_reduction dinamically
    (date, maxinc, doharvest, tranrf, yld, rld_s, re_reduction, re_q, re_w, h, init_pb, frr_s, md95_s) = sp_int.getSimplaceValuesExtended(sim,gram_per_cm, h1, h2, h3, h4, area=area) # |\label{l7_5_simplace:MaxInc_simplace}|
    
    # if not using re_reduction from simplace set args.elongationrestriction == 0
    if elongationrestriction == 0:
        re_reduction = [1.]*(layers+1)
    scale_elongation.data = np.array(re_reduction)
    
    maxinc = round(maxinc,2)
    inc = 0.0
    print("Simulating: ",date," MaxIncr:", round(maxinc,2), " Tranrf:", round(tranrf,2), " Yield:", round(yld,2), " Step:", s)
    # run CPlantBox if there's any root increment
    if(maxinc > 0): # |\label{l7_5_simplace:RunCPB}|
        # simulate root system
        rs.simulate(dt,maxinc,se,True)
        
        # get simulated root length increment
        l = np.sum(rs.getParameter('length'))
        inc =  l - ol
        ol = l
        print('increase was '+str(round(inc, 2)))
        if inc > maxinc+0.1:
            msg = 'Warning: Root increment exceeds max increment by '+str(round(inc-maxinc,1))+'cm'
            print(msg); warns.append(msg)
        
        # calculate RLD for each soil layer using pb.SegmentAnalyser(rs)
        vRLD = sp_int.calculateRLD(rs, soildepth, layers, area) # |\label{l7_5_simplace:RLD_CPB}|
    
    # update vRLD in simplace
    sp_int.setSimplaceRoots(sim, vRLD, maxinc-inc, gram_per_cm, init_pb_switch, frr_s, md95_s, re_reduction, re_q, re_w, h) # |\label{l7_5_simplace:RLD_update_simplace}|
    
    # get pb-related variables
    if s == 0:
      out_pb = getRootOutputs.asDataFrame(rs, date)
    else:
      out_pb = out_pb.append(getRootOutputs.asDataFrame(rs, date), ignore_index=True)    
    
    # Initialize PB in next timestep
    if init_pb: init_pb_switch = True
    
    # check if harvest occurs before end of sims
    if doharvest: print("Harvested"); break # |\label{l7_5_simplace:LoopEnd}|

# write pb outputs |\label{l7_5_simplace:OutStart}|
rs.write(od+"/milena/PyPlantBox/lintul5/Lintul5Slim_PlantBox.vtp")
out_pb.to_csv(od + '/milena/PyPlantBox/lintul5/PlantBox_outputs.csv', index=False)

# plot on screen?
if plot:
    vp.plot_roots(rs, "type")

# check RootAges and shutdown simplace instance
warns = checkRuns.checkRootAge(rs, warns)
sim.closeProject()
sim.shutDown() # |\label{l7_5_simplace:OutEnd}|