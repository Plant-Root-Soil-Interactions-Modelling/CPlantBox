#!/usr/bin/python3

'''
#--------------------------#
#--- SIMPLACE-CPlantBox ---#
#--------------------------#

@date: 13-Apr-2022
@author: Murilo Vianna <mvianna@uni-bonn.de>

Before starting, please make sure you have SIMPLACE and CPlantBox installed.
A brief guidance is given in SIMPLACE_install_bonnares.sh and CPlantBox_install_bonnares.sh files
This implementation was developed based on previous SIMPLACE-CRootBox version: 
https://svn.simplace.net:8443/projects/simplace_run/repository/show/trunk/simulation/sabine/CRootBoxRuns
https://doi.org/10.3389/fpls.2022.865188

IMPORTANT: Before running this scripts, make sure the dependency modules are loaded in bonnaHPC:
# module load Python/3.8.2-GCCcore-9.3.0
# module load VTK/8.2.0-foss-2020a-Python-3.8.2
# module load CMake
'''

#--- import libs
import sys
import os
import argparse as ap
import numpy as np

#--- get home dir from ~
home_dir = os.path.expanduser('~')

#--- append CPlantBox module paths
sys.path.append(home_dir+'/workspace/CPlantBox')
sys.path.append(home_dir+'/workspace/CPlantBox/src/')
sys.path.append(home_dir+'/workspace/simplace_wrapper/python_simplace/trunk/PyPlantBox/misc')

#--- load CPlantBox modules
import plantbox as rb
import visualisation.vtk_plot as vp 

# import py_rootbox as rb
# from rb_tools import v2a, a2v
import simplace
from SimplaceRootbox.simplace import lintulslim_interact as sp_int
# from SimplaceRootbox.util import plotstep


# Get command line arguments

parser = ap.ArgumentParser(description=__doc__, formatter_class=ap.RawDescriptionHelpFormatter)
parser.add_argument('-s','--solutionid',type=int, default = 0, help='0: rootbox<->simplace; 1: rootbox//simplace')
parser.add_argument('-r','--rainscale',type=float, default = 1, help="scale rain from 0 to 1")
parser.add_argument('-c','--crop',type=int, default = 0, help="crop - 0:anagallis, 1:wheat")
parser.add_argument('-e','--elongationrestriction',type=int, default = 1, help="elongation restriction by soil/water - 0:no, 1:yes")
parser.add_argument('-n','--numberofsteps',type=int, default = 1, help="number of steps per day")
parser.add_argument('-p','--plot', type=int, default = 0, help='0: no plot; 1:plot during execution')
args = parser.parse_args()

#------------------------------#
#--- CRootBox configuration ---#
#------------------------------#
gram_per_cm = .000035
soildepth = 120
layers = 40
area = 6*12.5 # cm * cm
accuracy = 0.1 # cm
#rem maxiter = 10
simtime = 600 # days
maxinc = 20; # maximal length increment (cm/day), TODO base this value on some fancy model 
#rem crops = ['other/anagallis','wheat'] # name of parameter files located in folder modelparameter
plot3D_onscreen = True # want to see 3D roots by the end?

steps = args.numberofsteps
dt = 1/steps
N = simtime * steps # steps

#------------------------------#
#--- SIMPLACE configuration ---#
#------------------------------#

# Simplace configuration
# solutiondir = 'mvianna/CPlantBox/solution/' use this path to run on simplace_run
solutiondir = 'PyPlantBox/solution/'
solutions = ['LintulSlimRB.sol.xml','LintulSlimRef.sol.xml'] ## now with only one homogenous soil layer
# rem runs = ['rootbox','reference'] # only tested for 'rootbox' solution
jd = home_dir+'/workspace/'

# wd = jd + 'simplace_run/simulation/' use this path to run on simplace_run
wd = jd + 'simplace_wrapper/python_simplace/trunk/'
od = jd + 'simplace_run/output/'

# Set up Simplace and create the simulation
sol = wd + solutiondir +solutions[args.solutionid]

sim = simplace.SimplaceInstance(jd, wd, od, wd, wd)
sim.setLogLevel('ERROR')
sim.openProject(sol)

par = {'startdate':"20.03.1991","vRainScale":args.rainscale,"projectid":str(args.rainscale)}
sim.createSimulation(par)
ids = sim.getSimulationIDs()

# Parameters for elongation restriction due to soilwater in simplace
h1 = 0
h2 = -1
h3 = -5
h4 = -150

#--- read CRootBox Parameters
#--- testing for default anagalis parameters:
#--- https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/blob/master/modelparameter/rootsystem/Anagallis_femina_Leitner_2010.xml
# name = wd + 'mvianna/CPlantBox/data/modelparameter/' + 'Anagallis_femina_Leitner_2010.xml' use this path to run on simplace_run
name = wd + 'PyPlantBox/data/modelparameter/' + 'wheat.xml'

print(solutions[args.solutionid]," - vRainScale:",args.rainscale,' - root:',name,' - steps:',steps, ' - elongation restriction:',args.elongationrestriction)

#--- Set scale elongaton according to tutorial example "example5b_scaleelongation.py"
scale_elongation = rb.EquidistantGrid1D(0, -soildepth, layers)
scale_elongation.data = np.ones((layers))
se = rb.ProportionalElongation()
se.setBaseLookUp(scale_elongation)

# Initialize root system
rs = rb.RootSystem()
rs.setSeed(0)
rs.readParameters(name) 
for p in rs.getRootRandomParameter():
    p.f_se = se  # set scale elongation function

rs.initialize()  

# initialize length and rld states
ol = 0
vRLD=np.zeros(layers)

#if args.plot==1:
#    (ylim, p,a,g,h) = plotstep.init_plot();

# Simulation loop
for s in range(0,simtime):    
    
    # simulate simplace step to get maxinc and re_reduction dinamically
    (date, maxinc, doharvest, tranrf, yld, rld_s, re_reduction) = sp_int.getSimplaceValuesExtended(sim,gram_per_cm, h1, h2, h3, h4, area=area)
    
    #--- if not using re_reduction from simplace set args.elongationrestriction == 0
    if args.elongationrestriction == 0:
        re_reduction = [1.]*(layers+1)    
    
    #--- update scale elongation
    scale_elongation.data = np.array(re_reduction) 
        
    print("Simulating: ",date," MaxIncr:", round(maxinc,2), " Tranrf:", round(tranrf,2), " Yield:", round(yld,2), " Step:", s)
    inc = 0.0
    
    # run CPlantBox if there's any root increment
    if(maxinc > 0):
        for j in range(0, steps): 
            rs.simulate(dt,maxinc*dt, se, True) # not working: it runs but the code breaks when rb.SegmentAnalyser(rs)
            #rs.simulate(dt,True)            
        
        # get simulated root length increment
        l = np.sum(rs.getParameter('length'))
        inc =  l - ol
        ol = l
        
        # calculate RLD for each soil layer using rb.SegmentAnalyser(rs)
        vRLD = sp_int.calculateRLD(rs, soildepth, layers, area)
    
    # update vRLD in simplace      
    sp_int.setSimplaceRoots(sim, vRLD, maxinc-inc, gram_per_cm)
    
    #if args.plot==1:
    #    (ylim, p, a, g, h) = plotstep.plot_step(vRLD, rld_s, ylim, p, a, g, h)     

    if doharvest == True:
        print("Harvested")        
        rs.write(od+"/PyPlantBox/lintul/LintulSlim_PlantBox.vtp")
        if plot3D_onscreen: vp.plot_roots(rs, "type")
        
        print("Harvested")
        break


sim.closeProject()
sim.shutDown()

