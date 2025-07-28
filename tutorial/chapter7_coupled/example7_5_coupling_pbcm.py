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
import pandas as pd
import datetime
from pathlib import Path

#--- get home dir from ~
home_dir = os.path.expanduser('~')

#--- append CPlantBox module paths
sys.path.append(home_dir+'/workspace/CPlantBox')
sys.path.append(home_dir+'/workspace/CPlantBox/src/') 
sys.path.append(home_dir+'/workspace/simplace_run/simulation/milena/PyPlantBox/misc')

#--- load CPlantBox modules and simplace
import plantbox as rb
#import vtk_plot as vp
import visualisation.vtk_plot as vp # for Updated
import simplace
from SimplacePlantbox.simplace import lintulslim_interact as sp_int
#from SimplacePlantbox.util import plotstep
from SimplacePlantbox.util import getRootOutputs
from SimplacePlantbox.util import checkRuns

#--- Get command line arguments
parser = ap.ArgumentParser(description=__doc__, formatter_class=ap.RawDescriptionHelpFormatter)
parser.add_argument('-r','--rainscale',type=float, default = 1, help="scale rain from 0 to 1")
parser.add_argument('-c','--crop',type=int, default = 0, help="crop - 0:anagallis, 1:wheat")
parser.add_argument('-e','--elongationrestriction',type=int, default = 1, help="elongation restriction by soil/water - 0:no, 1:yes")
parser.add_argument('-n','--numberofsteps',type=int, default = 1, help="number of steps per day")
parser.add_argument('-p','--plot', type=int, default = 0, help='0: no plot; 1:plot during execution; 2: only at the end of sims')
parser.add_argument('-v','--verb', type=bool, default = True, help='verbose: True/False')
args = parser.parse_args()

#--------------------------------#
#--- Simulation configuration ---#
#--------------------------------#
gram_per_cm = .000035       # specific root length density used to calculate the maximum root increment in a timestep [g/cm]
area = 6*12.5               # plant area [cm * cm]
simtime = 600               # maximum simulation time-steps [days]
dt = 1                      # simulation timestep [days]

# Feddes parameters for root elongation restriction due to soil water potential
# Eq. 10 of https://doi.org/10.3389/fpls.2022.865188 (MvG parameters are in the SoilProperties.csv)
h1 = 0
h2 = -1
h3 = -5
h4 = -150

#------------------------------#
#--- SIMPLACE configuration ---#
#------------------------------#

#--- ModelSolution
solutiondir = 'milena/PyPlantBox/solution/Lintul5SlimPB.sol.xml'

#--- Paths
jd  = home_dir+'/workspace/'
wd  = jd + 'simplace_run/simulation/'
od  = jd + 'simplace_run/output/'
sol = wd + solutiondir

#--- read soil information from the ModelSolution
#--- note: we are using custom functions as reading xml was not easily done due to WIKI and comments
sol_lines = sp_int.ReadFile(sol)
vSoilFile = sp_int.get_SolutionVariable(sol_lines,'variables','vSoilFile')
vSoilName = sp_int.get_SolutionVariable(sol_lines,'variables','vSoilName') 
vSoilFile_path = sp_int.get_SolutionVariable(sol_lines,'interface','vSoilFile')
layerthickness = float(sp_int.get_SolutionVariable(sol_lines,'transform','layerthickness')) * 100 # cm

#--- read soil depth and layers for vertical grid from the soil file used in the simplace solution
vSoilFile_path = str(vSoilFile_path).replace('${_WORKDIR_}/', wd).replace('${vSoilFile}',vSoilFile)
SoilFile = pd.read_csv(vSoilFile_path, sep=';')
soildepth = max(SoilFile.loc[SoilFile['soilname'] == vSoilName,'depth']) * 100 # cm
layers = soildepth / layerthickness # number of soil layers [#]

soildepth = int(soildepth)
layers = int(layers)

#--- Instantiate simplace for the ModelSolution
sim = simplace.SimplaceInstance(jd, wd, od, wd, wd)
sim.setLogLevel('ERROR')
sim.openProject(sol)
par = {'startdate':"01.01.2022","vRainScale":args.rainscale,"projectid":str(args.rainscale)}
sim.createSimulation(par)
ids = sim.getSimulationIDs()

#--- read CRootBox Parameters
rootparname = wd + 'milena/PyPlantBox/data/modelparameter/sbarley_cPlantBox_Sabine_modified.xml'

#--- Set scale elongaton according to tutorial example "example5b_scaleelongation.py"
scale_elongation = rb.EquidistantGrid1D(0, -soildepth, layers)
scale_elongation.data = np.ones((layers))
se = rb.ProportionalElongation()
se.setBaseLookUp(scale_elongation)

#--- Update root parameters
rs = rb.RootSystem()
rs.setSeed(0)
rs.readParameters(rootparname)
for p in rs.getRootRandomParameter():
    p.f_se = se  # set scale elongation function

#--- intialize root system
rs.initialize()  

#--- initialize length and rld states
ol = 0
vRLD=np.zeros(layers)

#--- switcher to initialize PB [only True when first root growth happens]
init_pb_switch = False

#--- logs
warns = [str(datetime.datetime.now())]

if args.verb: 
    print(" - vRainScale:",args.rainscale,' - root:',rootparname, ' - elongation restriction:',args.elongationrestriction)

#--- Simulation loop
for s in range(0,simtime):    
    
    #--- simulate simplace step to get maxinc and re_reduction dinamically
    (date, maxinc, doharvest, tranrf, yld, rld_s, re_reduction, re_q, re_w, h, init_pb, frr_s, md95_s) = sp_int.getSimplaceValuesExtended(sim,gram_per_cm, h1, h2, h3, h4, area=area)
    
    #--- if not using re_reduction from simplace set args.elongationrestriction == 0
    if args.elongationrestriction == 0:
        re_reduction = [1.]*(layers+1)    
    scale_elongation.data = np.array(re_reduction)
    
    maxinc = round(maxinc,2)
    inc = 0.0       
    if args.verb: print("Simulating: ",date," MaxIncr:", round(maxinc,2), " Tranrf:", round(tranrf,2), " Yield:", round(yld,2), " Step:", s)
    #--- run CPlantBox if there's any root increment
    if(maxinc > 0):        
        #--- simulate root system
        rs.simulate(dt,maxinc,se,args.verb)
        
        #--- get simulated root length increment
        l = np.sum(rs.getParameter('length'))
        inc =  l - ol
        ol = l
        if args.verb: print('increase was '+str(round(inc, 2)))        
        if inc > maxinc+0.1:
            msg = 'Warning: Root increment exceeds max increment by '+str(round(inc-maxinc,1))+'cm'
            print(msg); warns.append(msg)           
        
        #--- calculate RLD for each soil layer using rb.SegmentAnalyser(rs)
        vRLD = sp_int.calculateRLD(rs, soildepth, layers, area)
    
    #--- update vRLD in simplace      
    sp_int.setSimplaceRoots(sim, vRLD, maxinc-inc, gram_per_cm, init_pb_switch, frr_s, md95_s, re_reduction, re_q, re_w, h)
    
    #--- get pb-related variables
    if s == 0:
      out_pb = getRootOutputs.asDataFrame(rs, date)
    else:
      out_pb = out_pb.append(getRootOutputs.asDataFrame(rs, date), ignore_index=True)    
    
    #--- Initialize PB in next timestep
    if init_pb: init_pb_switch = True
    
    #--- check if harvest occurs before end of sims
    if doharvest: print("Harvested"); break

#--- write pb outputs
rs.write(od+"/milena/PyPlantBox/lintul5/Lintul5Slim_PlantBox.vtp")
out_pb.to_csv(od + '/milena/PyPlantBox/lintul5/PlantBox_outputs.csv', index=False)

#--- plot on screen?        
if args.plot == 2: vp.plot_roots(rs, "type")

#--- check RootAges
warns = checkRuns.checkRootAge(rs, warns)
warns_fn = Path(od).joinpath("milena/PyPlantBox/lintul5/warnings.out")
if len(warns) > 1:
    np.savetxt(fname=warns_fn.absolute(), X=np.array(warns), fmt="%s")
else:
    if warns_fn.exists():
        warns_fn.unlink()

sim.closeProject()
sim.shutDown()