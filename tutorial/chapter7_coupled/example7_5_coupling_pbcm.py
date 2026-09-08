"""coupling pbcm"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp 

# load simplace wrapper package
import simplace # https://simplace.net/doc/python_wrapper.htm

# load custom modules to exchange outputs between simplace-pb states
from simplace_cpb.PyPlantBox.misc.SimplaceRootbox.simplace import lintulslim_interact as sp_int

# Simulation configuration
gram_per_cm = 0.000035  # specific root length density used to calculate the maximum root increment in a timestep (g cm-1)
area = 6 * 12.5  # estimated area of a wheat plant for RLD scaling (cm2)
sim_time = 600  # maximum simulation time-n_steps (days)
dt = 1  # simulation dt (days)
plot = True  # plot outputs at the end

# Feddes parameters for root elongation restriction due to soil water potential
# Eq. 10 of https://doi.org/10.3389/fpls.2022.865188
h1 = 0
h2 = -1
h3 = -5
h4 = -150

# SIMPLACE initialization |\label{l7_5_simplace:InitStart}|
home_dir = str(Path("~").expanduser())
jd = home_dir + '/workspace/'
wd = jd + 'CPlantBox/tutorial/chapter7_coupled/simplace_cpb/'
od = jd + 'simplace_run/output/'
sol = wd + 'PyPlantBox/solution/LintulSlimRB.sol.xml' # simplace model solution file: https://simplace.net/doc/index.html?solution.htm

# instantiate simplace simulation
sim = simplace.SimplaceInstance(jd, wd, od, wd, wd)
sim.setLogLevel('ERROR')
sim.openProject(sol)

par = {'startdate':"20.03.1991"} # overwrite startdate in the model solution
sim.createSimulation(par)
ids = sim.getSimulationIDs()

# CPlantBox initialization
# read root parameters
root_parfile = wd + '/PyPlantBox/data/modelparameter/wheat.xml'

# Soil discretization matching the 1D-profile used by simplace input file
soildepth = 120
layers = 40

# Set scale elongaton according to tutorial:
# tutorial/chapter3_responses/example3_1_carbon.py
scale_elongation = pb.EquidistantGrid1D(0, -soildepth, layers)
scale_elongation.data = np.ones((layers))
se = pb.ProportionalElongation()
se.setBaseLookUp(scale_elongation)

# Initialize root system
plant = pb.Plant()
plant.setSeed(0)
plant.readParameters(root_parfile) 
for p in plant.getOrganRandomParameter(pb.root):
    p.f_se = se  # set scale elongation function

# intialize root system
plant.initialize()

# initialize length and rld states
ol = 0
vRLD=np.zeros(layers) # |\label{l7_5_simplace:InitEnd}|

# Simulation loop
for s in range(0,sim_time): # |\label{l7_5_simplace:LoopStart}|
    
    # simulate simplace step to get maxinc and re_reduction dinamically
    (date, maxinc, doharvest, tranrf, yld, rld_s, re_reduction) = sp_int.getSimplaceValuesExtended(sim,gram_per_cm, h1, h2, h3, h4, area=area) # |\label{l7_5_simplace:MaxInc_simplace}|

    # update scale elongation
    scale_elongation.data = np.array(re_reduction)

    maxinc = round(maxinc, 2)
    inc = 0.0
    print("Simulating: ", date, " MaxIncr:", round(maxinc, 2), " Tranrf:", round(tranrf, 2), " Yield:", round(yld, 2), " Step:", s)
    
    # run CPlantBox if there's any root biomass allocated from simplace
    if(maxinc > 0): # |\label{l7_5_simplace:RunCPB}|
        
        plant.simulate(dt,maxinc*dt, se, True) # not working: it runs but the code breaks when pb.SegmentAnalyser(rs)        
        
        # get simulated root length increment
        l = np.sum(plant.getParameter('length'))
        inc =  l - ol
        ol = l
        
        # calculate RLD for each soil layer using pb.SegmentAnalyser(rs)
        vRLD = sp_int.calculateRLD(plant, soildepth, layers, area) # |\label{l7_5_simplace:RLD_CPB}|
    
    # update vRLD in simplace
    sp_int.setSimplaceRoots(sim, vRLD, maxinc-inc, gram_per_cm) # |\label{l7_5_simplace:RLD_SIMPLACE}|

    if doharvest == True:
        # exit loop if harvest happens 
        break # |\label{l7_5_simplace:LoopEnd}|

plant.write(od + "PyPlantBox/lintul/LintulSlim_PlantBox.vtp")
if plot:

    # read and plot some simplace outputs
    simplace_out = pd.read_csv(od + "PyPlantBox/lintul/slim_rootbox.csv", sep=";") # |\label{l7_5_simplace:ReadOuts}|
    simplace_out["CURRENT.DATE"] = pd.to_datetime(simplace_out["CURRENT.DATE"], format="%d.%m.%Y")
    
    simplace_out.plot(x="CURRENT.DATE", y=["Yield", "AGBm"], ylabel = "Dry Biomass [g/m2]")
    simplace_out.plot(x="CURRENT.DATE", y=["EVAP", "TRAN", "LAI"], ylabel = "Transpiration/Evaporation [mm/day], LAI [m2/m2]")

    # plot final root architecture
    vp.plot_roots(plant, "type")
    plt.show()

sim.closeProject()
sim.shutDown()