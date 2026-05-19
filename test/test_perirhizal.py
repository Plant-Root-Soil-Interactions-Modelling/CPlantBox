
# this is a test file by Erik Kopp to test the alternative implementation of the perirhizal resistances and the perirhizal diffusion


import sys; sys.path.append("../src/functional")
sys.path.append("../../dumux-rosi/build-cmake/cpp/python_binding/")
sys.path.append("../../../dumuxtest/dumux/dumux-rosi/build-cmake/cpp/python_binding/")
#from plantbox import Perirhizal
from plantbox.functional.Perirhizal import PerirhizalPython
from numpy import linalg as LA
import plantbox.functional.van_genuchten as vg
from rosi.richards_flat import RichardsFlatWrapper as RichardsWrapper  # Python part, macroscopic soil model
#from richards import RichardsWrapper  # Python part
from rosi.richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial
from rosi.rosi_richardsnc_cyl import RichardsNCCylFoam # C++ part (Dumux binding), macroscopic soil model
#from rosi_richardsnc_cyl import RichardsNCCylFoam as RichardsNC_cyl  # C++ part (Dumux binding)
from rosi_richards22c import RichardsNCSPILU as RichardsNCSP #test

import Perirhizal
import pandas as pd
import numpy as np
import time
import math
import matplotlib.pyplot as plt
from plantbox.visualisation import figure_style

# run the dumux implementation of root water and nitrate uptake an then compare it to the alpha omega model

n_tests = 1 #try everything here for this many random parameter sets
do_computation = True #should the computation be run or take the data from a saved file

# general parameters

max_time = 3 #d
n_times = 300 # number of time intervals
times = np.linspace(0,max_time,n_times)[1:]
r_prhiz = 0.3 # cm
r_root = 0.02 # cm


rho = r_prhiz / r_root
NC = 101 # number of spatial discretisations
n_sp = NC - 1

length = 1 #cm?

dt = max_time / n_times

#space for the oupputs
watercontent_dumux = np.zeros((n_tests, n_times, n_sp+1))
soluteconcentration_dumux = np.zeros((n_tests, n_times, n_sp+1))
soluteconcentration_steadystate = np.zeros((n_tests, n_times, n_sp+1))
watercontent_steadyrate = np.zeros((n_tests, n_times, n_sp+1))
soluteconcentration_steadyrate = np.zeros((n_tests, n_times, n_sp+1))

soilVG = [0.078, 0.43, 0.036, 1.56, 24.96*10]  # hydrus loam soil #test:conductivity times 10 TODO remove test

lb = 0.5
points = np.logspace(np.log(r_root) / np.log(lb), np.log(r_prhiz) / np.log(lb),
                                  NC, base = lb) 
CC = np.array([(points[i] + points[i+1])/2 for i in range(len(points)-1)])

def load_constants(cyl):

    cyl.molarMassC = 12.011
    cyl.molarMassN = 14.0
    cyl.yr_per_d = 1/365 # [yr/d]
    cyl.m3_per_cm3 = 1e-6; # m3/cm3
    cyl.cm3_per_m3 = 1e6; # cm3/m3
    cyl.cm3_per_L = 1e3; # cm3/dm3

    cyl.kgC_per_mol = (1/1000) * cyl.molarMassC
    cyl.kgN_per_mol = (1/1000) * cyl.molarMassN
    
    yr_per_d = 1/365 # [yr/d]
    m3_per_cm3 = 1e-6; # m3/cm3
    cm3_per_m3 = 1e6; # cm3/m3

    molMarssW = 18 # g/mol
    rhoW = 1 #g/cm3
    rhoWM = rhoW/molMarssW #[g/cm3]*[g/mol] = mol/cm3
    cyl.rhoWM = rhoWM
    mlFr=1/rhoWM

    cyl.molarMassC = 12.011
    cyl.molarMassN = 14.0
    cyl.mg_per_molC = cyl.molarMassC * 1000.
    cyl.mg_per_molN = cyl.molarMassN * 1000.
    yr_per_d = 1/365 # [yr/d]
    m3_per_cm3 = 1e-6; # m3/cm3
    cm3_per_m3 = 1e6; # cm3/m3

    molMarssW = 18 # g/mol
    molMassMulC = 30 #g/mol
    rhoW = 1 #g/cm3
    rhoWM = rhoW/molMarssW #[g/cm3]*[g/mol] = mol/cm3
    cyl.mlFr=1/rhoWM #only water phase is relevant right now




    return cyl

def getSoilTextureAndShape(s, soilVG):
    """ soil shape and texture data
        to adapt according to the soil represented
    """
    
    constants = {'solidDensity': 2650, # kg/m3
    'solidMolarMass': 60.08e-3 #kg/mol 
    }
    paramSet = pd.DataFrame()
    
    solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
    solidMolarMass = 60.08e-3 # [kg/mol]
    # theta_r, theta_s, alpha, n, Ks
    #soilVG = [paramSet['soilVG_theta_r'],paramSet['soilVG_theta_s'],paramSet['soilVG_alpha'],paramSet['soilVG_n'],paramSet['soilVG_Ks']]
    soilTextureAndShape = {
                           "solidDensity":solidDensity,
                        'solidMolarMass': solidMolarMass,
                           'soilVG':soilVG}
    #common transformations
    
    soilTexture =soilTextureAndShape
    
    s.solidDensity = soilTexture['solidDensity'] #[kg/m^3 solid]
    s.solidMolarMass = soilTexture['solidMolarMass']# [kg/mol]
    s.soil =  soilTexture['soilVG']
    
    s.vg_soil = vg.Parameters(s.soil)
    # [mol / m3 solid] =[kg/m^3 solid] / [kg/mol]
    s.solidMolDensity = s.solidDensity/s.solidMolarMass
    # [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
    s.bulkDensity_m3 = s.solidMolDensity*(1.- s.vg_soil.theta_S)
    # (kg/m3) * (volume soil / total volume) * (g/kg) * (m3/cm3)
    s.bulkMassDensity_gpercm3 = s.solidDensity*(1.- s.vg_soil.theta_S)*1000*1e-6
    
    s.masspercm3 = s.solidDensity * 1.e-6
    return soilTexture

def run_perirhizal_test():
    
    #space for the outputs
    water_dumux = np.zeros((n_times, n_sp+1))
    solute_dumux = np.zeros((n_times, n_sp+1))
    water_ss = np.zeros((n_times, n_sp+1))
    solute_ss = np.zeros((n_times, n_sp+1))
    water_sr = np.zeros((n_times, n_sp+1))
    solute_sr = np.zeros((n_times, n_sp+1))
    
    mean_water = np.zeros((n_times))
    mean_solutes = np.zeros((n_times))
    
    # determine some random parameters
    initial_waterpotential = -3000 + np.random.rand() * 50 #cm3/cm3 #or choose an initial pressure head?
    initial_soluteconcentration = 2e-5*(1.0+np.random.rand()) #g/cm3 #TODO: lookup realistic concnetration, maybe 10 times the Michaelis Menten half saturation?
    print("initial_soluteconcentration",initial_soluteconcentration)
    
    # root conductivity and solute uptake parameters, it is chosen to be constant throughout the entire simulation time
    root_conductivity = 1e-4 #1/d
    inner_kr = root_conductivity * r_root * 2 * 3.14
    Vmax = 4.0e-11 * 62 * 1 * (2*3.14*r_root) * (24*3600) #mol/(cm2 s) * (g/mol) * cm * cm * (s/d) -> g / d #TODO: remove test multiplication by 10
    Km = 1.5e-7 * 62   #mol/cm3 -> g/cm3
    
    #DS_W = 1.902e-5 #cm2/s
    #Ds = DS_W / 10000#m2/s
    Ds = 1.902e-5 #cm2/s
    Ds = Ds * 24 * 3600 #cm2/d #TODO: remove the division by 3 after testing

    # the xylem matrix potential varies over time (keep it low so that there is little to no outflow of water)
    rx_t = lambda t : -7000+200*np.sin(t) #cm

    # load the perirhizal model
    peri = PerirhizalPython() 
    sp = vg.Parameters(soilVG)
    peri.set_soil(sp)
    #no lookup tables are used here as there arent many simulations

    

    simtimes = np.linspace(max_time/n_times,max_time,n_times).tolist() # days #TODO remove

    # initialise the dumux model
    s = RichardsWrapper(RichardsNCCylFoam())

    s.initialize()
    s.createGrid1d(np.linspace(r_root, r_prhiz, NC), length = length/100)  # [cm]
    s.setVGParameters([soilVG])

    # theta = 0.378, benchmark is set be nearly fully saturated, so we don't care too much about the specific values
    s.setHomogeneousIC(initial_waterpotential)  # cm pressure head

    s.setTopBC("constantFluxCyl",0.0)  #  [cm/day] "noFlux")#
    s.setBotBC("constantFluxCyl",0.0) # "noFlux")#
    s.setParameter("Soil.BC.Top.SType", "3")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Top.CValue", "0.0")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Bot.SType", "3")  # michaelisMenten=8 (SType = Solute Type)
    s.setParameter("Soil.BC.Bot.CValue", "0.0")

    s.setParameter("Soil.IC.C", str(initial_soluteconcentration))  # g / cm3  # TODO specialised setter?

    #s.setParameter("Component.MolarMass", "1.8e-2")
    s.setParameter("Component.MolarMass", "62.0")    

    s.setParameter("Component.LiquidDiffusionCoefficient", str(Ds / 1.e4 / (24*3600)))  # m2 s-1
    s.initializeProblem(maxDt = 0.01)

    cellVolumes = s.getCellSurfacesCyl() * length # cm3
    print("cellVolumes",cellVolumes)
    print(s)

    s.ddt = 1.e-4  # days

    simtimes.insert(0, 0)
    dt_ = np.diff(simtimes)
    
    for r, dt in enumerate(dt_):
        

        time = simtimes[r] 
        print('time',time)
            
        #if time >= 5:
        #    s.setSoluteTopBC([1], [0.])
        print("*****", "#", r, "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

        Wvolbefore = cellVolumes * s.getWaterContent() # cm3
        Smassbefore = s.getSolution(1) * Wvolbefore # g
        
        #model root water uptake 
        rx = rx_t(r*dt)
        current_rs_potential = s.getSolutionHead()
        water_dumux[r,1:] = vg.water_content(current_rs_potential,peri.sp)
        current_rs_potential = current_rs_potential[0]
        root_wateruptake = inner_kr * (rx - current_rs_potential)
        water_dumux[r,0] = root_wateruptake
        print(water_dumux)
        #s.setParameter( "Soil.BC.Bot.Value", str(root_wateruptake))
        
        
        #model root solute uptake
        #note: solutes are given in g
        current_rs_concentration = s.getSolution(1)
        solute_dumux[r,1:] = current_rs_concentration
        current_rs_concentration = current_rs_concentration[0]
        root_soluteuptake = - Vmax * current_rs_concentration / (Km + current_rs_concentration)
        solute_dumux[r,0] = root_soluteuptake
        print(solute_dumux)
        #s.setParameter( "Soil.BC.Bot.SValue", str(root_soluteuptake))
        #s.initializeProblem(maxDt = 0.01)
        
        s.setSource({0: root_wateruptake * length})
        s.setSource({0: root_soluteuptake * length}, eq_idx = 1)
        
        s.solve(dt, saveInnerFluxes_ = True)
        Wvolafter = cellVolumes*s.getWaterContent() # cm3
        Smassafter = s.getSolution(1) * Wvolafter # g
        
            
        rootSoilFluxes = s.getInnerFlow(0, length) * dt # cm3
        rootSoilFluxesS = s.getInnerFlow(1, length) * dt # g
        soilSoilFluxes = s.getOuterFlow(0, length) * dt # cm3
        soilSoilFluxesS = s.getOuterFlow(1, length) * dt # g
        
        # TODO: currently, setSource not properly implemented for richards and richards 2c
        # so left out of the mass balance.
        # scvSources = s.getSource(0) * cellVolumes * dt # cm3
        # scvSourcesS = s.getSource(1) * cellVolumes * dt # kg
        
    # #cyl.initializeProblem(maxDt=cyl.maxDt )
    # #boundary conditions
    # #water
    # cyl.setParameter( "Soil.BC.Bot.Type", str(int(3)))
    # cyl.setParameter( "Soil.BC.Top.Type", str(int(3)))
    # cyl.setParameter( "Soil.BC.Bot.Value", str(0.0)) #will be prescribed at each timestep
    # cyl.setParameter( "Soil.BC.Top.Value", str(0.0))
    
    # #solutes
    # s.setParameter( "Soil.BC.Bot.SType", str(int(3))) # TODO change this to 8 as Michaelis mnetne
    # s.setParameter( "Soil.BC.Top.SType", str(int(3)))
    # s.setParameter( "Soil.BC.Bot.CValue", str(0.0)) #Michaelis Menten parameters are set elsewhere
    # s.setParameter( "Soil.BC.Top.CValue", str(0.0))
    # print("test")
    
    # # run the dumux model
    # #get water and solute next to the root
    # current_rs_potential = s.getSolutionHead()[0] #todo rmeove
    # print(current_rs_potential)
    # current_rs_concentration = s.getSolution(1)[0] * s.rhoWM # todo remove
    # print(current_rs_concentration[0])
    # for i in range(n_times):
        # rx = rx_t(i*dt)
        # #model root water uptake 
        # root_wateruptake = root_conductivity * (rx - current_rs_potential)
        # water_dumux[i,0] = root_wateruptake
        # #s.setParameter( "Soil.BC.Bot.Value", str(root_wateruptake))
        # #s.setSource( 0, str(root_wateruptake)) #TODO activat ehtis
        # print("test")
        # #model root solute uptake 
        # root_soluteuptake = Vmax * current_rs_concentration * (Km + current_rs_concentration)
        # solute_dumux[i,0] = root_wateruptake
        # #s.setParameter( "Soil.BC.Bot.Value", str(root_wateruptake)) #todo set BC solute
        # print("test5")
        # s.solve(dt)
        
        # # save outputs of dumux
        # print("test6")
        # water_dumux[i,1:] = cyl.getWaterContent()
        # solute_dumux[i,1:] = cyl.getSolution(0)
        
        # mean_water = np.mean(water_dumux[i,1:])
        # mean_solutes = np.mean(solutes_dumux[i,1:])
        
        
        # current_rs_potential = cyl.getSolutionHead()
        # current_rs_potential = current_rs_potential[0]
        # current_rs_concentration = solute_dumux[i,1]

    

    
        # run the alpah omega model on waterflow (steady rate)
        # mean water content of the rhizosphere
        total_water = 0
        total_solute = 0
        for j in range(NC-1):
            total_water = total_water + water_dumux[r,j+1]*(points[j+1]**2 - points[j]**2)
            total_solute = total_solute + water_dumux[r,j+1]*solute_dumux[r,j+1]*(points[j+1]**2 - points[j]**2)
        mean_water = total_water / (points[NC-1]**2 - points[0]**2)
        total_solute = total_solute / (points[NC-1]**2 - points[0]**2)
        mean_solute = total_solute / mean_water
        
        #translate the mean water content to a mean matrix potential
        sx = vg.pressure_head(mean_water, peri.sp)
        
        #start the steady rate solver
        h_sr = peri.soil_root_interface_(rx, sx, inner_kr, rho, peri.sp)
        
        #compute both matrix flux potentials:
        Phi_root = np.array(vg.fast_mfp[sp](h_sr))
        Phi_soil = np.array(vg.fast_mfp[sp](sx))
        
        #compute the spatial watercontents
        Phi_A, Phi_C = peri.determine_mfp_function(Phi_root, Phi_soil, rho)
        #outer Phi
        Phi_out = Phi_A+Phi_C
        
        for j in range(NC-1):
            r_rel = CC[j] / r_prhiz
            Phi = Phi_A*(r_rel**2-np.log(r_rel**2))+Phi_C
            water_sr[r,j+1] = vg.water_content(vg.fast_imfp[sp](Phi),peri.sp)

        # run the alpha omega model on solute flow (both steady state and steady rate)
        waterflow = root_wateruptake
        result_solutes_ss = peri.soil_root_solutes_ss_([Phi_root], [Phi_out], [mean_solute], [Vmax], [Km], Ds, [waterflow], peri.sp)
        result_solutes_sr = peri.soil_root_solutes_sr_([Phi_root], [Phi_soil], [rho], [mean_solute], [Vmax], [Km], Ds, [waterflow], peri.sp)
        result_solutes_ss = result_solutes_ss[0]
        result_solutes_sr = result_solutes_sr[0]
        print("steadystate", result_solutes_ss, "steadyrate", result_solutes_sr)
        
        solute_ss[r,1] = result_solutes_ss
        solute_ss[r,0] = Vmax * result_solutes_ss / (Km + result_solutes_ss)
        solute_sr[r,1] = result_solutes_sr
        solute_sr[r,0] = Vmax * result_solutes_sr / (Km + result_solutes_sr)
        
        #F0 = peri.lookup_table_solutes((Phi_root,0)) # for the ratio of concentration next to the root to somewhere in the perirhizal zone
        F0 = peri.integral_overDiffusion_(Phi_root,peri.sp)
        D_tilde = 1/Ds/math.pow(sp.theta_S-sp.theta_R,13/3)*(sp.theta_S*sp.theta_S)
        for j in range(NC-1):
            
            r_rel = CC[j] / r_prhiz
            Phi_current = Phi_A*(r_rel**2-np.log(r_rel**2))+Phi_C
            #F = peri.lookup_table_solutes((Phi_current,0))-F0
            F = peri.integral_overDiffusion_(Phi_current,peri.sp)-F0
            #print("Ds", Ds, "Dtilde",D_tilde,"F",F,"Dtilde*F",D_tilde*F)
            F_tilde=math.exp(D_tilde*F)
            solute_ss[r,j+1] = result_solutes_ss * F_tilde + (1-F_tilde) * solute_ss[r,0] / waterflow
            solute_sr[r,j+1] = result_solutes_sr * F_tilde + (1-F_tilde) * solute_sr[r,0] / waterflow
    return water_dumux, solute_dumux, solute_ss, water_sr, solute_sr

if do_computation:
    # save everything in the np arrays
    for i in range(n_tests):
        water_dumux, solute_dumux, solute_ss, water_sr, solute_sr = run_perirhizal_test()
        watercontent_dumux[i,:,:]=water_dumux[:,:]
        soluteconcentration_dumux[i,:,:]=solute_dumux[:,:]
        soluteconcentration_steadystate[i,:,:]=solute_ss[:,:]
        watercontent_steadyrate[i,:,:]=water_sr[:,:]
        soluteconcentration_steadyrate[i,:,:]=solute_sr[:,:]
    
    np.savez("test_perirhizal.npz", watercontent_dumux=watercontent_dumux, soluteconcentration_dumux=soluteconcentration_dumux, soluteconcentration_steadystate=soluteconcentration_steadystate, watercontent_steadyrate=watercontent_steadyrate, soluteconcentration_steadyrate=soluteconcentration_steadyrate)

else:
    simulation_results = np.load("test_perirhizal.npz")
    watercontent_dumux = simulation_results["watercontent_dumux"]
    soluteconcentration_dumux = simulation_results["soluteconcentration_dumux"]
    soluteconcentration_steadystate = simulation_results["soluteconcentration_steadystate"]
    watercontent_steadyrate = simulation_results["watercontent_steadyrate"]
    soluteconcentration_steadyrate = simulation_results["soluteconcentration_steadyrate"]


# compare both for the differint means of water / solute content
run = 0
timestep = [1,3,5,7,9]
for i in range(5):
    timestep[i] = int(n_times * timestep[i] / 10)

#plot water and nitrogen in that one voxel
#fig, ax1 = figure_style.subplots12(1, 5)
fig, ax1 = figure_style.subplots12(nrows=5, ncols=1)





linestyle_dumux = "solid"
linestyle_steadystate = "dotted"
linestyle_steadyrate = "dashed"
    
#print(water_dumux)    
#print(CC)
for i in range(5):
    ax2 = ax1[i].twinx()
    
    water_dumux = watercontent_dumux[run, timestep[i], 1:]
    water_perirhizal = watercontent_steadyrate[run, timestep[i], 1:]
    solute_dumux = soluteconcentration_dumux[run, timestep[i], 1:]
    solute_steadystate = soluteconcentration_steadystate[run, timestep[i], 1:]
    solute_steadyrate = soluteconcentration_steadyrate[run, timestep[i], 1:]

    ax1[i].plot(CC, water_dumux, "b", linestyle = linestyle_dumux, label = "water_dumux")
    ax1[i].plot(CC, water_perirhizal, "b", linestyle = linestyle_steadyrate, label = "water_perirhizal")
    ax2.plot(CC, solute_dumux, "m", linestyle = linestyle_dumux, label = "solute_dumux")
    ax2.plot(CC, solute_steadystate, "m", linestyle = linestyle_steadystate, label = "solute_steadystate")
    ax2.plot(CC, solute_steadyrate, "m", linestyle = linestyle_steadyrate, label = "solute_steadyrate")

    ax1[i].set_xlabel("distance root [cm]")
    ax1[i].set_ylabel("water")
    ax2.set_ylabel("nitrogen")
    ax1[i].legend(["watercontent cm3/cm3"], loc="upper left")
    ax2.legend(["nitrogen concentration mol/cm3"], loc="upper right")

    ax1[i].legend(loc="upper left")
    ax2.legend(loc="upper right")
#np.save("input/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_)))) 
plt.show()

# for i in range(5):
    # for j in range(2):
        # ax2 = ax1[i,j].twinx()
    
        # water_dumux = watercontent_dumux[run, timestep[i], 1:]
        # water_perirhizal = watercontent_steadyrate[run, timestep[i], 1:]
        # solute_dumux = soluteconcentration_dumux[run, timestep[i], 1:]
        # solute_steadystate = soluteconcentration_steadystate[run, timestep[i], 1:]
        # solute_steadyrate = soluteconcentration_steadyrate[run, timestep[i], 1:]

        # ax1[i,j].plot(CC, water_dumux, "b", linestyle = linestyle_dumux, label = "water_dumux")
        # ax1[i,j].plot(CC, water_perirhizal, "b", linestyle = linestyle_steadyrate, label = "water_perirhizal")
        # ax2.plot(CC, solute_dumux, "m", linestyle = linestyle_dumux, label = "solute_dumux")
        # if j==0:
            # ax2.plot(CC, solute_steadystate, "m", linestyle = linestyle_steadystate, label = "solute_steadystate")
        # else:
            # ax2.plot(CC, solute_steadyrate, "m", linestyle = linestyle_steadyrate, label = "solute_steadyrate")

        # ax1[i,j].set_xlabel("distance root [cm]")
        # ax1[i,j].set_ylabel("water")
        # ax2.set_ylabel("nitrogen concentration")
        # ax1[i,j].legend(["watercontent cm3/cm3"], loc="upper left")
        # ax2.legend(["nitrogen concentration mol/cm3"], loc="upper right")

        # ax1[i,j].legend(loc="upper left")
        # ax2.legend(loc="upper right")
# #np.save("input/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_)))) 
# plt.show()


# #---old stuff below


# # generate random inputs for one perirhizal segment

# """
        # rx             xylem matric potential [cm]
        # sx             bulk soil matric potential [cm]
        # inner_kr       root radius times hydraulic conductivity [cm/day] 
        # rho            geometry factor [1] (outer_radius / inner_radius)
        # c_sol          concentration of solutes in the bulk soil
        # Vmax           Michaelis Menten Kinetics maximal solute uptake rate
        # Km             Michaelis Menten Kinetics half saturation
        # sp             soil parameter: van Genuchten parameter set (type vg.Parameters)
        # sp_theta_r
        # sp_theta_s
        # sp_n
        # sp_alpha
        # sp_Ksat
# """
# labels = ["rx","sx","inner_kr","rho","c_sol","Vmax","Km","hsr","hsr_lookup","hsr_global","waterflow","sol_c","sol_c_ss","sol_c_sr","sol_U","sol_U_ss","sol_U_sr"]

# inputs = ["rx","sx","inner_kr","rho","c_sol","Vmax","Km"]
# outputs = ["hsr","hsr_lookup","hsr_global","sol_c","sol_c_ss","sol_c_sr","sol_U","sol_U_ss","sol_U_sr"]
# outputs = ["rx","sx","waterflow","sol_c","sol_c_ss","sol_c_sr","sol_U","sol_U_ss","sol_U_sr"]

# #kr is assumed to lie between .5 and 5e-6 cm/hPa/d = 0.5 to 5e-6 1/d 
# #r is assumed to be 0.03
# #multiply this number by 10 in order to not fall under the threshhold at which the soil hydraulic potential is chosen as a default

# Intervals = {
        # "rx": [-14000,-10],
        # "sx": [-1000,-10],
        # "inner_kr": [1.5e-7,1.5e-6], 
        # "rho": [10,199.0],
        # "c_sol": [1.0e-9,1.0e-5],
        # "Vmax": [1.0e-11,1.0e-10],
        # "Km": [1.0e-8,1.0e-6],
        # }
# Intervals = pd.DataFrame(Intervals, index = ["min","max"])

# ntests = 10 
# tests = pd.DataFrame(index = range(ntests),columns = labels)

# for one_label in Intervals.columns:
    # min_val = Intervals.loc["min",one_label]
    # max_val = Intervals.loc["max",one_label]

    # testvalues = (max_val - min_val) * np.random.rand(ntests) + min_val * np.ones(ntests)
    # tests[one_label] = testvalues

# peri = PerirhizalPython()
# #peri = PerirhizalPython(Perirhizal)


# rx = np.array(tests.loc[:,"rx"]) #* (-1)
# sx = np.array(tests.loc[:,"sx"]) #* (-1)
# inner_kr = np.array(tests.loc[:,"inner_kr"])
# rho = np.array(tests.loc[:,"rho"])
# c_sol = np.array(tests.loc[:,"c_sol"])
# Vmax = np.array(tests.loc[:,"Vmax"])
# Km = np.array(tests.loc[:,"Km"])
# Ds = 1.0e1
# r_root = 0.03

# hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]
# sp = vg.Parameters(hydrus_loam)
# peri.set_soil(sp)

# #create the lookup table
# #t = time.time()
# #peri.create_lookup("results/hydrus_loam", sp)
# #elapsed_time = time.time() - t
# #print("Creating the lookup table took about ", elapsed_time) #

# peri.open_lookup("results/hydrus_loam")


# #create the global lookup table
# #t = time.time()
# #peri.create_lookup_global("results/"+peri.water_filename, sp)
# #elapsed_time = time.time() - t
# #print("Creating the global lookup table took about ", elapsed_time) #

# peri.open_global_lookup("results/"+peri.water_filename)

# # create lookup tables for the solute flow
# Ds_=np.logspace(np.log10(1.0e0), np.log10(1.0e1), 2)
# #peri.create_integralDiffusion_lookup("results/hydrus_loam_ss_solutes", sp)
# #peri.create_integralconcentration_lookup("results/hydrus_loam_sr_solutes", Ds_, sp)

# # open lookup tables for the solute flow
# peri.open_lookup_solutes("results/hydrus_loam_ss_solutes")
# peri.open_lookup_sr_solutes("results/hydrus_loam_sr_solutes")


# #numerically solve each root segment, this is the reference solution 
# hsr = np.array([peri.soil_root_interface_(rx[i], sx[i], inner_kr[i], rho[i], peri.sp) for i in range(0, len(rx))])
# tests.loc[:,"hsr"] = hsr
# print("norm of the root soil matrix potentials: ", LA.norm(hsr))

# #the standard implementation
# hsr1= peri.soil_root_interface_potentials_table(rx, sx, inner_kr, rho)
# tests.loc[:,"hsr_lookup"] = hsr1
# print("Norm of the difference to basic lookup table:", LA.norm(hsr-hsr1))

# #the alternative global implementation
# hsr3 = peri.soil_root_interface_potentials_table_global(rx, sx, inner_kr, rho)
# tests.loc[:,"hsr_global"] = hsr3
# print("Norm of the difference to global lookup table:", LA.norm(hsr-hsr3))

# #waterflow = 2*3.14*np.multiply((hsr-rx),inner_kr)/r_root
# waterflow = 2*3.14*np.multiply((sx-rx),inner_kr)#/r_root
# Phi_root = np.array([vg.fast_mfp[sp](hsr[i]) for i in range(len(hsr))])
# Phi_soil = np.array([vg.fast_mfp[sp](sx[i]) for i in range(len(sx))])
# tests.loc[:,"waterflow"] = waterflow
# Ds = 1.0e1

# #base solutes
# tests.loc[:,"sol_c"] = c_sol
# tests.loc[:,"sol_U"] = np.array([Vmax[i]*c_sol[i]/(Km[i]+c_sol[i])*1e6 for i in range(len(c_sol))])

# # steady state solutes
# solutes_ss = peri.soil_root_solutes_ss_(Phi_root, Phi_soil, c_sol, Vmax, Km, Ds, waterflow, peri.sp)
# tests.loc[:,"sol_c_ss"] = solutes_ss
# tests.loc[:,"sol_U_ss"] = np.array([Vmax[i]*solutes_ss[i]/(Km[i]+solutes_ss[i])*1e6 for i in range(len(solutes_ss))])

# # steady rate solutes
# solutes_sr = peri.soil_root_solutes_sr_(Phi_root, Phi_soil, rho, c_sol, Vmax, Km, Ds, waterflow, peri.sp)
# tests.loc[:,"sol_c_sr"] = solutes_sr
# tests.loc[:,"sol_U_sr"] = np.array([Vmax[i]*solutes_sr[i]/(Km[i]+solutes_sr[i])*1e6 for i in range(len(solutes_sr))])




# print("All inputs:")
# print(tests[inputs])
# print("All outputs:")
# print(tests[outputs])




