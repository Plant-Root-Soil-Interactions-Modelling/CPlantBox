
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
#from rosi_richards22c import RichardsNCSPILU as RichardsNCSP #test

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

max_time = 40 #d
n_times = 400 # number of time intervals
times = np.linspace(0,max_time,n_times)[1:]
r_prhiz = 0.6 # cm
r_root = 0.02 # cm


rho = r_prhiz / r_root
NC = 41 # number of spatial discretisations
n_sp = NC - 1

length = 1 #cm?

dt = max_time / n_times

#space for the oupputs
watercontent_dumux = np.zeros((n_tests, n_times, n_sp+1))
soluteconcentration_dumux = np.zeros((n_tests, n_times, n_sp+1))
soluteconcentration_steadystate = np.zeros((n_tests, n_times, n_sp+1))
watercontent_steadyrate = np.zeros((n_tests, n_times, n_sp+1))
watercontent_steadyrate2 = np.zeros((n_tests, n_times, n_sp+1))
soluteconcentration_steadyrate = np.zeros((n_tests, n_times, n_sp+1))

matrixpotential_perirhizal2 = np.zeros((n_tests, n_times, n_sp+1))

soilVG = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam soil 

soilVG = [0.08, 0.43, 0.04, 1.6, 50] #loam from benchmark #da war am Ende noch eine 0.5, keine Ahnung warum, vermutlich eine vorherige Programmversion

lb = 0.5
points = np.logspace(np.log(r_root) / np.log(lb), np.log(r_prhiz) / np.log(lb),
                                  NC, base = lb) 
CC = np.array([(points[i] + points[i+1])/2 for i in range(len(points)-1)])



def run_perirhizal_test():
    
    #space for the outputs
    water_dumux = np.zeros((n_times, n_sp+1))
    water2_dumux = np.zeros((n_times, n_sp+1))
    solute_dumux = np.zeros((n_times, n_sp+1))
    water_ss = np.zeros((n_times, n_sp+1))
    solute_ss = np.zeros((n_times, n_sp+1))
    water_sr = np.zeros((n_times, n_sp+1))
    water_test = np.zeros((n_times, n_sp+1))
    solute_sr = np.zeros((n_times, n_sp+1))
    
    mp_perirhizal2 = np.zeros((n_times, n_sp+1))
    
    mean_water = np.zeros((n_times))
    mean_solutes = np.zeros((n_times))
    
    # determine some random parameters
    initial_waterpotential = -100 #+ np.random.rand() * 50 #cm3/cm3 #or choose an initial pressure head?
    initial_soluteconcentration = 2e-5*(1.0+np.random.rand()) #g/cm3 #TODO: lookup realistic concnetration, maybe 10 times the Michaelis Menten half saturation?
    print("initial_soluteconcentration",initial_soluteconcentration)
    
    # root conductivity and solute uptake parameters, it is chosen to be constant throughout the entire simulation time
    root_conductivity = 1e-4 #1/d
    inner_kr = root_conductivity * r_root * 2 * 3.14
    waterdemand = -0.1 #cm/d
    Vmax = 4.0e-11 * 62 * 1 * (2*3.14*r_root) * (24*3600)  #mol/(cm2 s) * (g/mol) * cm * cm * (s/d) -> g / d 
    Km = 1.5e-7 * 62   #mol/cm3 -> g/cm3
    
    #DS_W = 1.902e-5 #cm2/s
    #Ds = DS_W / 10000#m2/s
    Ds = 1.902e-5 #cm2/s
    Ds = Ds * 24 * 3600 #cm2/d 

    # the xylem matrix potential varies over time (keep it low so that there is little to no outflow of water)
    rx_t = lambda t : -1000+200*np.sin(t) #cm
    rx_t = lambda t : -14000+0*np.sin(t) #cm
    rx_t = lambda t : -1000-200*t #cm

    # load the perirhizal model
    peri = PerirhizalPython() 
    sp = vg.Parameters(soilVG)
    peri.set_soil(sp)
    #no lookup tables are used here as there arent many simulations

    

    simtimes = np.linspace(max_time/n_times,max_time,n_times).tolist() # days #TODO remove

    # initialise the dumux model
    s = RichardsWrapper(RichardsNCCylFoam())

    s.initialize()
    s.createGrid1d(np.linspace(r_root, r_prhiz, NC), length = length/100)  # [m] -> [cm]
    s.setVGParameters([soilVG])

    # theta = 0.378, benchmark is set be nearly fully saturated, so we don't care too much about the specific values
    s.setHomogeneousIC(initial_waterpotential)  # cm pressure head

    s.setTopBC("constantFluxCyl",0.0)  #  [cm/day] "noFlux")#
    #s.setBotBC("constantFluxCyl",0.0) # "noFlux")#
    s.setBotBC("constantFluxCyl",waterdemand) # "noFlux")# Flux in cm/d
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
        
        #alternative: adapt root matrix potential for a steady uptake:
        rx = max(waterdemand / root_conductivity + current_rs_potential, -14750) #note: waterdemand is assumed to be negative
        
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
        
        #s.setSource({0: root_wateruptake * length})
        s.setSource({0: root_soluteuptake * length}, eq_idx = 1)
        
        s.solve(dt, saveInnerFluxes_ = True)
        Wvolafter = cellVolumes*s.getWaterContent() # cm3
        Smassafter = s.getSolution(1) * Wvolafter # g
        
            
        rootSoilFluxes = s.getInnerFlow(0, length) * dt # cm3
        rootSoilFluxesS = s.getInnerFlow(1, length) * dt # g
        soilSoilFluxes = s.getOuterFlow(0, length) * dt # cm3
        soilSoilFluxesS = s.getOuterFlow(1, length) * dt # g
        

    

    
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
        
        #compute Phi_A based on the water demand
        #Phi_A = waterdemand * r_root / (2*(1-rho**2))
        #Phi_C = 
        #compute Phi_C based on the result
        
        #outer Phi
        Phi_out = Phi_A+Phi_C
        print("Phiroot", Phi_root, "Phi_soil", Phi_soil, "Phi_A", Phi_A, "Phi_C", Phi_C)
        r_rel = CC[0] / r_prhiz
        print("Phi_out", Phi_out, "Phi_in", Phi_A*(r_rel**2-np.log(r_rel**2))+Phi_C)
        
        for j in range(NC-1):
            r_rel = CC[j] / r_prhiz
            Phi = Phi_A*(r_rel**2-np.log(r_rel**2))+Phi_C
            #Phi_soil = Phi_out # TODO: remove this
            #if water_dumux[r,j+1]<soilVG[0]+0.02: # I model this as a threshhold for water stress
            #    Phi = Phi_soil+abs(waterdemand)*r_root*((r_rel*rho)**2/(2*(1-rho**2))-rho**2/(1-rho**2)*(np.log(r_rel**2)+0.5))
            #else:
            #    Phi = Phi_soil*r_root*((r_rel*rho)**2-1+2*rho**2*np.log(r_root/r_rel))/(rho**2-1+2*rho**np.log(1/(rho)))
            if Phi<=0: #water stress
                water_sr[r,j+1] = soilVG[0]+0.02
            else:
                water_sr[r,j+1] = vg.water_content(vg.fast_imfp[sp](Phi),peri.sp)
            Phi = Phi_soil+abs(waterdemand)*r_root*((r_rel*rho)**2/(2*(1-rho**2))+rho**2/(1-rho**2)*(np.log(1/r_rel)-0.5)) #waterdemand is usually negative
            Phi = Phi_out+abs(waterdemand)*r_root*((r_rel*rho)**2/(2*(1-rho**2))+rho**2/(1-rho**2)*(np.log(1/r_rel)-0.5)) #waterdemand is usually negative
            
            mp_perirhizal2[r,j+1]=vg.fast_imfp[sp](Phi)
            #alternative steady rate approximation            
            if Phi<=0: #water stress
                water2_dumux[r,j+1] = soilVG[0]+0.001
            else:
                water2_dumux[r,j+1] = vg.water_content(vg.fast_imfp[sp](Phi),peri.sp)

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
            #Phi_current = waterdemand*r_root*((r_rel*rho)**2/(2*(1-rho**2))-rho**2/(1-rho**2)*(np.log(r_rel**2)+0.5))+Phi_soil
            #F = peri.lookup_table_solutes((Phi_current,0))-F0
            F = peri.integral_overDiffusion_(Phi_current,peri.sp)-F0
            #print("Ds", Ds, "Dtilde",D_tilde,"F",F,"Dtilde*F",D_tilde*F)
            F_tilde=math.exp(D_tilde*F)
            solute_ss[r,j+1] = result_solutes_ss * F_tilde + (1-F_tilde) * solute_ss[r,0] / waterflow
            solute_sr[r,j+1] = result_solutes_sr * F_tilde + (1-F_tilde) * solute_sr[r,0] / waterflow
    return water_dumux, solute_dumux, solute_ss, water_sr, solute_sr, water2_dumux, mp_perirhizal2

if do_computation:
    # save everything in the np arrays
    for i in range(n_tests):
        water_dumux, solute_dumux, solute_ss, water_sr, solute_sr, water2_dumux, mp_perirhizal2 = run_perirhizal_test()
        watercontent_dumux[i,:,:]=water_dumux[:,:]
        soluteconcentration_dumux[i,:,:]=solute_dumux[:,:]
        soluteconcentration_steadystate[i,:,:]=solute_ss[:,:]
        watercontent_steadyrate[i,:,:]=water_sr[:,:]
        watercontent_steadyrate2[i,:,:]=water2_dumux[:,:] #TODO this is confusing
        soluteconcentration_steadyrate[i,:,:]=solute_sr[:,:]
        matrixpotential_perirhizal2[i,:,:]=mp_perirhizal2[:,:]
    
    np.savez("test_perirhizal.npz", watercontent_dumux=watercontent_dumux, watercontent_steadyrate2=watercontent_steadyrate2, matrixpotential_perirhizal2=matrixpotential_perirhizal2, soluteconcentration_dumux=soluteconcentration_dumux, soluteconcentration_steadystate=soluteconcentration_steadystate, watercontent_steadyrate=watercontent_steadyrate, soluteconcentration_steadyrate=soluteconcentration_steadyrate)

else:
    simulation_results = np.load("test_perirhizal.npz")
    watercontent_dumux = simulation_results["watercontent_dumux"]
    soluteconcentration_dumux = simulation_results["soluteconcentration_dumux"]
    soluteconcentration_steadystate = simulation_results["soluteconcentration_steadystate"]
    watercontent_steadyrate = simulation_results["watercontent_steadyrate"]
    soluteconcentration_steadyrate = simulation_results["soluteconcentration_steadyrate"]
    watercontent_steadyrate2 = simulation_results["watercontent_steadyrate2"]
    matrixpotential_perirhizal2 = simulation_results["matrixpotential_perirhizal2"]
    


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
#timestep = [100,300,5,7,9]
for i in range(5):
    ax2 = ax1[i].twinx()
    print(watercontent_dumux)
    water_dumux = watercontent_dumux[run, timestep[i], 1:]
    water_perirhizal = watercontent_steadyrate[run, timestep[i], 1:]
    water_steadyrate2 = watercontent_steadyrate2[run, timestep[i], 1:]
    solute_dumux = soluteconcentration_dumux[run, timestep[i], 1:]
    solute_steadystate = soluteconcentration_steadystate[run, timestep[i], 1:]
    solute_steadyrate = soluteconcentration_steadyrate[run, timestep[i], 1:]

    ax1[i].plot(CC, water_dumux, "b", linestyle = linestyle_dumux, label = "water_dumux")
    #ax1[i].plot(CC, water_perirhizal, "b", linestyle = linestyle_steadyrate, label = "water_perirhizal")
    ax1[i].plot(CC, water_steadyrate2, "b", linestyle = linestyle_steadyrate, label = "water_perirhizal2")
    #ax2.plot(CC, solute_dumux, "m", linestyle = linestyle_dumux, label = "solute_dumux")
    #ax2.plot(CC, solute_steadystate, "m", linestyle = linestyle_steadystate, label = "solute_steadystate")
    #ax2.plot(CC, solute_steadyrate, "m", linestyle = linestyle_steadyrate, label = "solute_steadyrate")

    ax1[i].set_xlabel("distance root [cm]")
    ax1[i].set_ylabel("water")
    ax2.set_ylabel("nitrogen")
    ax1[i].legend(["watercontent cm3/cm3"], loc="upper left")
    ax2.legend(["nitrogen concentration mol/cm3"], loc="upper right")

    ax1[i].legend(loc="upper left")
    ax2.legend(loc="upper right")
#np.save("input/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_)))) 
plt.show()

fig, ax1 = figure_style.subplots12(nrows=5, ncols=1)



peri = PerirhizalPython() 
sp = vg.Parameters(soilVG)
peri.set_soil(sp)

linestyle_dumux = "solid"
linestyle_steadystate = "dotted"
linestyle_steadyrate = "dashed"
    
#print(water_dumux)    
#print(CC)
#timestep = [100,300,5,7,9]
for i in range(5):
    ax2 = ax1[i].twinx()
    
    water_dumux = watercontent_dumux[run, timestep[i], 1:]
    water_perirhizal = watercontent_steadyrate[run, timestep[i], 1:]
    water_steadyrate2 = watercontent_steadyrate2[run, timestep[i], 1:]
    solute_dumux = soluteconcentration_dumux[run, timestep[i], 1:]
    solute_steadystate = soluteconcentration_steadystate[run, timestep[i], 1:]
    solute_steadyrate = soluteconcentration_steadyrate[run, timestep[i], 1:]
    
    mp_steadyrate2 = matrixpotential_perirhizal2[run, timestep[i], 1:]
    
    
    water_dumux = vg.pressure_head(water_dumux, peri.sp)
    water_perirhizal = vg.pressure_head(water_perirhizal, peri.sp)
    #water_perirhizal2 = vg.pressure_head(water_steadyrate2, peri.sp)

    ax1[i].plot(CC, water_dumux, "b", linestyle = linestyle_dumux, label = "water_dumux")
    #ax1[i].plot(CC, water_perirhizal, "b", linestyle = linestyle_steadyrate, label = "water_perirhizal")
    ax1[i].plot(CC, mp_steadyrate2, "b", linestyle = linestyle_steadyrate, label = "water_perirhizal")
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


