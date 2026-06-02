
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

from scipy.optimize import fsolve, root_scalar

# run the dumux implementation of root water and nitrate uptake an then compare it to the alpha omega model

n_tests = 1 #try everything here for this many random parameter sets
do_computation = True #should the computation be run or take the data from a saved file

# general parameters

max_time = 1 #d
n_times = 50 # number of time intervals
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

matrixpotential_perirhizal1 = np.zeros((n_tests, n_times, n_sp+1))
matrixpotential_perirhizal2 = np.zeros((n_tests, n_times, n_sp+1))

soluteconcentration_Tiina = np.zeros((n_tests, n_times, n_sp+1))

soilVG = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam soil 

soilVG = [0.08, 0.43, 0.04, 1.6, 50] #loam from benchmark #da war am Ende noch eine 0.5, keine Ahnung warum, vermutlich eine vorherige Programmversion

lb = 0.5
points = np.logspace(np.log(r_root) / np.log(lb), np.log(r_prhiz) / np.log(lb),
                                  NC, base = lb) 
CC = np.array([(points[i] + points[i+1])/2 for i in range(len(points)-1)])

volumes = np.array([(points[i+1]**2 - points[i]**2)*3.14 for i in range(len(points)-1)])
initial_waterpotential = -300
initial_soluteconcentration = 2e-5#mol/cm3
outer_conc = initial_soluteconcentration

def run_perirhizal_test():
    
    #space for the outputs
    water_dumux = np.zeros((n_times, n_sp+1))
    water_sr2 = np.zeros((n_times, n_sp+1))
    solute_dumux = np.zeros((n_times, n_sp+1))
    water_ss = np.zeros((n_times, n_sp+1))
    solute_ss = np.zeros((n_times, n_sp+1))
    water_sr = np.zeros((n_times, n_sp+1))
    water_test = np.zeros((n_times, n_sp+1))
    solute_sr = np.zeros((n_times, n_sp+1))
    
    mp_perirhizal1 = np.zeros((n_times, n_sp+1))
    mp_perirhizal2 = np.zeros((n_times, n_sp+1))
    
    mean_water = np.zeros((n_times))
    mean_solutes = np.zeros((n_times))
    
    solutes_Tiina = np.zeros((n_times, n_sp+1))
    
    # determine some random parameters
    #initial_waterpotential = -700 #+ np.random.rand() * 50 #cm3/cm3 #or choose an initial pressure head?
    #initial_soluteconcentration = 2e-5#*(1.0+np.random.rand()) #g/cm3 #TODO: lookup realistic concnetration, maybe 10 times the Michaelis Menten half saturation?
    print("initial_soluteconcentration",initial_soluteconcentration)
    
    # root conductivity and solute uptake parameters, it is chosen to be constant throughout the entire simulation time
    root_conductivity = 1e-4 #1/d
    inner_kr = root_conductivity * r_root * 2 * 3.14
    waterdemand = -0.1 #cm/d
    radial_waterdemand = 2*3.14*r_root * waterdemand #cm2/d
    Vmax = 4.0e-11 * 1 * (2*3.14*r_root) * (24*3600)  #mol/(cm2 s) * cm * cm * (s/d) -> mol / d 
    Vmax_per_area = Vmax / (1 * (2*3.14*r_root)) #mol / d /cm2 = mol/(cm2d)
    Km = 1.5e-7   #mol/cm3
    
    #DS_W = 1.902e-5 #cm2/s
    #Ds = DS_W / 10000#m2/s
    Ds = 1.902e-5 / 2#cm2/s # division by 2: from NO3 to H2PO4
    Ds = Ds * 24 * 3600 #cm2/d 
    
    outer_waterpotential = initial_waterpotential #cm
    
    outer_kr = root_conductivity * r_prhiz * 2 * 3.14 #TODO: this is just for testing purposes as I do not want a Dirichlet BC as it would be mixed BC
    outer_conc = initial_soluteconcentration

    # the xylem matrix potential varies over time (keep it low so that there is little to no outflow of water)
    rx_t = lambda t : -1000+200*np.sin(t) #cm
    rx_t = lambda t : -14000+0*np.sin(t) #cm
    rx_t = lambda t : -1000-200*t #cm

    # load the perirhizal model
    peri = PerirhizalPython() 
    sp = vg.Parameters(soilVG)
    peri.set_soil(sp)
    #no lookup tables are used here as there arent many simulations

    outer_watercontent = vg.water_content(outer_waterpotential, peri.sp)
    

    simtimes = np.linspace(max_time/n_times,max_time,n_times).tolist() # days #TODO remove

    # initialise the dumux model
    s = RichardsWrapper(RichardsNCCylFoam())

    s.initialize()
    s.createGrid1d(np.linspace(r_root, r_prhiz, NC), length = length/100)  # [m] -> [cm]
    s.setVGParameters([soilVG])

    # theta = 0.378, benchmark is set be nearly fully saturated, so we don't care too much about the specific values
    s.setHomogeneousIC(initial_waterpotential)  # cm pressure head

    s.setTopBC("constantFluxCyl",0.0)  #  [cm/day] "noFlux")#
    #s.setTopBC("constantPressure",-200)  #  [cm/day] "noFlux")#
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
        
        dampening = 0.01 #slow down the inflow from outside

        time = simtimes[r] 
        print('time',time)
            
        #if time >= 5:
        #    s.setSoluteTopBC([1], [0.])
        print("*****", "#", r, "external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")

        Wvolbefore = cellVolumes * s.getWaterContent() # cm3
        Smassbefore = s.getSolution(1) * Wvolbefore # g
        
        #model root water uptake 
        rx = rx_t(r*dt)
        
        
        current_potential = s.getSolutionHead()
        water_dumux[r,1:] = vg.water_content(current_potential,peri.sp)
        current_rs_potential = current_potential[0]
        current_outer_potential = current_potential[-1]
        current_outer_watercontent = vg.water_content(current_outer_potential, peri.sp)
        root_wateruptake = inner_kr * (rx - current_rs_potential)
        outer_watersource = outer_kr * (outer_waterpotential - current_outer_potential)
        outer_watersource = (outer_watercontent - current_outer_watercontent) /dt * dampening
        water_dumux[r,0] = root_wateruptake
        print(water_dumux)
        #s.setParameter( "Soil.BC.Bot.Value", str(root_wateruptake))
        
        #alternative: adapt root matrix potential for a steady uptake:
        rx = max(waterdemand / root_conductivity + current_rs_potential, -14750) #note: waterdemand is assumed to be negative
        #rx = current_rs_potential
        
        #model root solute uptake
        #note: solutes are given in g
        current_concentration = s.getSolution(1)
        solute_dumux[r,1:] = current_concentration
        current_rs_concentration = current_concentration[0]
        current_outer_concentration = current_concentration[-1]
        root_soluteuptake = - Vmax * max(current_rs_concentration,0) / (Km + current_rs_concentration) # mol / d
        outer_solutesource = ( outer_conc - current_outer_concentration ) * current_outer_watercontent / dt * dampening
        solute_dumux[r,0] = root_soluteuptake
        print(solute_dumux)
        #s.setParameter( "Soil.BC.Bot.SValue", str(root_soluteuptake))
        #s.initializeProblem(maxDt = 0.01)
        
        s.setSource({NC-2: 1 * outer_watersource * volumes[-1] *length})
        s.setSource({0: root_soluteuptake * 62, NC-2: 1* outer_solutesource * volumes[-1] * length * 62}, eq_idx = 1) #factor 1000 is important for dumux as in the richards wrapper it is divided by 1000
        
        s.solve(dt, saveInnerFluxes_ = True)
        Wvolafter = cellVolumes*s.getWaterContent() # cm3
        Smassafter = s.getSolution(1) * Wvolafter # mol
        
            
        rootSoilFluxes = s.getInnerFlow(0, length) * dt # cm3
        rootSoilFluxesS = s.getInnerFlow(1, length) * dt # mol
        soilSoilFluxes = s.getOuterFlow(0, length) * dt # cm3
        soilSoilFluxesS = s.getOuterFlow(1, length) * dt # mol
        

    

    
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
        
        #alternative computation using the waterdemand
        Phi_A = abs(waterdemand)*r_root*((rho)**2/(2*(1-rho**2)))
        Phi_C = Phi_soil - Phi_A
        
        #if waterstress, then Phi_A = - rho**2*Phi_C
        if Phi_A*(1/rho**2-np.log(1/rho**2))+Phi_C <=0:
            Phi_A = Phi_soil*rho**2/(rho**2-1)
            Phi_C = -Phi_A*(1/rho**2-np.log(1/rho**2))
        
        Phi_A_orig = Phi_A
        Phi_C_orig = Phi_C
        #compute Phi_A based on the water demand
        #Phi_A = waterdemand * r_root / (2*(1-rho**2))
        #Phi_C = 
        #compute Phi_C based on the result
        
        #outer Phi
        Phi_out = Phi_A+Phi_C
        print("Phiroot", Phi_root, "Phi_soil", Phi_soil, "Phi_A", Phi_A, "Phi_C", Phi_C)
        r_rel = CC[0] / r_prhiz
        print("Phi_out", Phi_out, "Phi_in", Phi_A*(r_rel**2-np.log(r_rel**2))+Phi_C)
        
        #determine the critical matrix potential
        h_out = [-200, -800, -3200, -6000,-10000,-13700, -14000] # cm
        r_rel = 1/rho
        MFP_root = lambda Phi : Phi+abs(waterdemand)*r_root*((r_rel*rho)**2/(2*(1-rho**2))+rho**2/(1-rho**2)*(np.log(1/r_rel)-0.5))
        MFP_root_stress = lambda Phi : Phi*((r_rel*rho)**2 - 1 + 2*rho**2*np.log(1/(r_rel*rho))/(rho**2 - 1 + 2*rho**2*np.log(1/rho)))
        MFP_mean_difference = lambda Phi : Phi+abs(waterdemand)*r_root*((0.53*rho)**2/(2*(1-rho**2))+rho**2/(1-rho**2)*(np.log(1/0.53)-0.5)) - Phi_soil
        print(MFP_root_stress(vg.fast_mfp[sp](h_out[0])),MFP_root_stress(vg.fast_mfp[sp](h_out[-1])))
        if MFP_mean_difference(vg.fast_mfp[sp](h_out[0]))*MFP_mean_difference(vg.fast_mfp[sp](h_out[-1]))<0:
            MFP_out = root_scalar(MFP_mean_difference, method="brentq", bracket=[ vg.fast_mfp[sp](h_out[0]), vg.fast_mfp[sp](h_out[-1])]).root
        else:
            MFP_out = 1.e-3
        print(MFP_out)
        #MFP_out = MFP_out /10000
        
        
        for j in range(NC-1):
            r_rel = CC[j] / r_prhiz
            Phi = Phi_A*(r_rel**2-np.log(r_rel**2))+Phi_C
            mp_perirhizal1[r,j+1]=vg.fast_imfp[sp](Phi) # TODO: adapt this to stress maybe?
            
            #Phi_soil = Phi_out # TODO: remove this
            #if water_dumux[r,j+1]<soilVG[0]+0.02: # I model this as a threshhold for water stress
            #    Phi = Phi_soil+abs(waterdemand)*r_root*((r_rel*rho)**2/(2*(1-rho**2))-rho**2/(1-rho**2)*(np.log(r_rel**2)+0.5))
            #else:
            #    Phi = Phi_soil*r_root*((r_rel*rho)**2-1+2*rho**2*np.log(r_root/r_rel))/(rho**2-1+2*rho**np.log(1/(rho)))
            if Phi<=0: #water stress
                water_sr[r,j+1] = soilVG[0]+0.02
            else:
                water_sr[r,j+1] = vg.water_content(mp_perirhizal1[r,j+1], peri.sp)
            Phi = Phi_soil+abs(waterdemand)*r_root*((r_rel*rho)**2/(2*(1-rho**2))+rho**2/(1-rho**2)*(np.log(1/r_rel)-0.5)) #waterdemand is usually negative
            
            Phi_outer = vg.fast_mfp[sp](h_out[1])
            Phi = Phi_outer+abs(waterdemand)*r_root*((r_rel*rho)**2/(2*(1-rho**2))+rho**2/(1-rho**2)*(np.log(1/r_rel)-0.5)) #waterdemand is usually negative
            
            
            MFP_root = lambda h_outer : vg.fast_mfp[sp](h_outer)+abs(waterdemand)*r_root*((r_rel*rho)**2/(2*(1-rho**2))+rho**2/(1-rho**2)*(np.log(1/r_rel)-0.5))
            
            #print(MFP_root(-14500),MFP_root(-10))
            if MFP_root(vg.fast_mfp[sp](-14500))*MFP_root(vg.fast_mfp[sp](-10))<0:
                Phi_outer = root_scalar(MFP_root, method="brentq", bracket=[-14500,-10]).root
            else:
                Phi_outer = 1.0e-3
            Phi = MFP_out+abs(waterdemand)*r_root*((r_rel*rho)**2/(2*(1-rho**2))+rho**2/(1-rho**2)*(np.log(1/r_rel)-0.5)) #waterdemand is usually negative
            Phi_A = abs(waterdemand)*r_root*rho**2/(2*(1-rho**2))
            Phi_C = MFP_out-Phi_A
            MFP_root = lambda Phi : Phi+abs(waterdemand)*r_root*((r_rel*rho)**2/(2*(1-rho**2))+rho**2/(1-rho**2)*(np.log(1/r_rel)-0.5))
            MFP_root_stress = lambda Phi : Phi*((r_rel*rho)**2 - 1 + 2*rho**2*np.log(1/(r_rel*rho)))/(rho**2 - 1 + 2*rho**2*np.log(1/rho))
            Phi = MFP_root(MFP_out)
            if Phi<=0:
                mp_perirhizal2[r,j+1]=vg.fast_imfp[sp](MFP_root_stress(MFP_out))
            else:
                mp_perirhizal2[r,j+1]=vg.fast_imfp[sp](MFP_root(MFP_out))
            
            #alternative steady rate approximation            
            if Phi<=0: #water stress
                water_sr2[r,j+1] = soilVG[0]+0.001
            else:
                water_sr2[r,j+1] = vg.water_content(vg.fast_imfp[sp](Phi),peri.sp)
        
        #use the original values for Phi_A and Phi_C
        Phi_A = Phi_A_orig
        Phi_C = Phi_C_orig
        Phi_root = Phi_A*(1/rho**2-np.log(1/rho**2))+Phi_C
        Phi_test = Phi_A*(0.1**2-np.log(0.1**2))+Phi_C
        Phi_soil = Phi_A*(0.53**2-np.log(0.53**2))+Phi_C
        Phi_out = Phi_A*(1.0**2-np.log(1.0**2))+Phi_C
        # run the alpha omega model on solute flow (both steady state and steady rate)
        waterflow = root_wateruptake
        result_solutes_ss = peri.soil_root_solutes_ss_([Phi_root], [Phi_test], [mean_solute], [Vmax_per_area], [Km], Ds, [-abs(radial_waterdemand)], peri.sp)
        result_solutes_sr = peri.soil_root_solutes_sr_([Phi_root], [Phi_soil], [rho], [mean_solute], [Vmax_per_area], [Km], Ds, [-abs(radial_waterdemand)], peri.sp)
        result_solutes_ss = result_solutes_ss[0]
        result_solutes_sr = result_solutes_sr[0]
        print("steadystate", result_solutes_ss, "steadyrate", result_solutes_sr)
        
        solute_ss[r,1] = result_solutes_ss
        solute_ss[r,0] = Vmax_per_area * result_solutes_ss / (Km + result_solutes_ss)
        solute_sr[r,1] = result_solutes_sr
        solute_sr[r,0] = Vmax_per_area * result_solutes_sr / (Km + result_solutes_sr)
        
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
            F_tilde=math.exp(-D_tilde*F)
            solute_ss[r,j+1] = result_solutes_ss * F_tilde - (1-F_tilde) * solute_ss[r,0] / abs(waterdemand)#waterdemand is assumed to be negative
            solute_sr[r,j+1] = result_solutes_sr * F_tilde - (1-F_tilde) * solute_sr[r,0] / abs(waterdemand)#waterdemand is assumed to be negative

        solutes_Tiina[r,0] = peri.solutesuptake_convdiff_([mean_water],[outer_conc], [Vmax], [Km], Ds, [waterflow], [r_root], [0.], [time], peri.sp)[0] 
    return water_dumux, solute_dumux, solute_ss, water_sr, solute_sr, water_sr2, mp_perirhizal1, mp_perirhizal2, solutes_Tiina

if do_computation:
    # save everything in the np arrays
    for i in range(n_tests):
        water_dumux, solute_dumux, solute_ss, water_sr, solute_sr, water_sr2, mp_perirhizal1, mp_perirhizal2, solutes_Tiina = run_perirhizal_test()
        watercontent_dumux[i,:,:]=water_dumux[:,:]
        soluteconcentration_dumux[i,:,:]=solute_dumux[:,:]
        soluteconcentration_steadystate[i,:,:]=solute_ss[:,:]
        watercontent_steadyrate[i,:,:]=water_sr[:,:]
        watercontent_steadyrate2[i,:,:]=water_sr2[:,:] #TODO this is confusing
        soluteconcentration_steadyrate[i,:,:]=solute_sr[:,:]
        matrixpotential_perirhizal1[i,:,:]=mp_perirhizal1[:,:]
        matrixpotential_perirhizal2[i,:,:]=mp_perirhizal2[:,:]
        soluteconcentration_Tiina[i,:,:]=solutes_Tiina[:,:]
    
    np.savez("test_perirhizal.npz", watercontent_dumux=watercontent_dumux, soluteconcentration_Tiina=soluteconcentration_Tiina, watercontent_steadyrate2=watercontent_steadyrate2, matrixpotential_perirhizal1=matrixpotential_perirhizal1, matrixpotential_perirhizal2=matrixpotential_perirhizal2, soluteconcentration_dumux=soluteconcentration_dumux, soluteconcentration_steadystate=soluteconcentration_steadystate, watercontent_steadyrate=watercontent_steadyrate, soluteconcentration_steadyrate=soluteconcentration_steadyrate)

else:
    simulation_results = np.load("test_perirhizal.npz")
    watercontent_dumux = simulation_results["watercontent_dumux"]
    soluteconcentration_dumux = simulation_results["soluteconcentration_dumux"]
    soluteconcentration_steadystate = simulation_results["soluteconcentration_steadystate"]
    watercontent_steadyrate = simulation_results["watercontent_steadyrate"]
    soluteconcentration_steadyrate = simulation_results["soluteconcentration_steadyrate"]
    watercontent_steadyrate2 = simulation_results["watercontent_steadyrate2"]
    matrixpotential_perirhizal1 = simulation_results["matrixpotential_perirhizal1"]
    matrixpotential_perirhizal2 = simulation_results["matrixpotential_perirhizal2"]
    soluteconcentration_Tiina = simulation_results["soluteconcentration_Tiina"]
    


# compare both for the differint means of water / solute content
run = 0
timestep = np.array([1,3,5,7,9])
timestep = np.array(np.linspace(1,7,num=5)) #TODO remove
print(timestep)
for i in range(5):
    timestep[i] = int(n_times * timestep[i] / 10)
timestep = timestep.astype(int)
#plot water and nitrogen in that one voxel
#fig, ax1 = figure_style.subplots12(1, 5)
fig, ax1 = figure_style.subplots12(nrows=5, ncols=1)

print(timestep)



linestyle_dumux = "solid"
linestyle_steadystate = "dotted"
linestyle_steadyrate = "dashed"
linestyle_special = "dashdot"
    
#print(water_dumux)    
#print(CC)
#timestep = [100,300,5,7,9]
for i in range(5):
    ax2 = ax1[i].twinx()
    #print("watercontent dumux", watercontent_dumux)
    water_dumux = watercontent_dumux[run, timestep[i], 1:]
    water_perirhizal = watercontent_steadyrate[run, timestep[i], 1:]
    water_steadyrate2 = watercontent_steadyrate2[run, timestep[i], 1:]
    solute_dumux = soluteconcentration_dumux[run, timestep[i], 1:]
    solute_steadystate = soluteconcentration_steadystate[run, timestep[i], 1:]
    solute_steadyrate = soluteconcentration_steadyrate[run, timestep[i], 1:]
    solute_Tiina = np.array([soluteconcentration_Tiina[run, timestep[i], 0] + 1*(outer_conc - soluteconcentration_Tiina[run, timestep[i], 0]) * cc / r_prhiz for cc in CC]) #np.linspace(soluteconcentration_steadyrate[run, timestep[i], 0], outer_conc, NC-1)

    ax1[i].plot(CC, water_dumux, "b", linestyle = linestyle_dumux, label = "water_dumux")
    ax1[i].plot(CC, water_perirhizal, "b", linestyle = linestyle_steadyrate, label = "water_perirhizal")
    #ax1[i].plot(CC, water_steadyrate2, "b", linestyle = linestyle_steadyrate, label = "water_perirhizal2")
    ax2.plot(CC, solute_dumux, "m", linestyle = linestyle_dumux, label = "solute_dumux")
    #ax2.plot(CC, solute_steadystate, "m", linestyle = linestyle_steadystate, label = "solute_steadystate")
    #ax2.plot(CC, solute_steadyrate, "m", linestyle = linestyle_steadyrate, label = "solute_steadyrate")
    ax2.plot(CC, solute_Tiina, "m", linestyle = linestyle_steadystate, label = "solute_Tiina")

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
    
    mp_steadyrate1 = matrixpotential_perirhizal1[run, timestep[i], 1:]
    mp_steadyrate2 = matrixpotential_perirhizal2[run, timestep[i], 1:]
    
    
    water_dumux = vg.pressure_head(water_dumux, peri.sp)
    water_perirhizal = vg.pressure_head(water_perirhizal, peri.sp)
    #water_perirhizal2 = vg.pressure_head(water_steadyrate2, peri.sp)

    ax1[i].plot(CC, water_dumux, "b", linestyle = linestyle_dumux, label = "water_dumux")
    #ax1[i].plot(CC, water_perirhizal, "b", linestyle = linestyle_steadyrate, label = "water_perirhizal") #this and the following line give the same results, as they should
    ax1[i].plot(CC, mp_steadyrate1, "b", linestyle = linestyle_steadyrate, label = "water_sr")
    #ax1[i].plot(CC, mp_steadyrate2, "b", linestyle = linestyle_steadystate, label = "water_sr_stress")
    ax2.plot(CC, solute_dumux, "m", linestyle = linestyle_dumux, label = "solute_dumux")
    #ax2.plot(CC, solute_steadystate, "m", linestyle = linestyle_steadystate, label = "solute_steadystate")
    #ax2.plot(CC, solute_steadyrate * solute_dumux[-1] / solute_steadyrate[-1], "m", linestyle = linestyle_steadyrate, label = "solute_steadyrate")
    ax2.plot(CC, solute_steadystate * solute_dumux[-1] / solute_steadystate[-1], "m", linestyle = linestyle_steadyrate, label = "solute_steadystate")
    #ax2.plot(CC, solute_Tiina, "m", linestyle = linestyle_steadystate, label = "solute_Tiina")
    ax2.plot(CC, solute_Tiina * solute_dumux[-1] / outer_conc, "m", linestyle = linestyle_special, label = "solute_Tiina_scaled")

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


