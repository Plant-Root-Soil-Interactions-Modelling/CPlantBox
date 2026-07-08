
# this is a test file by Erik Kopp to test the alternative implementation of the perirhizal resistances (water) and solute flow

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time

import sys; sys.path.append("../.."); sys.path.append("../../src/")
import plantbox as pb
#from plantbox import Perirhizal
from plantbox.functional.Perirhizal import PerirhizalPython
import plantbox.functional.van_genuchten as vg
from plantbox.visualisation import figure_style
#import Perirhizal


# dumux rosi imports
#sys.path.append("../../dumux-rosi/build-cmake/cpp/python_binding/")
#sys.path.append("../../../dumuxtest/dumux/dumux-rosi/build-cmake/cpp/python_binding/")
from rosi.richards_flat import RichardsFlatWrapper as RichardsWrapper  # Python part, macroscopic soil model
#from richards import RichardsWrapper  # Python part
from rosi.richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial
from rosi.rosi_richardsnc_cyl import RichardsNCCylFoam # C++ part (Dumux binding), macroscopic soil model
#from rosi_richardsnc_cyl import RichardsNCCylFoam as RichardsNC_cyl  # C++ part (Dumux binding)
#from rosi_richards22c import RichardsNCSPILU as RichardsNCSP #test



#numerics
import math
from scipy.optimize import fsolve, root_scalar
from numpy import linalg as LA

# run the dumux implementation of root water and nitrate uptake, later compare it to the analytical approximation

n_tests = 1 #try everything here for this many random parameter sets
do_computation = True #should the computation be run or take the data from a saved file

# general parameters
max_time = 0.1 # d
n_times = 10 # number of time intervals
r_prhiz = 0.6 # perirhizal radius[cm]
r_root = 0.02 # root radius [cm]
NC = 40 # number of spatial discretisations
length = 1 #default length of the segment, will not change the outcpme as all variables are constant in this direction [cm]

#two scenarios will be computed: one without inflow, another with a Dirichlet BC
n_scenarios = 2

#initial conditions
initial_waterpotential = -200
initial_soluteconcentration = 2e-5#mol/cm3

#space for the outputs
# number of tests, number of scenarios, number of timesteps, (rootuptake, inflow, discretisation)
#water outputs
#simulations in dumux
watercontent_dumux = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #watercontent in cm3/cm3, (radial) uptake and inflow in cm2/d 
waterpotential_dumux = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #waterpotential in cm, (radial) uptake and inflow in cm2/d 
#analytical approximations
watercontent_sr = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #watercontent in cm3/cm3, (radial) uptake and inflow in cm2/d 
waterpotential_sr = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #waterpotential in cm, (radial) uptake and inflow in cm2/d 

#solute outputs
#simulations in dumux
solutes_dumux = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, (radial) uptake and inflow in mol/(cm d) 
#analytical approximations
#steady state approximations
solutes_dumux_ss = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, (radial) uptake and inflow in mol/(cm d) 
#steady rate approximations
solutes_dumux_sr = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, (radial) uptake and inflow in mol/(cm d) 
#general steady rate with far field approximation
solutes_dumux_ff = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, (radial) uptake and inflow in mol/(cm d) 



#discretisation
lb = 0.5
points = np.logspace(np.log(r_root) / np.log(lb), np.log(r_prhiz) / np.log(lb), NC+1, base = lb) 
CC = np.array([(points[i] + points[i+1])/2 for i in range(NC)])
volumes = np.array([(points[i+1]**2 - points[i]**2)*3.14 for i in range(NC)])




def run_perirhizal_test(max_time, n_times, r_prhiz, r_root, NC, points, CC, volumes, length, n_scenarios, initial_waterpotential, initial_soluteconcentration):
    
    #space for the outputs
    watercontent_dumux = np.zeros((n_scenarios, n_times, NC+2)) #watercontent in cm3/cm3, (radial) uptake and inflow in cm/d 
    waterpotential_dumux = np.zeros((n_scenarios, n_times, NC+2)) #waterpotential in cm, (radial) uptake and inflow in cm/d 
    watercontent_sr = np.zeros((n_scenarios, n_times, NC+2)) #watercontent in cm3/cm3, (radial) uptake and inflow in cm/d 
    waterpotential_sr = np.zeros((n_scenarios, n_times, NC+2)) #waterpotential in cm, (radial) uptake and inflow in cm/d 
    solutes_dumux = np.zeros((n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, (radial) uptake and inflow in mol/(cm2d) 
    solutes_dumux_sr = np.zeros((n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, (radial) uptake and inflow in mol/(cm2d) 
    solutes_dumux_ss = np.zeros((n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, (radial) uptake and inflow in mol/(cm2d) 
    solutes_dumux_ff = np.zeros((n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, (radial) uptake and inflow in mol/(cm2d)
    
    simtimes = np.linspace(0,max_time,n_times+1)[1:]
    dt = max_time / n_times
    rho = r_prhiz / r_root
    
    #soil parameters
    soilVG = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam soil 
    
    #discretisation
    #lb = 0.5
    #points = np.logspace(np.log(r_root) / np.log(lb), np.log(r_prhiz) / np.log(lb), NC+1, base = lb) 
    #CC = np.array([(points[i] + points[i+1])/2 for i in range(NC)])
    #volumes = np.array([(points[i+1]**2 - points[i]**2)*3.14 for i in range(NC)])
    
    
    # root conductivity and solute uptake parameters, constant throughout the entire simulation time
    molarMassWater = 18 #g/mol
    molarMassSolute = 62 #g/mol
    root_conductivity = 1e-4 #1/d
    inner_kr = root_conductivity * r_root * 2 * 3.14
    waterdemand = -0.05 #cm/d
    radial_waterdemand = 2*3.14*r_root * waterdemand #cm2/d
    Vmax = 4.0e-11 * 1 * (2*3.14*r_root) * (24*3600)  #mol/(cm2 s) * cm * cm * (s/d) -> mol / d 
    Vmax_per_area = Vmax / (1 * (2*3.14*r_root)) #mol / d /cm2 = mol/(cm2d)
    Km = 1.5e-7   #mol/cm3
    
    #diffusion coefficient of nitrate
    Ds = 1.902e-5 * 24 * 3600  #cm2/s -> cm2/d
    

    # load the perirhizal model
    peri = PerirhizalPython() 
    sp = vg.Parameters(soilVG)
    peri.set_soil(sp)
    #no lookup tables are used here as there arent many simulations
    
    


    # initialise the dumux models for the scenarios
    s_nf = RichardsWrapper(RichardsNCCylFoam()) #no flux outer BC
    s_g = RichardsWrapper(RichardsNCCylFoam()) #Cauchy outer BC

    for s in [s_nf, s_g]:
        s.initialize()
        s.createGrid1d(points, length = length/100)  # [m] -> [cm]
        s.setVGParameters([soilVG])
        s.setHomogeneousIC(initial_waterpotential)  # cm pressure head
   
        s.setBotBC("constantFluxCyl",waterdemand) # "noFlux")# Flux in cm/d
        if s == s_nf:
            s.setTopBC("constantFluxCyl",0.0)  #  [cm/day] "noFlux")#default, will be changed for one scenario
            s.setParameter("Soil.BC.Top.SType", "3")  # constantFluxCyl=3 (SType = Solute Type)
            s.setParameter("Soil.BC.Top.CValue", "0.0") 
        else:
            s.setTopBC("constantPressure",initial_waterpotential)
            s.setParameter("Soil.BC.Top.SType", "10")  #advective flow
            s.setParameter("Soil.BC.Top.CValue", str(initial_soluteconcentration*molarMassSolute)) 
        s.setParameter("Soil.BC.Bot.SType", "8")  # michaelisMenten=8 (SType = Solute Type)
        s.setParameter("Soil.BC.Bot.CValue", "0.0") #should not matter
        #s.setParameter("RootSystem.Uptake.Vmax", str(Vmax_per_area*molarMassSolute)) #mol/(cm2d) -> g/(cm2 d)
        s.setParameter("RootSystem.Uptake.Vmax", str(Vmax*molarMassSolute)) #mol/d -> g/d #TODO: Vmax or Vmax per area?
        s.setParameter("RootSystem.Uptake.Km", str(Km*molarMassSolute)) # mol/cm3 -> g/cm3
        "RootSystem.Uptake.Vmax"
        "RootSystem.Uptake.Km"
        s.setParameter("Soil.IC.C", str(initial_soluteconcentration*molarMassSolute))  # g / cm3  # TODO specialised setter?
        s.setParameter("Component.MolarMass", str(molarMassWater/1000)) #g/mol -> kg/mol water
        s.setParameter("1.Component.MolarMass", str(molarMassSolute/1000)) #g/mol -> kg/mol nitrate
        s.setParameter("1.Component.LiquidDiffusionCoefficient", str(Ds / 1.e4 / (24*3600))) #cm^2/s -> m^2/s
        s.initializeProblem(maxDt = 0.01)
        s.ddt = 1.e-4  # days
    
    #in the general case the outer conditions are kept constant
    #s_g.setTopBC("constantPressure",initial_waterpotential)  
    #s_g.setParameter("Soil.BC.Top.SType", "1")  #Dirichlet BC
    #s_g.setParameter("Soil.BC.Top.SType", "10")  #advective flow
    #s_g.setParameter("Soil.BC.Top.CValue", str(initial_soluteconcentration*molarMassSolute)) 


    cellVolumes = s_g.getCellSurfacesCyl() * length # cm3
    
    #variables for the no flux case
    mean_watercontent_nf = 0.1 #cm3/cm3
    mean_waterpotential_nf = -100 #cm
    mean_solutecontent_nf = 0 #mol/cm3
    Phi_soil_nf = 1.0 #mfp of the mean soil
    Phi_root_nf = 1.0 #mfp next to the root
    #mfp parameters for the perirhizal model
    Phi_outer_nf = 1.0 #mfp at the outer perirhizal radius
    rootuptake_w_nf = 1.0 #water uptake of the root, cm / d
    inflow_w_nf = 0.0 #water uptake from outside the perirhizal zone, cm/d
    rootuptake_s_nf = 1.0 #water uptake of the root, mol / cm2d
    inflow_s_nf = 0.0 #water uptake from outside the perirhizal zone, mol / cm2d
    
    #variables for the general case
    mean_watercontent_g = 0.1 #cm3/cm3
    mean_waterpotential_g = -100 #cm
    mean_solutecontent_g = 0 #mol/cm3
    Phi_soil_g = 1.0 #mfp of the mean soil
    Phi_root_g = 1.0 #mfp next to the root
    #mfp parameters for the perirhizal model
    Phi_outer_g = 1.0 #mfp at the outer perirhizal radius
    rootuptake_w_g = 1.0 #water uptake of the root, cm / d
    inflow_w_g = 0.0 #water uptake from outside the perirhizal zone, cm/d
    rootuptake_s_g = 1.0 #water uptake of the root, mol / cm2d
    inflow_s_g = 0.0 #water uptake from outside the perirhizal zone, mol / cm2d
    
    
    for r, time in enumerate(simtimes):
        
        dt = time
        if r>0:
            dt = time - simtimes[r-1]
        
        print('time',time)
        print('no flux BC')
        print("*****", "#", r, "external time step", dt, " d, simulation time", s_nf.simTime, "d, internal time step", s_nf.ddt, "d")
        print('general')
        print("*****", "#", r, "external time step", dt, " d, simulation time", s_g.simTime, "d, internal time step", s_g.ddt, "d")
        
        #one timestep
        s_nf.solve(dt, saveInnerFluxes_ = True)
        s_g.solve(dt, saveInnerFluxes_ = True)
        
        #watercontent and solute content, discretised
        watercontent_nf = s_nf.getWaterContent() # cm3
        waterpotential_nf = np.array([vg.pressure_head(watercontent_nf[i],peri.sp) for i in range(NC)]) # cm
        watercontent_g = s_g.getWaterContent() # cm3
        waterpotential_g = np.array([vg.pressure_head(watercontent_g[i],peri.sp) for i in range(NC)]) # cm
        solutes_nf = s_nf.getSolution(1) / molarMassSolute # mol/cm3
        solutes_g = s_g.getSolution(1) / molarMassSolute # mol/cm3
            
        #inflow and outflow    
        rootuptake_w_nf = s_nf.getInnerFlow(0, length) /(length*(2*np.pi*r_root)) # cm /d
        rootuptake_s_nf = s_nf.getInnerFlow(1, length) / molarMassSolute /(length*(2*np.pi*r_root)) # mol / cm2d
        inflow_w_nf = s_nf.getOuterFlow(0, length) /(length*(2*np.pi*r_prhiz)) # cm /d
        inflow_s_nf = s_nf.getOuterFlow(1, length) / molarMassSolute /(length*(2*np.pi*r_prhiz)) # mol / cm2d
        
        rootuptake_w_g = s_g.getInnerFlow(0, length) /(length*(2*np.pi*r_root)) # cm /d
        rootuptake_s_g = s_g.getInnerFlow(1, length) / molarMassSolute /(length*(2*np.pi*r_root)) # mol / cm2d
        inflow_w_g = s_g.getOuterFlow(0, length) /(length*(2*np.pi*r_prhiz)) # cm /d
        inflow_s_g = s_g.getOuterFlow(1, length) / molarMassSolute /(length*(2*np.pi*r_prhiz)) # mol / cm2d
        
        rootuptake_w_nf = rootuptake_w_nf[0]
        rootuptake_s_nf = rootuptake_s_nf[0]
        inflow_w_nf = inflow_w_nf[0]
        inflow_s_nf = inflow_s_nf[0]
        rootuptake_w_g = rootuptake_w_g[0]
        rootuptake_s_g = rootuptake_s_g[0]
        inflow_w_g = abs(inflow_w_g[0])
        inflow_s_g = abs(inflow_s_g[0]) #TODO take all the flows to the root as negative?
        
        #store the dumux outputs
        #scenario no flux outer BC
        watercontent_dumux[0,r,0]=rootuptake_w_nf
        watercontent_dumux[0,r,1]=inflow_w_nf
        watercontent_dumux[0,r,2:]=watercontent_nf
        waterpotential_dumux[0,r,0]=rootuptake_w_nf
        waterpotential_dumux[0,r,1]=inflow_w_nf
        waterpotential_dumux[0,r,2:]=waterpotential_nf
        solutes_dumux[0,r,0]=rootuptake_s_nf
        solutes_dumux[0,r,1]=inflow_s_nf
        solutes_dumux[0,r,2:]=solutes_nf
        #general scenario
        watercontent_dumux[1,r,0]=rootuptake_w_g
        watercontent_dumux[1,r,1]=inflow_w_g
        watercontent_dumux[1,r,2:]=watercontent_g
        waterpotential_dumux[1,r,0]=rootuptake_w_g
        waterpotential_dumux[1,r,1]=inflow_w_g
        waterpotential_dumux[1,r,2:]=waterpotential_g
        solutes_dumux[1,r,0]=rootuptake_s_g
        solutes_dumux[1,r,1]=inflow_s_g
        solutes_dumux[1,r,2:]=solutes_g
        
        #determine means
        #no flux outer BC
        mean_watercontent_nf = np.average(watercontent_nf, weights=volumes)
        mean_waterpotential_nf = vg.pressure_head(mean_watercontent_nf, peri.sp)
        mean_soluteconcent_nf = np.average(solutes_nf, weights=np.multiply(mean_watercontent_nf, volumes))
        #general outer BC
        mean_watercontent_g = np.average(watercontent_g, weights=volumes)
        mean_waterpotential_g = vg.pressure_head(mean_watercontent_g, peri.sp)
        mean_soluteconcent_g = np.average(solutes_g, weights=np.multiply(mean_watercontent_g, volumes))
        
        #determine coefficients for the analytical approximation
        #equation [4] in Schroeder2008 doi:10.2136/vzj2007.0114
        #Phi(r)=Phi_outer + (q_root*r_root-q_out*r_prhiz)*(rho**2)/(1-rho**2)*(((r/r_prhiz)**2-1)/2-ln(r/r_prhiz))+q_out*r_prhiz*ln(r/r_prhiz)
        #Assume Phi(0.53r)=Phi(mean water potential)
        #both q_root and q_out are assumed to have positive signs if the water flows to the root
        #no flux outer BC
        Phi_soil_nf = vg.fast_mfp[peri.sp](mean_waterpotential_nf) #mfp of the mean soil
        Phi_outer_nf = Phi_soil_nf - (rootuptake_w_nf*r_root-inflow_w_nf*r_prhiz)*(rho**2)/(1-rho**2)*(((0.53)**2-1)/2-np.log(0.53))-(inflow_w_nf*r_prhiz)*r_prhiz*np.log(0.53) #mfp at the outer perirhizal radius
        Phi_nf = lambda r: Phi_outer_nf + (rootuptake_w_nf*r_root-inflow_w_nf*r_prhiz)*(rho**2)/(1-rho**2)*(((r/r_prhiz)**2-1)/2-np.log(r/r_prhiz))+(inflow_w_nf*r_prhiz)*r_prhiz*np.log(r/r_prhiz)#mfp function depending on radius
        Phi_root_nf = Phi_nf(r_root)#mfp next to the root 
        #general outer BC: Dirichlet BC to a fixed potential
        Phi_soil_g = vg.fast_mfp[peri.sp](mean_waterpotential_g) #mfp of the mean soil
        Phi_outer_g = Phi_soil_g - (rootuptake_w_g*r_root-inflow_w_g*r_prhiz)*(rho**2)/(1-rho**2)*(((0.53)**2-1)/2-np.log(0.53))-(inflow_w_g*r_prhiz)*np.log(0.53) #mfp at the outer perirhizal radius
        Phi_g = lambda r: Phi_outer_g + (rootuptake_w_g*r_root-inflow_w_g*r_prhiz)*(rho**2)/(1-rho**2)*(((r/r_prhiz)**2-1)/2-np.log(r/r_prhiz))+(inflow_w_g*r_prhiz)*np.log(r/r_prhiz)#mfp function depending on radius
        Phi_root_g = Phi_nf(r_root)#mfp next to the root 
        print("inflow", inflow_w_g)
    
        #write the steady rate approximations for the water as outputs
        waterpotential_sr[0,r,0] = rootuptake_w_nf
        waterpotential_sr[0,r,1] = inflow_w_nf
        waterpotential_sr[0,r,2:] = np.array([vg.fast_imfp[peri.sp](Phi_nf(CC[i])) for i in range(NC)])
        watercontent_sr[0,r,0] = rootuptake_w_nf
        watercontent_sr[0,r,1] = inflow_w_nf
        watercontent_sr[0,r,2:] = np.array([vg.water_content(waterpotential_sr[0,r,2+i],peri.sp) for i in range(NC)])
        
        waterpotential_sr[1,r,0] = rootuptake_w_g
        waterpotential_sr[1,r,1] = inflow_w_g
        waterpotential_sr[1,r,2:] = np.array([vg.fast_imfp[peri.sp](Phi_g(CC[i])) for i in range(NC)])
        watercontent_sr[1,r,0] = rootuptake_w_g
        watercontent_sr[1,r,1] = inflow_w_g
        watercontent_sr[1,r,2:] = np.array([vg.water_content(waterpotential_sr[1,r,2+i],peri.sp) for i in range(NC)])
        
        
        #compute the analytical approximations for the solute uptake
        #case of dumux no flux outer BC
        #use the outer concentration for the steady state approximation here, even if that is not necessarily known in the macroscopic model
        #(the outer concentration can be computed as a mean of all surrounding voxels?
        #TODO: check weather Vmax or Vmax per area
        #result_solutes_ss_nf = peri.soil_root_solutes_ss_([Phi_root_nf], [Phi_outer_nf], [solutes_nf[-1]], [Vmax_per_area], [Km], Ds, [radial_waterdemand], peri.sp)
        #result_solutes_sr_nf = peri.soil_root_solutes_sr_([Phi_root_nf], [Phi_outer_nf], [rho], [mean_soluteconcent_nf], [Vmax_per_area], [Km], Ds, [radial_waterdemand], peri.sp)
        
        result_solutes_ss_nf = peri.soil_root_solutes_steadyrate_simplified_([Phi_root_nf], [Phi_soil_nf], [r_root], [r_prhiz], [mean_soluteconcent_nf], [Vmax_per_area], [Km], Ds, [abs(waterdemand)], peri.sp, n_approx = 5)
        result_solutes_sr_nf = peri.soil_root_solutes_steadyrate_simplified_([Phi_root_nf], [Phi_soil_nf], [r_root], [r_prhiz], [mean_soluteconcent_nf], [Vmax_per_area], [Km], Ds, [abs(radial_waterdemand)], peri.sp, n_approx = 1)
        
        #safe the results
        #(here there is no computed inflow, so index 1 doesn't do anything)
        result_solutes_ss_nf = result_solutes_ss_nf[0]
        solutes_dumux_ss[0,r,0]=-Vmax_per_area * result_solutes_ss_nf / (Km + result_solutes_ss_nf)
        result_solutes_sr_nf = result_solutes_sr_nf[0]
        solutes_dumux_sr[0,r,0]=-Vmax_per_area * result_solutes_sr_nf / (Km + result_solutes_sr_nf)
        
        F0_nf = peri.integral_AdvectionDiffusion_(Phi_root_nf,peri.sp)
        F0_g = peri.integral_AdvectionDiffusion_(Phi_root_g,peri.sp)
        D_tilde = 1/Ds/math.pow(sp.theta_S-sp.theta_R,13/3)*(sp.theta_S*sp.theta_S)
        
        for j in range(NC):
            r_current = CC[j]
            F_nf = peri.integral_AdvectionDiffusion_(Phi_nf(r_current),peri.sp)-F0_nf
            F_g = peri.integral_AdvectionDiffusion_(Phi_g(r_current),peri.sp)-F0_g
            F_tilde_nf=math.exp(D_tilde*F_nf)
            F_tilde_g=math.exp(D_tilde*F_g)
            solutes_dumux_ss[0,r,2+j] = result_solutes_ss_nf * F_tilde_nf + (1-F_tilde_nf) * solutes_dumux_ss[0,r,0] / waterdemand#waterdemand is assumed to be negative #TODO: look at the signs
            solutes_dumux_sr[0,r,2+j] = result_solutes_sr_nf * F_tilde_nf + (1-F_tilde_nf) * solutes_dumux_sr[0,r,0] / radial_waterdemand #an uptake of both water and solute is assumed
        
        #case of general steady rate water uptake
        #for the steady state take again the outer concentration
        #TODO: check weather Vmax or Vmax per area
        rsc, Uptake, ss_uptake, sr_uptake, _ = peri.soil_root_solutes_sr([Phi_outer_g], [rootuptake_w_g*2*np.pi*r_root], [inflow_w_g*2*np.pi*r_prhiz], [r_root], [r_prhiz], [mean_soluteconcent_g], [initial_soluteconcentration], [Vmax_per_area], [Km], [Ds], peri.sp, mode = "ss")
        _, _, soluteconcentration = peri.watersolutes_disc(Phi_outer_g, rootuptake_w_g*r_root, inflow_w_g*r_prhiz, r_root, r_prhiz, CC, rsc, Ds, ss_uptake, sr_uptake, peri.sp)
        solutes_dumux_ss[1,r,0] = -Uptake[0]
        solutes_dumux_ss[1,r,1] = -ss_uptake[0]
        solutes_dumux_ss[1,r,2:] = soluteconcentration[:]
        solutes_dumux_ss[1,r,0] = -Vmax_per_area * soluteconcentration[0] / (Km + soluteconcentration[0]) #TODO: check weather Vmax or Vmax per area
        #steady rate no flux outer BC
        rsc, Uptake, ss_uptake, sr_uptake, _ = peri.soil_root_solutes_sr([Phi_outer_g], [rootuptake_w_g*2*np.pi*r_root], [inflow_w_g*2*np.pi*r_prhiz], [r_root], [r_prhiz], [mean_soluteconcent_g], [initial_soluteconcentration], [Vmax_per_area], [Km], [Ds], peri.sp, mode = "sr")
        _, _, soluteconcentration = peri.watersolutes_disc(Phi_outer_g, rootuptake_w_g*r_root, inflow_w_g*r_prhiz, r_root, r_prhiz, CC, rsc, Ds, ss_uptake, sr_uptake, peri.sp)
        solutes_dumux_sr[1,r,0] = -Uptake[0]
        solutes_dumux_sr[1,r,1] = -ss_uptake[0]
        solutes_dumux_sr[1,r,2:] = soluteconcentration[:]
        solutes_dumux_sr[1,r,0] = -Vmax_per_area * soluteconcentration[0] / (Km + soluteconcentration[0]) #TODO: check weather Vmax or Vmax per area
        #steady rate solute uptake with the farfield approximation 
        rsc, Uptake, ss_uptake, sr_uptake, _ = peri.soil_root_solutes_sr([Phi_outer_g], [rootuptake_w_g*2*np.pi*r_root], [inflow_w_g*2*np.pi*r_prhiz], [r_root], [r_prhiz], [mean_soluteconcent_g], [initial_soluteconcentration], [Vmax_per_area], [Km], [Ds], peri.sp, mode = "ff")
        _, waterpotential, soluteconcentration = peri.watersolutes_disc(Phi_outer_g, 2*np.pi * rootuptake_w_g*r_root, 2*np.pi * inflow_w_g*r_prhiz, r_root, r_prhiz, CC, rsc, Ds, ss_uptake, sr_uptake, peri.sp)
        waterpotential_sr[1,r,2:] = waterpotential #TODO: this is a test #vg.fast_mfp[peri.sp](waterpotential_dumux[1,r,-1])
        solutes_dumux_ff[1,r,0] = -Uptake[0]
        solutes_dumux_ff[1,r,1] = -ss_uptake[0]
        solutes_dumux_ff[1,r,2:] = soluteconcentration[:]
        solutes_dumux_ff[1,r,0] = -Vmax_per_area * soluteconcentration[0] / (Km + soluteconcentration[0]) #TODO: check weather Vmax or Vmax per area
        
        
    return watercontent_dumux, waterpotential_dumux, watercontent_sr, waterpotential_sr, solutes_dumux, solutes_dumux_sr, solutes_dumux_ss, solutes_dumux_ff
    

if do_computation:
    # save everything in the np arrays
    for i in range(n_tests):
        watercontent_dumux[i,:,:,:], waterpotential_dumux[i,:,:,:], watercontent_sr[i,:,:,:], waterpotential_sr[i,:,:,:], solutes_dumux[i,:,:,:], solutes_dumux_sr[i,:,:,:], solutes_dumux_ss[i,:,:,:], solutes_dumux_ff[i,:,:,:] = run_perirhizal_test(max_time, n_times, r_prhiz, r_root, NC, points, CC, volumes, length, n_scenarios, initial_waterpotential, initial_soluteconcentration)
    
    np.savez("test_perirhizal.npz", 
    watercontent_dumux=watercontent_dumux, 
    waterpotential_dumux=waterpotential_dumux, 
    watercontent_sr=watercontent_sr, 
    waterpotential_sr=waterpotential_sr, 
    solutes_dumux=solutes_dumux, 
    solutes_dumux_sr=solutes_dumux_sr, 
    solutes_dumux_ss=solutes_dumux_ss, 
    solutes_dumux_ff=solutes_dumux_ff)

else:
    simulation_results = np.load("test_perirhizal.npz")
    watercontent_dumux = simulation_results["watercontent_dumux"]
    waterpotential_dumux = simulation_results["waterpotential_dumux"]
    watercontent_sr = simulation_results["watercontent_sr"]
    waterpotential_sr = simulation_results["waterpotential_sr"]
    solutes_dumux = simulation_results["solutes_dumux"]
    solutes_dumux_sr = simulation_results["solutes_dumux_sr"]
    solutes_dumux_ss = simulation_results["solutes_dumux_ss"]
    solutes_dumux_ff = simulation_results["solutes_dumux_ff"]
    
    


# compare both for the differint means of water / solute content
run = 0
#timestep = np.array([1,3,5,7,9])
timestep = np.array(np.linspace(1,9,num=5)) 
for i in range(5):
    timestep[i] = int(n_times * timestep[i] / 10)
timestep = timestep.astype(int)


linestyle_dumux = "solid"
linestyle_steadystate = "dotted"
linestyle_steadyrate = "dashed"
linestyle_special = "dashdot"
    
#peri = PerirhizalPython() 
#sp = vg.Parameters(soilVG)
#peri.set_soil(sp)

fig, ax1 = figure_style.subplots12(nrows=5, ncols=2)
# dumux(both)
# left: sr no flux, ss for the no flux
# right: ss, sr, farfield for sr waterflow



for i in range(5):
    ax2_0 = ax1[i,0].twinx()
    ax2_1 = ax1[i,1].twinx()
    
    #load data
    water_dumux_nf = waterpotential_dumux[run, 0, timestep[i], 2:]
    water_dumux_g = waterpotential_dumux[run, 1, timestep[i], 2:]
    water_steadyrate_nf = waterpotential_sr[run, 0, timestep[i], 2:]
    water_steadyrate_g = waterpotential_sr[run, 1, timestep[i], 2:]
    solute_dumux_nf = solutes_dumux[run, 0, timestep[i], 2:]
    solute_sr_nf = solutes_dumux_sr[run, 0, timestep[i], 2:]
    solute_ss_nf = solutes_dumux_ss[run, 0, timestep[i], 2:]
    solute_dumux_g = solutes_dumux[run, 1, timestep[i], 2:]
    solute_sr_g = solutes_dumux_sr[run, 1, timestep[i], 2:]
    solute_ss_g = solutes_dumux_ss[run, 1, timestep[i], 2:]
    solute_ff_g = solutes_dumux_ff[run, 1, timestep[i], 2:]
    
    #left plot: no flux outer BC
    ax1[i,0].plot(CC, water_dumux_nf, "b", linestyle = linestyle_dumux, label = "water_dumux")
    ax1[i,0].plot(CC, water_steadyrate_nf, "b", linestyle = linestyle_steadyrate, label = "water_sr")
    ax2_0.plot(CC, solute_dumux_nf, "m", linestyle = linestyle_dumux, label = "solute_dumux")
    ax2_0.plot(CC, solute_ss_nf, "m", linestyle = linestyle_steadystate, label = "solute_ss")
    ax2_0.plot(CC, solute_sr_nf, "m", linestyle = linestyle_steadyrate, label = "solute_sr")
    
    #right plot: Dirichlet (initial conditions) outer BC
    ax1[i,1].plot(CC, water_dumux_g, "b", linestyle = linestyle_dumux, label = "water_dumux")
    ax1[i,1].plot(CC, water_steadyrate_g, "b", linestyle = linestyle_steadyrate, label = "water_sr")
    ax2_1.plot(CC, solute_dumux_g, "m", linestyle = linestyle_dumux, label = "solute_dumux")
    ax2_1.plot(CC, solute_ss_g, "m", linestyle = linestyle_steadystate, label = "solute_ss")
    ax2_1.plot(CC, solute_sr_g, "m", linestyle = linestyle_steadyrate, label = "solute_sr")
    ax2_1.plot(CC, solute_ff_g, "m", linestyle = linestyle_special, label = "solute_ff")
    #print(solute_ss_g)  
    #print(solute_sr_g)
    #print(solutes_dumux[run, 1, timestep[i], :])    
ax1[i,0].set_xlabel("distance root [cm]")
ax1[i,0].set_ylabel("water")
ax2_0.set_ylabel("nitrogen")
ax1[i,0].legend(["watercontent cm3/cm3"], loc="upper left")
ax2_0.legend(["nitrogen concentration mol/cm3"], loc="upper right")

ax1[i,0].legend(loc="upper left")
ax2_0.legend(loc="upper right")

ax1[i,1].set_xlabel("distance root [cm]")
ax1[i,1].set_ylabel("water")
ax2_1.set_ylabel("nitrogen")
ax1[i,1].legend(["watercontent cm3/cm3"], loc="upper left")
ax2_1.legend(["nitrogen concentration mol/cm3"], loc="upper right")

ax1[i,1].legend(loc="upper left")
ax2_1.legend(loc="upper right")

#np.save("input/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_)))) 
plt.show()
    


fig, ax1 = figure_style.subplots12(nrows=1, ncols=2)
# dumux(both)
# left: sr no flux, ss for the no flux
# right: ss, sr, farfield for sr waterflow

suptake_dumux_nf = solutes_dumux[run, 0, 1:, 0]
suptake_dumux = solutes_dumux[run, 1, 1:, 0]

solute_sr_nf = solutes_dumux_sr[run, 0, 1:, 0]
solute_ss_nf = solutes_dumux_ss[run, 0, 1:, 0]
solute_sr = solutes_dumux_sr[run, 1, 1:, 0]
solute_ss = solutes_dumux_ss[run, 1, 1:, 0]
solute_ff = solutes_dumux_ff[run, 1, 1:, 0]

#TODO: test segment here, remove
solute_sr_nf_conc = solutes_dumux_sr[run, 0, 1:, 2]
solute_ss_nf_conc = solutes_dumux_ss[run, 0, 1:, 2]
solute_sr_conc = solutes_dumux_sr[run, 1, 1:, 2]
solute_ss_conc = solutes_dumux_ss[run, 1, 1:, 2]
solute_ff_conc = solutes_dumux_ff[run, 1, 1:, 2]

Vmax = 4.0e-11 * 1 * (2*3.14*r_root) * (24*3600)  #mol/(cm2 s) * cm * cm * (s/d) -> mol / d 
Vmax_per_area = Vmax / (1 * (2*3.14*r_root)) #mol / d /cm2 = mol/(cm2d)
Km = 1.5e-7   #mol/cm3

for i in range(len(solute_sr_nf_conc)):
    solute_sr_nf[i]=-Vmax_per_area*solute_sr_nf_conc[i]/(Km+solute_sr_nf_conc[i])
    solute_ss_nf[i]=-Vmax_per_area*solute_ss_nf_conc[i]/(Km+solute_ss_nf_conc[i])
    solute_sr[i]=-Vmax_per_area*solute_sr_conc[i]/(Km+solute_sr_conc[i])
    solute_ss[i]=-Vmax_per_area*solute_ss_conc[i]/(Km+solute_ss_conc[i])
    solute_ff[i]=-Vmax_per_area*solute_ff_conc[i]/(Km+solute_ff_conc[i])

#print(suptake_dumux)
#print(solute_ss)
#print(solutes_dumux_ss[run, 1, 1:, 2])
ax1[0].plot(suptake_dumux_nf, suptake_dumux_nf, "m", linestyle = linestyle_dumux, label = "dumux")
ax1[0].plot(suptake_dumux_nf, abs(solute_sr_nf), "m", linestyle = linestyle_steadyrate, label = "steady rate nf")
ax1[0].plot(suptake_dumux_nf, abs(solute_ss_nf), "m", linestyle = linestyle_steadystate, label = "steady state nf")
ax1[0].scatter(suptake_dumux_nf, abs(solute_ss_nf), marker = "*")
ax1[1].plot(suptake_dumux, suptake_dumux, "m", linestyle = linestyle_dumux, label = "dumux")
ax1[1].scatter(suptake_dumux, abs(suptake_dumux), marker = "*")
ax1[1].plot(suptake_dumux, abs(solute_sr), "m", linestyle = linestyle_steadyrate, label = "steady rate")
ax1[1].scatter(suptake_dumux, abs(solute_sr), marker = "*")
ax1[1].plot(suptake_dumux, abs(solute_ss), "m", linestyle = linestyle_steadystate, label = "steady state")
ax1[1].scatter(suptake_dumux, abs(solute_ss), marker = "*")
ax1[1].plot(suptake_dumux, abs(solute_ff), "m", linestyle = linestyle_special, label = "far field approximation")
ax1[1].scatter(suptake_dumux, abs(solute_ff), marker = "*")

ax1[0].set_xlabel("dumux solute uptake")
ax1[0].set_ylabel("analytical approximation")
ax1[0].legend(["solute utake mol/cm2d"], loc="upper left")

ax1[1].set_xlabel("dumux solute uptake")
ax1[1].set_ylabel("analytical approximation")
ax1[1].legend(["solute utake mol/cm2d"], loc="upper left")


#np.save("input/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_)))) 
plt.show()
