
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
showuniform = True #should the uniform implementation be shown
showTiina = True #should Tiinas implementation be shown

# general parameters
max_time = 1.0  # d
n_times = 100+1 # number of time slots, -1 for intervals
r_prhiz = 0.6 # perirhizal radius[cm], computed for a RLD above 1cm/cm3
r_root = 0.02 # root radius [cm]
NC = 10 # number of spatial discretisations
length = 1 #default length of the segment, will not change the outcpme as all variables are assumed constant in this direction [cm]

#three scenarios will be computed: one without inflow, one with an advective flow (Dirichlet BC), another with a Dirichlet BC
n_scenarios = 3

#initial conditions
initial_waterpotential = -100
initial_soluteconcentration = 1.6e-6#mol/cm3, 103mg/L of NO3 (one of the Tereno measurements in 2015, TODO: look for another source) leads to slightly above 1.6

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
solutes_dumux = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
#analytical approximations
#steady state approximations
solutes_ss = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
#simple steady rate approximations
solutes_sr_simp = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
#steady rate approximations
solutes_sr = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
#general steady rate with far field approximation, not used here as there is no scenario for it
#solutes_ff = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
#dirichlet outer BC
solutes_d = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
#assume a uniform solute concentration
solutes_u = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
#Tiina Roose
solutes_TR = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 


#discretisation
lb = 0.5
points = np.logspace(np.log(r_root) / np.log(lb), np.log(r_prhiz) / np.log(lb), NC+1, base = lb) 
CC = np.array([(points[i] + points[i+1])/2 for i in range(NC)])
volumes = np.array([(points[i+1]**2 - points[i]**2)*3.14 for i in range(NC)])

soilVG = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam soil 
sp = vg.Parameters(soilVG)

def run_perirhizal_test(max_time, n_times, r_prhiz, r_root, NC, points, CC, volumes, length, n_scenarios, initial_waterpotential, initial_soluteconcentration):
    
    #space for the outputs
    watercontent_dumux = np.zeros((n_scenarios, n_times, NC+2)) #watercontent in cm3/cm3, (radial) uptake and inflow in cm/d 
    waterpotential_dumux = np.zeros((n_scenarios, n_times, NC+2)) #waterpotential in cm, (radial) uptake and inflow in cm/d 
    watercontent_sr = np.zeros((n_scenarios, n_times, NC+2)) #watercontent in cm3/cm3, (radial) uptake and inflow in cm/d 
    waterpotential_sr = np.zeros((n_scenarios, n_times, NC+2)) #waterpotential in cm, (radial) uptake and inflow in cm/d 
    solutes_dumux = np.zeros((n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, (radial) uptake and inflow in mol/(cm2d) 
    solutes_ss = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
    solutes_sr_simp = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
    solutes_sr = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
    solutes_ff = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
    solutes_d = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
    solutes_u = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
    solutes_TR = np.zeros((n_tests, n_scenarios, n_times, NC+2)) #solute concentration in mol/cm3, uptake and inflow in mol/(cm2d) 
    
    simtimes = np.linspace(0,max_time,n_times)#[1:]
    dt = max_time / n_times
    rho = r_prhiz / r_root
    
    #soil parameters
    soilVG = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam soil 
    
    # root conductivity and solute uptake parameters, constant throughout the entire simulation time
    molarMassWater = 18 #g/mol
    molarMassSolute = 62 #g/mol, NO3
    root_conductivity = 1e-4 #1/d
    inner_kr = root_conductivity * r_root * 2 * 3.14
    waterdemand = -0.05 #cm/d
    radial_waterdemand = 2*3.14*r_root * waterdemand #cm2/d
    Vmax = 4.0e-11 * (2*3.14*r_root) * (24*3600)  #mol/(cm2 s) * cm * cm * (s/d) -> mol / (cm d) 
    Vmax_per_area = Vmax / (2*3.14*r_root) #mol / (cm d) /cm2 = mol/(cm2d)
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
    s_af = RichardsWrapper(RichardsNCCylFoam()) #dirichlet BC for water, advective flow for solutes
    s_d = RichardsWrapper(RichardsNCCylFoam()) #Dirichlet outer BC

    for s in [s_nf, s_af, s_d]:
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
            if s == s_af:    
                s.setParameter("Soil.BC.Top.Type", "11") #dummy Dirichlet BC so dumux does not get mad
                s.setParameter("Soil.BC.Top.Value", str(initial_waterpotential))
                s.setParameter("Soil.BC.Top.SType", "10") #advective flow
                s.setParameter("Soil.BC.Top.CValue", str(0))
            else:
                s.setParameter("Soil.BC.Top.Type", "11") #dummy Dirichlet BC so dumux does not get mad
                s.setParameter("Soil.BC.Top.Value", str(initial_waterpotential))
                s.setParameter("Soil.BC.Top.SType", "11") #dummy Dirichlet BC so dumux does not get mad
                s.setParameter("Soil.BC.Top.CValue", str(initial_soluteconcentration*molarMassSolute)) 
        s.setParameter("Soil.BC.Bot.SType", "8")  # michaelisMenten=8 (SType = Solute Type)
        s.setParameter("Soil.BC.Bot.CValue", "0.0") #should not matter
        s.setParameter("RootSystem.Uptake.Vmax", str(Vmax_per_area*molarMassSolute)) #mol/(cm2d) -> g/(cm2 d)
        #s.setParameter("RootSystem.Uptake.Vmax", str(Vmax*molarMassSolute)) #mol/d -> g/d #TODO: Vmax or Vmax per area?
        s.setParameter("RootSystem.Uptake.Km", str(Km*molarMassSolute)) # mol/cm3 -> g/cm3
        s.setParameter("Soil.IC.C", str(initial_soluteconcentration*molarMassSolute))  # g / cm3
        s.setParameter("Component.MolarMass", str(molarMassWater/1000)) #g/mol -> kg/mol water
        s.setParameter("1.Component.MolarMass", str(molarMassSolute/1000)) #g/mol -> kg/mol nitrate
        s.setParameter("1.Component.LiquidDiffusionCoefficient", str(Ds / 1.e4 / (24*3600))) #cm^2/s -> m^2/s
        s.initializeProblem(maxDt = 0.01)
        s.ddt = 1.e-4  # days
    
    cellVolumes = s_d.getCellSurfacesCyl() * length # cm3
    
    #initial solute concentration for the analytical approximations
    #simplified steady rate no flux (1d lookup table)
    solutes_sr_simp = np.zeros((n_scenarios, n_times, NC+2))
    solutes_sr_simp[0,0,1]= initial_soluteconcentration
    solutes_sr_simp[1,0,1]= initial_soluteconcentration
    solutes_sr_simp[2,0,1]= initial_soluteconcentration
    
    #steady state solute flow
    solutes_ss = np.zeros((n_scenarios, n_times, NC+2))
    solutes_ss[1,0,1]= initial_soluteconcentration
    solutes_ss[2,0,1]= initial_soluteconcentration
    
    #steady rate no flux outer BC solute flow in the general water flow
    solutes_sr = np.zeros((n_scenarios, n_times, NC+2))
    solutes_sr[0,0,1]= initial_soluteconcentration
    solutes_sr[1,0,1]= initial_soluteconcentration
    solutes_sr[2,0,1]= initial_soluteconcentration
    
    #Dirichlet BC
    solutes_d = np.zeros((n_scenarios, n_times, NC+2))
    solutes_d[2,0,1]= initial_soluteconcentration
    
    #uniform concentration (= no analytical approximation)
    solutes_u = np.zeros((n_scenarios, n_times, NC+2))
    solutes_u[0,0,1]= initial_soluteconcentration
    solutes_u[1,0,1]= initial_soluteconcentration
    solutes_u[2,0,1]= initial_soluteconcentration
    
    #Tiina Roose approximation
    solutes_TR = np.zeros((n_scenarios, n_times, NC+2))
    
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
    rootuptake_s_nf = 0.0 #solute uptake of the root, mol / cm2d
    inflow_s_nf = 0.0 #solute uptake from outside the perirhizal zone, mol / cm2d
    
    #variables for the general case
    mean_watercontent_af = 0.1 #cm3/cm3
    mean_waterpotential_af = -100 #cm
    mean_solutecontent_af = 0 #mol/cm3
    Phi_soil_af = 1.0 #mfp of the mean soil
    Phi_root_af = 1.0 #mfp next to the root
    #mfp parameters for the perirhizal model
    Phi_outer_af = 1.0 #mfp at the outer perirhizal radius
    rootuptake_w_af = 1.0 #water uptake of the root, cm / d
    inflow_w_af = 0.0 #water uptake from outside the perirhizal zone, cm/d
    rootuptake_s_af = 0.0 #solute uptake of the root, mol / cm2d
    inflow_s_af = 0.0 #solute uptake from outside the perirhizal zone, mol / cm2d
    
    #variables for the Dirichlet case
    mean_watercontent_d = 0.1 #cm3/cm3
    mean_waterpotential_d = -100 #cm
    mean_solutecontent_d = 0 #mol/cm3
    Phi_soil_d = 1.0 #mfp of the mean soil
    Phi_root_d = 1.0 #mfp next to the root
    #mfp parameters for the perirhizal model
    Phi_outer_d = 1.0 #mfp at the outer perirhizal radius
    rootuptake_w_d = 1.0 #water uptake of the root, cm / d
    inflow_w_d = 0.0 #water uptake from outside the perirhizal zone, cm/d
    rootuptake_s_d = 0.0 #solute uptake of the root, mol / cm2d
    inflow_s_d = 0.0 #solute uptake from outside the perirhizal zone, mol / cm2d
    
    
    for r in range(1,len(simtimes)):
        
        dt = simtimes[r]
        if r>0:
            dt = simtimes[r] - simtimes[r-1]
        
        print('time',simtimes[r])
        print('no flux BC')
        print("*****", "#", r, "external time step", dt, " d, simulation time", s_nf.simTime, "d, internal time step", s_nf.ddt, "d")
        print('advective flow')
        print("*****", "#", r, "external time step", dt, " d, simulation time", s_af.simTime, "d, internal time step", s_af.ddt, "d")
        print('Dirichlet BC')
        print("*****", "#", r, "external time step", dt, " d, simulation time", s_d.simTime, "d, internal time step", s_d.ddt, "d")
        
        #one timestep
        s_nf.solve(dt, saveInnerFluxes_ = True)
        s_af.solve(dt, saveInnerFluxes_ = True)
        s_d.solve(dt, saveInnerFluxes_ = True) #TODO: this has numerical problems
        
        #watercontent and solute content, discretised
        watercontent_nf = s_nf.getWaterContent() # cm3
        waterpotential_nf = np.array([vg.pressure_head(watercontent_nf[i],peri.sp) for i in range(NC)]) # cm
        watercontent_af = s_af.getWaterContent() # cm3
        waterpotential_af = np.array([vg.pressure_head(watercontent_af[i],peri.sp) for i in range(NC)]) # cm
        watercontent_d = s_d.getWaterContent() # cm3
        waterpotential_d = np.array([vg.pressure_head(watercontent_d[i],peri.sp) for i in range(NC)]) # cm
        solutecontents_nf = s_nf.getSolution(1) / molarMassSolute # mol/cm3
        solutecontents_af = s_af.getSolution(1) / molarMassSolute # mol/cm3
        solutecontents_d = s_d.getSolution(1) / molarMassSolute # mol/cm3
            
        #inflow and outflow    
        rootuptake_w_nf = s_nf.getInnerFlow(0, length) /(length*(2*np.pi*r_root)) # cm /d
        rootuptake_s_nf = s_nf.getInnerFlow(1, length) / molarMassSolute /(length*(2*np.pi*r_root)) # mol / (cm2d)
        inflow_w_nf = s_nf.getOuterFlow(0, length) /(length*(2*np.pi*r_prhiz)) # cm /d
        inflow_s_nf = s_nf.getOuterFlow(1, length) / molarMassSolute /(length*(2*np.pi*r_prhiz)) # mol / (cm2d)
        
        rootuptake_w_af = s_af.getInnerFlow(0, length) /(length*(2*np.pi*r_root)) # cm /d
        rootuptake_s_af = s_af.getInnerFlow(1, length) / molarMassSolute /(length*(2*np.pi*r_root)) # mol / (cm2d)
        inflow_w_af = s_af.getOuterFlow(0, length) /(length*(2*np.pi*r_prhiz)) # cm /d
        inflow_s_af = s_af.getOuterFlow(1, length) / molarMassSolute /(length*(2*np.pi*r_prhiz)) # mol / (cm2d)
        
        rootuptake_w_d = s_d.getInnerFlow(0, length) /(length*(2*np.pi*r_root)) # cm /d
        rootuptake_s_d = s_d.getInnerFlow(1, length) / molarMassSolute /(length*(2*np.pi*r_root)) # mol / (cm2d)
        inflow_w_d = s_d.getOuterFlow(0, length) /(length*(2*np.pi*r_prhiz)) # cm /d
        inflow_s_d = s_d.getOuterFlow(1, length) / molarMassSolute /(length*(2*np.pi*r_prhiz)) # mol / (cm2d)
        
        rootuptake_w_nf = rootuptake_w_nf[0]
        rootuptake_s_nf = rootuptake_s_nf[0]
        inflow_w_nf = abs(inflow_w_nf[0])
        inflow_s_nf = inflow_s_nf[0]
        rootuptake_w_af = rootuptake_w_af[0]
        rootuptake_s_af = rootuptake_s_af[0]
        inflow_w_af = abs(inflow_w_af[0])
        inflow_s_af = inflow_s_af[0] 
        rootuptake_w_d = rootuptake_w_d[0]
        rootuptake_s_d = rootuptake_s_d[0]
        inflow_w_d = abs(inflow_w_d[0])
        inflow_s_d = inflow_s_d[0] 
        
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
        solutes_dumux[0,r,2:]=solutecontents_nf
        #advective flow scenario
        watercontent_dumux[1,r,0]=rootuptake_w_af
        watercontent_dumux[1,r,1]=inflow_w_af
        watercontent_dumux[1,r,2:]=watercontent_af
        waterpotential_dumux[1,r,0]=rootuptake_w_af
        waterpotential_dumux[1,r,1]=inflow_w_af
        waterpotential_dumux[1,r,2:]=waterpotential_af
        solutes_dumux[1,r,0]=rootuptake_s_af
        solutes_dumux[1,r,1]=inflow_s_af
        solutes_dumux[1,r,2:]=solutecontents_af
        #Dirichlet BC
        watercontent_dumux[2,r,0]=rootuptake_w_d
        watercontent_dumux[2,r,1]=inflow_w_d
        watercontent_dumux[2,r,2:]=watercontent_d
        waterpotential_dumux[2,r,0]=rootuptake_w_d
        waterpotential_dumux[2,r,1]=inflow_w_d
        waterpotential_dumux[2,r,2:]=waterpotential_d
        solutes_dumux[2,r,0]=rootuptake_s_d
        solutes_dumux[2,r,1]=inflow_s_d
        solutes_dumux[2,r,2:]=solutecontents_d
        
        #determine means 
        #for the explicit Euler method of the timesteps
        f_root = (2*np.pi*r_root)/(np.pi*(r_prhiz**2-r_root**2))
        f_prhiz = Ds* pow(0.25,13/3)/(0.4**2)*(2*np.pi*r_prhiz)/(np.pi*(r_prhiz**2-r_root**2)) #only have diffusion on the outside for now, assume water content of 0.25, saturated watercontent of 0.4
        IC = initial_soluteconcentration
        #no flux outer BC
        mean_watercontent_nf = np.average(watercontent_nf, weights=volumes)
        mean_waterpotential_nf = vg.pressure_head(mean_watercontent_nf, peri.sp)
        mean_soluteconcent_nf = np.average(solutecontents_nf, weights=np.multiply(mean_watercontent_nf, volumes))
        mean_soluteconcent_sr_simp_nf = solutes_sr_simp[0,r-1,1]+(f_root*solutes_sr_simp[0,r-1,0])*dt
        solutes_sr_simp[0,r,1] = mean_soluteconcent_sr_simp_nf
        #advective flow
        mean_watercontent_af = np.average(watercontent_af, weights=volumes)
        mean_waterpotential_af = vg.pressure_head(mean_watercontent_af, peri.sp)
        mean_soluteconcent_af = np.average(solutecontents_af, weights=np.multiply(mean_watercontent_af, volumes))
        mean_soluteconcent_sr_simp_af = solutes_sr_simp[1,r-1,1]+(f_root*solutes_sr_simp[1,r-1,0])*dt
        solutes_sr_simp[1,r,1] = mean_soluteconcent_sr_simp_af
        #Dirichlet BC
        mean_watercontent_d = np.average(watercontent_d, weights=volumes)
        mean_waterpotential_d = vg.pressure_head(mean_watercontent_d, peri.sp)
        mean_soluteconcent_d = np.average(solutecontents_d, weights=np.multiply(mean_watercontent_d, volumes))
        mean_soluteconcent_sr_simp_d = solutes_sr_simp[2,r-1,1]+(f_root*solutes_sr_simp[2,r-1,0]+f_prhiz*(IC-solutes_sr_simp[2,r-1,-1]))*dt
        solutes_sr_simp[2,r,1] = mean_soluteconcent_sr_simp_d
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
        #advective flow, Dirichlet for water transport
        Phi_soil_af = vg.fast_mfp[peri.sp](mean_waterpotential_af) #mfp of the mean soil
        Phi_outer_af = Phi_soil_af - (rootuptake_w_af*r_root-inflow_w_af*r_prhiz)*(rho**2)/(1-rho**2)*(((0.53)**2-1)/2-np.log(0.53))-(inflow_w_af*r_prhiz)*np.log(0.53) #mfp at the outer perirhizal radius
        Phi_af = lambda r: Phi_outer_af + (rootuptake_w_af*r_root-inflow_w_af*r_prhiz)*(rho**2)/(1-rho**2)*(((r/r_prhiz)**2-1)/2-np.log(r/r_prhiz))+(inflow_w_af*r_prhiz)*np.log(r/r_prhiz)#mfp function depending on radius
        Phi_root_af = Phi_af(r_root)#mfp next to the root 
        #Dirichlet BC
        Phi_soil_d = vg.fast_mfp[peri.sp](mean_waterpotential_d) #mfp of the mean soil
        Phi_outer_d = Phi_soil_d - (rootuptake_w_d*r_root-inflow_w_d*r_prhiz)*(rho**2)/(1-rho**2)*(((0.53)**2-1)/2-np.log(0.53))-(inflow_w_d*r_prhiz)*np.log(0.53) #mfp at the outer perirhizal radius
        Phi_d = lambda r: Phi_outer_d + (rootuptake_w_d*r_root-inflow_w_d*r_prhiz)*(rho**2)/(1-rho**2)*(((r/r_prhiz)**2-1)/2-np.log(r/r_prhiz))+(inflow_w_d*r_prhiz)*np.log(r/r_prhiz)#mfp function depending on radius
        Phi_root_d = Phi_d(r_root)#mfp next to the root 
    
        #write the steady rate approximations for the water as outputs
        waterpotential_sr[0,r,0] = rootuptake_w_nf
        waterpotential_sr[0,r,1] = inflow_w_nf
        waterpotential_sr[0,r,2:] = np.array([vg.fast_imfp[peri.sp](Phi_nf(CC[i])) for i in range(NC)])
        watercontent_sr[0,r,0] = rootuptake_w_nf
        watercontent_sr[0,r,1] = inflow_w_nf
        watercontent_sr[0,r,2:] = np.array([vg.water_content(waterpotential_sr[0,r,2+i],peri.sp) for i in range(NC)])
        
        waterpotential_sr[1,r,0] = rootuptake_w_af
        waterpotential_sr[1,r,1] = inflow_w_af
        waterpotential_sr[1,r,2:] = np.array([vg.fast_imfp[peri.sp](Phi_af(CC[i])) for i in range(NC)])
        watercontent_sr[1,r,0] = rootuptake_w_af
        watercontent_sr[1,r,1] = inflow_w_af
        watercontent_sr[1,r,2:] = np.array([vg.water_content(waterpotential_sr[1,r,2+i],peri.sp) for i in range(NC)])
        
        waterpotential_sr[2,r,0] = rootuptake_w_d
        waterpotential_sr[2,r,1] = inflow_w_d
        waterpotential_sr[2,r,2:] = np.array([vg.fast_imfp[peri.sp](Phi_d(CC[i])) for i in range(NC)])
        watercontent_sr[2,r,0] = rootuptake_w_d
        watercontent_sr[2,r,1] = inflow_w_d
        watercontent_sr[2,r,2:] = np.array([vg.water_content(waterpotential_sr[2,r,2+i],peri.sp) for i in range(NC)])
        
        
        
        #compute the analytical approximations for the solute uptake
        #case of dumux no flux outer BC
        result_solutes_sr_nf = peri.soil_root_solutes_steadyrate_simplified_([Phi_root_nf], [Phi_soil_nf], [r_root], [r_prhiz], [mean_soluteconcent_sr_simp_nf], [Vmax_per_area], [Km], Ds, [waterdemand], peri.sp, n_approx = 5)
        result_solutes_sr_af = peri.soil_root_solutes_steadyrate_simplified_([Phi_root_af], [Phi_soil_af], [r_root], [r_prhiz], [mean_soluteconcent_sr_simp_af], [Vmax_per_area], [Km], Ds, [waterdemand], peri.sp, n_approx = 5)
        result_solutes_sr_d = peri.soil_root_solutes_steadyrate_simplified_([Phi_root_d], [Phi_soil_d], [r_root], [r_prhiz], [mean_soluteconcent_sr_simp_d], [Vmax_per_area], [Km], Ds, [waterdemand], peri.sp, n_approx = 5)
        
        #safe the results
        #(here there is no computed inflow, so index 1 doesn't do anything)
        result_solutes_sr_nf = result_solutes_sr_nf[0]
        solutes_sr_simp[0,r,0]=-Vmax_per_area * result_solutes_sr_nf / (Km + result_solutes_sr_nf)
        result_solutes_sr_af = result_solutes_sr_af[0]
        solutes_sr_simp[1,r,0]=-Vmax_per_area * result_solutes_sr_af / (Km + result_solutes_sr_af)
        result_solutes_sr_d = result_solutes_sr_d[0]
        solutes_sr_simp[2,r,0]=-Vmax_per_area * result_solutes_sr_d / (Km + result_solutes_sr_d)
        
        F0_nf = peri.integral_AdvectionDiffusion_(Phi_nf(CC[0]),peri.sp)
        F0_af = peri.integral_AdvectionDiffusion_(Phi_af(CC[0]),peri.sp)
        F0_d = peri.integral_AdvectionDiffusion_(Phi_d(CC[0]),peri.sp)
        D_tilde = 1/Ds/math.pow(sp.theta_S-sp.theta_R,13/3)*(sp.theta_S*sp.theta_S)
        
        # the ratio of waterflow and soluteflow is assumed to remain constant throughout the perirhizal zone
        for j in range(NC):
            r_current = CC[j]
            F_nf = peri.integral_AdvectionDiffusion_(Phi_nf(r_current),peri.sp)-F0_nf
            F_af = peri.integral_AdvectionDiffusion_(Phi_af(r_current),peri.sp)-F0_af
            F_d = peri.integral_AdvectionDiffusion_(Phi_d(r_current),peri.sp)-F0_d
            F_tilde_nf=math.exp(-D_tilde*F_nf)
            F_tilde_af=math.exp(-D_tilde*F_af)
            F_tilde_d=math.exp(-D_tilde*F_d)
            solutes_sr_simp[0,r,2+j] = result_solutes_sr_nf / F_tilde_nf - (1-1/F_tilde_nf) * solutes_sr_simp[0,r,0] / (waterdemand)  #rewrite this using water solutes disc
            solutes_sr_simp[1,r,2+j] = result_solutes_sr_af / F_tilde_af - (1-1/F_tilde_af) * solutes_sr_simp[1,r,0] / (waterdemand) 
            solutes_sr_simp[2,r,2+j] = result_solutes_sr_d / F_tilde_d - (1-1/F_tilde_d) * solutes_sr_simp[2,r,0] / (waterdemand) 
        
        #Tiina Roose approximation
        E = 0 #minimal solute uptake
        DS_TR = Ds * math.pow(mean_watercontent_nf,10/3)/(sp.theta_S**2)
        waterflow_TR = 2*np.pi*r_root*waterdemand
        waterflow_TR = waterdemand
        rsc = peri.solutesuptake_convdiff_([mean_watercontent_nf], [initial_soluteconcentration], [Vmax_per_area], [Km], DS_TR, [waterflow_TR], [r_root], [E], [simtimes[r]], sp)
        solutes_TR[0,r,0]=-Vmax_per_area*rsc[0]/(Km + rsc[0])
        solutes_TR[0,r,2]=rsc[0]
        
        DS_TR = Ds * math.pow(mean_watercontent_af,10/3)/(sp.theta_S**2)
        rsc = peri.solutesuptake_convdiff_([mean_watercontent_af], [initial_soluteconcentration], [Vmax_per_area], [Km], DS_TR, [waterflow_TR], [r_root], [E], [simtimes[r]], sp)
        solutes_TR[1,r,0]=-Vmax_per_area*rsc[0]/(Km + rsc[0])
        solutes_TR[1,r,2]=rsc[0]
        
        DS_TR = Ds * math.pow(mean_watercontent_d,10/3)/(sp.theta_S**2)
        rsc = peri.solutesuptake_convdiff_([mean_watercontent_d], [initial_soluteconcentration], [Vmax_per_area], [Km], DS_TR, [waterflow_TR], [r_root], [E], [simtimes[r]], sp)
        solutes_TR[2,r,0]=-Vmax_per_area*rsc[0]/(Km + rsc[0])
        solutes_TR[2,r,2]=rsc[0]
        
        # case of general steady rate water uptake
        # safe the means, they are computed via the explicit Euler timestepping scheme
        
        mean_soluteconcent_ss_af = solutes_ss[1,r-1,1]+f_root*solutes_ss[1,r-1,0]*dt
        mean_soluteconcent_ss_d = solutes_ss[2,r-1,1]+(f_root*solutes_ss[2,r-1,0]+f_prhiz*(IC-solutes_ss[2,r-1,-1]))*dt
        mean_soluteconcent_sr_nf = solutes_sr[0,r-1,1]+f_root*solutes_sr[0,r-1,0]*dt
        mean_soluteconcent_sr_af = solutes_sr[1,r-1,1]+f_root*solutes_sr[1,r-1,0]*dt
        mean_soluteconcent_sr_d = solutes_sr[2,r-1,1]+f_root*solutes_sr[2,r-1,0]*dt
        mean_soluteconcent_d_d = solutes_d[2,r-1,1]+(f_root*solutes_d[2,r-1,0]+f_prhiz*(IC-solutes_d[2,r-1,-1]))*dt
        mean_soluteconcent_u_nf = solutes_u[0,r-1,1]+f_root*solutes_u[0,r-1,0]*dt
        mean_soluteconcent_u_af = solutes_u[1,r-1,1]+f_root*solutes_u[1,r-1,0]*dt
        mean_soluteconcent_u_d = solutes_u[2,r-1,1]+(f_root*solutes_u[2,r-1,0]+f_prhiz*(IC-solutes_u[2,r-1,-1]))*dt
        
        solutes_ss[1,r,1] = mean_soluteconcent_ss_af
        solutes_ss[2,r,1] = mean_soluteconcent_ss_d 
        solutes_sr[0,r,1] = mean_soluteconcent_sr_nf 
        solutes_sr[1,r,1] = mean_soluteconcent_sr_af 
        solutes_sr[2,r,1] = mean_soluteconcent_sr_d 
        solutes_d[2,r,1] = mean_soluteconcent_d_d 
        solutes_u[0,r,1] = mean_soluteconcent_u_nf
        solutes_u[1,r,1] = mean_soluteconcent_u_af 
        solutes_u[2,r,1] = mean_soluteconcent_u_d 
        
        #for the steady state take again the outer concentration
        #general steady state
        rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_af], [rootuptake_w_af*2*np.pi*r_root], [inflow_w_af*2*np.pi*r_prhiz], [CC[0]], [CC[-1]], [mean_soluteconcent_af], [initial_soluteconcentration], [Vmax], [Km], [Ds], peri.sp, mode = "ss")
        _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_af, 2*np.pi * rootuptake_w_af*r_root, 2*np.pi * inflow_w_af*r_prhiz, CC[0], CC[-1], CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        solutes_ss[1,r,0] = -(Uptake[0]+r_root**2*quadratic_flow[0]) / (2 * np.pi * r_root)
        #solutes_ss[1,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        solutes_ss[1,r,2:] = soluteconcentration[:]
        rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_d], [rootuptake_w_d*2*np.pi*r_root], [inflow_w_d*2*np.pi*r_prhiz], [CC[0]], [CC[-1]], [mean_soluteconcent_d], [initial_soluteconcentration], [Vmax], [Km], [Ds], peri.sp, mode = "ss")
        _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_d, 2*np.pi * rootuptake_w_d*r_root, 2*np.pi * inflow_w_d*r_prhiz, CC[0], CC[-1], CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        solutes_ss[2,r,0] = -(Uptake[0]+r_root**2*quadratic_flow[0]) / (2 * np.pi * r_root)
        #solutes_ss[2,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        solutes_ss[2,r,2:] = soluteconcentration[:]
        #general steady rate no flux outer BC
        rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_nf], [rootuptake_w_nf*2*np.pi*r_root], [inflow_w_nf*2*np.pi*r_prhiz], [CC[0]], [CC[-1]], [mean_soluteconcent_nf], [initial_soluteconcentration], [Vmax], [Km], [Ds], peri.sp, mode = "sr")
        _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_nf, 2*np.pi * rootuptake_w_nf*r_root, 2*np.pi * inflow_w_nf*r_prhiz, CC[0], CC[-1], CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        #print("rsc_sr_nf", rsc, Uptake, quadratic_flow, c_noflux)
        solutes_sr[0,r,0] = -(Uptake[0]+r_root**2*quadratic_flow[0]) / (2 * np.pi * r_root)
        #solutes_sr[0,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        solutes_sr[0,r,2:] = soluteconcentration[:]
        rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_af], [rootuptake_w_af*2*np.pi*r_root], [inflow_w_af*2*np.pi*r_prhiz], [CC[0]], [CC[-1]], [mean_soluteconcent_af], [initial_soluteconcentration], [Vmax], [Km], [Ds], peri.sp, mode = "sr")
        _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_af, 2*np.pi * rootuptake_w_af*r_root, 2*np.pi * inflow_w_af*r_prhiz, CC[0], CC[-1], CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        #print("rsc_sr_af", rsc, Uptake, quadratic_flow, c_noflux)
        solutes_sr[1,r,0] = -(Uptake[0]+r_root**2*quadratic_flow[0]) / (2 * np.pi * r_root)
        #solutes_sr[1,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        solutes_sr[1,r,2:] = soluteconcentration[:]
        rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_d], [rootuptake_w_d*2*np.pi*r_root], [inflow_w_d*2*np.pi*r_prhiz], [CC[0]], [CC[-1]], [mean_soluteconcent_d], [initial_soluteconcentration], [Vmax], [Km], [Ds], peri.sp, mode = "sr")
        _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_d, 2*np.pi * rootuptake_w_d*r_root, 2*np.pi * inflow_w_d*r_prhiz, CC[0], CC[-1], CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        #print("rsc_sr_d", rsc, Uptake, quadratic_flow, c_noflux)
        solutes_sr[2,r,0] = -(Uptake[0]+r_root**2*quadratic_flow[0]) / (2 * np.pi * r_root)
        #solutes_sr[2,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        solutes_sr[2,r,2:] = soluteconcentration[:]
        #Dirichlet BC
        rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_d], [rootuptake_w_d*2*np.pi*r_root], [inflow_w_d*2*np.pi*r_prhiz], [CC[0]], [CC[-1]], [mean_soluteconcent_d], [initial_soluteconcentration], [Vmax], [Km], [Ds], peri.sp, mode = "dirichlet")
        _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_d, 2*np.pi * rootuptake_w_d*r_root, 2*np.pi * inflow_w_d*r_prhiz, CC[0], CC[-1], CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        solutes_d[2,r,0] = -(Uptake[0]+r_root**2*quadratic_flow[0]) / (2 * np.pi * r_root)
        #solutes_d[2,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        solutes_d[2,r,2:] = soluteconcentration[:]
        #uniform concentration
        solutes_u[0,r,0] = -Vmax_per_area*mean_soluteconcent_u_nf/(Km + mean_soluteconcent_u_nf)
        solutes_u[1,r,0] = -Vmax_per_area*mean_soluteconcent_u_af/(Km + mean_soluteconcent_u_af)
        solutes_u[2,r,0] = -Vmax_per_area*mean_soluteconcent_u_d/(Km + mean_soluteconcent_u_d)
        solutes_u[0,r,2:] = np.ones(NC) * mean_soluteconcent_u_nf
        solutes_u[1,r,2:] = np.ones(NC) * mean_soluteconcent_u_af
        solutes_u[2,r,2:] = np.ones(NC) * mean_soluteconcent_u_d
        
        
    return watercontent_dumux, waterpotential_dumux, watercontent_sr, waterpotential_sr, solutes_dumux, solutes_sr_simp, solutes_ss, solutes_sr, solutes_d, solutes_u, solutes_TR

if do_computation:
    # save everything in the np arrays
    for i in range(n_tests):
        watercontent_dumux[i,:,:,:], waterpotential_dumux[i,:,:,:], watercontent_sr[i,:,:,:], waterpotential_sr[i,:,:,:], solutes_dumux[i,:,:,:], solutes_sr_simp[i,:,:,:], solutes_ss[i,:,:,:], solutes_sr[i,:,:,:], solutes_d[i,:,:,:], solutes_u[i,:,:,:], solutes_TR[i,:,:,:] = run_perirhizal_test(max_time, n_times, r_prhiz, r_root, NC, points, CC, volumes, length, n_scenarios, initial_waterpotential, initial_soluteconcentration)
    np.savez("test_perirhizal.npz", 
    watercontent_dumux=watercontent_dumux, 
    waterpotential_dumux=waterpotential_dumux, 
    watercontent_sr=watercontent_sr, 
    waterpotential_sr=waterpotential_sr, 
    solutes_dumux=solutes_dumux, 
    solutes_sr_simp=solutes_sr_simp, 
    solutes_ss=solutes_ss, 
    solutes_sr=solutes_sr,
    solutes_d=solutes_d,
    solutes_u=solutes_u,
    solutes_TR=solutes_TR)

else:
    simulation_results = np.load("test_perirhizal.npz")
    watercontent_dumux = simulation_results["watercontent_dumux"]
    waterpotential_dumux = simulation_results["waterpotential_dumux"]
    watercontent_sr = simulation_results["watercontent_sr"]
    waterpotential_sr = simulation_results["waterpotential_sr"]
    solutes_dumux = simulation_results["solutes_dumux"]
    solutes_sr_simp = simulation_results["solutes_sr_simp"]
    solutes_ss = simulation_results["solutes_ss"]
    solutes_sr = simulation_results["solutes_sr"]
    solutes_d = simulation_results["solutes_d"]
    solutes_u = simulation_results["solutes_u"]
    solutes_TR = simulation_results["solutes_TR"]
    
    


# compare both for the differint means of water / solute content
run = 0
timestep = np.array(np.linspace(1,9,num=5)) 
for i in range(5):
    timestep[i] = int(n_times * timestep[i] / 10)
timestep = timestep.astype(int)


linestyle_dumux = "solid"
linestyle_steadystate = "dotted"
linestyle_steadyrate = "dashed"
linestyle_special = "dashdot"
    

fig, ax1 = figure_style.subplots12(nrows=5, ncols=3)
# dumux(all 3)
# left: sr no flux
# middle: dirichlet BC water, advective flow solutes
# right: ss, sr, farfield for sr waterflow



for i in range(5):
    ax2_0 = ax1[i,0].twinx()
    ax2_1 = ax1[i,1].twinx()
    ax2_2 = ax1[i,2].twinx()
    
    #load data
    water_dumux_nf = waterpotential_dumux[run, 0, timestep[i], 2:]
    water_dumux_af = waterpotential_dumux[run, 1, timestep[i], 2:]
    water_dumux_d = waterpotential_dumux[run, 2, timestep[i], 2:]
    water_steadyrate_nf = waterpotential_sr[run, 0, timestep[i], 2:]
    water_steadyrate_af = waterpotential_sr[run, 1, timestep[i], 2:]
    water_steadyrate_d = waterpotential_sr[run, 2, timestep[i], 2:]
    solute_dumux_nf = solutes_dumux[run, 0, timestep[i], 2:]
    solute_dumux_af = solutes_dumux[run, 1, timestep[i], 2:]
    solute_dumux_d = solutes_dumux[run, 2, timestep[i], 2:]
    
    solutes_sr_simp_nf = solutes_sr_simp[run, 0, timestep[i], 2:]
    solutes_sr_simp_af = solutes_sr_simp[run, 1, timestep[i], 2:]
    solutes_sr_simp_d = solutes_sr_simp[run, 2, timestep[i], 2:]
    
    solutes_ss_af = solutes_ss[run, 1, timestep[i], 2:]
    solutes_ss_d = solutes_ss[run, 2, timestep[i], 2:]
    
    solutes_sr_nf = solutes_sr[run, 0, timestep[i], 2:]
    solutes_sr_af = solutes_sr[run, 1, timestep[i], 2:]
    solutes_sr_d = solutes_sr[run, 2, timestep[i], 2:]
    
    solutes_d_d = solutes_d[run, 2, timestep[i], 2:]
    
    solutes_u_nf = solutes_u[run, 0, timestep[i], 2:]
    solutes_u_af = solutes_u[run, 1, timestep[i], 2:]
    solutes_u_d = solutes_u[run, 2, timestep[i], 2:]
    
    
    #left plot: no flux outer BC
    ax1[i,0].plot(CC, water_dumux_nf, "b", linestyle = linestyle_dumux, label = "water_dumux")
    ax1[i,0].plot(CC, water_steadyrate_nf, "b", linestyle = linestyle_steadyrate, label = "water_sr")
    ax2_0.plot(CC, solute_dumux_nf, "m", linestyle = linestyle_dumux, label = "solute_dumux")
    ax2_0.plot(CC, solutes_sr_simp_nf, "m", linestyle = linestyle_steadyrate, label = "solute_sr_simp")
    ax2_0.plot(CC, solutes_sr_nf, "g", linestyle = linestyle_steadyrate, label = "solute_sr")
    ax2_0.plot(CC, solutes_u_nf, "r", linestyle = linestyle_dumux, label = "solute_u")
    
    #middle plot: water dirichlet, advective flow solutes
    ax1[i,1].plot(CC, water_dumux_af, "b", linestyle = linestyle_dumux, label = "water_dumux")
    ax1[i,1].plot(CC, water_steadyrate_af, "b", linestyle = linestyle_steadyrate, label = "water_sr")
    ax2_1.plot(CC, solute_dumux_af, "m", linestyle = linestyle_dumux, label = "solute_dumux")
    ax2_1.plot(CC, solutes_sr_simp_af, "m", linestyle = linestyle_steadyrate, label = "solute_sr_simp")
    ax2_1.plot(CC, solutes_sr_af, "g", linestyle = linestyle_steadyrate, label = "solute_sr")
    ax2_1.plot(CC, solutes_ss_af, "y", linestyle = linestyle_steadystate, label = "solute_ss")
    ax2_1.plot(CC, solutes_u_af, "r", linestyle = linestyle_dumux, label = "solute_u")
    
    #right plot: Dirichlet (initial conditions) outer BC
    ax1[i,2].plot(CC, water_dumux_d, "b", linestyle = linestyle_dumux, label = "water_dumux")
    ax1[i,2].plot(CC, water_steadyrate_d, "b", linestyle = linestyle_steadyrate, label = "water_sr")
    ax2_2.plot(CC, solute_dumux_d, "m", linestyle = linestyle_dumux, label = "solute_dumux")
    ax2_2.plot(CC, solutes_sr_simp_d, "m", linestyle = linestyle_steadyrate, label = "solute_sr_simp")
    ax2_2.plot(CC, solutes_sr_d, "g", linestyle = linestyle_steadyrate, label = "solute_sr")
    ax2_2.plot(CC, solutes_ss_d, "y", linestyle = linestyle_steadystate, label = "solute_ss")
    ax2_2.plot(CC, solutes_d_d, "c", linestyle = linestyle_special, label = "solute_d")
    ax2_2.plot(CC, solutes_u_d, "r", linestyle = linestyle_dumux, label = "solute_u")
 
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

ax1[i,2].set_xlabel("distance root [cm]")
ax1[i,2].set_ylabel("water")
ax2_2.set_ylabel("nitrogen")
ax1[i,2].legend(["watercontent cm3/cm3"], loc="upper left")
ax2_2.legend(["nitrogen concentration mol/cm3"], loc="upper right")

ax1[i,2].legend(loc="upper left")
ax2_2.legend(loc="upper right")

#np.save("input/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_)))) 
figure = plt.gcf() # get current figure
figure.set_size_inches(8, 6)
plt.savefig("concentrations", dpi = 100)
plt.show()
    


fig, ax1 = figure_style.subplots12(nrows=1, ncols=3)
# dumux(both)
# left: sr no flux, ss for the no flux
# right: ss, sr, farfield for sr waterflow

#load solute uptake
suptake_dumux_nf = solutes_dumux[run, 0, 1:, 0]
suptake_dumux_af = solutes_dumux[run, 1, 1:, 0]
suptake_dumux_d = solutes_dumux[run, 2, 1:, 0]

suptake_sr_simp_nf = solutes_sr_simp[run, 0, 1:, 0]
suptake_sr_simp_af = solutes_sr_simp[run, 1, 1:, 0]
suptake_sr_simp_d = solutes_sr_simp[run, 2, 1:, 0]

suptake_ss_af = solutes_ss[run, 1, 1:, 0]
suptake_ss_d = solutes_ss[run, 2, 1:, 0]

suptake_sr_nf = solutes_sr[run, 0, 1:, 0]
suptake_sr_af = solutes_sr[run, 1, 1:, 0]
suptake_sr_d = solutes_sr[run, 2, 1:, 0]

suptake_d_d = solutes_d[run, 2, 1:, 0]

suptake_u_nf = solutes_u[run, 0, 1:, 0]
suptake_u_af = solutes_u[run, 1, 1:, 0]
suptake_u_d = solutes_u[run, 2, 1:, 0]

suptake_TR_nf = solutes_TR[run, 0, 1:, 0]
suptake_TR_af = solutes_TR[run, 1, 1:, 0]
suptake_TR_d = solutes_TR[run, 2, 1:, 0]


ax1[0].plot(suptake_dumux_nf, suptake_dumux_nf, "m", linestyle = linestyle_dumux, label = "dumux")
#ax1[0].scatter(suptake_dumux_nf, abs(suptake_dumux_nf), "m", marker = "*")
ax1[0].plot(suptake_dumux_nf, abs(suptake_sr_simp_nf), "b", linestyle = linestyle_dumux, label = "steady rate simp")
ax1[0].plot(suptake_dumux_nf, abs(suptake_sr_nf), "g", linestyle = linestyle_dumux, label = "steady rate")
ax1[0].plot(suptake_dumux_nf, abs(suptake_u_nf), "r", linestyle = linestyle_dumux, label = "uniform")
ax1[0].plot(suptake_dumux_nf, abs(suptake_TR_nf), "k", linestyle = linestyle_dumux, label = "TR")

ax1[1].plot(suptake_dumux_af, suptake_dumux_af, "m", linestyle = linestyle_dumux, label = "dumux")
#ax1[1].scatter(suptake_dumux_af, abs(suptake_dumux_af), marker = "*")
ax1[1].plot(suptake_dumux_af, abs(suptake_sr_simp_af), "b", linestyle = linestyle_dumux, label = "steady rate simp")
ax1[1].plot(suptake_dumux_af, abs(suptake_sr_af), "g", linestyle = linestyle_dumux, label = "steady rate")
ax1[1].plot(suptake_dumux_af, abs(suptake_ss_af), "y", linestyle = linestyle_dumux, label = "steady state")
ax1[1].plot(suptake_dumux_af, abs(suptake_u_af), "r", linestyle = linestyle_dumux, label = "uniform")
ax1[1].plot(suptake_dumux_af, abs(suptake_TR_af), "k", linestyle = linestyle_dumux, label = "TR")

ax1[2].plot(suptake_dumux_d, suptake_dumux_d, "m", linestyle = linestyle_dumux, label = "dumux")
#ax1[2].scatter(suptake_dumux_d, abs(suptake_dumux_d), marker = "*")
ax1[2].plot(suptake_dumux_d, abs(suptake_sr_simp_d), "b", linestyle = linestyle_dumux, label = "steady rate simp")
ax1[2].plot(suptake_dumux_d, abs(suptake_sr_d), "g", linestyle = linestyle_dumux, label = "steady rate")
ax1[2].plot(suptake_dumux_d, abs(suptake_ss_d), "y", linestyle = linestyle_dumux, label = "steady state")
ax1[2].plot(suptake_dumux_d, abs(suptake_d_d), "c", linestyle = linestyle_dumux, label = "steady state")
ax1[2].plot(suptake_dumux_d, abs(suptake_u_d), "r", linestyle = linestyle_dumux, label = "uniform")
ax1[2].plot(suptake_dumux_d, abs(suptake_TR_d), "k", linestyle = linestyle_dumux, label = "TR")

ax1[0].set_xlabel("dumux solute uptake")
ax1[0].set_ylabel("analytical approximation")
ax1[0].legend(["solute utake mol/cm2d"], loc="upper left")

ax1[1].set_xlabel("dumux solute uptake")
ax1[1].set_ylabel("analytical approximation")
ax1[1].legend(["solute utake mol/cm2d"], loc="upper left")

ax1[2].set_xlabel("dumux solute uptake")
ax1[2].set_ylabel("analytical approximation")
ax1[2].legend(["solute utake mol/cm2d"], loc="upper left")

#np.save("input/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_))))
figure = plt.gcf() # get current figure
figure.set_size_inches(8, 6)
plt.savefig("relative_uptake", dpi = 100) 
plt.show()


fig, ax1 = figure_style.subplots12(nrows=1, ncols=3)
# dumux(both)
# left: sr no flux, ss for the no flux
# right: ss, sr, farfield for sr waterflow

simtimes = np.linspace(0,max_time,n_times)

#load solute uptake
suptake_dumux_nf = solutes_dumux[run, 0, 1:, 0]
#doublecheck if the uptake is correct
beginning_solutes = initial_soluteconcentration * vg.water_content(initial_waterpotential,sp) * (np.pi * (r_prhiz**2 - r_root**2) * length)
total_uptake_nf = [np.sum([solutes_dumux[0, 0, i, 2+j] *  watercontent_dumux[0, 0, i, 2+j] * np.pi * (points[j+1]**2-points[j]**2) * length for j in range(NC)]) for i in range(n_times)]
total_uptake_nf = np.array([(beginning_solutes - total_uptake_nf[i])/(2*np.pi*r_root*length)/max_time*(n_times-1) for i in range(n_times)])
suptake_dumux_af = solutes_dumux[run, 1, 1:, 0]
suptake_dumux_d = solutes_dumux[run, 2, 1:, 0]

suptake_sr_simp_nf = solutes_sr_simp[run, 0, 1:, 0]
suptake_sr_simp_af = solutes_sr_simp[run, 1, 1:, 0]
suptake_sr_simp_d = solutes_sr_simp[run, 2, 1:, 0]

suptake_ss_af = solutes_ss[run, 1, 1:, 0]
suptake_ss_d = solutes_ss[run, 2, 1:, 0]

suptake_sr_nf = solutes_sr[run, 0, 1:, 0]
suptake_sr_af = solutes_sr[run, 1, 1:, 0]
suptake_sr_d = solutes_sr[run, 2, 1:, 0]

suptake_d_d = solutes_d[run, 2, 1:, 0]

suptake_u_nf = solutes_u[run, 0, 1:, 0]
suptake_u_af = solutes_u[run, 1, 1:, 0]
suptake_u_d = solutes_u[run, 2, 1:, 0]

suptake_dumux_nf = np.array([ sum(suptake_dumux_nf[:i]) for i in range(n_times)])
suptake_dumux_af = np.array([ sum(suptake_dumux_af[:i]) for i in range(n_times)])
suptake_dumux_d = np.array([ sum(suptake_dumux_d[:i]) for i in range(n_times)])

suptake_sr_simp_nf = np.array([ sum(suptake_sr_simp_nf[:i]) for i in range(n_times)])
suptake_sr_simp_af = np.array([ sum(suptake_sr_simp_af[:i]) for i in range(n_times)])
suptake_sr_simp_d = np.array([ sum(suptake_sr_simp_d[:i]) for i in range(n_times)])

suptake_ss_af = np.array([ sum(suptake_ss_af[:i]) for i in range(n_times)])
suptake_ss_d = np.array([ sum(suptake_ss_d[:i]) for i in range(n_times)])

suptake_sr_nf = np.array([ sum(suptake_sr_nf[:i]) for i in range(n_times)])
suptake_sr_af = np.array([ sum(suptake_sr_af[:i]) for i in range(n_times)])
suptake_sr_d = np.array([ sum(suptake_sr_d[:i]) for i in range(n_times)])

suptake_d_d = np.array([ sum(suptake_d_d[:i]) for i in range(n_times)])

suptake_u_nf = np.array([ sum(suptake_u_nf[:i]) for i in range(n_times)])
suptake_u_af = np.array([ sum(suptake_u_af[:i]) for i in range(n_times)])
suptake_u_d = np.array([ sum(suptake_u_d[:i]) for i in range(n_times)])

suptake_TR_nf = np.array([ sum(suptake_TR_nf[:i]) for i in range(n_times)])
suptake_TR_af = np.array([ sum(suptake_TR_af[:i]) for i in range(n_times)])
suptake_TR_d = np.array([ sum(suptake_TR_d[:i]) for i in range(n_times)])



ax1[0].plot(simtimes, abs(suptake_dumux_nf), "m", linestyle = linestyle_dumux, label = "dumux")
#ax1[0].plot(simtimes[1:], abs(total_uptake_nf[1:]), "g", linestyle = linestyle_dumux, label = "dumux 2 uptake") #should produce the same line
ax1[0].plot(simtimes, abs(suptake_sr_simp_nf), "b", linestyle = linestyle_dumux, label = "steady rate simp")
ax1[0].plot(simtimes, abs(suptake_sr_nf), "g", linestyle = linestyle_dumux, label = "steady rate")


ax1[1].plot(simtimes, abs(suptake_dumux_af), "m", linestyle = linestyle_dumux, label = "dumux")
ax1[1].plot(simtimes, abs(suptake_sr_simp_af), "b", linestyle = linestyle_dumux, label = "steady rate simp")
ax1[1].plot(simtimes, abs(suptake_sr_af), "g", linestyle = linestyle_dumux, label = "steady rate")
ax1[1].plot(simtimes, abs(suptake_ss_af), "y", linestyle = linestyle_dumux, label = "steady state")

ax1[2].plot(simtimes, abs(suptake_dumux_d), "m", linestyle = linestyle_dumux, label = "dumux")
ax1[2].plot(simtimes, abs(suptake_sr_simp_d), "b", linestyle = linestyle_dumux, label = "steady rate simp")
ax1[2].plot(simtimes, abs(suptake_sr_d), "g", linestyle = linestyle_dumux, label = "steady rate nf BC")
ax1[2].plot(simtimes, abs(suptake_ss_d), "y", linestyle = linestyle_dumux, label = "steady state")
ax1[2].plot(simtimes, abs(suptake_d_d), "c", linestyle = linestyle_dumux, label = "steady rate D BC")

if showuniform:
    ax1[0].plot(simtimes, abs(suptake_u_nf), "r", linestyle = linestyle_dumux, label = "uniform")
    ax1[1].plot(simtimes, abs(suptake_u_af), "r", linestyle = linestyle_dumux, label = "uniform")
    ax1[2].plot(simtimes, abs(suptake_u_d), "r", linestyle = linestyle_dumux, label = "uniform")

if showTiina:
    ax1[0].plot(simtimes, abs(suptake_TR_nf), "k", linestyle = linestyle_dumux, label = "TR")    
    ax1[1].plot(simtimes, abs(suptake_TR_af), "k", linestyle = linestyle_dumux, label = "TR")
    ax1[2].plot(simtimes, abs(suptake_TR_d), "k", linestyle = linestyle_dumux, label = "TR")
#ax1[0].set_xlabel("dumux solute uptake")
#ax1[0].set_ylabel("analytical approximation")
#ax1[0].legend(["solute utake mol/cm2d"], loc="upper left")

#ax1[1].set_xlabel("dumux solute uptake")
#ax1[1].set_ylabel("analytical approximation")
#ax1[1].legend(["solute utake mol/cm2d"], loc="upper left")

#ax1[2].set_xlabel("dumux solute uptake")
#ax1[2].set_ylabel("analytical approximation")
#ax1[2].legend(["solute utake mol/cm2d"], loc="upper left")

ax1[0].legend(loc="upper left")
ax1[1].legend(loc="upper left")
ax1[2].legend(loc="upper left")

#np.save("input/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_)))) 
figure = plt.gcf() # get current figure
figure.set_size_inches(8, 6)
plt.savefig("cumulative_uptake", dpi = 100)
plt.show()
