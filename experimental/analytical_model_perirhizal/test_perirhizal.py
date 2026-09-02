
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


# dumux rosi imports
from rosi.richards_flat import RichardsFlatWrapper as RichardsWrapper  # Python part, macroscopic soil model
#from richards import RichardsWrapper  # Python part
from rosi.richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial
from rosi.rosi_richardsnc_cyl import RichardsNCCylFoam # C++ part (Dumux binding), macroscopic soil model


#numerics
import math
from scipy.optimize import fsolve, root_scalar
from numpy import linalg as LA


# run the dumux implementation of root water and nitrate uptake, later compare it to the analytical approximation


do_computation = True #should the computation be run or take the data from a saved file
showuniform = False #should the uniform implementation be shown
showTiina = True #should Tiinas implementation be shown

n_methods = 13 #how many different methods are used?
n_scenarios = 2 #two scenarios, no flux and Dirichlet BC
NC = 10 #number of discretisaitons of the axialsymmetric numerical model. More is a problem for dumux.
#0: mean, 1: root uptake (total radius), 2: influx into domain (total radius), 3 to (NC+2): discretisation


#methods:
wc_num = 0# numerical waterpotential
wp_num = 1# numerical watercontent
wp_ana_sr = 2# analytical steady rate waterpotential approximation
wc_ana_sr = 3# analytical steady rate watercontent approximation
s_num = 4# numerical solute
s_TR = 5# Tiina roose approximation solute
s_sr_simp = 6# simplified steady rate approximation
s_sr_nf = 7# steady rate no outer Flux approximation
s_ss = 8# steady state approximation
s_sr_d = 9# steady rate Dirichlet outer boundary condition approximation
s_u = 10# steady rate Dirichlet outer boundary condition approximation
s_sr = 11 # simplified steady rate no outer BC 
s_sr_lookup = 12# steady rate no flux outer BC with lookup

all_watercontents = [wc_num, wc_ana_sr]
all_waterpotentials = [wp_num, wp_ana_sr]
all_solutes = [s_num,s_TR,s_sr_simp,s_sr_nf,s_ss,s_sr_d,s_u,s_sr,s_sr_lookup]


#different scenarios
nf_idx = 0 #no flux is the first scenario
d_idx = 1 #Dirichlet B is the second scenario



# general parameters
#fraction_time = 100
max_time = 0.1 # d
n_times = 10+1 # number of time slots, -1 for intervals
r_prhiz = 0.6 # perirhizal radius[cm], computed for a RLD above 1cm/cm3
r_root = 0.02 # root radius [cm]
length = 1 #default length of the segment, will not change the outcpme as all variables are assumed constant in this direction [cm]
b = 0 #buffer power

bool_phosphate = False
if bool_phosphate:
    #ads = 0.05 L/(g,hr), kdes = 0.05 1/hr
    #b = rho_ Kd according to richards1p2cproblem.hh, but Kd = kdes/kads that does not make sense: increase in kdes should lower the amount of adsorption
    #try b = rho_b k_ads/k_des with a molar mass of 100g/mol, or density 1000g/L of water? probably the water
    #take rho_b = 1.5
    b = 0.16
f_b = 1
#initial conditions
initial_waterpotential = -300
initial_soluteconcentration = 1.6e-6#mol/cm3, 103mg/L of NO3 (one of the Tereno measurements in 2015, TODO: look for another source) leads to slightly above 1.6

#space for the outputs
values = np.zeros((n_methods, n_scenarios, n_times, NC+3))


#discretisation
lb = 0.5
points = np.logspace(np.log(r_root) / np.log(lb), np.log(r_prhiz) / np.log(lb), NC+1, base = lb) 
CC = np.array([(points[i] + points[i+1])/2 for i in range(NC)])
volumes = np.array([(points[i+1]**2 - points[i]**2)*3.14 for i in range(NC)])

soilVG = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam soil 
sp = vg.Parameters(soilVG)

initial_watercontent = vg.water_content(initial_waterpotential,sp) 

 
simtimes = np.linspace(0,max_time,n_times)#[1:]
dt = max_time / n_times
rho = r_prhiz / r_root

    
facto_w = 1#18*1000#?
    
# root conductivity and solute uptake parameters, constant throughout the entire simulation time
molarMassWater = 18 #g/mol
molarMassSolute = 62 #g/mol, NO3
root_conductivity = 1e-4 #1/d
inner_kr = root_conductivity * r_root * 2 * 3.14
waterdemand = -0.05/facto_w #cm/d
radial_waterdemand = 2*3.14*r_root * waterdemand #cm2/d
factor = 1#test
factor2 = 1/1000/100*2.512
Vmax = 4.0e-11 * (2*3.14*r_root) * (24*3600)*molarMassSolute /factor  #mol/(cm2 s) * cm * cm * (s/d) -> mol / (cm d) 
Vmax_per_area = Vmax / (2*3.14*r_root) #g / (cm d) /cm = g/(cm2d)
Vmax_per_area2 = 4.0e-11 * molarMassSolute * (24*3600) # g/(cm2d)
Km = 1.5e-7*molarMassSolute   #mol/cm3 -> g/cm3
    
total_volume = np.pi*(r_prhiz**2-r_root**2)*length    
    
initial_solutemass = initial_soluteconcentration * molarMassSolute
#diffusion coefficient of nitrate
Ds = 1.902e-5 * 24 * 3600  #cm2/s -> cm2/d
if bool_phosphate:
    Ds = 0.6e-5 * 24 * 3600  #cm2/s -> cm2/d
    

# load the perirhizal model
peri = PerirhizalPython() 
sp = vg.Parameters(soilVG)
peri.set_soil(sp)
    #no lookup tables are used here as there arent many simulations
    
    #peri.create_solute_sr_lookup("lookup_nitrate_sr", sp)
    #peri.open_lookup_solutes("lookup_nitrate_sr")

rho_b = 1.5 #bulk soil density

#alternative for lookup
peri2 = PerirhizalPython() 
sp = vg.Parameters(soilVG)
peri2.set_soil(sp)
#peri2.create_integralAdvectionDiffusion_lookup("processing/lookup_nitrate_simp", sp)
#peri2.open_lookup_solutes_simplified("processing/lookup_nitrate_simp")   
#peri2.create_solute_sr_lookup("processing/lookup_nitrate_sr", sp)
peri2.open_lookup_solutes("processing/lookup_nitrate_sr") 
    
Ds0 = peri.Ds0_ref 
scaling = Ds / Ds0#math.sqrt(Ds / Ds0)
print("scaling",scaling)


if do_computation:
    # initialise the dumux models for the scenarios
    s_nf = RichardsWrapper(RichardsNCCylFoam()) #no flux outer BC
    s_d = RichardsWrapper(RichardsNCCylFoam()) #Dirichlet outer BC
    
        #s_nf = RichardsWrapper(RichardsNCSP()) #no flux outer BC
    
    for s in [s_nf, s_d]:
        s.initialize()
        s.createGrid1d(points, length = length/100)  # [m] -> [cm]
        s.setVGParameters([soilVG])
        s.setHomogeneousIC(initial_waterpotential)  # cm pressure head
       
        s.setBotBC("constantFluxCyl",waterdemand*facto_w) # "noFlux")# Flux in cm/d
        if s == s_nf:
            s.setTopBC("constantFluxCyl",0.0)  #  [cm/day] "noFlux")#default, will be changed for one scenario
            s.setParameter("Soil.BC.Top.SType", "3")  # constantFluxCyl=3 (SType = Solute Type)
            s.setParameter("Soil.BC.Top.CValue", "0.0") 
        else:
            s.setParameter("Soil.BC.Top.Type", "11") #dummy Dirichlet BC so dumux does not get mad
            s.setParameter("Soil.BC.Top.Value", str(initial_waterpotential))
            s.setParameter("Soil.BC.Top.SType", "11") #dummy Dirichlet BC so dumux does not get mad
            s.setParameter("Soil.BC.Top.CValue", str(initial_soluteconcentration*molarMassSolute*initial_watercontent/(initial_watercontent+b*f_b))) 
        s.setParameter("Soil.BC.Bot.SType", "8")  # michaelisMenten=8 (SType = Solute Type)
        s.setParameter("Soil.BC.Bot.CValue", "0.0") #should not matter
        #s.setParameter("RootSystem.Uptake.Vmax", str(Vmax_per_area/1000/100*2.512*factor)) #mol/(cm2d) -> g/(cm2 d) #I have tested this scaling, yet I do not understand it #1/1000/100*2.512
        s.setParameter("RootSystem.Uptake.Vmax", str(Vmax_per_area2)) #mol/(cm2d) -> g/(cm2 d) #I have tested this scaling, yet I do not understand it #1/1000/100*2.512
        #s.setParameter("RootSystem.Uptake.Vmax", str(Vmax_per_area2)) #mol/(cm2d) -> g/(cm2 d) #I have tested this scaling, yet I do not understand it #1/1000/100*2.512
            #s.setParameter("RootSystem.Uptake.Vmax", str(Vmax*molarMassSolute)) #mol/d -> g/d #TODO: Vmax or Vmax per area?
        s.setParameter("RootSystem.Uptake.Km", str(Km)) # mol/cm3 -> g/cm3
        s.setParameter("Soil.IC.C", str(initial_soluteconcentration*molarMassSolute*initial_watercontent/(initial_watercontent+b*f_b)))  # g / L
        s.setParameter("Component.MolarMass", str(molarMassWater/1000)) #g/mol -> kg/mol water
        s.setParameter("1.Component.MolarMass", str(molarMassSolute/1000)) #g/mol -> kg/mol nitrate
        s.setParameter("1.Component.LiquidDiffusionCoefficient", str(Ds / 1.e4 / (24*3600) )) #cm^2/d -> m^2/s
        s.setParameter("Component.BufferPower", str(b*f_b))
        s.initializeProblem(maxDt = 0.01)
        s.ddt = 1.e-4  # days
        
        #s_nf.setParameter("Soil.BC.Bot.SType", s_nf.dumux_str(8))
        
    cellVolumes = s_d.getCellSurfacesCyl() * length # cm3
    
    #values[,,0,0] #beginning water potential
    #np.zeros(n_methods, n_scenarios, n_times, NC+3)
    
    values[all_waterpotentials,:,0,0] = np.ones((len(all_waterpotentials),n_scenarios))*initial_waterpotential
    values[all_waterpotentials,:,0,3:] = np.ones((len(all_waterpotentials),n_scenarios,NC))*initial_waterpotential
    values[all_watercontents,:,0,0] = np.ones((len(all_watercontents),n_scenarios))*initial_watercontent
    values[all_watercontents,:,0,3:] = np.ones((len(all_watercontents),n_scenarios,NC))*initial_watercontent
    values[all_solutes,:,0,0] = np.ones((len(all_solutes),n_scenarios))*initial_solutemass
    values[all_solutes,:,0,3:] = np.ones((len(all_solutes),n_scenarios,NC))*initial_solutemass
    values[s_ss,nf_idx,0,3:] = np.ones((NC))*initial_solutemass*initial_watercontent/(initial_watercontent+b) #test buffer
    
        
        
    for r in range(1,len(simtimes)):
        
        dt = simtimes[r]
        if r>0:
            dt = simtimes[r] - simtimes[r-1]
        
        #one timestep
        
        print('time',simtimes[r])
        print('no flux BC')
        print("*****", "#", r, "external time step", dt, " d, simulation time", s_nf.simTime, "d, internal time step", s_nf.ddt, "d")
        s_nf.solve(dt, saveInnerFluxes_ = True)
        #print('advective flow')
        #print("*****", "#", r, "external time step", dt, " d, simulation time", s_af.simTime, "d, internal time step", s_af.ddt, "d")
        #s_af.solve(dt, saveInnerFluxes_ = True)
        print('Dirichlet BC')
        print("*****", "#", r, "external time step", dt, " d, simulation time", s_d.simTime, "d, internal time step", s_d.ddt, "d")
        s_d.solve(dt, saveInnerFluxes_ = True) #TODO: this has numerical problems
        
        #watercontent and solute content, discretised
        watercontent_nf = s_nf.getWaterContent() # cm3
        mean_watercontent_nf = np.average(watercontent_nf, weights=volumes)
        mean_waterpotential_nf = vg.pressure_head(mean_watercontent_nf,peri.sp)
        
        values[wc_num,nf_idx,r,0] = np.array(mean_watercontent_nf) # cm3
        values[wp_num,nf_idx,r,0] = np.array(mean_waterpotential_nf) # cm3
        values[wc_num,nf_idx,r,3:] = np.array(watercontent_nf) # cm3
        values[wp_num,nf_idx,r,3:] = np.array([vg.pressure_head(watercontent_nf[i],peri.sp) for i in range(NC)]) # cm
        
        
        watercontent_d = s_d.getWaterContent() # cm3
        mean_watercontent_d = np.average(watercontent_d, weights=volumes)
        mean_waterpotential_d = vg.pressure_head(mean_watercontent_d,peri.sp)
        
        values[wc_num,d_idx,r,0] = np.array(mean_watercontent_d) # cm3
        values[wp_num,d_idx,r,0] = np.array(mean_waterpotential_d) # cm3
        values[wc_num,d_idx,r,3:] = np.array(watercontent_d) # cm3
        values[wp_num,d_idx,r,3:] = np.array([vg.pressure_head(watercontent_d[i],peri.sp) for i in range(NC)]) # cm
        
        solutecontents_nf = s_nf.getConcentration(1) # mol/cm3
        solutecontents_d = s_d.getConcentration(1)  # mol/cm3
        values[s_num,nf_idx,r,3:] = solutecontents_nf
        mean_solutecontent_nf = np.average(solutecontents_nf, weights=np.multiply(watercontent_nf,volumes))
        values[s_num,nf_idx,r,0] = mean_solutecontent_nf
        
        values[s_num,d_idx,r,3:] = solutecontents_d
        mean_solutecontent_d = np.average(solutecontents_d, weights=np.multiply(watercontent_d,volumes))
        values[s_num,d_idx,r,0] = mean_solutecontent_d
        
            
        #inflow and outflow    
        rootuptake_w_nf = s_nf.getInnerFlow(0, length) /(length*(2*np.pi*r_root)) # cm /d
        rootuptake_s_nf = s_nf.getInnerFlow(1, length) /(length*(2*np.pi*r_root)) # g / (cm2d) #/ molarMassSolute
        inflow_w_nf = s_nf.getOuterFlow(0, length) /(length*(2*np.pi*r_prhiz)) # cm /d
        inflow_s_nf = s_nf.getOuterFlow(1, length)  /(length*(2*np.pi*r_prhiz)) # mol / (cm2d) #/ molarMassSolute
        print("solutions",solutecontents_nf, Vmax,Vmax*solutecontents_nf[0]/(Km+solutecontents_nf[0]),rootuptake_s_nf[0])  
        
        rootuptake_w_d = s_d.getInnerFlow(0, length) /(length*(2*np.pi*r_root)) # cm /d
        rootuptake_s_d = s_d.getInnerFlow(1, length)  /(length*(2*np.pi*r_root)) # mol / (cm2d)#/ molarMassSolute
        inflow_w_d = s_d.getOuterFlow(0, length) /(length*(2*np.pi*r_prhiz)) # cm /d
        inflow_s_d = s_d.getOuterFlow(1, length)  /(length*(2*np.pi*r_prhiz)) # mol / (cm2d)#/ molarMassSolute
        print("rootwateruptake", rootuptake_w_d, inflow_w_d)
        rootuptake_w_nf = rootuptake_w_nf[0]
        rootuptake_s_nf = rootuptake_s_nf[0]
        inflow_w_nf = abs(inflow_w_nf[0])
        inflow_s_nf = inflow_s_nf[0]

        rootuptake_w_d = rootuptake_w_d[0]
        rootuptake_s_d = rootuptake_s_d[0]
        inflow_w_d = abs(inflow_w_d[0])
        inflow_s_d = inflow_s_d[0] 
        
        
        #no flux outer BC
        values[wc_num,nf_idx,r,1]=rootuptake_w_nf
        values[wc_num,nf_idx,r,2]=inflow_w_nf
        values[wp_num,nf_idx,r,1]=rootuptake_w_nf
        values[wp_num,nf_idx,r,2]=inflow_w_nf
        values[s_num,nf_idx,r,1]=rootuptake_s_nf
        values[s_num,nf_idx,r,2]=inflow_s_nf
        
        #Dirichlet BC
        values[wc_num,d_idx,r,1]=rootuptake_w_d
        values[wc_num,d_idx,r,2]=inflow_w_d
        values[wp_num,d_idx,r,1]=rootuptake_w_d
        values[wp_num,d_idx,r,2]=inflow_w_d
        values[s_num,d_idx,r,1]=rootuptake_s_d
        values[s_num,d_idx,r,2]=inflow_s_d
        
        #determine means 
        #for the explicit Euler method of the timesteps
        f_root = (2*np.pi*r_root)/(np.pi*(r_prhiz**2-r_root**2))
        f_prhiz = Ds* pow(0.25,13/3)/(0.4**2)*(2*np.pi*r_prhiz)/(np.pi*(r_prhiz**2-r_root**2)) #only have diffusion on the outside for now, assume water content of 0.25, saturated watercontent of 0.4
        IC = initial_soluteconcentration*molarMassSolute
        
        change_water_nf = 2*np.pi * (rootuptake_w_nf*r_root - inflow_w_nf* r_prhiz) / (np.pi * (r_prhiz**2 - r_root**2)) * dt
        #change_water_af = 2*np.pi * (waterdemand*r_root + inflow_w_af* r_prhiz) / (np.pi * (r_prhiz**2 - r_root**2)) * dt
        change_water_d = 2*np.pi * (rootuptake_w_d*r_root - inflow_w_d* r_prhiz) / (np.pi * (r_prhiz**2 - r_root**2)) * dt
        
        

        #equation [4] in Schroeder2008 doi:10.2136/vzj2007.0114
        #Phi(r)=Phi_outer + (q_root*r_root-q_out*r_prhiz)*(rho**2)/(1-rho**2)*(((r/r_prhiz)**2-1)/2-ln(r/r_prhiz))+q_out*r_prhiz*ln(r/r_prhiz)
        #Assume Phi(0.53r)=Phi(mean water potential)
        #both q_root and q_out are assumed to have positive signs if the water flows to the root
        #no flux outer BC
        Phi_soil_nf = vg.fast_mfp[peri.sp](mean_waterpotential_nf) #mfp of the mean soil
        Phi_outer_nf = Phi_soil_nf - (rootuptake_w_nf*r_root-inflow_w_nf*r_prhiz)*(rho**2)/(1-rho**2)*(((0.53)**2-1)/2-np.log(0.53))-(inflow_w_nf*r_prhiz)*r_prhiz*np.log(0.53) #mfp at the outer perirhizal radius
        Phi_nf = lambda r: Phi_outer_nf + (rootuptake_w_nf*r_root-inflow_w_nf*r_prhiz)*(rho**2)/(1-rho**2)*(((r/r_prhiz)**2-1)/2-np.log(r/r_prhiz))+(inflow_w_nf*r_prhiz)*r_prhiz*np.log(r/r_prhiz)#mfp function depending on radius
        Phi_root_nf = Phi_nf(r_root)#mfp next to the root 
        #Phi_0_nf, Phi_1_nf, Phi_2_nf, Phi_nf = peri.determine_mfp_function_influx(Phi_soil_nf, rootuptake_w_nf*2*np.pi*r_root, inflow_w_nf, rho)
        #Phi_outer_nf=Phi_0_nf
        #Dirichlet BC
        Phi_soil_d = vg.fast_mfp[peri.sp](mean_waterpotential_d) #mfp of the mean soil
        Phi_outer_d = Phi_soil_d - (rootuptake_w_d*r_root-inflow_w_d*r_prhiz)*(rho**2)/(1-rho**2)*(((0.53)**2-1)/2-np.log(0.53))-(inflow_w_d*r_prhiz)*np.log(0.53) #mfp at the outer perirhizal radius
        Phi_d = lambda r: Phi_outer_d + (rootuptake_w_d*r_root-inflow_w_d*r_prhiz)*(rho**2)/(1-rho**2)*(((r/r_prhiz)**2-1)/2-np.log(r/r_prhiz))+(inflow_w_d*r_prhiz)*np.log(r/r_prhiz)#mfp function depending on radius
        Phi_root_d = Phi_d(r_root)#mfp next to the root 
        #Phi_0_d, Phi_1_d, Phi_2_d, Phi_d = peri.determine_mfp_function_influx(Phi_soil_d, rootuptake_w_d*2*np.pi*r_root, inflow_w_d, rho)
        #Phi_outer_d=Phi_0_d
    
        
        

        def mean_watercontent_simp(Phi_root, Phi_soil, r_root, r_prhiz, CC, volumes):
            #this function determines the mean watercontent for a given steady rate parameterisation of the perirhizal zone
            #Phi_root and Phi_soil are the parameters corresponding to the matrix flux potential at r_root and 0.53*r_prhiz, CC is the discretisation
            
            #determine the "real" parameterisation
            Phi_A, Phi_C = peri.determine_mfp_function(Phi_root, Phi_soil, r_prhiz/r_root)
            Phi = lambda s : Phi_A * (s**2 - 2 * np.log(s)) + Phi_C
            
            #determine the water contents
            watercontents = np.array([vg.water_content(vg.fast_imfp[peri.sp](Phi(pos/r_prhiz)),peri.sp) for pos in CC])

            #determine the mean watercontent
            mean_watercontent = np.average(watercontents, weights=volumes)
            return mean_watercontent
        

        #compute the analytical approximations for the solute uptake
        #case of dumux no flux outer BC
        Phi_root_nf = vg.fast_mfp[peri.sp](values[wp_num,nf_idx,r,3])
        Phi_outer_nf = Phi_nf(r_prhiz)
        Phi_root_d = vg.fast_mfp[peri.sp](values[wp_num,d_idx,r,3])
        Phi_outer_d = Phi_d(r_prhiz)

        # determine Phi_soil such that the mean watercontent is accurate
        #case nf
        Phi_soil_min = Phi_nf(r_root)
        Phi_soil_max = Phi_outer_nf
        for ind_bis in range(10):
            Phi_soil_current = (Phi_soil_min+Phi_soil_max)/2
            test_watercontent = mean_watercontent_simp(Phi_root_nf, Phi_soil_current, r_root, r_prhiz, CC, volumes)
            if test_watercontent < mean_watercontent_nf:
                Phi_soil_min = Phi_soil_current
            else:
                Phi_soil_max = Phi_soil_current
        Phi_soil_nf = Phi_soil_current
        print("results bisection nf", Phi_nf(r_root), Phi_outer_nf, Phi_soil_nf)

        Phi_A_nf, Phi_C_nf = peri.determine_mfp_function(Phi_root_nf, Phi_soil_nf, r_prhiz/r_root)
        Phi_nf = lambda radius: Phi_A_nf * ((radius/r_prhiz)**2 - 2 * np.log(radius/r_prhiz)) + Phi_C_nf#mfp function depending on radius

        #case d
        Phi_soil_min = Phi_d(r_root)
        Phi_soil_max = Phi_outer_d
        for ind_bis in range(10):
            Phi_soil_current = (Phi_soil_min+Phi_soil_max)/2
            test_watercontent = mean_watercontent_simp(Phi_root_d, Phi_soil_current, r_root, r_prhiz, CC, volumes)
            if test_watercontent < mean_watercontent_d:
                Phi_soil_min = Phi_soil_current
            else:
                Phi_soil_max = Phi_soil_current
        Phi_soil_d = Phi_soil_current
        print("results bisection d", Phi_d(r_root), Phi_outer_d, Phi_soil_d)

        Phi_A_d, Phi_C_d = peri.determine_mfp_function(Phi_root_d, Phi_soil_d, r_prhiz/r_root)
        Phi_d = lambda radius: Phi_A_d * ((radius/r_prhiz)**2 - 2 * np.log(radius/r_prhiz)) + Phi_C_d#mfp function depending on radius


        values[wp_ana_sr,nf_idx,r,1]=rootuptake_w_nf
        values[wp_ana_sr,nf_idx,r,2]=inflow_w_nf
        values[wp_ana_sr,nf_idx,r,3:]=np.array([vg.fast_imfp[peri.sp](Phi_nf(CC[i])) for i in range(NC)])
        #print("sign inflow", inflow_w_nf)
        values[wc_ana_sr,nf_idx,r,1]=rootuptake_w_nf
        values[wc_ana_sr,nf_idx,r,2]=inflow_w_nf
        values[wc_ana_sr,nf_idx,r,3:]=np.array([vg.water_content(values[wp_ana_sr,nf_idx,r,3+i],peri.sp) for i in range(NC)])
        
        
        values[wp_ana_sr,d_idx,r,1]=rootuptake_w_d
        values[wp_ana_sr,d_idx,r,2]=inflow_w_d
        values[wp_ana_sr,d_idx,r,3:]=np.array([vg.fast_imfp[peri.sp](Phi_d(CC[i])) for i in range(NC)])
        
        values[wc_ana_sr,d_idx,r,1]=rootuptake_w_d
        values[wc_ana_sr,d_idx,r,2]=inflow_w_d
        values[wc_ana_sr,d_idx,r,3:]=np.array([vg.water_content(values[wp_ana_sr,nf_idx,r,3+i],peri.sp) for i in range(NC)])
        
        watercontent_ana_nf = values[wc_ana_sr,nf_idx,r,3:]
        mean_watercontent_ana_nf = np.average(watercontent_ana_nf, weights=volumes)
        mean_waterpotential_ana_nf = vg.pressure_head(mean_watercontent_ana_nf, peri.sp)

        watercontent_ana_d = values[wc_ana_sr,d_idx,r,3:]
        mean_watercontent_ana_d = np.average(watercontent_ana_d, weights=volumes)
        mean_waterpotential_ana_d = vg.pressure_head(mean_watercontent_ana_d, peri.sp)

        #no flux outer BC
        mean_waterpotential_ana_nf = vg.pressure_head(mean_watercontent_ana_nf, peri.sp)
        mean_soluteconcent_nf = np.average(solutecontents_nf, weights=np.multiply(watercontent_ana_nf, volumes))
        #mean_soluteconcent_sr_simp_nf = solutes_sr_simp[0,r-1,1]+(f_root*solutes_sr_simp[0,r-1,0])*dt/mean_watercontent_nf
        mean_soluteconcent_sr_simp_nf = values[s_sr_simp,nf_idx,r-1,0]+(f_root*values[s_sr_simp,nf_idx,r-1,1])*dt/mean_watercontent_ana_nf # the Dirichlet BC doesnt really work with the simplified sr as it always overestimates the outer concentration

        #solutes_sr_simp[0,r,1] = mean_soluteconcent_sr_simp_nf * mean_watercontent_nf/(mean_watercontent_nf+change_water_nf)
        mean_soluteconcent_sr_simp_nf = mean_soluteconcent_sr_simp_nf * mean_watercontent_ana_nf/(mean_watercontent_ana_nf+change_water_nf)
        values[s_sr_simp,nf_idx,r,0] = mean_soluteconcent_sr_simp_nf

        print("jump solute test", values[s_sr_simp,nf_idx,r-1,0], f_root, values[s_sr_simp,nf_idx,r-1,1], dt, mean_watercontent_ana_nf, change_water_nf, values[s_sr_simp,nf_idx,r,0])
        
        #Dirichlet BC
        mean_waterpotential_ana_d = vg.pressure_head(mean_watercontent_ana_d, peri.sp)
        mean_soluteconcent_d = np.average(solutecontents_d, weights=np.multiply(watercontent_ana_d, volumes))
        mean_soluteconcent_sr_simp_d = values[s_sr_simp,d_idx,r-1,0]+(f_root*values[s_sr_simp,d_idx,r-1,1]+f_prhiz*max(IC-values[s_sr_simp,d_idx,r-1,-1],0))*dt/mean_watercontent_ana_d # the Dirichlet BC doesnt really work with the simplified sr as it always overestimates the outer concentration
        # solutes_sr_simp[2,r,1] = mean_soluteconcent_sr_simp_d * mean_watercontent_d/(mean_watercontent_d+change_water_d)
        mean_soluteconcent_sr_simp_d = mean_soluteconcent_sr_simp_d * mean_watercontent_ana_d/(mean_watercontent_ana_d+change_water_d)
        values[s_sr_simp,d_idx,r,0] = mean_soluteconcent_sr_simp_d
        #determine coefficients for the analytical approximation
        
        print("mean watercontents", mean_watercontent_nf, mean_watercontent_d)

        #result_solutes_sr_nf = peri.soil_root_solutes_steadyrate_simplified_([vg.fast_mfp[peri.sp](values[wp_num,nf_idx,r,3])], [vg.fast_mfp[peri.sp](values[wp_num,nf_idx,r,-1])], [r_root], [r_prhiz], [mean_soluteconcent_sr_simp_nf], [-Vmax_per_area], [Km], Ds, [waterdemand], peri.sp, n_approx = 10)
        result_solutes_sr_simp_nf = peri.soil_root_solutes_steadyrate_simplified_([Phi_root_nf], [Phi_soil_nf], [r_root], [r_prhiz], [mean_soluteconcent_sr_simp_nf], [-Vmax_per_area], [Km], Ds, [waterdemand], peri.sp, n_approx = 5)
        #result_solutes_sr_nf = peri.soil_root_solutes_steadyrate_simplified_([Phi_root_nf], [Phi_outer_nf], [r_root], [r_prhiz], [mean_soluteconcent_sr_simp_nf], [-Vmax_per_area], [Km], Ds, [waterdemand], peri.sp, n_approx = 10)
        #result_solutes_sr_d = peri.soil_root_solutes_steadyrate_simplified_([vg.fast_mfp[peri.sp](values[wp_num,d_idx,r,3])], [vg.fast_mfp[peri.sp](values[wp_num,d_idx,r,-1])], [r_root], [r_prhiz], [mean_soluteconcent_sr_simp_d], [-Vmax_per_area], [Km], Ds, [waterdemand], peri.sp, n_approx = 10)
        result_solutes_sr_simp_d = peri.soil_root_solutes_steadyrate_simplified_([Phi_root_d], [Phi_soil_d], [r_root], [r_prhiz], [mean_soluteconcent_sr_simp_d], [-Vmax_per_area], [Km], Ds, [waterdemand], peri.sp, n_approx = 5)
        #result_solutes_sr_d = peri.soil_root_solutes_steadyrate_simplified_([Phi_root_d], [Phi_outer_d], [r_root], [r_prhiz], [mean_soluteconcent_sr_simp_d], [-Vmax_per_area], [Km], Ds, [waterdemand], peri.sp, n_approx = 10)
        
        print("special_simp",Phi_root_nf,Phi_soil_nf,Phi_root_d,Phi_soil_d)
        
        #safe the results
        #(here there is no computed inflow, so index 2 doesn't do anything)
        result_solutes_sr_simp_nf = result_solutes_sr_simp_nf[0]
        #solutes_sr_simp[0,r,0]=-Vmax_per_area * result_solutes_sr_nf / (Km + result_solutes_sr_nf)
        values[s_sr_simp,nf_idx,r,1]=-Vmax_per_area * result_solutes_sr_simp_nf / (Km + result_solutes_sr_simp_nf) 
        result_solutes_sr_simp_d = result_solutes_sr_simp_d[0]
        values[s_sr_simp,d_idx,r,1]=-Vmax_per_area * result_solutes_sr_simp_d / (Km + result_solutes_sr_simp_d)#TODO: track influx from the outside
        
        print("simp_sr_uptakes",mean_soluteconcent_sr_simp_nf,mean_soluteconcent_sr_simp_d,result_solutes_sr_simp_nf,result_solutes_sr_simp_d,values[s_sr_simp,nf_idx,r,1],values[s_sr_simp,d_idx,r,1])
        
        F0_nf = peri.integral_AdvectionDiffusion_(Phi_nf(r_root),peri.sp)
        #F0_af = peri.integral_AdvectionDiffusion_(Phi_af(r_root),peri.sp)
        F0_d = peri.integral_AdvectionDiffusion_(Phi_d(r_root),peri.sp)
        D_tilde = 1/Ds/math.pow(sp.theta_S-sp.theta_R,13/3)*(sp.theta_S*sp.theta_S)
        
        #print("Comparison waterfunction", Phi_soil_nf, Phi_nf(CC[-1]))
        
        # the ratio of waterflow and soluteflow is assumed to remain constant throughout the perirhizal zone
        for j in range(NC):
            r_current = CC[j]
            F_nf = peri.integral_AdvectionDiffusion_(Phi_nf(r_current),peri.sp)-F0_nf
            F_d = peri.integral_AdvectionDiffusion_(Phi_d(r_current),peri.sp)-F0_d
            F_tilde_nf=math.exp(-F_nf*D_tilde)
            F_tilde_d=math.exp(-F_d*D_tilde)
            print("test F", F_tilde_nf)
            #solutes_sr_simp[0,r,2+j] = result_solutes_sr_nf * F_tilde_nf + (1-F_tilde_nf) * solutes_sr_simp[0,r,0] / (waterdemand)  #rewrite this using water solutes disc
            values[s_sr_simp,nf_idx,r,3+j] = result_solutes_sr_simp_nf * F_tilde_nf + (1-F_tilde_nf) * values[s_sr_simp,nf_idx,r,1] / (waterdemand)  #rewrite this using water solutes disc
            #solutes_sr_simp[1,r,2+j] = result_solutes_sr_af * F_tilde_af + (1-F_tilde_af) * solutes_sr_simp[1,r,0] / (waterdemand) 
            #solutes_sr_simp[2,r,2+j] = result_solutes_sr_d * F_tilde_d + (1-F_tilde_d) * solutes_sr_simp[2,r,0] / (waterdemand) 
            values[s_sr_simp,d_idx,r,3+j] = result_solutes_sr_simp_d * F_tilde_d + (1-F_tilde_d) * values[s_sr_simp,d_idx,r,1] / (waterdemand) 
            #uptake = -Vmax_per_area * result_solutes_sr_nf / (result_solutes_sr_nf + Km)
            #print("outer c2", 1/F_tilde_nf, uptake, result_solutes_sr_nf, waterdemand, result_solutes_sr_nf /F_tilde_nf   - (1-1/F_tilde_nf) *  uptake / waterdemand)
        
        
        #Tiina Roose approximation
        E = 0 #minimal solute uptake
        DS_TR = Ds * mean_watercontent_nf#math.pow(mean_watercontent_nf,10/3)/(sp.theta_S**2)
        waterflow_TR = 2*np.pi*r_root*waterdemand
        waterflow_TR = waterdemand
        rsc = peri.solutesuptake_convdiff_([initial_watercontent], [initial_soluteconcentration*molarMassSolute], [Vmax_per_area], [Km], DS_TR, [waterflow_TR], [CC[0]], [E], [simtimes[r]], sp)
        #solutes_TR[0,r,0]=-Vmax_per_area*rsc[0]/(Km + rsc[0])
        values[s_TR,d_idx,r,1]=-Vmax_per_area*rsc[0]/(Km + rsc[0])
        #solutes_TR[0,r,2]=rsc[0]
        values[s_TR,nf_idx,r,3]=rsc[0]
        
        #DS_TR = Ds * mean_watercontent_af#math.pow(mean_watercontent_af,10/3)/(sp.theta_S**2)
        #rsc = peri.solutesuptake_convdiff_([initial_watercontent], [initial_soluteconcentration*molarMassSolute], [Vmax_per_area], [Km], DS_TR, [waterflow_TR], [CC[0]], [E], [simtimes[r]], sp)
        #solutes_TR[1,r,0]=-Vmax_per_area*rsc[0]/(Km + rsc[0])
        #solutes_TR[1,r,2]=rsc[0]
        
        DS_TR = Ds * mean_watercontent_d#math.pow(mean_watercontent_d,10/3)/(sp.theta_S**2)
        rsc = peri.solutesuptake_convdiff_([initial_watercontent], [initial_soluteconcentration*molarMassSolute], [Vmax_per_area], [Km], DS_TR, [waterflow_TR], [CC[0]], [E], [simtimes[r]], sp)
        #solutes_TR[2,r,0]=-Vmax_per_area*rsc[0]/(Km + rsc[0])
        values[s_TR,d_idx,r,1]=-Vmax_per_area*rsc[0]/(Km + rsc[0])
        #solutes_TR[2,r,2]=rsc[0]
        values[s_TR,d_idx,r,2]=rsc[0]
        
        # case of general steady rate water uptake
        # safe the means, they are computed via the explicit Euler timestepping scheme
        
        #change_water_nf = 2*np.pi * (waterdemand*r_root + inflow_w_nf* r_prhiz) / (np.pi * (r_prhiz**2 - r_root**2)) * dt
        #change_water_af = 2*np.pi * (waterdemand*r_root + inflow_w_af* r_prhiz) / (np.pi * (r_prhiz**2 - r_root**2)) * dt
        #change_water_d = 2*np.pi * (waterdemand*r_root + inflow_w_d* r_prhiz) / (np.pi * (r_prhiz**2 - r_root**2)) * dt
        
        
        #values[s_sr_simp,d_idx,r,1]
        #values[s_u,d_idx,r-1
        
        mean_soluteconcent_ss_nf = (values[s_ss,nf_idx,r-1,0]*mean_watercontent_nf+f_root*values[s_ss,nf_idx,r-1,1]*dt)/(mean_watercontent_nf+change_water_nf)
        #mean_soluteconcent_ss_af = (solutes_ss[1,r-1,1]*mean_watercontent_af+f_root*solutes_ss[1,r-1,0]*dt)/(mean_watercontent_af+change_water_af)
        mean_soluteconcent_ss_d = (values[s_ss,d_idx,r-1,0]*mean_watercontent_d+(f_root*values[s_ss,d_idx,r-1,1]+f_prhiz*(IC-values[s_ss,d_idx,r-1,-1]))*dt)/(mean_watercontent_d+change_water_d)
        mean_soluteconcent_sr_nf = (values[s_sr,nf_idx,r-1,0]*mean_watercontent_nf+f_root*values[s_sr,nf_idx,r-1,1]*dt)/(mean_watercontent_nf+change_water_nf)
        mean_soluteconcent_sr_lookup_nf = (values[s_sr_lookup,nf_idx,r-1,0]*mean_watercontent_nf+f_root*values[s_sr_lookup,nf_idx,r-1,1]*dt)/(mean_watercontent_nf+change_water_nf)
        #mean_soluteconcent_sr_af = (solutes_sr[1,r-1,1]*mean_watercontent_af+f_root*solutes_sr[1,r-1,0]*dt)/(mean_watercontent_af+change_water_af)
        mean_soluteconcent_sr_d = (values[s_sr,d_idx,r-1,0]*mean_watercontent_d+(f_root*values[s_sr,d_idx,r-1,1]+f_prhiz*(IC-values[s_sr,d_idx,r-1,-1]))*dt)/(mean_watercontent_d+change_water_d)
        mean_soluteconcent_sr_lookup_d = (values[s_sr_lookup,d_idx,r-1,0]*mean_watercontent_d+(f_root*values[s_sr_lookup,d_idx,r-1,1]+f_prhiz*(IC-values[s_sr_lookup,d_idx,r-1,-1]))*dt)/(mean_watercontent_d+change_water_d)
        mean_soluteconcent_d_d = (values[s_sr_d,d_idx,r-1,0]*mean_watercontent_d+(f_root*values[s_sr_d,d_idx,r-1,1]+f_prhiz*(IC-values[s_sr_d,d_idx,r-1,-1]))*dt)/(mean_watercontent_d+change_water_d)
        mean_soluteconcent_u_nf = (values[s_u,nf_idx,r-1,0]*mean_watercontent_nf+f_root*values[s_u,nf_idx,r-1,1]*dt)/(mean_watercontent_nf+change_water_nf)
        #mean_soluteconcent_u_af = (solutes_u[1,r-1,1]*mean_watercontent_af+f_root*solutes_u[1,r-1,0]*dt)/(mean_watercontent_af+change_water_af)
        mean_soluteconcent_u_d = (values[s_u,d_idx,r-1,0]*mean_watercontent_d+(f_root*values[s_u,d_idx,r-1,1]+f_prhiz*(IC-values[s_u,d_idx,r-1,-1]))*dt)/(mean_watercontent_d+change_water_d)
        
        print("compare mean", mean_soluteconcent_ss_nf, mean_soluteconcent_ss_d, change_water_nf, change_water_d)
        
        values[s_ss,nf_idx,r,0] = mean_soluteconcent_ss_nf
        #solutes_ss[1,r,1] = mean_soluteconcent_ss_af
        values[s_ss,d_idx,r,0] = mean_soluteconcent_ss_d 
        values[s_sr,nf_idx,r,0] = mean_soluteconcent_sr_nf 
        values[s_sr_lookup,nf_idx,r,0] = mean_soluteconcent_sr_lookup_nf 
        #solutes_sr[1,r,1] = mean_soluteconcent_sr_af 
        values[s_sr,d_idx,r,0] = mean_soluteconcent_sr_d 
        values[s_sr_lookup,d_idx,r,0] = mean_soluteconcent_sr_lookup_d 
        values[s_sr_d,d_idx,r,0] = mean_soluteconcent_d_d 
        values[s_u,nf_idx,r,0] = mean_soluteconcent_u_nf
        #solutes_u[1,r,1] = mean_soluteconcent_u_af 
        values[s_u,d_idx,r,0] = mean_soluteconcent_u_d 
        
        CC_tilde=CC.copy()
        CC_tilde[0:-1]=CC[1:]
        CC_tilde[-1]=r_prhiz
        print("perirhiz",r_prhiz,CC_tilde)
        
        #for the steady state take again the outer concentration
        #general steady state
        rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_nf], [rootuptake_w_nf*2*np.pi*r_root], [-inflow_w_nf*2*np.pi*r_prhiz], [r_root], [r_prhiz], [mean_soluteconcent_ss_nf*mean_watercontent_nf/(mean_watercontent_nf+b)], [0], [Vmax], [Km], [Ds], peri.sp, mode = "ss")
        _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_nf, 2*np.pi * rootuptake_w_nf*r_root, -2*np.pi * inflow_w_nf*r_prhiz, r_root , r_prhiz , CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        values[s_ss,nf_idx,r,1] = -(Uptake[0]+(r_root / scaling)**2*quadratic_flow[0]) / (2 * np.pi *r_root )
        #solutes_ss[0,r,0] = -Vmax_per_area*rsc[0] / (Km+rsc[0])
        #solutes_ss[0,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        values[s_ss,nf_idx,r,3:] = soluteconcentration[:]*(mean_watercontent_nf+b)/mean_watercontent_nf
        
        # rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_af], [-rootuptake_w_af*2*np.pi*r_root], [inflow_w_af*2*np.pi*r_prhiz], [r_root], [r_prhiz], [mean_soluteconcent_ss_af], [0], [Vmax], [Km], [Ds], peri.sp, mode = "ss")
        # _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_af, -2*np.pi * rootuptake_w_af*r_root, 2*np.pi * inflow_w_af*r_prhiz, r_root , r_prhiz , CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        # solutes_ss[1,r,0] = -(Uptake[0]+(r_root / scaling)**2*quadratic_flow[0]) / (2 * np.pi *r_root )
        # #solutes_ss[1,r,0] = -Vmax_per_area*rsc[0] / (Km+rsc[0])
        # #solutes_ss[1,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        # solutes_ss[1,r,2:] = soluteconcentration[:]
        # #solutes_ss[1,r,0] = -Vmax_per_area*soluteconcentration[0]/(Km + soluteconcentration[0])
        
        rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_d], [rootuptake_w_d*2*np.pi*r_root], [-inflow_w_d*2*np.pi*r_prhiz], [r_root], [r_prhiz], [mean_soluteconcent_ss_d*mean_watercontent_d/(mean_watercontent_d+b)], [0], [Vmax], [Km], [Ds], peri.sp, mode = "ss")
        _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_d, 2*np.pi * rootuptake_w_d*r_root, -2*np.pi * inflow_w_d*r_prhiz, r_root , r_prhiz , CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        values[s_ss,d_idx,r,1] = -(Uptake[0]+(r_root / scaling)**2*quadratic_flow[0]) / (2 * np.pi * r_root )
        #solutes_ss[2,r,0] = -Vmax_per_area*rsc[0] / (Km+rsc[0])
        #solutes_ss[2,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        values[s_ss,d_idx,r,3:] = soluteconcentration[:]*(mean_watercontent_d+b)/mean_watercontent_d
        
        #general steady rate no flux outer BC
        rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_nf], [rootuptake_w_nf*2*np.pi*r_root], [-inflow_w_nf*2*np.pi*r_prhiz], [r_root], [r_prhiz], [mean_soluteconcent_sr_nf], [0], [Vmax], [Km], [Ds], peri.sp, mode = "sr")
        _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_nf, 2*np.pi * rootuptake_w_nf*r_root, -2*np.pi * inflow_w_nf*r_prhiz, r_root , r_prhiz , CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        #print("rsc_sr_nf", rsc, Uptake, quadratic_flow, c_noflux)
        values[s_sr,nf_idx,r,1] = -(Uptake[0]+(r_root / scaling)**2*quadratic_flow[0]) / (2 * np.pi * r_root )
        #solutes_sr[0,r,0] = -Vmax_per_area*rsc[0] / (Km+rsc[0])
        #solutes_sr[0,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        values[s_sr,nf_idx,r,3:] = soluteconcentration[:]
        
        # rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_af], [-rootuptake_w_af*2*np.pi*r_root], [inflow_w_af*2*np.pi*r_prhiz], [r_root], [r_prhiz], [mean_soluteconcent_sr_af], [0], [Vmax], [Km], [Ds], peri.sp, mode = "sr")
        # _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_af, -2*np.pi * rootuptake_w_af*r_root, 2*np.pi * inflow_w_af*r_prhiz, r_root , r_prhiz , CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        # #print("rsc_sr_af", rsc, Uptake, quadratic_flow, c_noflux)
        # solutes_sr[1,r,0] = -(Uptake[0]+(r_root / scaling)**2*quadratic_flow[0]) / (2 * np.pi * r_root )
        # #solutes_sr[1,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        # solutes_sr[1,r,2:] = soluteconcentration[:]
        
        rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_d], [rootuptake_w_d*2*np.pi*r_root], [-inflow_w_d*2*np.pi*r_prhiz], [r_root], [r_prhiz], [mean_soluteconcent_sr_d], [0], [Vmax], [Km], [Ds], peri.sp, mode = "sr")
        _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_d, 2*np.pi * rootuptake_w_d*r_root, -2*np.pi * inflow_w_d*r_prhiz, r_root , r_prhiz , CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        #print("outerBC_sr", Uptake, quadratic_flow, r_prhiz, CC[-1], Uptake[0]+r_prhiz**2*quadratic_flow[0], mean_soluteconcent_sr_d, rsc)
        
        values[s_sr,d_idx,r,1] = -(Uptake[0]+(CC[0] / scaling)**2*quadratic_flow[0]) / (2 * np.pi * r_root )
        values[s_sr,d_idx,r,1] = -Vmax /(2 * np.pi * r_root) *rsc[0] / (Km+rsc[0])
        print("rsc_sr_d", rsc, Uptake, quadratic_flow, c_noflux, values[s_sr,d_idx,r,1], -(Uptake[0]+(r_root / scaling)**2*quadratic_flow[0])  / (2 * np.pi * r_root ))
        #solutes_sr[2,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        values[s_sr,d_idx,r,3:] = soluteconcentration[:]
        
        #general steady rate no flux outer BC with lookup table
        #rsc, Uptake, quadratic_flow, c_noflux = peri2.soil_root_solutes_sr([Phi_outer_nf], [-rootuptake_w_nf*2*np.pi*r_root], [inflow_w_nf*2*np.pi*r_prhiz], [r_root], [r_prhiz], [mean_soluteconcent_sr_lookup_nf], [0], [Vmax], [Km], [Ds], peri.sp, mode = "sr")
        #_, _, soluteconcentration, soluteconcentration_mean = peri2.watersolutes_disc(Phi_outer_nf, -2*np.pi * rootuptake_w_nf*r_root, 2*np.pi * inflow_w_nf*r_prhiz, r_root , r_prhiz , CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        #print("rsc_sr_nf", rsc, Uptake, quadratic_flow, c_noflux)
        values[s_sr_lookup,nf_idx,r,1] = -(Uptake[0]+(r_root / scaling)**2*quadratic_flow[0]) / (2 * np.pi * r_root )
        #solutes_sr[0,r,0] = -Vmax_per_area*rsc[0] / (Km+rsc[0])
        #solutes_sr[0,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        values[s_sr_lookup,nf_idx,r,3:] = soluteconcentration[:]
        
        # rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_af], [-rootuptake_w_af*2*np.pi*r_root], [inflow_w_af*2*np.pi*r_prhiz], [r_root], [r_prhiz], [mean_soluteconcent_sr_af], [0], [Vmax], [Km], [Ds], peri.sp, mode = "sr")
        # _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_af, -2*np.pi * rootuptake_w_af*r_root, 2*np.pi * inflow_w_af*r_prhiz, r_root , r_prhiz , CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        # #print("rsc_sr_af", rsc, Uptake, quadratic_flow, c_noflux)
        # solutes_sr[1,r,0] = -(Uptake[0]+(r_root / scaling)**2*quadratic_flow[0]) / (2 * np.pi * r_root )
        # #solutes_sr[1,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        # solutes_sr[1,r,2:] = soluteconcentration[:]
        
        #rsc, Uptake, quadratic_flow, c_noflux = peri2.soil_root_solutes_sr([Phi_outer_d], [-rootuptake_w_d*2*np.pi*r_root], [inflow_w_d*2*np.pi*r_prhiz], [r_root], [r_prhiz], [mean_soluteconcent_sr_lookup_d], [0], [Vmax], [Km], [Ds], peri.sp, mode = "sr")
        #_, _, soluteconcentration, soluteconcentration_mean = peri2.watersolutes_disc(Phi_outer_d, -2*np.pi * rootuptake_w_d*r_root, 2*np.pi * inflow_w_d*r_prhiz, r_root , r_prhiz , CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        #print("outerBC_sr", Uptake, quadratic_flow, r_prhiz, CC[-1], Uptake[0]+r_prhiz**2*quadratic_flow[0], mean_soluteconcent_sr_d, rsc)
        
        values[s_sr_lookup,d_idx,r,1] = -(Uptake[0]+(CC[0] / scaling)**2*quadratic_flow[0]) / (2 * np.pi * r_root )
        values[s_sr_lookup,d_idx,r,1] = -Vmax /(2 * np.pi * r_root) *rsc[0] / (Km+rsc[0])
        print("rsc_sr_d", rsc, Uptake, quadratic_flow, c_noflux, values[s_sr_lookup,d_idx,r,1], -(Uptake[0]+(r_root / scaling)**2*quadratic_flow[0])  / (2 * np.pi * r_root ))
        #solutes_sr[2,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        values[s_sr_lookup,d_idx,r,3:] = soluteconcentration[:]
        
        #Dirichlet BC
        rsc, Uptake, quadratic_flow, c_noflux = peri.soil_root_solutes_sr([Phi_outer_d], [rootuptake_w_d*2*np.pi*r_root], [-inflow_w_d*2*np.pi*r_prhiz], [r_root], [r_prhiz], [mean_soluteconcent_d_d], [initial_soluteconcentration*molarMassSolute], [Vmax], [Km], [Ds], peri.sp, mode = "dirichlet")
        _, _, soluteconcentration, soluteconcentration_mean = peri.watersolutes_disc(Phi_outer_d, 2*np.pi * rootuptake_w_d*r_root, -2*np.pi * inflow_w_d*r_prhiz, r_root , r_prhiz , CC, Ds, Uptake, quadratic_flow, c_noflux, peri.sp)
        values[s_sr_d,d_idx,r,1] = -(Uptake[0]+(r_root / scaling)**2*quadratic_flow[0]) / (2 * np.pi * r_root )
        #solutes_d[2,r,1] = -(Uptake[0] + r_prhiz**2 * quadratic_flow[0])
        values[s_sr_d,d_idx,r,3:] = soluteconcentration[:]
        #print("outerBC_d", Uptake, quadratic_flow, r_prhiz, CC[-1], Uptake[0]+r_prhiz**2*quadratic_flow[0], rsc)
        #uniform concentration
        values[s_u,nf_idx,r,1] = -Vmax_per_area*mean_soluteconcent_u_nf/(Km + mean_soluteconcent_u_nf)
        #solutes_u[1,r,0] = -Vmax_per_area*mean_soluteconcent_u_af/(Km + mean_soluteconcent_u_af)
        values[s_u,d_idx,r,1] = -Vmax_per_area*mean_soluteconcent_u_d/(Km + mean_soluteconcent_u_d)
        values[s_u,nf_idx,r,3:] = np.ones(NC) * mean_soluteconcent_u_nf
        #solutes_u[1,r,2:] = np.ones(NC) * mean_soluteconcent_u_af
        values[s_u,d_idx,r,3:] = np.ones(NC) * mean_soluteconcent_u_d
            
            
        #return watercontent_dumux, waterpotential_dumux, watercontent_sr, waterpotential_sr, solutes_dumux, solutes_sr_simp, solutes_ss, solutes_sr, solutes_d, solutes_u, solutes_TR

    np.savez("test_perirhizal.npz", values = values)
    

else:
    simulation_results = np.load("test_perirhizal.npz")
    values = simulation_results["values"]
    
    

# maybe for future use: look at the full radial uptake rather than per surface as that does not get impeded by changing radii



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
    

fig, ax1 = figure_style.subplots12(nrows=5, ncols=2)
# dumux(all 3)
# left: sr no flux
# middle: dirichlet BC water, advective flow solutes
# right: ss, sr, farfield for sr waterflow



for i in range(5):
    ax2_0 = ax1[i,0].twinx()
    ax2_1 = ax1[i,1].twinx()
    #ax2_2 = ax1[i,2].twinx()
    
    #load data
    
    water_dumux_nf = values[wp_num,nf_idx, timestep[i], 3:]
    water_dumux_d = values[wp_num,d_idx, timestep[i], 3:]
    water_steadyrate_nf = values[wp_ana_sr,nf_idx, timestep[i], 3:]
    water_steadyrate_d = values[wp_ana_sr,d_idx, timestep[i], 3:]
    solute_dumux_nf = values[s_num,nf_idx, timestep[i], 3:]
    solute_dumux_d = values[s_num,d_idx, timestep[i], 3:]
    
    solutes_sr_simp_nf = values[s_sr_simp,nf_idx, timestep[i], 3:]
    solutes_sr_simp_d = values[s_sr_simp,nf_idx, timestep[i], 3:]
    
    solutes_ss_nf = values[s_ss,nf_idx, timestep[i], 3:]
    solutes_ss_d = values[s_ss,d_idx, timestep[i], 3:]
    
    solutes_sr_nf = values[s_sr,nf_idx, timestep[i], 3:]
    solutes_sr_d = values[s_sr,d_idx, timestep[i], 3:]
    
    solutes_d_d = values[s_sr_d,d_idx, timestep[i], 3:]
    
    solutes_sr_lookup_nf = values[s_sr_lookup,nf_idx, timestep[i], 3:]
    solutes_sr_lookup_d = values[s_sr_lookup,d_idx, timestep[i], 3:]
    
    solutes_u_nf = values[s_u,nf_idx, timestep[i], 3:]
    solutes_u_d = values[s_u,d_idx, timestep[i], 3:]
    
    
    #left plot: no flux outer BC
    ax1[i,0].plot(CC, water_dumux_nf, "b", linestyle = linestyle_dumux, label = "water_dumux")
    ax1[i,0].plot(CC, water_steadyrate_nf, "b", linestyle = linestyle_steadyrate, label = "water_sr")
    ax2_0.plot(CC, solute_dumux_nf, "m", linestyle = linestyle_dumux, label = "solute_dumux")
    ax2_0.plot(CC, solutes_sr_simp_nf, "m", linestyle = linestyle_steadyrate, label = "solute_sr_simp")
    ax2_0.plot(CC, solutes_sr_nf, "g", linestyle = linestyle_steadyrate, label = "solute_sr")
    ax2_0.plot(CC, solutes_sr_lookup_nf, "k", linestyle = linestyle_steadyrate, label = "solute_sr_lookup")
    ax2_0.plot(CC, solutes_ss_nf, "y", linestyle = linestyle_steadystate, label = "solute_ss")
    ax2_0.plot(CC, solutes_u_nf, "r", linestyle = linestyle_dumux, label = "solute_u")
    
    
    #right plot: Dirichlet (initial conditions) outer BC
    ax1[i,1].plot(CC, water_dumux_d, "b", linestyle = linestyle_dumux, label = "water_dumux")
    ax1[i,1].plot(CC, water_steadyrate_d, "b", linestyle = linestyle_steadyrate, label = "water_sr")
    ax2_1.plot(CC, solute_dumux_d, "m", linestyle = linestyle_dumux, label = "solute_dumux")
    ax2_1.plot(CC, solutes_sr_simp_d, "m", linestyle = linestyle_steadyrate, label = "solute_sr_simp")
    ax2_1.plot(CC, solutes_sr_d, "g", linestyle = linestyle_steadyrate, label = "solute_sr")
    ax2_1.plot(CC, solutes_sr_lookup_d, "k", linestyle = linestyle_steadyrate, label = "solute_sr_lookup")
    ax2_1.plot(CC, solutes_ss_d, "y", linestyle = linestyle_steadystate, label = "solute_ss")
    ax2_1.plot(CC, solutes_d_d, "c", linestyle = linestyle_special, label = "solute_d")
    ax2_1.plot(CC, solutes_u_d, "r", linestyle = linestyle_dumux, label = "solute_u")
 
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

# ax1[i,2].set_xlabel("distance root [cm]")
# ax1[i,2].set_ylabel("water")
# ax2_2.set_ylabel("nitrogen")
# ax1[i,2].legend(["watercontent cm3/cm3"], loc="upper left")
# ax2_2.legend(["nitrogen concentration mol/cm3"], loc="upper right")

ax1[i,1].legend(loc="upper left")
ax2_1.legend(loc="upper right")

#np.save("input/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_)))) 
figure = plt.gcf() # get current figure
figure.set_size_inches(8, 6)
plt.savefig("concentrations", dpi = 100)
plt.show()
    


fig, ax1 = figure_style.subplots12(nrows=1, ncols=2)
# dumux(both)
# left: sr no flux, ss for the no flux
# right: ss, sr, farfield for sr waterflow

#load solute uptake
suptake_dumux_nf = values[s_num,nf_idx,1:,1]# solutes_dumux[run, 0, 1:, 0]
suptake_dumux_d = values[s_num,d_idx,1:,1]# solutes_dumux[run, 2, 1:, 0]

suptake_sr_simp_nf = values[s_sr_simp,nf_idx,1:,1]# solutes_sr_simp[run, 0, 1:, 0]
suptake_sr_simp_d = values[s_sr_simp,d_idx,1:,1]# solutes_sr_simp[run, 2, 1:, 0]

suptake_ss_nf = values[s_ss,nf_idx,1:,1]# solutes_ss[run, 0, 1:, 0]
suptake_ss_d = values[s_ss,d_idx,1:,1]# solutes_ss[run, 2, 1:, 0]

suptake_sr_nf = values[s_sr,nf_idx,1:,1]# solutes_sr[run, 0, 1:, 0]
suptake_sr_d = values[s_sr,d_idx,1:,1]# solutes_sr[run, 2, 1:, 0]

suptake_d_d = values[s_sr_d,d_idx,1:,1]# solutes_d[run, 2, 1:, 0]

suptake_u_nf = values[s_u,nf_idx,1:,1]# solutes_u[run, 0, 1:, 0]
suptake_u_d = values[s_u,d_idx,1:,1]# solutes_u[run, 2, 1:, 0]

suptake_TR_nf = values[s_TR,nf_idx,1:,1]# solutes_TR[run, 0, 1:, 0]
suptake_TR_d = values[s_TR,d_idx,1:,1]# solutes_TR[run, 2, 1:, 0]



ax1[0].plot(suptake_dumux_nf, suptake_dumux_nf, "m", linestyle = linestyle_dumux, label = "dumux")
#ax1[0].scatter(suptake_dumux_nf, abs(suptake_dumux_nf), "m", marker = "*")
ax1[0].plot(suptake_dumux_nf, abs(suptake_sr_simp_nf), "b", linestyle = linestyle_dumux, label = "steady rate simp")
ax1[0].plot(suptake_dumux_nf, abs(suptake_sr_nf), "g", linestyle = linestyle_dumux, label = "steady rate")
ax1[0].plot(suptake_dumux_nf, abs(suptake_ss_nf), "y", linestyle = linestyle_dumux, label = "steady state")
ax1[0].plot(suptake_dumux_nf, abs(suptake_u_nf), "r", linestyle = linestyle_dumux, label = "uniform")
ax1[0].plot(suptake_dumux_nf, abs(suptake_TR_nf), "k", linestyle = linestyle_dumux, label = "TR")

ax1[1].plot(suptake_dumux_d, suptake_dumux_d, "m", linestyle = linestyle_dumux, label = "dumux")
#ax1[2].scatter(suptake_dumux_d, abs(suptake_dumux_d), marker = "*")
ax1[1].plot(suptake_dumux_d, abs(suptake_sr_simp_d), "b", linestyle = linestyle_dumux, label = "steady rate simp")
ax1[1].plot(suptake_dumux_d, abs(suptake_sr_d), "g", linestyle = linestyle_dumux, label = "steady rate")
ax1[1].plot(suptake_dumux_d, abs(suptake_ss_d), "y", linestyle = linestyle_dumux, label = "steady state")
ax1[1].plot(suptake_dumux_d, abs(suptake_d_d), "c", linestyle = linestyle_dumux, label = "steady state")
ax1[1].plot(suptake_dumux_d, abs(suptake_u_d), "r", linestyle = linestyle_dumux, label = "uniform")
ax1[1].plot(suptake_dumux_d, abs(suptake_TR_d), "k", linestyle = linestyle_dumux, label = "TR")

ax1[0].set_xlabel("dumux solute uptake")
ax1[0].set_ylabel("analytical approximation")
ax1[0].legend(["solute utake mol/cm2d"], loc="upper left")

# ax1[1].set_xlabel("dumux solute uptake")
# ax1[1].set_ylabel("analytical approximation")
# ax1[1].legend(["solute utake mol/cm2d"], loc="upper left")

ax1[1].set_xlabel("dumux solute uptake")
ax1[1].set_ylabel("analytical approximation")
ax1[1].legend(["solute utake mol/cm2d"], loc="upper left")

#np.save("input/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_))))
figure = plt.gcf() # get current figure
figure.set_size_inches(8, 6)
plt.savefig("relative_uptake", dpi = 100) 
plt.show()


fig, ax1 = figure_style.subplots12(nrows=1, ncols=2)
# dumux(both)
# left: sr no flux, ss for the no flux
# right: ss, sr, farfield for sr waterflow

#determine mean watercontent
initial_watercontent = vg.water_content(initial_waterpotential,sp)
mean_watercontent_nf = np.average(values[wc_num,nf_idx, -1, 3:], weights=volumes)

simtimes = np.linspace(0,max_time,n_times)

#load solute uptake
suptake_dumux_nf = values[s_num,nf_idx,1:,1]# solutes_dumux[run, 0, 1:, 0]
#doublecheck if the uptake is correct
beginning_solutes = initial_soluteconcentration * vg.water_content(initial_waterpotential,sp) * (np.pi * (r_prhiz**2 - r_root**2) * length)
#total_uptake_nf = [np.sum([solutes_dumux[0, 0, i, 2+j] *  watercontent_dumux[0, 0, i, 2+j] * np.pi * (points[j+1]**2-points[j]**2) * length for j in range(NC)]) for i in range(n_times)]
total_uptake_nf = [np.sum([values[s_num,nf_idx, i, 3+j] *  values[wc_num,nf_idx, i, 3+j] * np.pi * (points[j+1]**2-points[j]**2) * length for j in range(NC)]) for i in range(n_times)]
total_uptake_nf = np.array([(beginning_solutes - total_uptake_nf[i])/(2*np.pi*r_root*length)/max_time*(n_times-1) for i in range(n_times)])
suptake_dumux_d = values[s_num, d_idx, 1:, 1]

suptake_sr_simp_nf = values[s_sr_simp, nf_idx, 1:, 1]# solutes_sr_simp[run, 0, 1:, 0]
suptake_sr_simp_d = values[s_sr_simp, d_idx, 1:, 1]# solutes_sr_simp[run, 2, 1:, 0]

suptake_ss_nf = values[s_ss, nf_idx, 1:, 1]# solutes_ss[run, 0, 1:, 0]
suptake_ss_d = values[s_ss, d_idx, 1:, 1]# solutes_ss[run, 2, 1:, 0]

suptake_sr_nf = values[s_sr, nf_idx, 1:, 1]# solutes_sr[run, 0, 1:, 0]
suptake_sr_d = values[s_sr, d_idx, 1:, 1]# solutes_sr[run, 2, 1:, 0]

suptake_d_d = values[s_sr_d, d_idx, 1:, 1]# solutes_d[run, 2, 1:, 0]

suptake_u_nf = values[s_u, nf_idx, 1:, 1]# solutes_u[run, 0, 1:, 0]
suptake_u_d = values[s_u, d_idx, 1:, 1]# solutes_u[run, 2, 1:, 0]

#print("watercontent",np.multiply(volumes,watercontent_dumux[run, 0, 0, 2:]))
#already set in the beginning?
#watercontent_dumux[run, 0, 0, 2:] = np.ones(len(watercontent_dumux[run, 0, 0, 2:]))*initial_watercontent
#mean_watercontents_nf = np.array([np.average(watercontent_dumux[run, 0, i, 2:], weights=volumes) for i in range(n_times)])
mean_watercontents_nf = np.array([np.average(values[wc_num,nf_idx, i, 3:], weights=volumes) for i in range(n_times)])
mean_watercontents_nf[0] = initial_watercontent
mean_watercontents_ana_nf = np.array([np.average(values[wc_ana_sr,nf_idx, i, 3:], weights=volumes) for i in range(n_times)])
mean_watercontents_ana_nf[0] = initial_watercontent
mean_solutecontents_nf = np.array([np.average(values[s_num,nf_idx, i, 3:], weights=np.multiply(volumes,values[wc_num,nf_idx, i, 3:])) for i in range(n_times)])
mean_solutecontents_nf[0] = initial_solutemass
#mean_solute_sr_nf = np.array([np.average(solutes_sr[run, 0, i, 2:], weights=np.multiply(volumes,values[wc_num,nf_idx, i, 3:])) for i in range(n_times)])
mean_solute_sr_simp_nf = np.array([np.average(values[s_sr_simp,nf_idx, i, 3:], weights=np.multiply(volumes,values[wc_ana_sr,nf_idx, i, 3:])) for i in range(n_times)])
mean_solute_sr_simp_nf[0] = initial_solutemass
mean_solute_sr_nf = np.array([np.average(values[s_sr,nf_idx, i, 3:], weights=np.multiply(volumes,values[wc_ana_sr,nf_idx, i, 3:])) for i in range(n_times)])
mean_solute_sr_nf[0] = initial_solutemass

suptake_dumux_nf = np.array([ abs(sum(suptake_dumux_nf[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt #+ np.multiply(mean_solutecontents_nf,mean_watercontents_nf) * np.pi*(r_prhiz**2-r_root**2)  #test passed
suptake_dumux_d = np.array([ abs(sum(suptake_dumux_d[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt

#suptake_sr_simp_nf = np.array([ abs(sum(suptake_sr_simp_nf[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt + np.multiply(mean_solute_sr_simp_nf,mean_watercontents_ana_nf) * np.pi*(r_prhiz**2-r_root**2) #test successful: constant, as it should be
suptake_sr_simp_nf = np.array([ abs(sum(suptake_sr_simp_nf[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt #+ np.multiply(values[s_sr_simp,nf_idx, :, 0],mean_watercontents_ana_nf) * np.pi*(r_prhiz**2-r_root**2) #test successful: constant, as it should be
suptake_sr_simp_d = np.array([ abs(sum(suptake_sr_simp_d[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt

print("where is the jump coming from?", np.array([ suptake_sr_simp_nf[i] for i in range(n_times)])*(2*np.pi*r_root)*dt, np.multiply(mean_solute_sr_simp_nf,mean_watercontents_ana_nf) * np.pi*(r_prhiz**2-r_root**2)) #jump is coming from the integral, i.e. the second part
print("where is the jump coming from2?", mean_solute_sr_simp_nf,mean_watercontents_ana_nf) #jump is coming from a sudden increase in solute at timestep 1
#that jump is coming from manually computing the mean solute concentration instead of using the precomputed one, i.e. the algorithm does not accurtely portray the solute content, may be because of the low discritisation?
print("How bad is the jump in the full calculation?", mean_solute_sr_nf)

suptake_ss_nf = np.array([ abs(sum(suptake_ss_nf[:i])) for i in range(n_times)]) *(2*np.pi*r_root)*dt
suptake_ss_d = np.array([ abs(sum(suptake_ss_d[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt

#suptake_sr_nf = np.array([ abs(sum(suptake_sr_nf[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt + np.multiply(mean_solute_sr_nf,mean_watercontents_ana_nf) * np.pi*(r_prhiz**2-r_root**2) #test passed
suptake_sr_nf = np.array([ abs(sum(suptake_sr_nf[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt #+ np.multiply(values[s_sr,nf_idx, :, 0],mean_watercontents_ana_nf) * np.pi*(r_prhiz**2-r_root**2) #test passed
suptake_sr_d = np.array([ abs(sum(suptake_sr_d[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt

suptake_d_d = np.array([ abs(sum(suptake_d_d[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt

suptake_u_nf = np.array([ abs(sum(suptake_u_nf[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt
suptake_u_d = np.array([ abs(sum(suptake_u_d[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt

suptake_TR_nf = np.array([ abs(sum(suptake_TR_nf[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt
suptake_TR_d = np.array([ abs(sum(suptake_TR_d[:i])) for i in range(n_times)])*(2*np.pi*r_root)*dt

#facto_w = 18*1000
#waterdemand = -0.05 / facto_w

#ax1[0].plot(simtimes, abs(mean_watercontents_nf)+abs(waterdemand)*np.array(range(n_times)) /(r_prhiz**2-r_root**2)*(2*r_root)*(n_times-1)/max_time, "c", linestyle = linestyle_dumux, label = "dumux_water")


ax1[0].plot(simtimes, abs(suptake_dumux_nf), "m", linestyle = linestyle_dumux, label = "dumux")
#ax1[0].plot(simtimes[1:], abs(total_uptake_nf[1:]), "g", linestyle = linestyle_dumux, label = "dumux 2 uptake") #should produce the same line
ax1[0].plot(simtimes, abs(suptake_sr_simp_nf), "b", linestyle = linestyle_dumux, label = "steady rate simp")
ax1[0].plot(simtimes, abs(suptake_sr_nf), "g", linestyle = linestyle_dumux, label = "steady rate")
ax1[0].plot(simtimes, abs(suptake_ss_nf), "y", linestyle = linestyle_dumux, label = "steady state")


ax1[1].plot(simtimes, abs(suptake_dumux_d), "m", linestyle = linestyle_dumux, label = "dumux")
ax1[1].plot(simtimes, abs(suptake_sr_simp_d), "b", linestyle = linestyle_dumux, label = "steady rate simp")
ax1[1].plot(simtimes, abs(suptake_sr_d), "g", linestyle = linestyle_dumux, label = "steady rate nf BC")
ax1[1].plot(simtimes, abs(suptake_ss_d), "y", linestyle = linestyle_dumux, label = "steady state")
ax1[1].plot(simtimes, abs(suptake_d_d), "c", linestyle = linestyle_dumux, label = "steady rate D BC")

if showuniform:
    ax1[0].plot(simtimes, abs(suptake_u_nf), "r", linestyle = linestyle_dumux, label = "uniform")
    ax1[1].plot(simtimes, abs(suptake_u_d), "r", linestyle = linestyle_dumux, label = "uniform")

if showTiina:
    ax1[0].plot(simtimes, abs(suptake_TR_nf), "k", linestyle = linestyle_dumux, label = "TR")    
    ax1[1].plot(simtimes, abs(suptake_TR_d), "k", linestyle = linestyle_dumux, label = "TR")


ax1[0].legend(loc="upper left")
ax1[1].legend(loc="upper left")

#np.save("input/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_)))) 
figure = plt.gcf() # get current figure
figure.set_size_inches(8, 6)
plt.savefig("cumulative_uptake", dpi = 100)
plt.show()



fig, ax1 = figure_style.subplots12(nrows=1, ncols=2)
# dumux(all 3)
# left: sr no flux
# middle: dirichlet BC water, advective flow solutes
# right: ss, sr, farfield for sr waterflow



for i in range(1):
    
    ax_0 = ax1[0]
    ax_1 = ax1[1]
    i = 4
    
    #load data
    
    water_dumux_nf = values[wp_num,nf_idx, timestep[i], 3:]
    water_dumux_d = values[wp_num,d_idx, timestep[i], 3:]
    water_steadyrate_nf = values[wp_ana_sr,nf_idx, timestep[i], 3:]
    water_steadyrate_d = values[wp_ana_sr,d_idx, timestep[i], 3:]
    solute_dumux_nf = values[s_num,nf_idx, timestep[i], 3:]
    solute_dumux_d = values[s_num,d_idx, timestep[i], 3:]
    
    solutes_sr_simp_nf = values[s_sr_simp,nf_idx, timestep[i], 3:]
    solutes_sr_simp_d = values[s_sr_simp,nf_idx, timestep[i], 3:]
    
    solutes_ss_nf = values[s_ss,nf_idx, timestep[i], 3:]
    solutes_ss_d = values[s_ss,d_idx, timestep[i], 3:]
    
    solutes_sr_nf = values[s_sr,nf_idx, timestep[i], 3:]
    solutes_sr_d = values[s_sr,d_idx, timestep[i], 3:]
    
    solutes_d_d = values[s_sr_d,d_idx, timestep[i], 3:]
    
    solutes_sr_lookup_nf = values[s_sr_lookup,nf_idx, timestep[i], 3:]
    solutes_sr_lookup_d = values[s_sr_lookup,d_idx, timestep[i], 3:]
    
    solutes_u_nf = values[s_u,nf_idx, timestep[i], 3:]
    solutes_u_d = values[s_u,d_idx, timestep[i], 3:]
    
    
    #left plot: no flux outer BC
    ax_0.plot(CC, water_dumux_nf, "b", linestyle = linestyle_dumux, label = "water_dumux")
    ax_0.plot(CC, water_steadyrate_nf, "b", linestyle = linestyle_steadyrate, label = "water_sr")
    ax_1.plot(CC, solute_dumux_nf, "m", linestyle = linestyle_dumux, label = "solute_dumux")
    ax_1.plot(CC, solutes_sr_simp_nf, "m", linestyle = linestyle_steadyrate, label = "solute_sr_simp")
    ax_1.plot(CC, solutes_sr_nf, "g", linestyle = linestyle_steadyrate, label = "solute_sr")
    ax_1.plot(CC, solutes_ss_nf, "y", linestyle = linestyle_steadyrate, label = "solute_ss")
    #ax_1.plot(CC, solutes_u_nf, "r", linestyle = linestyle_dumux, label = "solute_u")
    
    print("Compare steady rate and simplified", solutes_sr_nf, solutes_sr_simp_nf)
    
    #right plot: Dirichlet (initial conditions) outer BC
    #ax1[i,2].plot(CC, water_dumux_d, "b", linestyle = linestyle_dumux, label = "water_dumux")
    #ax1[i,2].plot(CC, water_steadyrate_d, "b", linestyle = linestyle_steadyrate, label = "water_sr")
    #ax2_2.plot(CC, solute_dumux_d, "m", linestyle = linestyle_dumux, label = "solute_dumux")
    #ax2_2.plot(CC, solutes_sr_simp_d, "m", linestyle = linestyle_steadyrate, label = "solute_sr_simp")
    #ax2_2.plot(CC, solutes_sr_d, "g", linestyle = linestyle_steadyrate, label = "solute_sr")
    #ax2_2.plot(CC, solutes_ss_d, "y", linestyle = linestyle_steadystate, label = "solute_ss")
    #ax2_2.plot(CC, solutes_d_d, "c", linestyle = linestyle_special, label = "solute_d")
    #ax2_2.plot(CC, solutes_u_d, "r", linestyle = linestyle_dumux, label = "solute_u")
 
ax_0.set_xlabel("distance root [cm]")
ax_0.set_ylabel("water")
ax_0.legend(["watercontent cm3/cm3"], loc="upper left")
ax_1.legend(["nitrogen concentration mol/cm3"], loc="upper right")

ax_0.legend(loc="upper left")
ax_1.legend(loc="upper right")



#np.save("input/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_)))) 
figure = plt.gcf() # get current figure
figure.set_size_inches(8, 6)
plt.savefig("concentrations", dpi = 100)
plt.show()

