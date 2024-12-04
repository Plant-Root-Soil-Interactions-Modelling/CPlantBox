#!/usr/bin/env python
# coding: utf-8


import os
import sys
sys.path.append("../.."); sys.path.append("../../src/")


# # Coupled carbon and water flow in CPlantBox (with a static soil)
# 
# ## Simulation of water and carbon movement 
# 
# 
# In the following we will show how to compute the coupled water and carbon flow in the plant. 
# 
# 
# 
# We consider a dynamic plant and a static soil. 
# To compute the carbon flux, we use the code developped by Lacointe et al. (2019).
# 
# **Reference**
# 
# A Lacointe and P. Minchin. A mechanistic model to predict distribution of carbon among multiple sinks. *Methods in molecular biology* (Clifton, N.J.) vol. 2014, 2019.
# 	

# The sucrose flow depends on several plant, soil and atmospheric variables. For clarity, the basic functions defining those variables were moved to the file "parametersSucroseFlow".



sys.path.append("../../modelparameter/functional")
import plantbox as pb
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import numpy as np
import visualisation.vtk_plot as vp # for quick 3d vizualisations
import matplotlib.pyplot as plt
from functional.phloem_flux import PhloemFluxPython  
from plant_photosynthesis.wheat_FcVB_Giraud2023adapted import *
from plant_hydraulics.wheat_Giraud2023adapted import *
from plant_sucrose.wheat_phloem_Giraud2023adapted import *
from climate.dummyWeather import *
import numpy as np


# ## 1. Define initial conditions




#we start with a small plant to have a lower computation time
simInit = 7 # [day] init simtime
simMax = 8
dt =2/24
depth = 60
weatherInit = weather(simInit)
simDuration = simInit


# plant system 
pl = pb.MappedPlant(seednum = 2) #seednum: gives the option of setting a random seed to make the simulations replicable
path = "../../modelparameter/structural/plant/"
name = "Triticum_aestivum_adapted_2023"
pl.readParameters(path + name + ".xml")

sdf = pb.SDF_PlantBox(np.inf, np.inf, depth )
pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


verbose = False
pl.initialize(verbose )
pl.simulate(simInit, verbose)


#for post-processing
Q_out = 0 #sucrose lost by the plant
AnSum = 0 #assimilation
filename = "phloemoutputs.txt" 

Q_Rmbu      = np.array([0.])
Q_Grbu      = np.array([0.])
Q_Exudbu    = np.array([0.])
Q_STbu    = np.array([0.])

Q_Rmall      = np.array([])
Q_Grall      = np.array([])
Q_Exudall    = np.array([])
lengthTotall    = np.array([])
time = np.array([])
lengthTotInit = sum(pl.segLength())
lengthTotBU = sum(pl.segLength())


# ## 2. Define static soil


min_b = [-3./2, -12./2, -61.]#distance between wheat plants
max_b = [3./2, 12./2, 0.]
cell_number = [6, 24, 61] #soil resolution
layers = depth; soilvolume = (depth / layers) * 3 * 12
k_soil = [] #conductivity of soil when in contact with roots
p_mean = weatherInit['p_mean'] #mean soil water potential
p_bot = p_mean + depth/2
p_top = p_mean - depth/2
sx = np.linspace(p_top, p_bot, depth) #soil water potential per voxel

picker = lambda x,y,z : max(int(np.floor(-z)),-1) 
pl.setSoilGrid(picker)  # maps segment


# ## 3. create object to compute carbon and water flux
# The PhloemFluxPython class containes the functionalities of PhotosynthesisPython as well as the sucrose-related functions.



#give initial guess of leaf water potential and internal CO2 partial pressure (to start computation loop)
r = PhloemFluxPython(pl,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5)


# ## 4. set other parameters and initial variable
# We present bellow some of the main sucrose-related parameters.



r = setPhotosynthesisParameters(r,weatherInit)

r = setKrKx_phloem(r) # conductivity of the sieve tube
r.setKrm2([[2e-5]]) #effect of the sucrose content on maintenance respiration 
r.setKrm1([[10e-2]]) #effect of structural sucrose content on maintenance respiration
r.setRhoSucrose([[0.51],[0.65],[0.56]])  #sucrose density per organ type (mmol/cm3)
r.setRmax_st([[14.4,9.0,6.0,14.4],[5.,5.],[15.]]) #maximum growth rate when water and carbon limitation is activated
r.KMfu = 0.11                                     #michaelis menten coefficient for usage of sucrose
r.beta_loading = 0.6 #feedback effect of sieve tube concentraiton on loading from mesophyll
r.Vmaxloading = 0.05 #mmol/d, max loading rate from mesophyll
r.Mloading = 0.2                                      #michaelis menten coefficient for loading of sucrose
r.Gr_Y = 0.8 # efficiency of sucrose usage for growth. if <1, we have growth respiration
r.CSTimin = 0.4 #minimum sucrose concentration below which no sucrose usage occures
r.Csoil = 1e-4 #mean soil concentration in sucrose


r.update_viscosity = True #update sucrose viscosity according to concentraiton ?
r.atol = 1e-12 #max absolute error for sucrose flow solver
r.rtol = 1e-8 #max relative error for sucrose flow solver


# ## 5. launch simulation
# In this simulation, we use the same time step for all the modules. The first time steps tend to require longer computation time. increasing the maximum errors allowed for the sucrose computation (r.atol, r.rtol) and the minium and maximum plant segment length (dxMin, dx) can help decrease the computaiton time.




while simDuration <= simMax: 
    
    Nt = len(r.plant.nodes) 
    weatherX = weather(simDuration) #update weather variables
    r.Qlight = weatherX["Qlight"] #
    r = setKrKx_xylem(weatherX["TairC"], weatherX["RH"], r) #update xylem conductivity data
    
    #compute plant water flow
    r.solve_photosynthesis(ea_=weatherX["ea"], es_ = weatherX["es"], 
                           sim_time_ = simDuration, sxx_=sx, cells_ = True,
                           verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"] )
    
    
    AnSum += np.sum(r.Ag4Phloem)*dt #total cumulative carbon assimilaiton
    errLeuning = sum(r.outputFlux) #should be 0 : no storage of water in the plant
    fluxes = np.array(r.outputFlux)
    fluxesSoil = r.soilFluxes(simDuration, r.psiXyl, sx, approx=False) #root water flux per soil voxel
    
    
    #simulation of phloem flow
    startphloem= simDuration
    endphloem = startphloem + dt
    stepphloem = 1
    r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , True, filename)
        
    #get ouput of sucrose flow computation    
    Q_ST    = np.array(r.Q_out[0:Nt])          #sieve tube sucrose content
    Q_meso  = np.array(r.Q_out[Nt:(Nt*2)])     #mesophyll sucrose content
    Q_Rm    = np.array(r.Q_out[(Nt*2):(Nt*3)]) #sucrose used for maintenance respiration
    Q_Exud  = np.array(r.Q_out[(Nt*3):(Nt*4)]) #sucrose used for exudation
    Q_Gr    = np.array(r.Q_out[(Nt*4):(Nt*5)]) #sucrose used for growth and growth respiration
    
    C_ST    = np.array(r.C_ST)           #sieve tube sucrose concentraiton
    volST   = np.array(r.vol_ST)         #sieve tube volume
    volMeso   = np.array(r.vol_Meso)      #mesophyll volume     
    C_meso  = Q_meso/volMeso              #sucrose concentration in mesophyll
    Q_out   = Q_Rm + Q_Exud + Q_Gr       #total sucrose lost/used by the plant
    error   = sum(Q_ST + Q_meso + Q_out )- AnSum  #balance residual (error)
    
    lengthTot = sum(r.plant.segLength()) #total plant length 
    
    #variation of sucrose content at the last time step (mmol)
    Q_ST_i        = Q_ST      - Q_STbu #in the sieve tubes
    Q_Rm_i        = Q_Rm      - Q_Rmbu #for maintenance
    Q_Gr_i        = Q_Gr      - Q_Grbu #for growth
    Q_Exud_i      = Q_Exud    - Q_Exudbu #for exudation
    Q_out_i       = Q_Rm_i    + Q_Exud_i      + Q_Gr_i #total usage
    
    #print some outputs
    print("\n\n\n\t\tat ", int(np.floor(simDuration)),"d", int((simDuration%1)*24),"h, PAR:",  round(r.Qlight *1e6),"mumol m-2 s-1")
    print("Error in sucrose balance:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(error, div0f(error,AnSum, 1.)))
    print("Error in water balance:\n\tabs (cm3/day) {:5.2e}".format(errLeuning))
    print("water fluxes (cm3/day):\n\ttranspiration {:5.2e}".format(sum(fluxesSoil.values())))
    print("assimilated sucrose (cm)\tAn {:5.2e}".format(AnSum)) 
    print("sucrose concentration in sieve tube (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e} at {:d} segs \tmax  {:5.2e}".format(np.mean(C_ST), min(C_ST), len(np.where(C_ST == min(C_ST) )[0]), max(C_ST)))        
    print('cumulated \tRm   {:.2e}\tGr   {:.2e}\tExud {:5.2e}'.format(sum(Q_Rm), sum(Q_Gr), sum(Q_Exud)))
    print("aggregated sink repartition at last time step (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm_i)/sum(Q_out_i)*100, 
         sum(Q_Gr_i)/sum(Q_out_i)*100,sum(Q_Exud_i)/sum(Q_out_i)*100))
    print("total aggregated sink repartition (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm)/sum(Q_out)*100, 
         sum(Q_Gr)/sum(Q_out)*100,sum(Q_Exud)/sum(Q_out)*100))
    print("growth rate (cm/day)\ttotal {:5.2e}\tlast time step {:5.2e}".format(lengthTot - lengthTotInit, lengthTot - lengthTotBU))      
    
    #plant growth based on Gr * Gr_Y
    r.plant.simulate(dt, verbose)
    simDuration += dt
    
    #for post processing
    Ntbu = Nt
    Nt = len(r.plant.nodes)
    lengthTotBU = lengthTot
    Q_STbu       =   np.concatenate((Q_ST, np.full(Nt - Ntbu, 0.)))
    Q_Rmbu       =   np.concatenate((Q_Rm, np.full(Nt - Ntbu, 0.)))
    Q_Grbu       =   np.concatenate((Q_Gr, np.full(Nt - Ntbu, 0.))) 
    Q_Exudbu     =   np.concatenate((Q_Exud, np.full(Nt - Ntbu, 0.))) 
    
    
    Q_Rmall    = np.append( Q_Rmall  ,sum(Q_Rm_i))
    Q_Grall    = np.append( Q_Grall  ,sum(Q_Gr_i))
    Q_Exudall  = np.append( Q_Exudall,sum(Q_Exud_i))
    lengthTotall  = np.append( lengthTotall,lengthTot)
    time       = np.append( time ,simDuration)
    


# ## 8. plot some results



fig, axs = plt.subplots(2,2)
axs[0,0].plot(time, Q_Rmall/dt)
axs[0,0].set(xlabel='day of growth', ylabel='total Rm rate (mmol/day)')
axs[1,0].plot(time, Q_Grall/dt, 'tab:red')
axs[1,0].set(xlabel='day of growth', ylabel='total Gr rate (mmol/day)')
axs[0,1].plot(time, Q_Exudall/dt , 'tab:brown')
axs[0,1].set(xlabel='day of growth', ylabel='total exudation\nrate (mmol/day)')
axs[1,1].plot(time, lengthTotall , 'tab:green')
axs[1,1].set(xlabel='day of growth', ylabel='total plant\nlength (cm)')
fig.tight_layout()
plt.show()


# ## Take away messages
# 
# * Basic idea how to use the class *PhloemFlow*
# * The plant growth follows the rate of carbon usage for growth
