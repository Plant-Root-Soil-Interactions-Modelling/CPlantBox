# This example reproduces the infiltration example from Schnepf et al. (2023, doi.org/10.1093/insilicoplants/diad005). 
# Water infiltrates into an initially dry soil from the soil surface. Only the vertical water movement is 
# considered. A constant Neumann boundary condition is set at the upper boundary and a free drainage 
# boundary condition is set at the lower boundary. There is an analytical solution for this simple example, 
# which can optionally be plotted for comparison. 


# This example solves the Richards equation with DuMux. The github repository "dumux-rosi" (https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi.git) needs to be cloned
# into the same path level as CPlantBox. 


# add paths to the folders containing CPlantBox and dumux-rosi
import sys; #|\label{l61:paths_a}|
sys.path.append("../../../dumux-rosi/python/modules");
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox/src")   #|\label{l61:paths_e}|

import matplotlib.pyplot as plt   #|\label{l61:bibs1_a}|
import numpy as np  # |\label{l61:bibs1_e}|

from rosi_richards import RichardsSPnum  # C++ part (Dumux binding) |\label{l61:bibs2_a}|
from richards import RichardsWrapper  # Python part  |\label{l61:bibs2_e}|
#from analytic_solution import *  # plots the analytical solutions to ax1, ax2, ax3 # optional |\label{l61:bib3}|

# Define van Genuchten parameters for sand, loam and clay  |\label{l61:genuchten_a}|
# theta_r (-), theta_s (-), alpha (1/cm), n (-), Ks (cm d-1)   
sand = [0.045, 0.43, 0.15, 3, 1000] 
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil=loam   # Select soil type for simulation   |\label{l61:genuchten_e}|
# simulation time, days 
sim_time=1;   #|\label{l61:simtime}|
dt = 720 / (24 * 3600)  # time step [days]    |\label{l61:timestep}|

# Solve the Richards equation using the Python wrapper of dumux-rosi 
s = RichardsWrapper(RichardsSPnum())   # |\label{l61:initialize_a}|
s.initialize()   # |\label{l61:initialize_e}|
s.setTopBC("atmospheric", 0.5, [[-1., 1.e10], [100., 100.]])  #  [cm/day] atmospheric is with surface run-off   |\label{l61:top_bc}|
s.setBotBC("freeDrainage")   # |\label{l61:bottom_bc}|
N = 199
s.createGrid([-5., -5., -200.], [5., 5., 0.], [1, 1, N])  # [cm] N   |\label{l61:grid}|
s.setHomogeneousIC(-400.)  # cm pressure head    |\label{l61:ic}|
s.setVGParameters([soil])   #|\label{l61:set_vg}|
s.initializeProblem()   #|\label{l61:initialise}|
s.setCriticalPressure(-15000)   #|\label{l61:criticalPsi}|
s.ddt = 1.e-5  # initial dumux time step [days]  |\label{l61:initialDT}|

top_ind = s.pick([0., 0., -0.5])   #|\label{l61:picker_a}|
bot_ind = s.pick([0., 0., -199.5])  #|\label{l61:picker_b}|
top_new, bot_new, soil_times = [], [], []  # |\label{l61:initialize_flux_vectors}|

N = int(np.ceil(sim_time / dt)) # |\label{l61:loop}|
for i in range(0, N):
    t = i * dt  # current simulation time   |\label{l61:current time_a}|
    soil_times.append(t)  # |\label{l61:current time_e}|
    s.solve(dt)  #|\label{l61:solve}|
    velocities = s.getVelocities_()  # |\label{l61:velocities_a}|
    top_new.append(velocities[top_ind])
    bot_new.append(velocities[bot_ind]) # |\label{l61:velocities_e}|  

top_new = np.array(top_new)   #|\label{l61:vectors_a}|  
bot_new = np.array(bot_new)    
soil_times = np.array(soil_times)   #|\label{l61:vectors_e}|  

# Extract and plot numerical solution  |\label{l61:plot_profile_a}|
points = s.getDofCoordinates() 
theta = s.getWaterContent()
plt.figure(0)
plt.plot(theta,points[:, 2],linewidth=2) 
plt.xlabel(r'$\theta$ (cm$^3$ cm$^{-3}$)', fontsize=18)
plt.ylabel('depth (cm)', fontsize=18)
plt.xticks(fontsize=14); plt.yticks(fontsize=14)
#plt.title('Infiltration front in loam after 1 day')
plt.show()   # |\label{l61:plot_profile_e}|

plt.figure(1)    #|\label{l61:plot_fluxes_a}|
plt.plot(soil_times, top_new[:, 2],label="surface flux")
plt.plot(soil_times, bot_new[:, 2], label="bottom flux")
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('Vertical water flux (cm/day)', fontsize=18)
plt.xticks(fontsize=14); plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.show()  # |\label{l61:plot_fluxes_e}|