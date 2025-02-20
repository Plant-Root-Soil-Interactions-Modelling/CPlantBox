# This example reproduces the infiltration example from Schnepf et al. (2023, doi.org/10.1093/insilicoplants/diad005). 
# Water infiltrates into an initially dry soil from the soil surface. Only the vertical water movement is 
# considered. A constant Neumann boundary condition is set at the upper boundary and a free drainage 
# boundary condition is set at the lower boundary. There is an analytical solution for this simple example, 
# which can optionally be plotted for comparison. 


# This example solves the Richards equation with DuMux. The github repository "dumux-rosi" (https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi.git) needs to be cloned
# into the same path level as CPlantBox. 


# add paths to the folders containing CPlantBox and dumux-rosi
import sys; 
sys.path.append("../../../dumux-rosi/python/modules");
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox/src")

import matplotlib.pyplot as plt 
import numpy as np

from rosi_richards import RichardsSPnum  # C++ part (Dumux binding) |\label{l61:bibs2}|
from richards import RichardsWrapper  # Python part
#from analytic_solution import *  # plots the analytical solutions to ax1, ax2, ax3 # optional |\label{l61:analytic}|

# Define van Genuchten parameters for sand, loam and clay 
# theta_r (-), theta_s (-), alpha (1/cm), n (-), Ks (cm d-1)   
sand = [0.045, 0.43, 0.15, 3, 1000] 
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil=loam   # Select soil type for simulation
# simulation time, days 
sim_time=1; 
dt = 360 / (24 * 3600)  # time step [day]

# Solve the Richards equation using the Python wrapper of dumux-rosi 
s = RichardsWrapper(RichardsSPnum())  
s.initialize()
s.setTopBC("atmospheric", 0.5, [[-1., 1.e10], [100., 100.]])  #  [cm/day] atmospheric is with surface run-off
s.setBotBC("freeDrainage")
N = 199
s.createGrid([-5., -5., -200.], [5., 5., 0.], [1, 1, N])  # [cm] N
s.setHomogeneousIC(-400.)  # cm pressure head
s.setVGParameters([soil])
s.initializeProblem()
s.setCriticalPressure(-15000)
s.ddt = 1.e-5  # initial dumux time step [days]

top_ind = s.pick([0., 0., -0.5])
bot_ind = s.pick([0., 0., -199.5])
top_new, bot_new, soil_times = [], [], []

N = int(np.ceil(sim_time / dt))
for i in range(0, N):
    t = i * dt  # current simulation time
    soil_times.append(t)
    s.solve(dt)
    velocities = s.getVelocities_()
    top_new.append(velocities[top_ind])
    bot_new.append(velocities[bot_ind])    

top_new = np.array(top_new)
bot_new = np.array(bot_new)    
soil_times = np.array(soil_times)

# Extract and plot numerical solution 
points = s.getDofCoordinates() 
theta = s.getWaterContent()
plt.figure(0)
plt.plot(theta,points[:, 2],linewidth=2) 
plt.xlabel(r'$\theta$ (cm$^3$ cm$^{-3}$)')
plt.ylabel('depth (cm)')
plt.title('Infiltration front in loam after 1 day')
plt.show()

plt.figure(1)
plt.plot(soil_times, top_new[:, 2])
plt.plot(soil_times, bot_new[:, 2])
plt.show()