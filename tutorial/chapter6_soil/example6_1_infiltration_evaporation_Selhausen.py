# This example builds on the infiltration and evaporation problem M2.1 and M2.2 from Schnepf et al. (2023, \url{doi.org/10.1093/insilicoplants/diad005}) and extends the simulation to a multi-layered soil profile at the Field Minirhizotron Facilities Selhausen (\url{https://www.fz-juelich.de/en/ibg/ibg-3/research-groups/modelling-terrestrial-systems/soil-root-systems-and-rhizosphere-processes/field-minirhizotron-facilities}). 
# The hydraulic properties of the soil profile are taken from Bauer et al. (2011, table 3, \url{https://doi.org/10.1007/s10533-011-9583-1}). 
# The Richards equations is solved to simulate water infiltration over 1 day into an initially dry soil followed by evaporation over 2 days. Only the vertical water movement is considered. A Neumann boundary condition is set at the upper boundary and a free drainage boundary condition is set at the lower boundary. 

# The code solves the Richards equation with DuMux. The github repository "dumux-rosi" (https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi.git) needs to be cloned
# into the same path level as CPlantBox.

# add paths to the folders containing CPlantBox and dumux-rosi
import sys
base_path= "./CPlantBox/dumux/"
sys.path.append(base_path+"dumux-rosi/python/modules")
sys.path.append(base_path+"dumux-rosi/build-cmake/cpp/python_binding/")
sys.path.append(base_path+"CPlantBox/src")

import matplotlib.pyplot as plt  #
import numpy as np  #

from rosi_richards import RichardsSPnum  # C++ part (Dumux binding) 
from richards import RichardsWrapper  # Python part  

# Define Mualem van Genuchten parameters for Selhausen soil profile according to Bauer et al. (2011, table 3, \url{https://doi.org/10.1007/s10533-011-9583-1}) |\label{l61ies:genuchten_a}|
# theta_r (-), theta_s (-), alpha (1/cm), n (-), Ks (cm d-1) 
l1 = [0.008, 0.389, 0.012, 1.97, 91.68] # 0-20 cm |\label{l61ies:hyprop_s}|
l2 = [0.008, 0.389, 0.023, 1.23, 63.36] # 20-33 cm
l3 = [0.008, 0.389, 0.01, 1.1, 10] # 33-57 cm
l4 = [0.008, 0.389, 0.01, 1.1, 10] # 57-120 cm
# Combine the hydraulic conductivity vectors from all soil layers to define soil type for simulation
soil = [l1,l2,l3,l4] # |\label{l61ies:hyprop_e}|
sim_time = 3.05  # 
dt = 0.05 #720 / (24 * 3600)  # time step [days]   

# Solve the Richards equation using the Python wrapper of dumux-rosi
s = RichardsWrapper(RichardsSPnum())  # 
s.initialize()  #
s.setTopBC("atmospheric", 0.5, [[0., 1.,1.,3.], [10., 10.,-0.1,-0.1]])  #  [cm/day] atmospheric is with surface run-off   |\label{l61ies:top_bc}|
s.setBotBC("freeDrainage")  # |\label{l61ies:bottom_bc}|
N = 119*10 # use a fine grid resolution of 1 mm per grid point in z direction |\label{l61ies:grid}|
s.createGrid([-5., -5., -120.], [5., 5., 0.], [1, 1, N])  # [cm] N   |\label{l61ies:grid}|
# define soil layers 
layers_ID=[4, 4, 3, 3, 2, 2, 1, 1]  # |\label{l61ies:layers_s}|
layers_pos= [-120., -57., -57., -33.,-33,-20,-20,0] # |\label{l61ies:layers_e}|
s.setLayersZ(layers_ID,layers_pos)
s.setHomogeneousIC(-400.)  # cm pressure head    |\label{l61ies:ic}|
s.setVGParameters(soil)  # |\label{l61ies:set_vg}|
s.initializeProblem()  # |\label{l61ies:initialise}|
s.setCriticalPressure(-15000)  #
s.ddt = 1.e-5  # initial dumux time step [days] 

top_ind = s.pick([0., 0., -0.5])  
bot_ind = s.pick([0., 0., -119.5]) #  |\label{l61ies:bot_ind}|
top_new, bot_new, soil_times = [], [], []  

N = int(np.ceil(sim_time / dt))  #
z_, x_, h_ = [], [], [] # initialite solution vectors
for i in range(0, N):
    t = i * dt  # current simulation time   
    soil_times.append(t)  
    s.solve(dt)
    if i/2==np.ceil(i/2): 
        print("***** external time step", dt, " d, simulation time", s.simTime, "d, internal time step", s.ddt, "d")
        
    # ---- TO DO: replace getVelocity function update code to new functions for getting fluxes
    #velocities = s.getVelocities_()  
    #top_new.append(velocities[top_ind])
    #bot_new.append(velocities[bot_ind])
    # extract numerical solution
    
    points = s.getDofCoordinates() # coordinates
    theta = s.getWaterContent() # volumetric water content
    h = s.getSolutionHead() # head
    z_.append(points[:, 2])
    x_.append(theta)
    h_.append(h) 

top_new = np.array(top_new)  
bot_new = np.array(bot_new)
soil_times = np.array(soil_times)  

# define output times
sel_idx=np.searchsorted(soil_times, [0., 0.2, 0.5, 1.,1.2, 2.,3.])

# Plot solutions |\label{l61ies:plt_prof}|
fig1, axs = plt.subplots(2,1, sharex='row',figsize=(12,12))
cols=["k-","r-","r--","r-.","b-.","b--","b-"]
ii=0
for i in sel_idx:
    axs[0].plot(x_[i], z_[i], cols[ii],label = f'{soil_times[i]:.2f} d')
    axs[1].plot(h_[i], z_[i], cols[ii],label = f'{soil_times[i]:.2f} d')
    ii=ii+1
axs[0].set_title("Infiltration 1 $d$ with 10 $cm$ $d^{-1}$ \n and evaporation 2 $d$ with 0.1 $cm$ $d^{-1}$ afterwards")
axs[0].set_xlabel('Volumetric water content [$cm^{3}$ $cm^{-3}$]')
axs[0].set_ylabel('Depth [cm]')
axs[1].set_xlabel('Head [cm]')
axs[1].set_ylabel('Depth [cm]')
axs[0].legend(loc="best")
axs[1].legend(loc="best")

plt.show()

#plt.figure(1)  
#plt.plot(soil_times, top_new[:, 2], label = "surface flux")
#plt.plot(soil_times, bot_new[:, 2], label = "bottom flux")
#plt.xlabel('Time (days)', fontsize = 18)
#plt.ylabel('Vertical water flux (cm/day)', fontsize = 18)
#plt.xticks(fontsize = 14); plt.yticks(fontsize = 14)
#plt.legend(fontsize = 14)
#plt.show() 