"""This example reproduces the infiltration example from Schnepf et al. (2023, doi.org/10.1093/insilicoplants/diad005).
Water infiltrates into an initially dry soil from the soil surface. Only the vertical water movement is
considered. A constant Neumann boundary condition is set at the upper boundary and a free drainage
boundary condition is set at the lower boundary. There is an analytical solution for this simple example,
which can optionally be plotted for comparison.

This example solves the Richards equation with DuMux. The github repository "dumux-rosi" (https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi.git) needs to be cloned
into the same path level as CPlantBox."""

import matplotlib.pyplot as plt  # |\label{l61i:bibs1_a}|
import numpy as np  # |\label{l61i:bibs1_e}|

from plantbox.visualisation import figure_style
from rosi.richards import RichardsWrapper  # Python part  |\label{l61i:bibs2_e}|
from rosi.rosi_richards import RichardsSPnum  # C++ part (Dumux binding) |\label{l61i:bibs2_a}|

# Define van Genuchten parameters for sand, loam and clay  |\label{l61i:genuchten_a}|
# theta_r (-), theta_s (-), alpha (cm-1), n (-), Ks (cm day-1)
sand = [0.045, 0.43, 0.15, 3, 1000]
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = loam  # Select soil type for simulation   |\label{l61i:genuchten_e}|

# simulation time, days
sim_time = 1  # |\label{l61i:sim_time}|
dt = 720 / (24 * 3600)  # time step (day)    |\label{l61i:timestep}|

# Solve the Richards equation using the Python wrapper of dumux-rosi
s = RichardsWrapper(RichardsSPnum())  # |\label{l61i:initialize_a}|
s.initialize()  # |\label{l61i:initialize_e}|
s.setTopBC("atmospheric", 0.5, [[-1.0, 1.0e10], [100.0, 100.0]])  #  [cm/day] atmospheric is with surface run-off   |\label{l61i:top_bc}|
s.setBotBC("freeDrainage")  # |\label{l61i:bottom_bc}|
n_steps = 199
s.createGrid([-5.0, -5.0, -200.0], [5.0, 5.0, 0.0], [1, 1, n_steps])  # [cm] n_steps   |\label{l61i:grid}|
s.setHomogeneousIC(-400.0)  # matric potential (cm)   |\label{l61i:ic}|
s.setVGParameters([soil])  # |\label{l61i:set_vg}|
s.initializeProblem()  # |\label{l61i:initialise}|
s.setCriticalPressure(-15000)  # |\label{l61i:criticalPsi}|
s.ddt = 1.0e-5  # initial dumux time step (days)  |\label{l61i:initialDT}|

top_ind = s.pick([0.0, 0.0, -0.5])  # |\label{l61i:picker_a}|
bot_ind = s.pick([0.0, 0.0, -199.5])  # |\label{l61i:picker_b}|
top_new, bot_new, times = [], [], []  # |\label{l61i:initialize_flux_vectors}|

n_steps = round(sim_time / dt)  # |\label{l61i:loop}|
for i in range(0, n_steps):
    t = i * dt  # current simulation time |\label{l61i:current time_a}|
    times.append(t)  # |\label{l61i:current time_e}|
    s.solve(dt)  # |\label{l61i:solve}|
    velocities = s.getVelocities_()  # |\label{l61i:velocities_a}|
    top_new.append(velocities[top_ind])
    bot_new.append(velocities[bot_ind])  # |\label{l61i:velocities_e}|

top_new = np.array(top_new)  # |\label{l61i:vectors_a}|
bot_new = np.array(bot_new)
times = np.array(times)  # |\label{l61i:vectors_e}|

# Extract and plot numerical solution  |\label{l61i:plot_profile_a}|
points = s.getDofCoordinates()
theta = s.getWaterContent()

fig, ax = figure_style.subplots11()
ax.plot(theta, points[:, 2])
ax.set_xlabel("Volumetric water content (cm$^3$ cm$^{-3}$)")
ax.set_ylabel("Depth (cm)")
# ax.set_title('Infiltration front in loam after one day')
plt.show()  # |\label{l61i:plot_profile_e}|

fig, ax = figure_style.subplots11()  # |\label{l61i:plot_fluxes_a}|
ax.plot(times, top_new[:, 2], label="surface flux")
ax.plot(times, bot_new[:, 2], label="bottom flux")
ax.set_xlabel("Time (day)")
ax.set_ylabel("Vertical water flux (cm day$^{-1}$)")
ax.legend()
plt.show()  # |\label{l61i:plot_fluxes_e}|
