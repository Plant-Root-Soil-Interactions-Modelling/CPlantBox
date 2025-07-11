# This example reproduces the evaporation example M2.2 from Schnepf et al. (2023, doi.org/10.1093/insilicoplants/diad005).
# Water evaporated from the surface of an initially moist soil. Only the vertical water movement is
# considered. Atmospheric boundary conditions are set at the upper boundary and a free drainage
# boundary condition is set at the lower boundary. There is an analytical solution for this simple example,
# which can optionally be plotted for comparison.

# This example solves the Richards equation with DuMux. The github repository "dumux-rosi" (https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi.git) needs to be cloned
# into the same path level as CPlantBox.

import sys; sys.path.append("../modules"); sys.path.append("../../build-cmake/cpp/python_binding/");
# add paths to the folders containing CPlantBox and dumux-rosi
import sys;  # |\label{l61:paths_a}|
sys.path.append("../../../dumux-rosi/python/modules");
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/");
sys.path.append("../../../CPlantBox/src")  # |\label{l61:paths_e}|

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
# from example6_2_evaporation_analytic_solution import *  # analytical solution

import matplotlib.pyplot as plt
import numpy as np

# Define Van Genuchten and other parameters
sand = [0.045, 0.43, 0.15, 3, 1000]  # |\label{l61:params_a}|
loam = [0.08, 0.43, 0.04, 1.6, 50]
clay = [0.1, 0.4, 0.01, 1.1, 10]
soil = loam  # Select soil type for simulation   |\label{l61:params_e}|
sim_time = 10;  # |\label{l61:simtime}|
N = 2000
dt = sim_time / N  # time step [days]    |\label{l61:timestep}|
ic = -200  # |\label{l61:ic}|
evap = -0.1  # cm/d       |\label{l61:evap}|

# Solve the Richards equation using the Python wrapper of dumux-rosi
s = RichardsWrapper(RichardsSP())  # |\label{l61:initialize_a}|
s.initialize()  # |\label{l61:initialize_e}|
s.setTopBC("atmospheric", 0.5, [[0., 1.e10], [evap, evap]])  #  [cm/day] atmospheric is with surface run-off   |\label{l61:top_bc}|
# s.setTopBC("flux", evap)  #  [cm/day]
# s.setTopBC("constantPressure", -10000)
s.setBotBC("freeDrainage")  # BC freeDrainage   |\label{l61:bot_bc}|
NZ = 1399
s.createGrid([-5., -5., -100.], [5., 5., 0.], [1, 1, NZ])  # [cm]   |\label{l61:grid}|
# vols = (100. / NZ) * np.ones((NZ,)) * 100.  # cm3
s.setVGParameters([soil])  # |\label{l61:set_vg}|
s.setHomogeneousIC(ic)  # cm pressure head  |\label{l61:set_ic}|
# s.setParameter("Problem.EnableGravity", "false")
s.initializeProblem()  # |\label{l61:initialise}|
s.setCriticalPressure(-10000)  # |\label{l61:criticalPsi}|
s.setRegularisation(1.e-6, 0.)  # |\label{l61:regularize}|
idx_top = s.pickCell([0.0, 0.0, 0.0])  # index to watch surface flux  #|\label{l61:picker}|
initial_water = s.getWaterVolume()  # |\label{l61:gettheta}|
print(initial_water)
s.ddt = 1.e-5  # initial Dumux time step [days]  |\label{l61:initialDT}|
maxDt = 1.  # maximal Dumux time step [days]   |\label{l61:maxDT}|

x_, y_ = [], []
for i in range(0, N):  # |\label{l61:loop}|
    s.solve(dt, maxDt)  # |\label{l61:solve}|
    f = s.getNeumann(idx_top)  # f = s.getSolutionHeadAt(idx_top)   |\label{l61:Neumann_a}|
    #   current_water = s.getWaterVolume()
    #   f = (initial_water - current_water) / dt / 1.e2
    #   print(current_water, f)
    #   initial_water = current_water
    x_.append(s.simTime)
    y_.append(f)  # |\label{l61:Neumann_e}|

# Extract and plot numerical solution
plt.plot(x_, y_, 'b')  # |\label{l61:plot_fluxes_a}|
plt.ylabel('$E_{act}$ (cm day$^{-1}$)', fontsize = 20)
plt.xlim(0, 10)
plt.title(soil, fontsize = 20)
plt.xlabel('$t$ (days)', fontsize = 20)
plt.ylabel('$E_{act}$ (cm day$^{-1}$)', fontsize = 20)
plt.tick_params(axis = 'both', which = 'major', labelsize = 16)
plt.show()  # |\label{l61:plot_fluxes_e}|
