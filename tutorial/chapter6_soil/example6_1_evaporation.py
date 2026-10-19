"""This example reproduces the evaporation example M2.2 from Schnepf et al. (2023, doi.org/10.1093/insilicoplants/diad005). Water evaporated from the surface of an initially moist soil. Only the vertical water movement is considered. Atmospheric boundary conditions are set at the upper boundary and a free drainage boundary condition is set at the lower boundary. There is an analytical solution for this simple example, which can optionally be plotted for comparison.
This example solves the Richards equation with DuMux. The github repository "dumux-rosi" (https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi.git) needs to be cloned into the same path level as CPlantBox."""

import matplotlib.pyplot as plt

from plantbox.visualisation import figure_style
from rosi.richards import RichardsWrapper  # Python part |\label{l61:paths_a}|
from rosi.rosi_richards import RichardsSP  # C++ part (Dumux binding) |\label{paths_e}|

# Define Van Genuchten and other parameters
soils = {
    "Sand": [0.045, 0.43, 0.15, 3, 1000],  # |\label{l61:params_a}|
    "Loam": [0.08, 0.43, 0.04, 1.6, 50],
    "Clay": [0.1, 0.4, 0.01, 1.1, 10],
}
soil = "Loam"  # Select soil type for simulation   |\label{l61:params_e}|
sim_time = 10  # |\label{l61:sim_time}|
n_steps = 10 * 24 * 60 // 10
dt = sim_time / n_steps  # time step (days)    |\label{l61:timestep}|
ic = -200  # |\label{l61:ic}|
evap = -0.1  # cm days-1       |\label{l61:evap}|

# Solve the Richards equation using the Python wrapper of dumux-rosi
s = RichardsWrapper(RichardsSP())  # |\label{l61:initialize_a}|
s.initialize()  # |\label{l61:initialize_e}|
s.setTopBC("atmospheric", 0.5, [[0.0, 1.0e10], [evap, evap]])  #  (cm day-1) atmospheric is with surface run-off   |\label{l61:top_bc}|
# s.setTopBC("flux", evap)  #  cm day-1
# s.setTopBC("constantPressure", -10000)
s.setBotBC("freeDrainage")  # BC freeDrainage   |\label{l61:bot_bc}|
NZ = 100  # 1399
s.createGrid([-5.0, -5.0, -100.0], [5.0, 5.0, 0.0], [1, 1, NZ])  # [cm]   |\label{l61:grid}|
# vols = (100. / NZ) * np.ones((NZ,)) * 100.  # cm3
s.setVGParameters([soils[soil]])  # |\label{l61:set_vg}|
s.setHomogeneousIC(ic)  # cm pressure head  |\label{l61:set_ic}|
# s.setParameter("Problem.EnableGravity", "false")
maxDt = 1.0  # maximal Dumux time step (day)   |\label{l61:maxDT}|
s.initializeProblem(maxDt)  # |\label{l61:initialise}|
s.setCriticalPressure(-10000)  # |\label{l61:criticalPsi}|
s.setRegularisation(1.0e-6, 0.0)  # |\label{l61:regularize}|
idx_top = s.pickCell([0.0, 0.0, 0.0])  # index to watch surface flux  #|\label{l61:picker}|
initial_water = s.getWaterVolume()  # |\label{l61:gettheta}|
s.ddt = 1.0e-5  # initial Dumux time step (days)  |\label{l61:initialDT}|

x_, y_ = [], []
for i in range(0, n_steps):  # |\label{l61:loop}|
    print(f" {i * dt:g} days")
    s.solve(dt)  # |\label{l61:solve}|
    f = s.getNeumann(idx_top)  # f = s.getSolutionHeadAt(idx_top)   |\label{l61:Neumann_a}|
    x_.append(s.simTime)
    y_.append(f)  # |\label{l61:Neumann_e}|

# Extract and plot numerical solution
fig, ax = figure_style.subplots11()
ax.plot(x_, y_, "b")  # |\label{l61:plot_fluxes_a}|
ax.set_xlim(0, 10)
ax.set_title(soil)
ax.set_xlabel("Time (day)")
ax.set_ylabel("Actual evaporation (cm day$^{-1}$)")
ax.tick_params(axis="both", which="major")
plt.show()  # |\label{l61:plot_fluxes_e}|
