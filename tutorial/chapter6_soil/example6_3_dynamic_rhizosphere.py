"""solute transport example - radially symmetric 1D model for nitrate uptake"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from plantbox.visualisation import figure_style
from rosi.richards_flat import RichardsFlatWrapper  # Python part of cylindrcial
from rosi.rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding), macroscopic soil model


def sinusoidal(t):
    """sinusoidal function (used for transpiration) (integral over one day is 1)"""
    return np.sin(2.0 * np.pi * np.array(t) - 0.5 * np.pi) + 1.0


def plot_profile(cc, h, c):  # |\label{l63:plot_profile_start}|
    """shows soil matric potential and concentration in the profile"""
    _, ax1 = figure_style.subplots11()
    ax1.plot(cc, h, color="tab:red")
    ax1.set_ylabel("Soil water potential (hPa)", color="tab:red")
    ax1.set_xlabel("Distance from root surface (cm)")
    ax1.tick_params(axis="y", labelcolor="tab:red")
    ax2 = ax1.twinx()
    ax2.plot(cc, c, color="tab:blue")
    ax2.set_ylabel("Nitrate concentration (g L$^{-1}$)", color="tab:blue")
    ax2.set_xlabel("Distance from root surface (cm)")
    ax2.tick_params(axis="y", labelcolor="tab:blue")
    plt.tight_layout()  # |\label{l63:plot_profile_end}|


def plot_history(area, w, c, n_steps):  # |\label{l63:plot_history_start}|
    """plots concentration per liquid phase and concentration per soil volume"""
    _, ax1 = figure_style.subplots11()
    c_ = np.array([np.sum(np.multiply(area, np.multiply(c[i], w[i]))) for i in range(0, n_steps)])  # nitrate concentration per soil volume
    ax1.plot(np.linspace(0, sim_time, n_steps), np.sum(np.multiply(area, c), axis=1), color="tab:red")
    ax1.set_ylabel("Liquid phase (g L$^{-1}$)", color="tab:red")
    ax1.set_xlabel("Time (day)")
    ax1.tick_params(axis="y", labelcolor="tab:red")
    ax2 = ax1.twinx()
    ax2.plot(np.linspace(0, sim_time, n_steps), c_, color="tab:blue")
    ax2.set_ylabel("Soil (kg m$^{-3}$)", color="tab:blue")
    ax2.set_xlabel("Time (day)")
    ax2.tick_params(axis="y", labelcolor="tab:blue")
    plt.tight_layout()  # |\label{l63:plot_history_end}|


# Soil |\label{l63:soil_start}|
s = RichardsFlatWrapper(RichardsNCCylFoam())  # water & single solute
logbase = 0.5
NC = 10  # spatial resolution (1D model)
a_in = 0.02  # cm
a_out = 0.6  # cm
points = np.logspace(np.log(a_in) / np.log(logbase), np.log(a_out) / np.log(logbase), NC, base=logbase)
soil = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam
s.initialize()
s.createGrid1d(points)
s.setVGParameters([soil])  # |\label{l63:soil_end}|

# Inital conditions |\label{l63:ic_start}|
s.setHomogeneousIC(-100.0)  # hPa homogeneous pressure head
nitrate_initial_values = 5.0e-3 / soil[1]  # concentration in the soil to concentration in the water phase (g L-1)
s.setICZ_solute(nitrate_initial_values)  # step-wise function, ascending order  # |\label{l63:ic_end}|

# Boundary conditions |\label{l63:bc_start}|
rs_uptake_Wmax = 1.0  # (cm day-1)
s.setOuterBC("fluxCyl", 0.0)  #  Neumann boundary condition (cm day-1)
s.setInnerBC("fluxCyl", 0.0)  # |\label{l63:bc_end}|

rs_uptake_Vmax = 2.7e-6  # (g cm-2 day-1), Roose and Kirk (2009) # |\label{l63:MM_start}|
rs_uptake_km = 3.1e-6  # (g cm-3), Roose and Kirk (2009)
s.setInnerBC_solute(8)  # Michaelis Menten uptake
s.setParameter("RootSystem.Uptake.Vmax", s.dumux_str(rs_uptake_Vmax))  # active uptake parameters
s.setParameter("RootSystem.Uptake.Km", s.dumux_str(rs_uptake_km))
s.setParameter("Flux.UpwindWeight", "1")  # |\label{l63:MM_end}|

# Initialze problem   # |\label{l63:init_start}|
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "True")
s.setParameter("Newton.EnableChop", "True")
s.setParameter("Component.MolarMass", "6.2e-2")  # (kg mol-1)
s.setParameter("Component.LiquidDiffusionCoefficient", "1.9e-9")  # (m2 s-1)
s.initializeProblem()
wilting_point = -10000
s.setCriticalPressure(wilting_point)
s.ddt = 1e-4  # initial Dumux time step (day)  |\label{l63:init_end}|

# Simulation loop |\label{l63:simt_start}|
sim_time = 3  # days
dt = 3600.0 / (24.0 * 3600)  # days
cc = np.array(s.getCellCenters())  # (cm)
points = np.array(s.getPoints())  # (cm) cell faces
area = np.pi * (points[1:] ** 2 - points[:-1] ** 2)  # (cm2) area of each cell
print("area per cell", area, "cm2", "distance from root surface", cc, "cm")
theta = np.array(s.getWaterContent())
volume0 = np.sum(np.multiply(theta, area))
print("\ndomain water volume", volume0, "cm3 cm-1 = ", volume0 / 1000, "l/cm")
print("water content to water volume", volume0, "cm3 cm-1")  # |\label{l63:simt_end}|

n_steps = int(np.ceil(sim_time / dt))  # |\label{l63:param_start}|
c, h, w = [], [], []  # results
cmin = 0.0  # |\label{l63:param_end}|
for i in range(0, n_steps):  # |\label{l63:loop_start}|
    t = i * dt  # current simulation time
    print(t, "days")
    if cmin < 0.0:
        raise Exception
    s.setInnerBC("fluxCyl", -rs_uptake_Wmax * sinusoidal(t))  # [cm/day]
    s.solve(dt)
    h.append(s.getSolutionHead())  # cm
    w.append(s.getWaterContent())  # 1
    c.append(s.getSolution(1))  # g/L
    cmin = cmin if cmin < min(s.getSolution(1)) else min(s.getSolution(1))  # |\label{l63:loop_end}|

theta = w[-1]  # |\label{l63:wc}|
volumef = np.sum(np.multiply(theta, area))  # |\label{l63:vol}|
print("domain water volume", volumef, "cm3 cm-1  = ", volumef / 1000.0, "l cm-1")  # |\label{l63:results}|
print("change in water volume", volumef - volume0, "cm3 cm-1 = ", 1.0e-3 * (volumef - volume0), "l cm-1")

area = np.squeeze(area, -1)
plot_history(area, w, c, n_steps)  # |\label{l63:plot_history}|
plot_profile(cc, h[-1], c[-1])  # |\label{l63:plot_profile}|
plt.show()
