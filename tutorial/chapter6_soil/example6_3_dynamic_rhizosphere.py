""" 
Solute transport example - radially symmetric 1D model for nitrate uptake
"""
import sys; sys.path.append("../"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src") # |\label{l63:lib_start}|
sys.path.append("../../../dumux-rosi/python/modules"); sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/");

import datetime
import pickle
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import matplotlib as mpl
import plantbox as pb  # CPlantBox
from functional.xylem_flux import *  # root system Python hybrid solver
from rosi_richardsnc_cyl import RichardsNCCylFoam   # C++ part (Dumux binding), macroscopic soil model
from richards_flat import RichardsFlatWrapper  # Python part of cylindrcial
import functional.van_genuchten as vg # |\label{l63:lib_end}|

SMALL_SIZE = 20
MEDIUM_SIZE = 20
BIGGER_SIZE = 20
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
mpl.rcParams['mathtext.default'] = 'regular'

def plot_profile(cc, h, c,  depth = -100.): # |\label{l63:plot_profile_start}|
    """ shows soil matric potential and concentration in the profile"""
    fig, ax1 = plt.subplots()
    color = 'tab:red'
    ax1.plot(cc, h, color = color)
    ax1.set_ylabel("Soil water potential [hPa]", color = color)
    ax1.set_xlabel("Distance from root surface [cm]")
    ax1.tick_params(axis = 'y', labelcolor = color)
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.plot(cc, c, color = color)
    ax2.set_ylabel("Nitrate concentration [g/L]", color = color)
    ax2.set_xlabel("Distance from root surface [cm]")
    ax2.tick_params(axis = 'y', labelcolor = color)
    plt.tight_layout() # |\label{l63:plot_profile_end}|


def plot_history(area, w, c, N): # |\label{l63:plot_history_start}|
    """ plots concentration per liquid phase and concentration per soil volume"""
    c_ = np.array([np.sum(np.multiply(area, np.multiply(c[i], w[i]))) for i in range(0, N)])  # nitrate concentration per soil volume
    fig, ax1 = plt.subplots()
    color = 'tab:red'
    ax1.plot(np.linspace(0, simtime, N), np.sum(np.multiply(area,c), axis = 1), color = color)
    ax1.set_ylabel("[g/L] liquid phase", color = color)
    ax1.set_xlabel("Time [day]")
    ax1.tick_params(axis = 'y', labelcolor = color)
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.plot(np.linspace(0, simtime, N), c_, color = color)
    ax2.set_ylabel("[kg/m$^3$] soil", color = color)
    ax2.set_xlabel("Time [day]")
    ax2.tick_params(axis = 'y', labelcolor = color)
    plt.tight_layout() # |\label{l63:plot_history_end}|


""" Soil """  # |\label{l63:soil_start}|
s = RichardsFlatWrapper(RichardsNCCylFoam())  # water & single solute
logbase = 0.5
NC = 10 # [1] spatial resolution (1D model)
a_in = 0.02 # cm
a_out = 0.6 # cm
points = np.logspace(np.log(a_in) / np.log(logbase), np.log(a_out) / np.log(logbase), NC, base = logbase)
soil = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam
s.initialize()
s.createGrid1d(points)
s.setVGParameters([soil]) # |\label{l63:soil_end}|

""" Inital conditions """  # |\label{l63:ic_start}|
s.setHomogeneousIC(-100.)  # hPa homogeneous pressure head
nitrate_initial_values = 5.e-3 / soil[1] # [g/L] concentration in the soil to concentration in the water phase
s.setICZ_solute(nitrate_initial_values )  # step-wise function, ascending order  # |\label{l63:ic_end}|

""" Boundary conditions """  # |\label{l63:bc_start}|
RS_Uptake_Wmax = 1. # [cm/day]
s.setOuterBC("fluxCyl", 0.)  #  [cm/day] Neumann boundary condition
s.setInnerBC("fluxCyl", 0.) # |\label{l63:bc_end}|

RS_Uptake_Vmax = 2.7e-6 # [g cm-2 day-1], Roose and Kirk (2009) # |\label{l63:MM_start}|
RS_Uptake_km =  3.1e-6 # [g cm-3], Roose and Kirk (2009)        
s.setInnerBC_solute(8) # Michaelis Menten uptake
s.setParameter("RootSystem.Uptake.Vmax", s.dumux_str(RS_Uptake_Vmax))   # active uptake parameters
s.setParameter("RootSystem.Uptake.Km", s.dumux_str(RS_Uptake_km))
s.setParameter("Flux.UpwindWeight", "1") # |\label{l63:MM_end}|

""" Initialze problem """  # |\label{l63:init_start}|
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "True")  #
s.setParameter("Newton.EnableChop", "True")
s.setParameter("Component.MolarMass", "6.2e-2")  # [kg/mol]
s.setParameter("Component.LiquidDiffusionCoefficient", "1.9e-9")  # [m2/s]
s.initializeProblem()
wilting_point = -10000
s.setCriticalPressure(wilting_point)
s.ddt = 1e-4  # [day] initial Dumux time step # |\label{l63:init_end}|

""" Simulation loop """  # |\label{l63:simt_start}|
simtime = 3
dt = 3600. / (24.*3600) 
cc = np.array(s.getCellCenters())  # [cm] 
points =  np.array(s.getPoints())  # [cm] cell faces
area = np.pi * (points[1:]**2 - points[:-1]**2)  # [cm2] area of each cell
print("area per cell", area, "cm2", "distance from root surface", cc, "cm")
theta =  np.array(s.getWaterContent()) 
volume0 =  np.sum(np.multiply(theta , area))
print("\ndomain water volume", volume0, "cm3/cm  = ", volume0 / 1000, "l/cm")  
print("water content to water volume",volume0 , "cm3/cm") # |\label{l63:simt_end}|

N = int(np.ceil(simtime / dt)) # |\label{l63:param_start}|
c, h, w = [], [], []  # results
cmin = 0.  # |\label{l63:param_end}|
for i in range(0, N):  # |\label{l63:loop_start}|
    t = i * dt  # current simulation time
    print(t, "days")
    if cmin < 0.:
        raise Exception
    s.setInnerBC("fluxCyl", -RS_Uptake_Wmax * sinusoidal(t))  # [cm/day]
    s.solve(dt)
    h.append(s.getSolutionHead())  # [cm]
    w.append(s.getWaterContent())  # [1]
    c.append(s.getSolution(1))  # [g/L]
    cmin = cmin if cmin < min(s.getSolution(1)) else min(s.getSolution(1)) # |\label{l63:loop_end}| 
   

theta = w[-1] # |\label{l63:wc}|
volumef = np.sum(np.multiply(theta , area)) # |\label{l63:vol}|
print("domain water volume",  volumef, "cm3/cm  = ", volumef / 1000., "l/cm")  # |\label{l63:results}|
print("change in water volume",  volumef -  volume0, "cm3/cm = ", 1.e-3 * (volumef -  volume0), "l/cm")

plot_history(area, w, c, N) # |\label{l63:plot_history}|
plot_profile(cc, h[-1], c[-1]) # |\label{l63:plot_profile}|
plt.show()