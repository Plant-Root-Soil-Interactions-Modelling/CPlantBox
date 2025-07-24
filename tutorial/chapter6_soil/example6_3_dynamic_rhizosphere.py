""" 
Solute transport example - nitrate in movement in soil 
"""
import sys; sys.path.append("../"); sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src")
sys.path.append("../../../dumux-rosi/python/modules"); sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/");

import datetime
import pickle
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd

import plantbox as pb  # CPlantBox
from rosi_richardsnc_cyl import RichardsNCCylFoam   # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial
import functional.van_genuchten as vg


def plot_results(h, c , times, net_inf, depth = -100.):
    """ creates a figures presenting soil water matric potential and 
    nitrate concentration over time """
    c = np.transpose(c)
    c = c[::-1,:]
    h = np.transpose(h)
    h = h[::-1,:]
    fig, ax = plt.subplots(3, 1, figsize = (18, 10), gridspec_kw = {'height_ratios': [1.5, 3, 3]})
    bar = ax[0].bar(times[::2], -10 * net_inf[::2], 0.8)  # cm -> mm
    ax[0].set_ylabel("net inf [mm/day]")
    ax[0].set_xlim(times[0] - 0.5, times[-1] + 0.5)
    divider = make_axes_locatable(ax[0])
    cax0 = divider.append_axes('right', size = '5%', pad = 0.05)
    cax0.axis('off')
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes('right', size = '5%', pad = 0.05)
    cmap_reversed = matplotlib.cm.get_cmap('jet_r')
    im = ax[1].imshow(h, cmap = cmap_reversed, aspect = 'auto', vmin = -1.e3, extent = [0 , sim_time, depth, 0.])
    cb = fig.colorbar(im, cax = cax, orientation = 'vertical')
    cb.ax.get_yaxis().labelpad = 30
    cb.set_label('soil matric potential [cm]', rotation = 270)
    ax[1].set_ylabel("depth [cm]")
    ax[1].set_xlabel("time [days]")
    divider = make_axes_locatable(ax[2])
    cax = divider.append_axes('right', size = '5%', pad = 0.05)
    cmap_ = matplotlib.cm.get_cmap('jet')
    im = ax[2].imshow(c, cmap = cmap_, aspect = 'auto', extent = [0 , sim_time, depth, 0.])
    cb = fig.colorbar(im, cax = cax, orientation = 'vertical')
    cb.ax.get_yaxis().labelpad = 30
    cb.set_label('nitrate concentration [g/L]', rotation = 270)
    ax[2].set_ylabel("depth [cm]")
    ax[2].set_xlabel("time [days]")
    print("range", np.min(h), np.max(h), "cm")
    print("range", np.min(c), np.max(c), "g/L")
    plt.tight_layout()
    plt.show()


def plot_profile(h, c, cc, depth = -100.):
    """ shows soil matric potential and concentration in the profile"""
    fig, ax1 = plt.subplots()
    color = 'tab:red'
    ax1.plot(cc, h, color = color)
    ax1.set_ylabel("soil matric potential [cm]", color = color)
    ax1.set_xlabel("distance from root surface [cm]")
    ax1.tick_params(axis = 'y', labelcolor = color)
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.plot(cc, c, ':', color = color)
    ax2.set_ylabel("nitrate concentration [g/L]", color = color)
    ax2.set_xlabel("distance from root surface [cm]")
    ax2.tick_params(axis = 'y', labelcolor = color)
    plt.tight_layout()
    plt.show()


def plot_history(w, c, N):
    """ plots concentration per liquid phase and concentration per soil volume"""
    c_ = np.array([np.sum(np.multiply(c[i], w[i])) for i in range(0, N)])  # nitrate concentration per soil volume
    fig, ax1 = plt.subplots()
    color = 'tab:red'
    ax1.plot(np.linspace(0, sim_time, N), np.sum(c, axis = 1), color = color)
    ax1.set_ylabel("[g/L] liquid phase", color = color)
    ax1.set_xlabel("Time [day]")
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.plot(np.linspace(0, sim_time, N), c_, color = color)
    ax2.set_ylabel("[kg/m$^3$] soil", color = color)
    ax2.set_xlabel("Time [day]")
    plt.tight_layout()
    plt.show()




""" Soil """  # |\label{l63:init_soil}|
s = RichardsNoMPIWrapper(RichardsNCCylFoam())  # water & single solute
logbase = 0.5
NC = 10 # [1] spatial resolution (1D model)
a_in = 0.02 # cm
a_out = 0.6
points = np.logspace(np.log(a_in) / np.log(logbase), np.log(a_out) / np.log(logbase), 
                             NC, base = logbase)

soil = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam
s.initialize()
s.createGrid1d(points)
s.setVGParameters([soil])

""" Inital conditions """  # |\label{l63:init_ic}|
s.setHomogeneousIC(-100.)  # cm homogeneous pressure head
# x = np.logspace(np.log(-10000) / np.log(logbase), np.log(-100) / np.log(logbase), 
#                            NC, base = logbase)
# s.setParameter("Soil.IC.P", s.dumux_str(x))# cm
nitrate_z = [0., -30., -30., -100.]  # top soil layer of 30 cm
nitrate_initial_values = np.array([5.e-3, 5.e-3, 1.e-3, 1.e-3]) / 0.43# / 1000  #  [kg/m3] -> [g/L]
s.setICZ_solute(nitrate_initial_values[::-1], nitrate_z[::-1])  # step-wise function, ascending order

""" Boundary conditions """  # |\label{l63:init_bc}|
s.setOuterBC("fluxCyl", 0.)  #  [cm/day] Neumann boundary condition
s.setInnerBC("fluxCyl", -1.)  # [cm/day]

RS_Uptake_Vmax = 2.7e-6 # [g cm-2 day-1]
RS_Uptake_km =  3.1e-6 # [g cm-3]
             
s.setInnerBC_solute(8) # Michaelis Menten uptake
s.setParameter("RootSystem.Uptake.Vmax", s.dumux_str(RS_Uptake_Vmax))   # active uptake parameters
s.setParameter("RootSystem.Uptake.Km", s.dumux_str(RS_Uptake_km))

s.setParameter("Flux.UpwindWeight", "1")

""" Initialze problem """  # |\label{l63:init}|
s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Newton.MaxAbsoluteResidual", "1.e-10")
s.setParameter("Newton.SatisfyResidualAndShiftCriterion", "True")  #

s.setParameter("Newton.EnableChop", "True")
s.setParameter("Component.MolarMass", "6.2e-2")  # [kg/mol]
s.setParameter("Component.LiquidDiffusionCoefficient", "1.9e-9")  # [m2/s]
s.initializeProblem()
wilting_point = -10000
s.setCriticalPressure(wilting_point)
s.ddt = 1e-4  # [day] initial Dumux time step

""" Simulation loop """  # |\label{l63:loop_init}|
sim_time = 3
dt = 3600. / (24.*3600)

cc = np.array(s.getCellCenters())  # [cm] 
points =  np.array(s.getPoints())  # [cm] cell faces
area = np.pi * (points[1:]**2 - points[:-1]**2)  # [cm2] area of each cell
print("area per cell", area, "cm2", "distance from root surface", cc, "cm")
theta =  np.array(s.getWaterContent())
volume0 =  np.sum(np.multiply(theta , area))
print("\ndomain water volume", volume0, "cm3/cm  = ", volume0 / 1000, "l/cm")  
print("water content to water volume",volume0 , "cm3/cm") 

N = int(np.ceil(sim_time / dt))
c, h, w = [], [], []  # results
cmin = 0.
print('sol',s.getSolution_(1).flatten())
for i in range(0, N):  # |\label{l63:loop_loop}|

    t = i * dt  # current simulation time
    print(t, "days")
    if cmin < 0.:
        raise Exception

    s.solve(dt)
    h.append(s.getSolutionHead_())  # [cm]
    w.append(s.getWaterContent())  # [1]
    c.append(s.getSolution_(1))  # [g/L]
    cmin = cmin if cmin < min(s.getSolution_(1).flatten()) else min(s.getSolution_(1).flatten())      

theta = s.getWaterContent()
volumef = np.sum(np.multiply(theta , area))
print("domain water volume",  volumef, "cm3/cm  = ", volumef / 1000., "l/cm")  # |\label{l63:results}|
print("change in water volume",  volumef -  volume0, "cm3/cm = ", 1.e-3 * (volumef -  volume0), "l/cm")

plot_history(w, c, N)
plot_profile(h[-1], c[-1],cc )

