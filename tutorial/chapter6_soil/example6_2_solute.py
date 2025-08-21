""" 
Solute transport example - nitrate in movement in soil 
"""
import sys; sys.path.append("../"); sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/python/modules"); sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/");

import datetime
import pickle
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import figure_style

import plantbox as pb  # CPlantBox
from rosi_richardsnc import RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from rosi_richards import RichardsSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
import functional.van_genuchten as vg


def plot_results(h, c , times, net_inf, fw, depth = -100.):
    """ creates a figures presenting soil water matric potential and 
    nitrate concentration over time """
    c = np.transpose(c)
    c = c[::-1,:]
    h = np.transpose(h)
    h = h[::-1,:]
    fig, ax = plt.subplots(3, 1, figsize = (18, 10), gridspec_kw = {'height_ratios': [1.5, 3, 3]})
    bar = ax[0].bar(times[::2], -10 * (net_inf[::2] - np.array(fw)), 0.8)  # cm -> mm
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


def plot_profile(h, c , depth = -100.):
    """ shows soil matric potential and concentration in the profile"""
    fig, ax1 = plt.subplots()
    color = 'tab:red'
    ax1.plot(h, np.linspace(depth, 0., h.shape[0]), color = color)
    ax1.set_xlabel("soil matric potential [cm]", color = color)
    ax1.set_ylabel("depth [cm]")
    ax1.tick_params(axis = 'x', labelcolor = color)
    ax2 = ax1.twiny()
    color = 'tab:blue'
    ax2.plot(c, np.linspace(depth, 0., c.shape[0]), ':', color = color)
    ax2.set_xlabel("nitrate concentration [g/L]", color = color)
    ax2.set_ylabel("depth [cm]")
    ax2.tick_params(axis = 'x', labelcolor = color)
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


def net_infiltration_csv(filename, start_date, end_date):
    """ Load net infiltration data from a CSV file and compute time as scalar 
    days since start_date"""
    df = pd.read_csv(filename)
    df['date'] = pd.to_datetime(df['date'])  # Convert 'date' column to datetime
    start_dt = pd.to_datetime(start_date)
    end_dt = pd.to_datetime(end_date)
    mask = (df['date'] >= start_dt) & (df['date'] <= end_dt)
    df_filtered = df.loc[mask]
    start_dt = pd.to_datetime(start_date)  # Reference start date
    times = (df_filtered['date'] - start_dt).dt.total_seconds() / (24 * 3600)
    netinf = df_filtered['NetInfiltration_mm_day']
    return times.to_numpy(), netinf.to_numpy() / 10  # mm/day > cm/day


""" Soil """  # |\label{l62:init_soil}|
s = RichardsWrapper(RichardsNCSP())  # water & single solute
min_b = [-35., -10., -100.]  # [cm]
max_b = [35., 10., 0.]  # [cm]
cell_number = [1, 1, 100]  # [1] spatial resolution (1D model)
area = (max_b[0] - min_b[0]) * (max_b[1] - min_b[1])  # [cm2]
vol = area * (max_b[2] - min_b[2])  # [cm3]
soil = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam
print("Area", area, "cm2,", "volume", vol, "cm3")
s.initialize()
s.createGrid(min_b, max_b, cell_number, False)
s.setVGParameters([soil])

""" Inital conditions """  # |\label{l62:init_ic}|
p_top = -400.  # initial matric potential [cm]
p_bot = -300.  # initial matric potential [cm]
s.setLinearIC(p_top, p_bot)  # [cm] pressure head
nitrate_z = [0., -30., -30., -100.]  # top soil layer of 30 cm
nitrate_initial_values = np.array([5.e-3, 5.e-3, 1.e-3, 1.e-3]) / 0.43 / 1000  #  [kg/m3] -> [g/L]
s.setICZ_solute(nitrate_initial_values[::-1], nitrate_z[::-1])  # step-wise function, ascending order

""" Boundary conditions """  # |\label{l62:init_bc}|
start_date_str = '2013-05-01'
end_date_str = '2013-08-01'
times, net_inf = net_infiltration_csv("RO_AKRW_003.2007-01-01T00_00_00.2015-01-01T00_00_00_net_infiltration.csv", start_date_str, end_date_str)
# plt.bar(times, net_inf, color = 'skyblue', edgecolor = 'k', width = 0.8);,
# print("net gain loss over the period: ", np.sum(net_inf), "l/m2")
# plt.show()
netinf_ = np.repeat(net_inf, 2)
times_ = np.repeat(times, 2)
s.setTopBC("atmospheric", 0.5, [times_[1:], netinf_[:-1] ])  # 0.5 is dummy value
s.setBotBC("freeDrainage")
# s.setBotBC("noflux")
s.setTopBC_solute(["constantFlux"], [0.], [0.])
# s.setTotBC_solute(["outflow"])
# s.setBotBC_solute(["constantFlux"], [0.])
s.setBotBC_solute(["outflow"])

""" Fertilizer """  # |\label{l62:init_source}|
fertilization_time = 31  # [day] fertilisation event
fertilization_amount = 80 * 1.e-5 * area / 1000  #  80 [kg/ha] = 80*1.e-5 [g/cm2]; -> [kg/day]

""" Initialze problem """  # |\label{l62:init}|
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

""" Simulation loop """  # |\label{l62:loop_init}|
start_date = datetime.datetime.strptime(start_date_str , '%Y-%m-%d')
end_date = datetime.datetime.strptime(end_date_str , '%Y-%m-%d')
timedelta_ = end_date - start_date
sim_time = timedelta_.days
dt = 3600. / (24.*3600)

theta = s.getWaterContent()
volume0 = s.getWaterVolume()
print("\ndomain water volume", volume0 , "cm3  = ", volume0 / 1000, "l")  # OK
print("water content to water volume", np.sum(theta) * area * 1, "cm3")  # OK
print("domain water volume", s.getWaterVolume() / area, "cm3/cm2  = ", s.getWaterVolume() / area * 10, "l/m2")  # OK
print("water content to water volume", np.sum(theta) * 1, "cm3/cm  = ", np.sum(theta) * 1 * 10, "l/m2")  # OK
print("sum net inf", 10 * np.sum(net_inf), "mm")
# plot_profile(s.getSolutionHead(), s.getSolution_(1))

N = int(np.ceil(sim_time / dt))
c, h, w = [], [], [] # results
fw = np.zeros(sim_time)

for i in range(0, N):  # |\label{l62:loop_loop}|

    t = i * dt  # current simulation time
    print(t, "days")

    if t >= fertilization_time and t < fertilization_time + 1:  # |\label{l62:fert_start}|
        ind = s.pick([0, 0, -0.5])
        soil_sol_fluxes = { ind: 0.7 * fertilization_amount}
        s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)
    elif t <= 1:
        ind = s.pick([0, 0, -0.5])
        soil_sol_fluxes = { ind: 0.3 * fertilization_amount}
        s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)
    else:
        ind = s.pick([0, 0, -0.5])
        soil_sol_fluxes = { ind: 0. * fertilization_amount}
        s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)  # |\label{l62:fert_end}|

    s.solve(dt, saveInnerFluxes_ = True)
    h.append(s.getSolutionHead_())  # [cm]
    w.append(s.getWaterContent())  # [1]
    c.append(s.getSolution_(1))  # [g/L]
    fw[int(t)] += s.getLowerBoundaryFluxesPerCell(0).sum() # [cm/day]  # |\label{l62:BCfluxes}| 

print("domain water volume", s.getWaterVolume(), "cm3  = ", s.getWaterVolume() / 1000., "l")  # |\label{l62:results}|
print("water content to water volume", np.sum(w[-1]) * area * 1, "cm3")
print("change in water volume", s.getWaterVolume() - volume0, "cm3 = ", 1.e-3 * (s.getWaterVolume() - volume0), "l")

plot_history(w, c, N)
plot_profile(h[-1], c[-1])
plot_results(h, c, times_, netinf_, fw, min_b[2])

