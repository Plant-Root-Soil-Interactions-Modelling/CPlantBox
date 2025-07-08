""" 
Solute transport example - nitrate in movement in soil 
"""
import sys; sys.path.append("../../../CPlantBox/"); sys.path.append("../../../CPlantBox/src")
sys.path.append("../../../dumux-rosi/python/modules"); sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/");

import datetime
import pickle
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd

import plantbox as pb  # CPlantBox
from rosi_richardsnc import RichardsNCSP  # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
import functional.van_genuchten as vg


def plot_results(h, c , times, net_inf, depth = -100.):
    """ creates a figures presenting soil water matric potential and nitrate concentration over time """
    c = np.transpose(c)
    c = c[::-1,:]
    h = np.transpose(h)
    h = h[::-1,:]

    fig, ax = plt.subplots(3, 1, figsize = (18, 10), gridspec_kw = {'height_ratios': [1.5, 3, 3]})
    bar = ax[0].bar(times[::2], np.array(net_inf[::2]) * 10, 0.1)  # cm -> mm
    ax[0].set_ylabel("net inf [mm/day]")
    ax[0].set_xlim(times[0], times[-1])
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
    im = ax[2].imshow(c, cmap = cmap_, aspect = 'auto', vmin = 0., vmax = 1.e-4, extent = [0 , sim_time, depth, 0.])  # vmax = 1.e-3, interpolation = 'bicubic', interpolation = 'nearest',
    cb = fig.colorbar(im, cax = cax, orientation = 'vertical')
    cb.ax.get_yaxis().labelpad = 30
    cb.set_label('nitrate concentration [g/cm$^3$]', rotation = 270)
    ax[2].set_ylabel("depth [cm]")
    ax[2].set_xlabel("time [days]")

    print("range", np.min(h), np.max(h), "cm")
    print("range", np.min(c), np.max(c), "g/cm3")
    plt.tight_layout()
    plt.show()


def plot_final_profile(h, c , times, net_inf, depth = -200.):
    """ shows the final profile"""
    fig, ax1 = plt.subplots()
    h = h[:, -1]
    color = 'tab:red'
    ax1.plot(h, np.linspace(0., depth, h.shape[0]), color = color)
    ax1.set_xlabel("soil matric potential [cm]", color = color)
    ax1.set_ylabel("depth [cm]")
    ax1.tick_params(axis = 'x', labelcolor = color)
    ax2 = ax1.twiny()
    c = c[:, -1]
    color = 'tab:blue'
    ax2.plot(c, np.linspace(0., depth, c.shape[0]), ':', color = color)
    ax2.set_xlabel("nitrate concentration [g/cm$^3$]", color = color)
    ax2.set_ylabel("depth [cm]")
    ax2.tick_params(axis = 'x', labelcolor = color)
    plt.tight_layout()
    plt.show()


def net_infiltration_table_pickle(filename, start_data, end_data):
    """ takes the net infiltration of filename 
    from start_date (format 1994-11-06 01:00:00) to end_data (format 1995-08-09 19:00:00) from a pickle file"""
    with open(filename, 'rb') as f:
        df = pickle.load(f)
    print(df)
    y = df["Net_infilteration"].loc[start_data: end_data].values / 10. * 24.  # mm -> cm, /hour -> /day
    y_, x_ = [], []
    for i in range(0, y.shape[0]):
        x_.extend([float(i) / 24, float(i + 1) / 24])  # hour -> day
        y_.extend([y[i], y[i]])
    return x_, y_


def add_nitrificatin_source(s, soil_sol_fluxes, nit_flux = 0.):
    """ adds a consant nitrate source @param nit_flux due to nitrification [g/day] 
    to the dictionary soil_sol_fluxes """
    z_ = np.linspace(-0.5, -29.5, 30)  # top 30 cm layers
    for z in z_:
        i = s.pick([0, 0, z])  # cell index
        soil_sol_fluxes[i] = nit_flux


if __name__ == '__main__':

    """ parameters """
    start_date_str = '1995-05-01 00:00:00'
    end_date_str = '1995-08-09 00:00:00'
    times, net_inf = net_infiltration_table_pickle("95.pkl", start_date_str, end_date_str)

    # start_date_str = '1996-05-01 00:00:00'
    # end_date_str = '1996-08-09 00:00:00'
    # times, net_inf = net_infiltration_table_pickle("96.pkl", start_date_str, end_date_str)

    # start_date_str = '1997-05-01 00:00:00'
    # end_date_str = '1997-08-09 00:00:00'
    # times, net_inf = net_infiltration_table_pickle("97.pkl", start_date_str, end_date_str)

    min_b = [-38., -8., -100.]
    max_b = [38., 8., 0.]
    area = (max_b[0] - min_b[0]) * (max_b[1] - min_b[1])  # [cm2]

    soil = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam
    vg.create_mfp_lookup(vg.Parameters(soil), -1.e5, 1000)

    # water potentials
    p_top = -300.  # initial matric potential [cm]
    p_bot = -100.  # initial matric potential [cm]

    # nitrate

    # Maize example
    # fertilization timings: 1. May - 1. November, fertilization 30% (at planting), 70% 1. June
    # fertilization amount: up to 60–80 kg/ha NO₃⁻-N -> 70kg/ha -> 70000g / 10000 m2 -> 7.e-4 g cm-2
    # nitrification range: 0.5-5 kg/ha/day;  1 kg/ha/day = 1.e-5 g/cm2/day

    fertilization_time = 31  # [day] fertilisation event
    fertilization_amount = 80 * 1.e-5  #  [g/cm2] = 80 [kg/ha]

    nitrification_rate = 0.1 * 1.e-6  # 1 mg/dm3/day = 1.e-6 / cm3 /day
    nitrate_z = [0., -30., -30., -100.]  # initial nitrate: top soil layer of 30 cm

    nitrate_initial_values = 1.e-2 * np.array([2.6e-4, 2.6e-4, 0.75 * 2.6e-4, 0.75 * 2.6e-4])  #  initial nitrate concentrations: kg / m3 (4.e-4)

    cell_number = [1, 1, 100]  # resolution (1D model)
    dt = 3600. / (24.*3600)

    """ start """
    start_date = datetime.datetime.strptime(start_date_str , '%Y-%m-%d %H:%M:%S')
    end_date = datetime.datetime.strptime(end_date_str , '%Y-%m-%d %H:%M:%S')
    timedelta_ = end_date - start_date
    sim_time = timedelta_.days

    s = RichardsWrapper(RichardsNCSP())  # water and one solute
    s.initialize()
    s.createGrid(min_b, max_b, cell_number, False)  # [cm]

    # IC
    s.setLinearIC(p_top, p_bot)  # [cm] pressure head, linearly interpolated
    s.setICZ_solute(nitrate_initial_values[::-1], nitrate_z[::-1])  # step-wise function, ascending order

    # BC
    s.setVGParameters([soil])
    lai_noroots = lambda x: 0.

    s.setTopBC("atmospheric", 0.5, [times, net_inf])  # 0.5 is dummy value
    s.setBotBC("freeDrainage")
    # s.setBotBC("noflux")  # to approximate deep drainage

    sol_times = np.array([0., 1.,
                          1., fertilization_time,
                          fertilization_time, fertilization_time + 1,
                          fertilization_time + 1, sim_time])  #

    sol_influx = np.array([0. * fertilization_amount, 0. * fertilization_amount,
                           0., 0.,
                           0.7 * fertilization_amount, 0.7 * fertilization_amount,
                           0., 0.])  # g/(cm2 day)

    # plt.plot(sol_times, sol_influx)  # quick check
    # plt.show()

    # s.setTopBC_solute(["managed"], [0.5], [sol_times, sol_influx])
    s.setTopBC_solute(["constantFlux"], [0.], [0.])
    s.setBotBC_solute(["outflow"])

    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    s.setParameter("Component.MolarMass", "6.2e-2")
    s.setParameter("Component.LiquidDiffusionCoefficient", "1.7e-9")
    s.initializeProblem()

    print()
    theta = s.getWaterContent()
    print("1 m2 / area", 1.e4 / area)
    print("domain water volume", s.getWaterVolume(), "cm3  = ", s.getWaterVolume() / 1000, "l")  # OK
    print("water content to water volume", np.sum(theta) * area, "cm3")  # OK
    print("domain water volume", s.getWaterVolume() / area, "cm3/cm2  = ", s.getWaterVolume() / area * 10, "l/m2")  # OK
    print("water content to water volume", np.sum(theta) * 1, "cm3/cm  = ", np.sum(theta) * 1 * 10, "l/m2")  # OK
    print("sum net inf", np.sum(net_inf) / 2 * 10, "mm  = ", np.sum(net_inf) / 2 * 1 * 10, "l/m2")  # (cm m2)/m2 = 10 l / m2, (each value is twice)
    print()

    wilting_point = -10000
    s.setCriticalPressure(wilting_point)
    s.ddt = 1e-4  # [day] initial Dumux time step
    c, h = [], []  # resulting solute concentration

    N = int(np.ceil(sim_time / dt))
    for i in range(0, N):
        t = i * dt  # current simulation time

        print(t, "days")

        if t >= fertilization_time and t < fertilization_time + 1:
            ind = s.pick([0, 0, -0.5])
            print("ferilizing")
            soil_sol_fluxes = { ind: fertilization_amount}
            s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)
        else:
            s.setSource({}, eq_idx = 1)

        # soil_sol_fluxes = {}  # empy dict
        # add_nitrificatin_source(s, soil_sol_fluxes, nit_flux = nitrification_rate)  # nitrification debendent on tillage practice
        # s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)  # richards.py [g/day]

        s.solve(dt)
        c.append(s.getSolution_(1))
        h.append(s.getSolutionHead_())

    theta = s.getWaterContent()
    print("domain water volume", s.getWaterVolume(), "cm3  = ", s.getWaterVolume() / 1000, "l")
    print("water content to water volume", np.sum(theta) * area, "cm3")

    plot_results(h, c, times, net_inf, min_b[2])
