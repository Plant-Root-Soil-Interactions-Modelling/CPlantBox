import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")

import plantbox as pb
import functional.xylem_flux
import viewer_conductivities

import numpy as np


def plot_depth_profile(analyser, ax, j:int):
    """ 
    analyser     pb.SegmentAnalyser class
    ax           matplotlib axis
    j            type of plot (0: "length", 1: "surface", 2: "volume")
    """
    type_str = ["length", "surface", "volume"]
    unit_str = ["(cm)", "(cm$^2$)", "(cm$^3$)"]
    ax.clear()
    n = int(np.ceil(-analyser.getMinBounds().z))
    z_ = np.linspace(-0.5, -n + 0.5, n)
    d = analyser.distribution(type_str[j], 0., float(-n), int(n), True)
    ax.plot(d, z_, "-*", label = "total")
    max_type = int(np.max(analyser.data["subType"]))
    for i in range(0, max_type + 1):
        ana = pb.SegmentAnalyser(analyser)  # copy
        ana.filter("subType", i)
        segn = len(ana.segments)
        if segn > 0:
            d = ana.distribution(type_str[j], 0., float(-n), int(n), True)
            ax.plot(d, z_, "-*", label = "type {:g}".format(i))
    ax.set_ylabel("Depth (cm)")
    ax.set_xlabel("Root system " + type_str[j] + " per 1 cm layer " + unit_str[j])
    ax.legend()


def plot_rootsystem_development(analyser, ax2, j):
    """ 
    analyser     pb.SegmentAnalyser class
    ax2           matplotlib axis
    j            type of plot (0: "length", 1: "surface", 2: "volume")
    """
    type_str = ["length", "surface", "volume"]
    unit_str = ["(cm)", "(cm$^2$)", "(cm$^3$)"]
    ax2.clear()
    radii = analyser.data["radius"]  # copy once
    if j == 0:
        weights = [analyser.getSegmentLength(i) for i in range(0, len(analyser.segments))]
    elif j == 1:
        weights = [2 * np.pi * radii[i] * analyser.getSegmentLength(i) for i in range(0, len(analyser.segments))]
    elif j == 2:
        weights = [np.pi * radii[i] * radii[i] * analyser.getSegmentLength(i) for i in range(0, len(analyser.segments))]
    cts = np.array(analyser.data["creationTime"])
    try:
        l_, t_ = np.histogram(cts, 100, weights = weights)
        ax2.plot(0.5 * (t_[1:] + t_[:-1]), np.cumsum(l_), "-", label = "total")
        max_type = int(np.max(analyser.data["subType"]))
        for i in range(0, max_type + 1):
            ana = pb.SegmentAnalyser(analyser)  # copy
            ana.filter("subType", i)
            n = len(ana.segments)
            if n > 0:
                radii = ana.data["radius"]  # copy once
                if j == 0:
                    weights = [ana.getSegmentLength(i) for i in range(0, len(ana.segments))]
                elif j == 1:
                    weights = [2 * np.pi * radii[i] * ana.getSegmentLength(i) for i in range(0, len(ana.segments))]
                elif j == 2:
                    weights = [np.pi * radii[i] * radii[i] * ana.getSegmentLength(i) for i in range(0, len(ana.segments))]
                cts = np.array(ana.data["creationTime"])
                l_, t_ = np.histogram(cts, 100, weights = weights)
                ax2.plot(0.5 * (t_[1:] + t_[:-1]), np.cumsum(l_), "-", label = "type {:g}".format(i))
        ax2.legend()
    except:
        pass
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Root system " + type_str[j] + " " + unit_str[j])


def plot_suf(data, ax3, j):
    """ 
    plots suf versus depth per root type 
    
    data                 DataModel (in viewer_data.poy)
    ax3                  matplotlib axis
    j                    scenario hard coded in viewer_conductivities.py
    """
    ax3.clear()
    if j == 0:
        viewer_conductivities.init_constant_scenario1(data.xylem_flux)
    elif j == 1:
        viewer_conductivities.init_constant_scenario2(data.xylem_flux)
    elif j == 2:
        viewer_conductivities.init_dynamic_scenario1(data.xylem_flux)
    elif j == 3:
        viewer_conductivities.init_dynamic_scenario2(data.xylem_flux)
    elif j == 4:
        viewer_conductivities.init_constant_scenario_wine(data.xylem_flux)

    print(data.max_ct)
    print(data.base_segs)
    krs, _ = data.xylem_flux.get_krs(data.max_ct, data.base_segs)
    suf = data.xylem_flux.get_suf(data.max_ct)
    data.analyser.addData("SUF", suf)
    n = int(np.ceil(-data.analyser.getMinBounds().z))
    z_ = np.linspace(-0.5, -n + 0.5, n)
    d = data.analyser.distribution("SUF", 0., float(-n), int(n), False)  # False!!!
    print("SUF total", np.min(d), np.max(d), np.mean(d), np.sum(d))
    ax3.plot(d, z_, "-*", label = "total")
    max_type = int(np.max(data.analyser.data["subType"]))
    for i in range(0, max_type + 1):
        ana = pb.SegmentAnalyser(data.analyser)  # copy
        ana.filter("subType", i)
        segn = len(ana.segments)
        if segn > 0:
            d = ana.distribution("SUF", 0., float(-n), int(n), False)
            ax3.plot(d, z_, "-*", label = "type {:g}".format(i))
            print("SUF", i, np.min(d), np.max(d), np.mean(d), np.sum(d))
    ax3.set_title("Root system krs {:g} cm$^2$/day".format(krs))
    ax3.set_ylabel("Depth (cm)")
    ax3.set_xlabel("Root system surface uptake fraction (SUF) per 1 cm layer (1)")
    ax3.legend()


def plot_krs(data, ax, j):
    """ 
    plots suf versus depth per root type 
    
    data                 DataModel (in viewer_data.poy)
    ax                   matplotlib axis
    j                    scenario hard coded in viewer_conductivities.py
    """
    ax.clear()
    if j == 0:
        viewer_conductivities.init_constant_scenario1(data.xylem_flux)
    elif j == 1:
        viewer_conductivities.init_constant_scenario2(data.xylem_flux)
    elif j == 2:
        viewer_conductivities.init_dynamic_scenario1(data.xylem_flux)
    elif j == 3:
        viewer_conductivities.init_dynamic_scenario2(data.xylem_flux)
    elif j == 4:
        viewer_conductivities.init_constant_scenario_wine(data.xylem_flux)

    t = int(np.ceil(data.max_ct))
    t_ = np.linspace(1, t, t)
    krs_ = []
    for t in t_:
        krs, _ = data.xylem_flux.get_krs(t, data.base_segs)
        krs_.append(krs)
    ax.plot(t_, krs_)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Root system conductivity Krs [cm$^2$/day]")

    # for Debugging
    kr, kx = data.xylem_flux.get_conductivities(-1.)  # (default) -1. means: current time is maximum of node creation times
    print("\nkr")
    print(kr)
    print("\n")
