import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")

import plantbox as pb
import functional.xylem_flux
import viewer_conductivities

import numpy as np
import colorsys
from matplotlib.lines import Line2D


def _subtype_color(value: int, max_label: int):
    """Match the Rhizomancer / VTK subtype palette used in the 3D viewer."""
    palette = [
        (1.0, 0.0, 0.0),   # type 0 / tap
        (0.0, 0.7, 0.0),
        (0.0, 0.55, 1.0),
        (1.0, 0.65, 0.0),
        (0.65, 0.0, 0.9),
        (0.0, 0.8, 0.8),
        (0.8, 0.8, 0.0),
    ]
    label = max(0, int(value))
    if label < len(palette):
        return palette[label]
    hue = (label - len(palette)) / max(1, (max_label + 1) - len(palette))
    return colorsys.hsv_to_rgb((hue + 0.12) % 1.0, 0.85, 1.0)


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
    ax.plot(d, z_, "-*", label = "total", color = "black")
    max_type = int(np.max(analyser.data["subType"]))
    for i in range(0, max_type + 1):
        ana = pb.SegmentAnalyser(analyser)  # copy
        ana.filter("subType", i)
        segn = len(ana.segments)
        if segn > 0:
            d = ana.distribution(type_str[j], 0., float(-n), int(n), True)
            ax.plot(d, z_, "-*", label = "type {:g}".format(i), color = _subtype_color(i, max_type))
    ax.set_ylabel("Depth (cm)")
    ax.set_xlabel("Root system " + type_str[j] + " per 1 cm layer " + unit_str[j])
    ax.legend()


def plot_rootsystem_development(analyser, ax2, j, include_unsupported_endpoints = False):
    """ 
    analyser     pb.SegmentAnalyser class
    ax2           matplotlib axis
    j            type of plot (0: "length", 1: "surface", 2: "volume")
    """
    type_str = ["length", "surface", "volume"]
    unit_str = ["(cm)", "(cm$^2$)", "(cm$^3$)"]
    ax2.clear()
    handles = []
    radii = analyser.data["radius"]  # copy once
    if j == 0:
        weights = [analyser.getSegmentLength(i) for i in range(0, len(analyser.segments))]
    elif j == 1:
        weights = [2 * np.pi * radii[i] * analyser.getSegmentLength(i) for i in range(0, len(analyser.segments))]
    elif j == 2:
        weights = [np.pi * radii[i] * radii[i] * analyser.getSegmentLength(i) for i in range(0, len(analyser.segments))]
    weights = np.array(weights)
    cts = np.array(analyser.data["creationTime"])
    # Time-resolved development plots should ignore geometry without supported timing,
    # otherwise ET <= 0 creates artificial bins and late jumps that are not observed data.
    valid = cts > 0
    try:
        global_endpoint_x = None
        if np.any(valid):
            l_, t_ = np.histogram(cts[valid], 100, weights = weights[valid])
            x_line = 0.5 * (t_[1:] + t_[:-1])
            y_line = np.cumsum(l_)
            line = ax2.plot(x_line, y_line, "-", label = "total", color = "black")[0]
            handles.append(line)
            if include_unsupported_endpoints and x_line.size > 0:
                global_endpoint_x = x_line[-1]
                ax2.plot(
                    [global_endpoint_x], [np.sum(weights)],
                    linestyle = "None",
                    marker = "o",
                    markersize = 7,
                    markerfacecolor = line.get_color(),
                    markeredgecolor = line.get_color(),
                    markeredgewidth = 1.5,
                )
        elif include_unsupported_endpoints and np.sum(weights) > 0:
            fallback_x = 1.0
            ax2.plot(
                [fallback_x], [np.sum(weights)],
                linestyle = "None",
                marker = "o",
                markersize = 7,
                markerfacecolor = "black",
                markeredgecolor = "black",
                markeredgewidth = 1.5,
            )
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
                weights = np.array(weights)
                cts = np.array(ana.data["creationTime"])
                valid = cts > 0
                if np.any(valid):
                    l_, t_ = np.histogram(cts[valid], 100, weights = weights[valid])
                    x_line = 0.5 * (t_[1:] + t_[:-1])
                    y_line = np.cumsum(l_)
                    line = ax2.plot(
                        x_line,
                        y_line,
                        "-",
                        label = "type {:g}".format(i),
                        color = _subtype_color(i, max_type)
                    )[0]
                    handles.append(line)
                    if include_unsupported_endpoints and x_line.size > 0:
                        ax2.plot(
                            [x_line[-1]], [np.sum(weights)],
                            linestyle = "None",
                            marker = "o",
                            markersize = 7,
                            markerfacecolor = line.get_color(),
                            markeredgecolor = line.get_color(),
                            markeredgewidth = 1.5,
                        )
                elif include_unsupported_endpoints and np.sum(weights) > 0:
                    fallback_x = global_endpoint_x if global_endpoint_x is not None else 1.0
                    subtype_color = _subtype_color(i, max_type)
                    ax2.plot(
                        [fallback_x], [np.sum(weights)],
                        linestyle = "None",
                        marker = "o",
                        markersize = 7,
                        markerfacecolor = subtype_color,
                        markeredgecolor = subtype_color,
                        markeredgewidth = 1.5,
                    )
        if include_unsupported_endpoints and handles:
            handles.append(Line2D(
                [0], [0],
                linestyle = "None",
                marker = "o",
                markersize = 7,
                markerfacecolor = "black",
                markeredgecolor = "black",
                markeredgewidth = 1.5,
                label = "Final value including ET<=0",
            ))
        if handles:
            ax2.legend(handles = handles)
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

    print(data.max_ct)
    print(data.base_segs)
    krs, _ = data.xylem_flux.get_krs(data.max_ct, data.base_segs)
    suf = data.xylem_flux.get_suf(data.max_ct)
    data.analyser.addData("SUF", suf)
    n = int(np.ceil(-data.analyser.getMinBounds().z))
    z_ = np.linspace(-0.5, -n + 0.5, n)
    d = data.analyser.distribution("SUF", 0., float(-n), int(n), False)  # False!!!
    print("SUF total", np.min(d), np.max(d), np.mean(d), np.sum(d))
    ax3.plot(d, z_, "-*", label = "total", color = "black")
    max_type = int(np.max(data.analyser.data["subType"]))
    for i in range(0, max_type + 1):
        ana = pb.SegmentAnalyser(data.analyser)  # copy
        ana.filter("subType", i)
        segn = len(ana.segments)
        if segn > 0:
            d = ana.distribution("SUF", 0., float(-n), int(n), False)
            ax3.plot(d, z_, "-*", label = "type {:g}".format(i), color = _subtype_color(i, max_type))
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
