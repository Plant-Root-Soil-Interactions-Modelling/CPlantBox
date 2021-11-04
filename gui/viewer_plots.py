import sys; sys.path.append("../src/python_modules/"); sys.path.append("../")

import plantbox as pb
import xylem_flux
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
    weights = [analyser.getSegmentLength(i) for i in range(0, len(analyser.segments))]
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
                weights = [ana.getSegmentLength(ii) for ii in range(0, n)]
                cts = np.array(ana.data["creationTime"])
                l_, t_ = np.histogram(cts, 100, weights = weights)
                ax2.plot(0.5 * (t_[1:] + t_[:-1]), np.cumsum(l_), "-", label = "type {:g}".format(i))
        ax2.legend()
    except:
        pass
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Root system " + type_str[j] + " " + unit_str[j])


def plot_suf(analyser, mapped_segments, max_ct, ax3, j):
    """ 
    plots suf versus depth per root type 
    
    analyser             pb.SegmentAnalyser class
    mapped_segments      pb.MappedSegments class
    max_ct               maximal creation time  
    ax3                  matplotlib axis
    j                    scenario hard coded in viewer_conductivities.py
    """
    ax3.clear()
    r = xylem_flux.XylemFluxPython(mapped_segments)
    if j == 0:
        viewer_conductivities.init_constant_scenario1(r)
    elif j == 1:
        viewer_conductivities.init_constant_scenario2(r)
    elif j == 2:
        viewer_conductivities.init_dynamic_scenario1(r)
    elif j == 3:
        viewer_conductivities.init_dynamic_scenario2(r)
    krs, _ = r.get_krs(max_ct)  # TODO move to label
    suf = r.get_suf(max_ct)
    analyser.addData("SUF", suf)
    n = int(np.ceil(-analyser.getMinBounds().z))
    z_ = np.linspace(-0.5, -n + 0.5, n)
    d = analyser.distribution("SUF", 0., float(-n), int(n), True)
    ax3.plot(d, z_, "-*", label = "total")
    max_type = int(np.max(analyser.data["subType"]))
    for i in range(0, max_type + 1):
        ana = pb.SegmentAnalyser(analyser)  # copy
        ana.filter("subType", i)
        segn = len(ana.segments)
        if segn > 0:
            d = ana.distribution("SUF", 0., float(-n), int(n), True)
            ax3.plot(d, z_, "-*", label = "type {:g}".format(i))
    ax3.set_title("Root system krs {:g}".format(krs))
    ax3.set_ylabel("Depth (cm)")
    ax3.set_xlabel("Root system surface uptake fraction (SUF) (1)")
    ax3.legend()

