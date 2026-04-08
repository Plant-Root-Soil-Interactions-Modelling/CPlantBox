"""
plots histograms of the outer radius of the perirhizal zone for various discretisations
"""

import matplotlib.pyplot as plt
import numpy as np
import scenario_setup as scenario
import vtk

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.functional.Perirhizal import *


def get_outer_radii(rootsystem, type_str, cell_number):
    """setup"""
    if rootsystem == "soybean":
        soil_, table_name, min_b, max_b, _, area, Kc = scenario.soybean(0)  # 0 = envirotype
        xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
        simtime = 42  # 87.5
    elif rootsystem == "maize":
        soil_, table_name, min_b, max_b, _, area, Kc = scenario.maize(0)  # 0 = envirotype
        xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
        simtime = 56  # 95
    elif rootsystem == "springbarley":
        soil_, table_name, min_b, max_b, _, area, Kc = scenario.springbarley(0)  # 0 = envirotype
        xml_name = "data/spring_barley_CF12.xml"  # root growth model parameter file
        simtime = 49  # 95
    else:
        print("get_outer_radii() unknown rootsystem name", rootsystem)
        raise

    s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type=1)
    r_ = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth

    """ simulation and post processing """
    r = r_.ms
    r.initialize(False)
    r.simulate(simtime, False)
    peri = PerirhizalPython(r)
    if type_str == "voronoi":
        outer_radii = peri.get_outer_radii_bounded_voronoi()
    else:
        outer_radii = peri.get_outer_radii(type_str)

    seg2cell = np.array(r.getSegmentMapper())  # r.seg2cell as list

    ana = pb.SegmentAnalyser(r.mappedSegments())
    ana.addData("outer_r", outer_radii)
    ana.addData("length", ana.getParameter("length"))
    # organic_layer = pb.SDF_Cuboid(pb.Vector3d([-1e6, -1e6, -organic]), pb.Vector3d([1e6, 1e6, max_b[2]]))
    topsoil_layer = pb.SDF_Cuboid(pb.Vector3d([-1e6, -1e6, -topsoil]), pb.Vector3d([1e6, 1e6, max_b[2]]))
    subsoil_layer = pb.SDF_Cuboid(pb.Vector3d([-1e6, -1e6, min_b[2]]), pb.Vector3d([1e6, 1e6, -topsoil]))
    # ana0 = pb.SegmentAnalyser(ana)
    # ana0.crop(organic_layer)
    # ana0.pack()
    ana1 = pb.SegmentAnalyser(ana)
    ana1.crop(topsoil_layer)
    ana1.pack()
    ana2 = pb.SegmentAnalyser(ana)
    ana2.crop(subsoil_layer)
    ana2.pack()
    # outer0 = ana0.data["outer_r"]
    outer1 = ana1.data["outer_r"]
    outer2 = ana2.data["outer_r"]
    length1 = ana1.data["length"]
    length2 = ana2.data["length"]
    print("outer1", np.nanmin(outer1), np.nanmax(outer1), "median", np.nanmedian(outer1), "mean_top", np.nanmean(outer1), np.nanstd(outer1))
    print("outer2", np.nanmin(outer2), np.nanmax(outer2), "median", np.nanmedian(outer2), "mean_sub", np.nanmean(outer2), np.nanstd(outer2))  #
    print(len(outer_radii), len(outer1) + len(outer2))  # len(outer0) +

    return outer1, outer2, seg2cell, length1, length2


topsoil = 30  # 10 * 2.54
subsoil = 150  # 30 * 2.54

titles = []
titles.append("3D - (1cm)$^3$")
titles.append("3D - (2cm)$^3$")
titles.append("3D - (4cm)$^3$")
titles.append("1D - Layer thickness 1 cm ")
titles.append("1D - Layer thickness 2 cm")
titles.append("1D - Layer thickness 4 cm")

""" maize """
cell_numbers = []
cell_numbers.append([76, 16, 150])
cell_numbers.append([38, 8, 75])
cell_numbers.append([19, 4, 38])
cell_numbers.append([1, 1, 150])
cell_numbers.append([1, 1, 75])
cell_numbers.append([1, 1, 38])

for split_type in ["surface"]:  # ["length", "surface", "volume"]:
    fig, axes = plt.subplots(3, 2, figsize=(20, 18))
    stats = [["min", "max", "median", "mean", "std"]]
    for i, cell_number in enumerate(cell_numbers):
        outer1, outer2, seg2cell, length1, length2 = get_outer_radii("maize", split_type, cell_number)

        # print("outer_r.shape", outer_r.shape)
        # print("seg2cell.shape", seg2cell.shape)
        # print("cell_ids", np.min(seg2cell), np.max(seg2cell))
        # layer0 = outer_r[seg2cell == 140]
        # print("layer0", len(layer0))
        # print(layer0)

        outer1, length1 = PerirhizalPython.to_range_(None, outer1, length1, 0.0, 2.0)
        outer2, length2 = PerirhizalPython.to_range_(None, outer2, length2, 0.0, 2.0)
        ax = axes[i % 3, i // 3]
        # ax.hist(outer_r, bins = 40, rwidth = 0.9)
        ax.hist([outer1, outer2], weights=[length1, length2], bins=40, rwidth=0.9, label=["topsoil", "subsoil"], stacked=True)
        # ax.set_xlim(0, 4)
        ax.set_ylim(0, 4000)
        ax.set_title(titles[i])
        stats.append(
            [
                np.min(outer1),
                np.max(outer1),
                np.median(outer1),
                np.mean(outer1),
                np.std(outer1),
                np.min(outer2),
                np.max(outer2),
                np.median(outer2),
                np.mean(outer2),
                np.std(outer2),
            ]
        )

        if i < 3:
            ax.set_ylabel("Root length [cm] ")
        if i == 2 or i == 5:
            ax.set_xlabel("Perirhizal outer radius [cm], 40 bins")

    print(stats[0])
    for i, cell_number in enumerate(cell_numbers):
        print("Summary", cell_number, ":")
        print(stats[i + 1])

    plt.tight_layout()
    plt.savefig("res_maize_" + split_type + ".png")
    plt.show()

# """ spring barley """
cell_numbers = []
cell_numbers.append([13, 3, 150])
cell_numbers.append([7, 2, 75])
cell_numbers.append([3, 1, 38])
cell_numbers.append([1, 1, 150])
cell_numbers.append([1, 1, 75])
cell_numbers.append([1, 1, 38])

for split_type in ["surface"]:  # ["length", "surface", "volume"]:
    fig, axes = plt.subplots(3, 2, figsize=(20, 18))
    stats = [["min", "max", "median", "mean", "std"]]
    for i, cell_number in enumerate(cell_numbers):
        outer1, outer2, seg2cell, length1, length2 = get_outer_radii("springbarley", split_type, cell_number)
        # outer_r = np.minimum(outer_r, 2)
        outer1, length1 = PerirhizalPython.to_range_(None, outer1, length1, 0.0, 2.0)
        outer2, length2 = PerirhizalPython.to_range_(None, outer2, length2, 0.0, 2.0)
        ax = axes[i % 3, i // 3]
        # ax.hist(outer_r, bins = 40, rwidth = 0.9)
        # ax.hist([outer1, outer2], bins = 40, rwidth = 0.9, label = ["topsoil", "subsoil"], stacked = True)
        ax.hist([outer1, outer2], weights=[length1, length2], bins=40, rwidth=0.9, label=["topsoil", "subsoil"], stacked=True)
        # ax.set_xlim(0, 2)
        ax.set_ylim(0, 170)
        ax.set_title(titles[i])
        stats.append(
            [
                np.min(outer1),
                np.max(outer1),
                np.median(outer1),
                np.mean(outer1),
                np.std(outer1),
                np.min(outer2),
                np.max(outer2),
                np.median(outer2),
                np.mean(outer2),
                np.std(outer2),
            ]
        )

        if i < 3:
            ax.set_ylabel("Root length [cm] ")
        if i == 2 or i == 5:
            ax.set_xlabel("Perirhizal outer radius [cm], 40 bins")

    print(stats[0])
    for i, cell_number in enumerate(cell_numbers):
        print("Summary", cell_number, ":")
        print(stats[i + 1])

    plt.tight_layout()
    plt.savefig("res_springbarley_" + split_type + ".png")
    plt.show()

print("fin")
