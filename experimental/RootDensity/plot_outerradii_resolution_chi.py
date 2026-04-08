"""
plots histograms of the outer radius of the perirhizal zone for various discretisations
"""

import matplotlib.pyplot as plt
import numpy as np
import scenario_setup as scenario
import vtk
from scipy.stats import chi2_contingency, chisquare

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
        outer_radii = peri.get_outer_radii_voronoi()
    else:
        outer_radii = peri.get_outer_radii(type_str)

    seg2cell = np.array(r.getSegmentMapper())  # r.seg2cell as list

    # outer_radii = peri.to_range_(outer_radii, 0., 2.)

    return outer_radii, seg2cell


root_system = "maize"  # maize, spring_barley

cell_numbers = []
cell_numbers.append([1, 1, 1])
if root_system == "maize":
    cell_numbers.append([76, 16, 150])
    cell_numbers.append([38, 8, 75])
    cell_numbers.append([19, 4, 38])
    cell_numbers.append([1, 1, 150])
    cell_numbers.append([1, 1, 75])
    cell_numbers.append([1, 1, 38])
elif root_system == "springbarley":
    cell_numbers.append([13, 3, 150])
    cell_numbers.append([7, 2, 75])
    cell_numbers.append([3, 1, 38])
    cell_numbers.append([1, 1, 150])
    cell_numbers.append([1, 1, 75])
    cell_numbers.append([1, 1, 38])

outer = [[], [], [], []]
stats = [[], [], [], []]  #     stats = [["min", "max", "median", "mean", "std"]]

for i, split_type in enumerate(["voronoi", "length", "surface", "volume"]):  # ["length", "surface", "volume"]:

    for j, cell_number in enumerate(cell_numbers):
        outer_r, seg2cell = get_outer_radii(root_system, split_type, cell_number)
        print("iteration", i, j)
        print("outer_r", outer_r.shape)
        stats[i].append([np.min(outer_r), np.max(outer_r), np.median(outer_r), np.mean(outer_r), np.std(outer_r)])
        outer[i].append(outer_r)

    # print(stats[0])
    # for i, cell_number in enumerate(cell_numbers):
    #     print("Summary", cell_number, ":")
    #     print(stats[i + 1])

p_ = np.zeros((4, len(cell_numbers) - 1))
mean_ = np.zeros((4, len(cell_numbers) - 1))

for i, split_type in enumerate(["voronoi", "length", "surface", "volume"]):  # ["length", "surface", "volume"]:

    for j, cell_number in enumerate(cell_numbers[1:]):

        print("\nIteration", i, j + 1)

        # print(outer[0][0].shape)
        # print(outer[i][j + 1].shape)
        observed = np.array([outer[0][0], outer[i][j + 1]])

        observed = np.minimum(observed, 3.99 * np.ones(observed.shape))  # limit to range(0,2)

        exp, bin_edges = np.histogram(observed[0, :], bins=20, range=(0, 4))
        obs, bin_edges = np.histogram(observed[1, :], bins=20, range=(0, 4))

        print("exp", exp.shape, np.sum(exp), exp)
        print("obs", obs.shape, np.sum(obs), obs)

        res = chisquare(f_obs=obs.transpose() / np.sum(obs) * 20, f_exp=exp.transpose() / np.sum(exp) * 20)
        print("Chi-Quadrat-Statistik:", res.statistic)
        print("P-Wert:", res.pvalue)
        p_[i, j] = res.pvalue
        mean_[i, j] = stats[i][j + 1][3]

print("summary")
print(p_)

print("\nmean")
# mean_ = np.array([[stats[i][j][3] for j in range(0, len(cell_numbers) - 1)] for i in range(0, 4)])
print("ref", stats[0][0][3])
print(mean_)
print("fin")

#
# Mean value comparison indicates as large as possible
# Chi square says resolution (2cm)^3 (but not sure i am doing it right)
#
