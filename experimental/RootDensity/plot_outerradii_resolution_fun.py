"""
TODO fitting to probabilty density function (PDF)
"""

import matplotlib.pyplot as plt
import numpy as np
import scenario_setup as scenario
import vtk
from scipy.optimize import curve_fit, minimize
from scipy.stats import chi2_contingency, chisquare, gamma

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.functional.Perirhizal import *


def fit_gamma(data):

    print("data2", data.shape, np.min(data), np.max(data), np.mean(data))

    # # Generate example data
    # data = gamma.rvs(a = 2, loc = 0, scale = 1, size = 1000)  # Example data from a gamma distribution

    def neg_log_likelihood(params):  # Define the negative log-likelihood function
        k, theta = params
        if k <= 0 or theta <= 0:
            return np.inf
        return -np.sum(gamma.logpdf(data, a=k, loc=0, scale=theta))

    # Initial parameter values
    initial_params = [np.mean(data), 1.0]  # (k, theta) mean = ktheta, variance = k theta^2; You can choose your own initial values
    print("initial_params", initial_params)

    # Minimize the negative log-likelihood function
    result = minimize(neg_log_likelihood, initial_params, method="L-BFGS-B")
    optimal_params = result.x

    # Print the estimated parameters
    k, theta = optimal_params
    # print("Estimated Parameters:")
    # print("k =", k)
    # print("θ =", theta)

    return k, theta


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


root_system = "springbarley"  # maize, spring_barley

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
outer_hist = [[], [], [], []]
stats = [[], [], [], []]  #     stats = [["min", "max", "median", "mean", "std"]]

for i, split_type in enumerate(["voronoi", "length", "surface", "volume"]):  # ["length", "surface", "volume"]:
    for j, cell_number in enumerate(cell_numbers):

        outer_r, seg2cell = get_outer_radii(root_system, split_type, cell_number)
        print("iteration", i, j)
        stats[i].append([np.min(outer_r), np.max(outer_r), np.median(outer_r), np.mean(outer_r), np.std(outer_r)])
        outer[i].append(outer_r)
        outer_r = np.minimum(outer_r, 2 * np.ones(outer_r.shape))
        outer_hist_, bin_edges = np.histogram(outer_r, bins=40, range=(0, 2))
        # print("outer_r", outer_r.shape, np.min(outer_r), np.max(outer_r))
        # print("outer_hist_", outer_hist_.shape, np.min(outer_hist_), np.max(outer_hist_))
        # print("outer_hist", np.sum(outer_hist_))
        outer_hist[i].append(np.array(outer_hist_) / np.sum(outer_hist_))

print("\nFITTING\n")

for i, split_type in enumerate(["voronoi", "length", "surface", "volume"]):  # ["length", "surface", "volume"]:
    for j, cell_number in enumerate(cell_numbers):

        print("\nIteration", i, j)
        data = outer_hist[i][j]
        print("data", data.shape, np.min(data), np.max(data), np.mean(data))
        k, theta = fit_gamma(data)
        print(k, theta)
