import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")
import numpy as np

import plantbox as pb
from estimate_params import *


def plot_baseroots(data, base_method, ax):
    """ 
    data     EstimateDataModel
    base_method    base method index
    """
    ax.clear()

    # tap root data
    age_, l_ = [], []
    for i, j_ in enumerate(data.tap_root_indices):
        for j in j_:
            age_.append(data.estimates[i][(j, "age")])
            l_.append(data.rsmls[i].properties["length"][j])
            if age_[-1] > 6 and l_[-1] < 5:
                print("short! measurement", i, "root", j)
    ax.plot(age_, l_, "r*", label = "tap roots")  # scatter plot

    # basal data
    age_, l_ = [], []
    for i, j_ in enumerate(data.basal_root_indices):
        for j in j_:
            age_.append(data.estimates[i][(j, "age")])
            l_.append(data.rsmls[i].properties["length"][j])
    if base_method < 2:  # basal = tap
        ax.plot(age_, l_, "r*")  # scatter plot
    else:
        ax.plot(age_, l_, "b*", label = "basal roots")  # scatter plot

    # plotting fit
    k = data.parameters[0].lmax
    r = data.parameters[0].r
    max_time = np.max(data.times)
    t_ = np.linspace(0, max_time, 200)
    l_ = negexp_length(t_, r, k)
    ax.plot(t_, l_, label = "tap root fit")
    ax.set_xlabel("Root age [day]")
    ax.set_ylabel("Measured root length [cm]")
    ax.legend()


def plot_laterals(data, base_method, calibration_method, ax):
    """
    """
    col_ = ["y*", "r*", "b*", "g*", "c*", "m*"] * 2
    ax.clear()

    order = 1
    indices = data.pick_order(order)
    c = np.array([len(x) for x in indices])

    while np.sum(c) > 0:
        age_ = []
        l_ = []
        for i, j_ in enumerate(indices):
            for j in j_:
                age_.append(data.estimates[i][(j, "age")])
                l_.append(data.rsmls[i].properties["length"][j])
        ax.plot(age_, l_, col_[order], label = "{:g} order lateras".format(order))
        order += 1
        indices = data.pick_order(order)  # update index set (TODO it must be always per order, but different target_types are possible for clustering and aggregation)
        c = np.array([len(x) for x in indices])
    # plotting fit
    names_ = ["", "1st order", "2nd order"]  # TODO clustering
    col_ = ["y", "r", "b"]
    for i in range(1, 3):
        k = data.parameters[i].lmax
        r = data.parameters[i].r
        max_time = np.max(data.times)
        t_ = np.linspace(0, max_time, 200)
        l_ = negexp_length(t_, r, k)
        ax.plot(t_, l_, col_[i], label = names_[i])
    ax.set_xlabel("Estimated root age [day]")
    ax.set_ylabel("Measured root length [cm]")
    ax.legend()
