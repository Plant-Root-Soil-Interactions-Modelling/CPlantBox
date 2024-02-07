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
            l_.append(data.rsmls[i].properties["length"][j])  # /float(data.rsmls[i].metadata.resolution)
            # if age_[-1] > 6 and l_[-1] < 5:
            #     print("short! measurement", i, "root", j)
    ax.plot(age_, l_, "r*", label = "Tap roots")  # scatter plot

    # # basal data
    # age_, l_ = [], []
    # for i, j_ in enumerate(data.basal_root_indices):
    #     for j in j_:
    #         age_.append(data.estimates[i][(j, "age")])
    #         l_.append(data.rsmls[i].properties["length"][j])  # /float(data.rsmls[i].metadata.resolution)
    # # l_ = np.sort(l_) #shortest basal root is the youngest
    # if base_method < 2:  # basal = tap
    #     ax.plot(age_, l_, "r*")  # scatter plot
    # else:
    #     ax.plot(age_, l_, "b*", label = "basal roots")  # scatter plot

    # plotting fit
    k = data.parameters[0].lmax
    r = data.parameters[0].r
    max_time = np.max(data.times)
    t_ = np.linspace(0, max_time, 200)
    l_ = negexp_length(t_, r, k)
    ax.plot(t_, l_, label = "Tap root fit")

    # # TODO hack
    # k = 100  # data.parameters[0].lmax
    # r = 1.129  # data.parameters[0].r
    # max_time = np.max(data.times)
    # t_ = np.linspace(0, max_time, 200)
    # l_ = negexp_length(t_, r, k)
    # ax.plot(t_, l_, label = "Tap root fit (lmax fixed)")

    ax.set_xlabel("Root age [day]")
    ax.set_ylabel("Measured root length [cm]")
    ax.legend()


def plot_laterals(data, base_method, calibration_method, ax, orders_ = [1, 2, 3]):
    """
    """
    names_ = ["", "1st order", "2nd order", "3nd order", "higher order",
              "higher order", "higher order", "higher order", "higher order", "higher order", "higher order"]
    col_ = ["y*", "r*", "b*", "g*", "c*", "m*"] * 2
    ax.clear()

    for order in orders_:

        indices = data.pick_order(order)
        c = np.array([len(x) for x in indices])
        print("plot_laterals(): number of order", order, "roots", c)
        age_ = []
        l_ = []
        for i, j_ in enumerate(indices):
            for j in j_:
                age_.append(data.estimates[i][(j, "age")])
                l_.append(data.rsmls[i].properties["length"][j])
        ax.plot(age_, l_, col_[order], label = names_[order])  # , alpha = i / len(data.times) * 0.8 + 0.2

    # plotting fit
    col_ = ["y", "r", "b", "g", "c", "m"] * 2
    for i in orders_:
        k = data.parameters[i].lmax
        r = data.parameters[i].r
        print("plot_laterals", k, r)
        max_time = np.max(age_)
        t_ = np.linspace(0, max_time, 200)
        l_ = negexp_length(t_, r, k)
        ax.plot(t_, l_, col_[i], label = names_[i] + " fit")

    ax.set_xlabel("Estimated root age [day]")
    ax.set_ylabel("Measured root length [cm]")
    ax.legend()

