import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")
import numpy as np

import plantbox as pb
from estimate_params import *


def plot_baseroots(data, index, ax):
    """ 
    data     EstimateDataModel
    """
    ax.clear()  # plotting raw data

    # tap root data
    length_tap, age_tap = [], []
    for i, j in enumerate(data.tap_root_indices):
        age_ = data.times[i]
        l_ = data.rsmls[i].properties["length"][j[0]]
        age_tap.append(age_)  # for fitting
        length_tap.append(l_)  # for fitting
        if i == 0:
            ax.plot(age_, l_, "r*", label = "tap roots")  # scatter plot
        else:
            ax.plot(age_, l_, "r*")  # scatter plot
    # print("tap roots")
    # print(data.tap_root_indices)
    # print(age_tap)
    # print(length_tap)

    # basal data
    length_basals, age_basals = [], []
    for i in range(0, len(data.times)):  # plants
        l_, age_ = [], []
        for j in data.basal_root_indices[i]:  # other base roots
            age_.append(data.times[i])
            l_.append(data.rsmls[i].properties["length"][j])
            age_basals.append(data.times[i])
            length_basals.append(data.rsmls[i].properties["length"][j])
        if j == 0:
            if index < 2:  # basal = tap
                ax.plot(age_, l_, "r*")  # scatter plot
            else:
                ax.plot(age_, l_, "b*", label = "basal roots")  # scatter plot
        else:
            if index < 2:  # basal = tap
                ax.plot(age_, l_, "r*")  # scatter plot
            else:
                ax.plot(age_, l_, "b*")  # scatter plot

    # fitting
    if index == 0:
        r, k, res = fit_taproot_rk([*length_tap, *length_basals], [*age_tap, *age_basals])
        print("r", r, "k", k, "res", res)
    elif index == 1:
        k = 100
        r, res = fit_taproot_r([*length_tap, *length_basals], [*age_tap, *age_basals], k)
        print("r", r, "k", k, "res", res)

    elif index == 2:  # TODO
        r, k, res_ = fit_taproot_rk(length_tap, age_tap)
        res, f = estiamte_emergance_order0(np.array(length_basals), np.array(age_basals), r, k)
        print(res)
        print(res.x[0])

    elif index == 3:  # TODO
        k = 100
        r, res_ = fit_taproot_r(length_tap, age_tap, k)
        res, f = estiamte_emergance_order0(np.array(length_basals), np.array(age_basals), r, k)
        print(res)
        print(res.x[0])

    data.set_rk(r, k)

    # plotting fit
    max_time = np.max(data.times) - 0.5
    t_ = np.linspace(0, max_time, 200)
    l_ = negexp_length(t_, r, k)
    ax.plot(t_, l_, label = "tap root fit")
    ax.set_xlabel("Root age [day]")
    ax.set_ylabel("Measured root length [cm]")
    ax.legend()

