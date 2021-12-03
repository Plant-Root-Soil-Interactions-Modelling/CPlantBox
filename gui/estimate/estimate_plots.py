import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")
import numpy as np

import plantbox as pb


def plot_baseroots(data, ax):
    """ 
    data     EstimateDataModel
    """
    for j in range(0, len(data.times)):  # plants

        bri = data.base_root_indices[j]

        l_, age_ = [], []
        i = bri[0]  # tap roots
        age_.append(data.times[j])
        l_.append(data.rsmls[j].properties["length"][i])
        ax.plot(age_, l_, "r*")  # scatter plot

        l_, age_ = [], []
        for i in bri[1:]:  # other base roots
            print(i)
            age_.append(data.times[j])
            l_.append(data.rsmls[j].properties["length"][i])
        ax.plot(age_, l_, "b*")  # scatter plot

    ax.set_xlabel("Root age [day]")
    ax.set_ylabel("Measured root length [cm]")

