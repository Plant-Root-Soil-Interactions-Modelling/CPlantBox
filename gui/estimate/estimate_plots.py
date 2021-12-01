import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")
import numpy as np

import plantbox as pb


def baseroots(data):
    """ 
    analyser     EstimateDataModel
    """
    col = ["r*", "g*", "b*", "m*", "c*"] * 10

    for j in range(0, len(data.times)):
        bri = data.base_root_indices[i]
        for i in bri:  # len(names)
            lengths = properties[i][j]["length"]
            ages = properties[i][j]["age"]
            for k in range(1, len(polylines[i][j])):  # skip tap root
                age_.append(data.times[j])
                l_.append(lengths[k])
        plt.plot(l_, age_, col[j])  # scatter plot
    plt.ylabel("Root age [day]")
    plt.xlabel("Measured root length [cm]")

