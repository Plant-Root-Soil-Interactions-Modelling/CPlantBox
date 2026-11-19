"""maize parameters quickly dervied form chatgpt
values from
Vidal et al. (2021), He et al. (2021), Fournier & Andrieu (2000, 2010, 2014), dos Santos et al. (2022)
"""

import numpy as np


#
# Leaf geometery
#
def get_leaf_sheath_length():
    leaf_sheath_length = np.array([0.0, 0.03, 0.04, 0.05, 0.07, 0.09, 0.12, 0.15, 0.18, 0.22, 0.25, 0.28, 0.30, 0.28, 0.25, 0.22, 0.18, 0.14, 0.10, 0.06])  # [m]
    return leaf_sheath_length * 100  # [m] -> [cm]


def get_leaf_blade_lenght():
    leaf_length = np.array(
        [
            0.0,  # one leaf more than internodes at the base
            0.35,
            0.35,
            0.35,  # 1–3: basal
            0.55,
            0.55,
            0.55,  # 4–6
            0.75,
            0.75,
            0.75,  # 7–9
            0.85,
            0.85,
            0.85,  # 10–12
            0.675,
            0.675,
            0.675,  # 13–15
            0.40,
            0.40,
            0.40,
            0.40,  # 16–19
        ]
    )  # [m]
    return leaf_length * 100  # [m] -> [cm]


def get_leaf_width():
    leaf_width = np.array(
        [
            0.0,  # one leaf more than internodes at the base
            0.045,
            0.045,
            0.045,  # 1–3
            0.065,
            0.065,
            0.065,  # 4–6
            0.085,
            0.085,
            0.085,  # 7–9
            0.10,
            0.10,
            0.10,  # 10–12
            0.08,
            0.08,
            0.08,  # 13–15
            0.055,
            0.055,
            0.055,
            0.055,  # 16–19
        ]
    )  # [m]
    return leaf_width * 100  # [m] -> [cm]


#
# Leaf emergence
#
def get_leaf_collar_emergence_gdd():
    leaf_collar_emergence_gdd = np.array([50, 108, 166, 224, 282, 340, 398, 456, 514, 572, 603, 634, 665, 696, 727, 758, 789, 820, 851])  # [°C day], base temperature 8°C
    return leaf_collar_emergence_gdd


def get_leaf_collar_emergence_days():
    leaf_collar_emergence_days = np.array([0, 5, 11, 17, 22, 28, 34, 40, 46, 51, 57, 60, 63, 67, 70, 73, 76, 79, 82, 85])  # [days after VE]
    return leaf_collar_emergence_days - 4 * np.ones(leaf_collar_emergence_days.shape)


#
# Stem geometry
#
def get_internode_lenght():
    internode_length = np.array(
        [
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,  # 1–5: basal
            0.085,
            0.085,
            0.085,
            0.085,  # 6–9: lower-middle
            0.16,
            0.16,
            0.16,
            0.16,  # 10–13: middle
            0.115,
            0.115,
            0.115,
            0.115,  # 14–17: upper
            0.02,
            0.02,  # 18–19: below tassel
        ]
    )  # [m]
    return internode_length * 100  # [m] -> [cm]


def get_internode_lenght2():
    internode_length2 = np.array(
        [0.005, 0.005, 0.005, 0.010, 0.020, 0.040, 0.060, 0.080, 0.100, 0.140, 0.160, 0.180, 0.180, 0.160, 0.140, 0.120, 0.100, 0.070, 0.040]  # 1–5: basal  # 6–9: lower-middle  # 10–13: middle  # 14–17: upper  # 18–19: below tassel
    )  # [m]
    return internode_length2 * 100  # [m] -> [cm]


#
# Crown root parameters
#
def get_crown_root_emergence():
    crown_root_emergence = np.array([12, 22, 29, 37, 45, 52])  # days after germination
    return crown_root_emergence


def get_crown_root_number():
    crown_root_number = np.array([5, 6, 8, 10, 12, 14])  # number of roots per whorl
    return crown_root_number


def get_crown_root_angles():
    crown_root_angles = np.array([25, 25, 30, 30, 35, 40])  # degrees below horizontal
    return crown_root_angles
