"""maize parameters quickly dervied form chatgpt
values from
Vidal et al. (2021), He et al. (2021)
"""

import numpy as np


def get_leaf_sheath_length():
    leaf_sheath_length = np.array(
        [0.0, 0.03, 0.04, 0.05, 0.07, 0.09, 0.12, 0.15, 0.18, 0.22, 0.25, 0.28, 0.30, 0.28, 0.25, 0.22, 0.18, 0.14, 0.10, 0.06]
    )  # [m]
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
