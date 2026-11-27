"""maize parameters quickly dervied form chatgpt
values from
Vidal et al. (2021), He et al. (2021), Fournier & Andrieu (2000, 2010, 2014), dos Santos et al. (2022)
"""

import matplotlib.pyplot as plt
import numpy as np

import plantbox as pb


class ExoGraviTropism(pb.Tropism):
    """User defined tropism depending on root age use either exotropism or gravitropism"""

    def __init__(self, plant, n, sigma, age):
        super(ExoGraviTropism, self).__init__(plant)
        self.exo = pb.Exotropism(plant, 0.0, 0.0)
        self.gravi = pb.Gravitropism(plant, 0.0, 0.0)
        self.setTropismParameter(n, sigma)
        self.age = age

    def tropismObjective(self, pos, old, a, b, dx, root):
        # age = root.getAge()
        age = root.calcAge(root.getLength())  # works also with large time steps
        if age < self.age:
            return self.exo.tropismObjective(pos, old, a, b, dx, root)
        else:
            return self.gravi.tropismObjective(pos, old, a, b, dx, root)


#
# Leaf geometery
#
def get_leaf_sheath_length():
    leaf_sheath_length = np.array([0.0, 0.03, 0.04, 0.05, 0.07, 0.09, 0.12, 0.15, 0.18, 0.22, 0.25, 0.28, 0.30, 0.28, 0.25, 0.22, 0.18, 0.10, 0.09, 0.07])  # [m] last
    return leaf_sheath_length * 100  # [m] -> [cm]


def get_leaf_blade_lenght():
    leaf_length = np.array(
        [
            0.0,  # one leaf more than internodes at the base
            0.5 * 0.35,
            0.5 * 0.35,
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
def get_leaf_emergence_gdd():
    leaf_emergence_gdd = np.array([50, 108, 166, 224, 282, 340, 398, 456, 514, 572, 603, 634, 665, 696, 727, 758, 789, 820, 851])  # [°C day], base temperature 8°C
    return leaf_emergence_gdd


def get_leaf_emergence_days():
    leaf_emergence_days = np.array([0, 5, 11, 17, 22, 28, 34, 40, 46, 51, 57, 60, 63, 67, 70, 73, 76, 79, 82, 85])  # [days after VE]
    internode_final_time_days = np.array(
        [0, 4.78481013, 6.17721519, 7.56962025, 8.96202532, 14.87974684, 20.79746835, 26.71518987, 32.63291139, 43.7721519, 54.91139241, 66.05063291, 77.18987342, 85.19620253, 93.20253165, 101.20886076, 109.21518987, 110.60759494, 112.0, 200]
    )  # quickfix
    return np.minimum(leaf_emergence_days - 2 * np.ones(leaf_emergence_days.shape), internode_final_time_days)
    # return leaf_emergence_days


def get_leaf_emergence_days_internodes():
    lmax = 197  # total stem length [cm]
    r = lmax / 110  # assuming lmax is reached after 110 days
    ln = [0.0]
    ln = ln + list(np.cumsum(get_internode_lenght2()))
    dt = np.array(ln) / r  # arrival times
    return dt


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
            0.03,
            0.03,  # 18–19: below tassel
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
    crown_root_angles = np.array([25, 30, 35, 40, 45, 50])  # degrees below horizontal
    return crown_root_angles


def get_crown_root_radii():
    basal_radius_min = np.array([0.40, 0.55, 0.75, 0.95, 1.20, 1.60])
    basal_radius_max = np.array([0.65, 0.85, 1.10, 1.35, 1.65, 2.50])
    basal_radius_mean = (basal_radius_min + basal_radius_max) / 2.0
    return basal_radius_mean / 10.0  # mm -> cm


if __name__ == "__main__":

    ln = get_internode_lenght()
    ln2 = get_internode_lenght2()
    plt.ylim(0, max(ln) * 1.1)
    plt.plot(ln, "*", label="parametrisation 1")
    plt.plot(ln2, "*", label="parametrisation 2")
    plt.xlabel("Internode number")
    plt.ylabel("Internode length [cm]")
    plt.legend()
    plt.show()

    bl = get_leaf_blade_lenght()
    sln = np.cumsum(ln)
    plt.plot(bl[1:], sln, "g*")
    plt.xlabel("Blade length [cm]")
    plt.ylabel("Cumulative internode length [cm]")
    plt.show()
