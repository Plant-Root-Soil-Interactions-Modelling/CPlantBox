"""
Helper to define age depent tabular values for root conductivities for XylemFluxPython (values are hard coded)

usage:  after calling:
init_conductivities(r :XylemFluxPython, age_dependent :bool)

use:
r.kr_f(age, type)
r.kx_f(age, type)
"""

import numpy as np


def convert_axial(kx):
    """converts axial conductivity with units of [m4 s-1 MPa-1] to [cm3 / day]"""
    rho = 1000  # [kg/m3]
    g = 9.81  # [m s-2]
    rho_g = rho * g  # [kg m-2 s-2] = [Pa m-1]
    kx_ = kx * rho_g  # m4 s-1 1e-6 Pa-1 * Pa m-1 = m3 s-1 1.e-6 = cm3 s-1
    kx_ = kx_ * 60.0 * 60.0 * 24  #  cm3 s-1 = 60*60*24 cm3 /day
    return kx_


def convert_radial(kr):
    """converts radia conductivity with units of [m s-1 MPa-1] to [1 / day]"""
    rho = 1000  # [kg/m3]
    g = 9.81  # [m s-2]
    rho_g = rho * g  # [kg m-2 s-2] = [Pa m-1]
    kr_ = kr * rho_g  # m s-1 1e-6 Pa-1 * Pa m-1 = s-1 1.e-6 = s-1 1.e-6
    kr_ = kr_ * 60.0 * 60.0 * 24 * 1.0e-6  # s-1 = 60*60*24 /day
    return kr_


def add_shoot_conductivity(params):
    """adds shoot conductivities to the tables, at subType 10 (only for constant scenarios!)"""
    params.set_kr_const(0.0, subType=10)
    params.set_kx_const(1.0, subType=10)


def init_constant_scenario1(params):
    """call to initialize age dependent or independent conductivities,
    initializes functions kr(age, type) and kx(age, type)"""
    kr_const = 1.73e-4  # [1/day]
    params.set_kr_const(kr_const)
    kx_const = 4.32e-2  # [cm3/day]
    params.set_kx_const(kx_const)
    add_shoot_conductivity(params)  # replace subType 10


def init_constant_scenario2(params):
    """call to initialize age dependent or independent conductivities,
    initializes functions kr(age, type) and kx(age, type)"""
    kr0 = 0.1 * 1.73e-4  # [1/day] # 0-order, mainly transport, smaller radial fluxes
    kr1 = 1.73e-4  # [1/day]
    kx0 = 10 * 4.32e-2  # [cm3/day]
    kx1 = 4.32e-2  # [cm3/day]
    params.set_kr_const(kr0, subType=[1, 4, 5])
    params.set_kr_const(kr1, subType=[2, 3])
    params.set_kx_const(kx0, subType=[1, 4, 5])
    params.set_kx_const(kx1, subType=[2, 3])
    add_shoot_conductivity(params)


def init_constant_scenario_wine(params):
    """call to initialize age dependent or independent conductivities,
    initializes functions kr(age, type) and kx(age, type)"""
    kx0 = 1000
    kx1 = 213.868458207658  # cm3 day-1
    kx2 = 53.8961351157623
    kx3 = 6.06479773916046
    kr0 = 0.0
    kr1 = 3.39e-05  # day-1
    kr2 = 5.93e-05
    kr3 = 3.39e-04
    params.set_kr_const(kr0, subType=1)
    params.set_kr_const(kr1, subType=2)
    params.set_kr_const(kr2, subType=3)
    params.set_kr_const(kr3, subType=4)
    params.set_kx_const(kx0, subType=1)
    params.set_kx_const(kx1, subType=2)
    params.set_kx_const(kx2, subType=3)
    params.set_kx_const(kx3, subType=4)
    add_shoot_conductivity(params)


def init_dynamic_scenario1(params):
    """call to initialize age dependent or independent conductivities,
    initializes functions kr(age, type) and kx(age, type)"""

    # tabular values for age and type depenent case [age, value] for type 0 (kr0), and type 1 (kr1)
    kr0_ = np.array(
        [
            [-1e4, 0.0],
            [-1.0e-9, 0.0],
            [0, 1.14e-03],
            [2, 1.09e-03],
            [4, 1.03e-03],
            [6, 9.83e-04],
            [8, 9.35e-04],
            [10, 8.90e-04],
            [12, 8.47e-04],
            [14, 8.06e-04],
            [16, 7.67e-04],
            [18, 7.30e-04],
            [20, 6.95e-04],
            [22, 6.62e-04],
            [24, 6.30e-04],
            [26, 5.99e-04],
            [28, 5.70e-04],
            [30, 5.43e-04],
            [32, 5.17e-04],
        ]
    )
    kr0_age = list(kr0_[:, 0])
    kr0_value = list(kr0_[:, 1])
    kr1_ = np.array(
        [
            [-1e4, 0.0],
            [-1.0e-9, 0.0],
            [0, 4.11e-03],
            [1, 3.89e-03],
            [2, 3.67e-03],
            [3, 3.47e-03],
            [4, 3.28e-03],
            [5, 3.10e-03],
            [6, 2.93e-03],
            [7, 2.77e-03],
            [8, 2.62e-03],
            [9, 2.48e-03],
            [10, 2.34e-03],
            [11, 2.21e-03],
            [12, 2.09e-03],
            [13, 1.98e-03],
            [14, 1.87e-03],
            [15, 1.77e-03],
            [16, 1.67e-03],
            [17, 1.58e-03],
        ]
    )
    kr1_age = list(kr1_[:, 0])
    kr1_value = list(kr1_[:, 1])
    kx0_ = np.array(
        [
            [0, 6.74e-02],
            [2, 7.48e-02],
            [4, 8.30e-02],
            [6, 9.21e-02],
            [8, 1.02e-01],
            [10, 1.13e-01],
            [12, 1.26e-01],
            [14, 1.40e-01],
            [16, 1.55e-01],
            [18, 1.72e-01],
            [20, 1.91e-01],
            [22, 2.12e-01],
            [24, 2.35e-01],
            [26, 2.61e-01],
            [28, 2.90e-01],
            [30, 3.21e-01],
            [32, 3.57e-01],
        ]
    )
    kx0_age = list(kx0_[:, 0])
    kx0_value = list(kx0_[:, 1])
    kx1_ = np.array(
        [
            [0, 4.07e-04],
            [1, 5.00e-04],
            [2, 6.15e-04],
            [3, 7.56e-04],
            [4, 9.30e-04],
            [5, 1.14e-03],
            [6, 1.41e-03],
            [7, 1.73e-03],
            [8, 2.12e-03],
            [9, 2.61e-03],
            [10, 3.21e-03],
            [11, 3.95e-03],
            [12, 4.86e-03],
            [13, 5.97e-03],
            [14, 7.34e-03],
            [15, 9.03e-03],
            [16, 1.11e-02],
            [17, 1.36e-02],
        ]
    )
    kx1_age = list(kx1_[:, 0])
    kx1_value = list(kx1_[:, 1])
    params.set_kr_age_dependent(kr0_age, kr0_value, subType=[1, 4, 5])
    params.set_kr_age_dependent(kr1_age, kr1_value, subType=[2, 3])
    params.set_kx_age_dependent(kx0_age, kx0_value, subType=[1, 4, 5])
    params.set_kx_age_dependent(kx1_age, kx1_value, subType=[2, 3])
    # add_shoot_conductivity for age dependent:
    params.set_kr_age_dependent([-1.0, 1.0e6], [0.0, 0.0], subType=10)
    params.set_kx_age_dependent([-1.0, 1.0e6], [1.0, 1.0], subType=10)


def init_dynamic_scenario2(params):
    """call to initialize age dependent or independent conductivities,
    initializes functions kr(age, type) and kx(age, type)"""

    # TODO
    init_constant_scenario1(params)  # for now, use constant values for dynamic scenario 2, to be replaced with tabular values


if __name__ == "__main__":

    print("Radial", convert_radial(4.0e-8))
    print("Radial", convert_radial(7.0e-8))
    print("Radial", convert_radial(4.0e-7))
    print()
    print("Axial E", convert_axial(2.52327153659882e-07))
    print("Axial E", convert_axial(6.35879572004218e-08))
    print("Axial E", convert_axial(7.15539431980838e-09))
