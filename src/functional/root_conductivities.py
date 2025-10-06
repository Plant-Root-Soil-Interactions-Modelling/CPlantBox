import numpy as np
from scipy import interpolate

from functional.PlantHydraulicParameters import PlantHydraulicParameters

"""
Helper to define age depent tabular values for root conductivities for XylemFluxPython (values are hard coded)

TODO this should be improved !

usage:  after calling:
init_conductivities(r :XylemFluxPython, age_dependent :bool)

use:
r.kr_f(age, type)
r.kx_f(age, type) 

"""


def convert_axial(kx):
    """ converts axial conductivity with units of [m4 s-1 MPa-1] to [cm3 / day] """
    rho = 1000  # [kg/m3]
    g = 9.81  # [m s-2]
    rho_g = rho * g  # [kg m-2 s-2] = [Pa m-1]
    kx_ = kx * rho_g  # m4 s-1 1e-6 Pa-1 * Pa m-1 = m3 s-1 1.e-6 = cm3 s-1
    kx_ = kx_ * 60.*60.*24  #  cm3 s-1 = 60*60*24 cm3 /day
    return kx_


def convert_radial(kr):
    """ converts radia conductivity with units of [m s-1 MPa-1] to [1 / day] """
    rho = 1000  # [kg/m3]
    g = 9.81  # [m s-2]
    rho_g = rho * g  # [kg m-2 s-2] = [Pa m-1]
    kr_ = kr * rho_g  # m s-1 1e-6 Pa-1 * Pa m-1 = s-1 1.e-6 = s-1 1.e-6
    kr_ = kr_ * 60.*60.*24 * 1.e-6  # s-1 = 60*60*24 /day
    return kr_


def init_conductivities(r, age_dependent:bool = False):
    """ call to initialize age dependent or independent conductivities, 
    initializes functions kr(age, type) and kx(age, type) """

    kx_const_ = 4.32e-2  # [cm3/day] # values for age indepenent case
    kr_const_ = 1.73e-4  # [1/day]

    # tabular values for age and type depenent case [age, value] for type 0 (kr0), and type 1 (kr1)
    kr0 = np.array([[-1e4, 1.e-9], [0., 1.e-9], [0, 1.14e-03], [2, 1.09e-03], [4, 1.03e-03], [6, 9.83e-04], [8, 9.35e-04], [10, 8.90e-04], [12, 8.47e-04], [14, 8.06e-04], [16, 7.67e-04], [18, 7.30e-04], [20, 6.95e-04], [22, 6.62e-04], [24, 6.30e-04], [26, 5.99e-04], [28, 5.70e-04], [30, 5.43e-04], [32, 5.17e-04]])
    kr1 = np.array([[-1e4, 1.e-9], [0., 1.e-9], [0, 4.11e-03], [1, 3.89e-03], [2, 3.67e-03], [3, 3.47e-03], [4, 3.28e-03], [5, 3.10e-03], [6, 2.93e-03], [7, 2.77e-03], [8, 2.62e-03], [9, 2.48e-03], [10, 2.34e-03], [11, 2.21e-03], [12, 2.09e-03], [13, 1.98e-03], [14, 1.87e-03], [15, 1.77e-03], [16, 1.67e-03], [17, 1.58e-03]])

    kx0 = np.array([[0, 6.74e-02], [2, 7.48e-02], [4, 8.30e-02], [6, 9.21e-02], [8, 1.02e-01], [10, 1.13e-01], [12, 1.26e-01], [14, 1.40e-01], [16, 1.55e-01], [18, 1.72e-01], [20, 1.91e-01], [22, 2.12e-01], [24, 2.35e-01], [26, 2.61e-01], [28, 2.90e-01], [30, 3.21e-01], [32, 3.57e-01]])
    kx1 = np.array([[0, 4.07e-04], [1, 5.00e-04], [2, 6.15e-04], [3, 7.56e-04], [4, 9.30e-04], [5, 1.14e-03], [6, 1.41e-03], [7, 1.73e-03], [8, 2.12e-03], [9, 2.61e-03], [10, 3.21e-03], [11, 3.95e-03], [12, 4.86e-03], [13, 5.97e-03], [14, 7.34e-03], [15, 9.03e-03], [16, 1.11e-02], [17, 1.36e-02]])

    if isinstance(r, PlantHydraulicParameters):
        if age_dependent:
            r.set_kr_age(kr0[:, 0], kr0[:, 1], subType = [0, 1])
            r.set_kr_age(kr1[:, 0], kr1[:, 1], subType = [2, 3, 4, 5])
            r.set_kx_age(kx0[:, 0], kx0[:, 1], subType = [0, 1])
            r.set_kx_age(kx1[:, 0], kx1[:, 1], subType = [2, 3, 4, 5])
        else:  # we set it as table to be able to make the rootsystem grow in a predefined way
            kr = np.array([[-1e4, 1.e-9], [-1.e-9, 1.e-9], [0., kr_const_], [1e4, kr_const_]])
            kx = np.array([[0, kx_const_], [1e4, kx_const_]])
            r.set_kr_age_dependent(kr[:, 0], kr[:, 1])
            r.set_kx_age_dependent(kx[:, 0], kx[:, 1])
    else:

        if age_dependent:
            r.setKrTables([kr0[:, 1], kr0[:, 1], kr1[:, 1], kr1[:, 1], kr1[:, 1], kr1[:, 1]],
                          [kr0[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr1[:, 0], kr1[:, 0]])  # [cm^3/day]
            r.setKxTables([kx0[:, 1], kx0[:, 1], kx1[:, 1], kx1[:, 1], kx1[:, 1], kx1[:, 1]],
                          [kx0[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx1[:, 0], kx1[:, 0]])  # [1/day]

        else:  # we set it as table to be able to make the rootsystem grow in a predefined way
            kr = np.array([[-1e4, 1.e-9], [-1.e-9, 1.e-9], [0., kr_const_], [1e4, kr_const_]])
            kx = np.array([[0, kx_const_], [1e4, kx_const_]])
            r.setKrTables([kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1]],
                          [kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0]])  # [cm^3/day]
            r.setKxTables([kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1]],
                          [kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0]])  # [1/day]


def init_conductivities_growth(r, age_dependent:bool = False, dt = 1):
    """ same as init_conductivities but with a 1 day slope if segments emerge
        @param dt     time span after emergence to reach full kr value
    """

    # values for age indepenent case
    kx_const_ = 4.32e-2  # [cm3/day]
    kr_const_ = 1.73e-4  # [1/day]

    # tabular values for age and type depenent case [age, value] for type 0 (kr0), and type 1 (kr1)
    kr0 = np.array([[-1e4, 1.e-9], [0., 1.e-9], [0, 1.14e-03], [2, 1.09e-03], [4, 1.03e-03], [6, 9.83e-04], [8, 9.35e-04], [10, 8.90e-04], [12, 8.47e-04], [14, 8.06e-04], [16, 7.67e-04], [18, 7.30e-04], [20, 6.95e-04], [22, 6.62e-04], [24, 6.30e-04], [26, 5.99e-04], [28, 5.70e-04], [30, 5.43e-04], [32, 5.17e-04]])
    kr1 = np.array([[-1e4, 1.e-9], [0., 1.e-9], [0, 4.11e-03], [1, 3.89e-03], [2, 3.67e-03], [3, 3.47e-03], [4, 3.28e-03], [5, 3.10e-03], [6, 2.93e-03], [7, 2.77e-03], [8, 2.62e-03], [9, 2.48e-03], [10, 2.34e-03], [11, 2.21e-03], [12, 2.09e-03], [13, 1.98e-03], [14, 1.87e-03], [15, 1.77e-03], [16, 1.67e-03], [17, 1.58e-03]])

    for i in range(2, kr0.shape[0]):
        kr0[i, 0] = kr0[i, 0] + dt
        kr1[i, 0] = kr1[i, 0] + dt

    kx0 = np.array([[0, 6.74e-02], [2, 7.48e-02], [4, 8.30e-02], [6, 9.21e-02], [8, 1.02e-01], [10, 1.13e-01], [12, 1.26e-01], [14, 1.40e-01], [16, 1.55e-01], [18, 1.72e-01], [20, 1.91e-01], [22, 2.12e-01], [24, 2.35e-01], [26, 2.61e-01], [28, 2.90e-01], [30, 3.21e-01], [32, 3.57e-01]])
    kx1 = np.array([[0, 4.07e-04], [1, 5.00e-04], [2, 6.15e-04], [3, 7.56e-04], [4, 9.30e-04], [5, 1.14e-03], [6, 1.41e-03], [7, 1.73e-03], [8, 2.12e-03], [9, 2.61e-03], [10, 3.21e-03], [11, 3.95e-03], [12, 4.86e-03], [13, 5.97e-03], [14, 7.34e-03], [15, 9.03e-03], [16, 1.11e-02], [17, 1.36e-02]])

    if age_dependent:
        r.setKrTables([kr0[:, 1], kr1[:, 1], kr1[:, 1], kr1[:, 1], kr1[:, 1], kr1[:, 1]],
                      [kr0[:, 0], kr1[:, 0], kr1[:, 0], kr1[:, 0], kr1[:, 0], kr1[:, 0]])  # [cm^3/day]
        r.setKxTables([kx0[:, 1], kx1[:, 1], kx1[:, 1], kx1[:, 1], kx1[:, 1], kx1[:, 1]],
                      [kx0[:, 0], kx1[:, 0], kx1[:, 0], kx1[:, 0], kx1[:, 0], kx1[:, 0]])  # [1/day]

    else:  # we set it as table to be able to make the rootsystem grow in a predefined way
        kr = np.array([[-1e4, 1.e-9], [0., 1.e-9], [dt, kr_const_], [1e4, kr_const_]])
        kx = np.array([[0, kx_const_], [1e4, kx_const_]])
        r.setKrTables([kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1]],
                      [kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0]])  # [cm^3/day]
        r.setKxTables([kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1]],
                      [kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0]])  # [1/day]


def init_conductivities_scenario_jan(r):
    """ Hydraulic conductivities - for Jans scenarios """

    # radial values # [cm^3/day]
    # (ageLr(1,i),LrRoot(1,i),i=1,nLr(1..3))
    kr0 = np.array([[0., 0.000181], [8., 0.000181], [10, 0.0000648], [18, 0.0000648], [25, 0.0000173], [300, 0.0000173]])
    kr1 = np.array([[0., 0.000181], [10., 0.000181], [16, 0.0000173], [300, 0.0000173]])
    kr2 = np.array([[0., 0.000181], [10., 0.000181], [16, 0.0000173], [300, 0.0000173]])

    # axial values # [1/day]
    # (ageKh(1,i),KhRoot(1,i),i=1,nKh(1..3))
    kx0 = np.array([[0., 0.0000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [300., 0.432]])
    kx1 = np.array([[0., 0.0000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])
    kx2 = np.array([[0., 0.0000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])

    # first entry twice, because rsml starts at type 1
    r.setKrTables([kr0[:, 1], kr1[:, 1], kr2[:, 1]], [kr0[:, 0], kr1[:, 0], kr2[:, 0]])  # values, age
    r.setKxTables([kx0[:, 1], kx1[:, 1], kx2[:, 1]], [kx0[:, 0], kx1[:, 0], kx2[:, 0]])  # values, age


def init_conductivities_scenario_jan_const(r):
        """ Hydraulic conductivities - for Jans scenarios, BUT with constanct kr (if we use a lookup table) """

        # kr_const = 6.48e-5
        kr_const = 1.73e-4  # in case of table look up, the values must agree
        kr = np.array([[0., kr_const], [1e4, kr_const]])
        r.setKrTables([kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1]],
                      [kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0]])  # [cm^3/day]
        # first entry twice, because rsml starts at type 1

        kx_const = 2.95e-2  # [cm3/day]
        kx = np.array([[0., kx_const], [1e4, kx_const]])

        r.setKxTables([kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1]],
                      [kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0]])  # values, age


def init_singleroot_contkrkx(r):
    """ Hydraulic conductivities - for Jans scenarios, BUT with constanct kr (if we use a lookup table) """
    kr_const = 2.843148e-5 / (2 * np.pi * 0.05 * 0.5)  # in case of table look up, the values must agree, = 0.00018100042

    kr = np.array([[0., kr_const], [1e4, kr_const]])
    r.setKrTables([kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1]],
                  [kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0]])  # [cm^3/day]

    # first entry twice, because rsml starts at type 1
    kx_const = 0.346 * 0.5  # [cm3/day] 1e6 *

    kx = np.array([[0., kx_const], [1e4, kx_const]])

    r.setKxTables([kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1]],
                  [kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0]])  # values, age

