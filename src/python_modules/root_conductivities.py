import numpy as np
from scipy import interpolate

"""
Helper to define age depent tabular values for root conductivities for XylemFluxPython (values are hard coded)

TODO this should be improved !

after calling:
    init_conductivities(r :XylemFluxPython, age_dependent :bool)

use:
r.kr_f(age, type)
r.kx_f(age, type) 

"""

# values for age indepenent case
kx_const_ = 4.32e-2  # [cm3/day]
kr_const_ = 1.73e-4  # [1/day]

# tabular values for age and type depenent case
kr0 = np.array([[-1e4, 1.e-9], [0., 1.e-9], [0, 1.14e-03], [2, 1.09e-03], [4, 1.03e-03], [6, 9.83e-04], [8, 9.35e-04], [10, 8.90e-04], [12, 8.47e-04], [14, 8.06e-04], [16, 7.67e-04], [18, 7.30e-04], [20, 6.95e-04], [22, 6.62e-04], [24, 6.30e-04], [26, 5.99e-04], [28, 5.70e-04], [30, 5.43e-04], [32, 5.17e-04]])
kr1 = np.array([[-1e4, 1.e-9], [0., 1.e-9], [0, 4.11e-03], [1, 3.89e-03], [2, 3.67e-03], [3, 3.47e-03], [4, 3.28e-03], [5, 3.10e-03], [6, 2.93e-03], [7, 2.77e-03], [8, 2.62e-03], [9, 2.48e-03], [10, 2.34e-03], [11, 2.21e-03], [12, 2.09e-03], [13, 1.98e-03], [14, 1.87e-03], [15, 1.77e-03], [16, 1.67e-03], [17, 1.58e-03]])

kx0 = np.array([[0, 6.74e-02], [2, 7.48e-02], [4, 8.30e-02], [6, 9.21e-02], [8, 1.02e-01], [10, 1.13e-01], [12, 1.26e-01], [14, 1.40e-01], [16, 1.55e-01], [18, 1.72e-01], [20, 1.91e-01], [22, 2.12e-01], [24, 2.35e-01], [26, 2.61e-01], [28, 2.90e-01], [30, 3.21e-01], [32, 3.57e-01]])
kx1 = np.array([[0, 4.07e-04], [1, 5.00e-04], [2, 6.15e-04], [3, 7.56e-04], [4, 9.30e-04], [5, 1.14e-03], [6, 1.41e-03], [7, 1.73e-03], [8, 2.12e-03], [9, 2.61e-03], [10, 3.21e-03], [11, 3.95e-03], [12, 4.86e-03], [13, 5.97e-03], [14, 7.34e-03], [15, 9.03e-03], [16, 1.11e-02], [17, 1.36e-02]])


def init_conductivities(r, age_dependent :bool = False):
    """ call to initialize age dependent or independent conductivities, 
    initializes functions kr(age, type) and kx(age, type) """

    if age_dependent:
        r.setKrTables([kr0[:, 1], kr1[:, 1], kr1[:, 1], kr1[:, 1], kr1[:, 1], kr1[:, 1]],
                      [kr0[:, 0], kr1[:, 0], kr1[:, 0], kr1[:, 0], kr1[:, 0], kr1[:, 0]])  # [cm^3/day]
        r.setKxTables([kx0[:, 1], kx1[:, 1], kx1[:, 1], kx1[:, 1], kx1[:, 1], kx1[:, 1]],
                      [kx0[:, 0], kx1[:, 0], kx1[:, 0], kx1[:, 0], kx1[:, 0], kx1[:, 0]])  # [1/day]

    else:  # we set it as table to be able to make the rootsystem grow in a predefined way
        kr = np.array([[-1e4, 1.e-9], [0., 1.e-9], [0., kr_const_], [1e4, kr_const_]])
        kx = np.array([[0, kx_const_], [1e4, kx_const_]])
        r.setKrTables([kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1], kr[:, 1]],
                      [kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0], kr[:, 0]])  # [cm^3/day]
        r.setKxTables([kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1], kx[:, 1]],
                      [kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0], kx[:, 0]])  # [1/day]
#         r.setKr([kr_const_])  # [cm^3/day] # todo check order of age (segs are sorted)
#         r.setKx([kx_const_])  # [1/day]

# kr = lambda age, type : kr_const_
# kx = lambda age, type : kx_const_
#
# krf_ = [interpolate.interp1d(kr0[:, 0], kr0[:, 1]), interpolate.interp1d(kr1[:, 0], kr1[:, 1]) ]
# kxf_ = [interpolate.interp1d(kx0[:, 0], kx0[:, 1]), interpolate.interp1d(kx1[:, 0], kx1[:, 1]) ]
#
#
# def kr_table_(age :float, type :int):
#     if type == 0:
#         return krf_[0](age, int)
#     else:
#         return krf_[1](age, int)
#
#
# def kx_table_(age :float, type :int):
#     if type == 0:
#         return kxf_[0](age, int)
#     else:
#         return kxf_[1](age, int)

