#
# Mualem - van Genuchten model, equations from van Genuchten, MT 1980
#
import numpy as np
from scipy import optimize
from scipy import integrate
from scipy import interpolate
import matplotlib.pyplot as plt


def pa2head(pa_, pnref = 1.e5, g = 9.81):
    """ pascal to pressure head, converts a numpy array """
    h = np.zeros(len(pa_))
    for i, p in enumerate(pa_):
        h[i] = (p - pnref) * 100. / 1000. / g
    return h


def head2pa(h_, pnref = 1.e5, g = 9.81):
    """ pressure to pascal, converts a numpy array """
    pa = np.zeros(len(h_))
    for i, h in enumerate(h_):
        pa[i] = pnref + h / 100.*1000.*g
    return pa


class Parameters:
    """ class containing the van genuchten parameters """

    def __init__(self, p):
        self.theta_R = p[0]
        self.theta_S = p[1]
        self.alpha = p[2]  # [1/cm]
        self.n = p[3]
        self.m = 1. - 1. / self.n
        self.Ksat = p[4]

    def __iter__(self):  # for conversion to list with list(soil)
        return (i for i in [self.theta_R, self.theta_S, self.alpha, self.n, self.Ksat])


def plot_retention_curve(param, label_ = ""):
    """ plots the retention curve """
    y_ = np.logspace(1., 4., 100)
    x_ = water_content(-y_, param)
    plt.plot(x_, y_, label = label_)
    plt.xlabel("water content [1]")
    plt.ylabel("- matric potential [cm]")
    plt.yscale('log')

def plot_hydraulic_conductivity(param, label_ = ""):
    """ plots the matric flux potential"""
    y_ = np.logspace(1., 4., 100)
    x_ = hydraulic_conductivity(-y_, param)
    plt.plot(y_, x_, label = label_)
    plt.ylabel("hydraulic_conductivity [cm/day]")
    plt.xlabel("- matric potential [cm]")
    plt.xscale('log')
    
def plot_matric_flux_potential(param, label_ = ""):
    """ plots the matric flux potential"""
    y_ = np.logspace(1., 4., 100)
    x_ =[ matric_flux_potential(-y__, param) for y__ in y_]
    plt.plot(y_, x_, label = label_)
    plt.ylabel("matric flux potential [cm2/day]")
    plt.xlabel("- matric potential [cm]")
    plt.xscale('log')

def pressure_head(theta, sp):
    """ returns pressure head at a given volumetric water content according to the van genuchten model """
    try:
        if isinstance(theta, (type(list()), type(np.array([])))):
            assert (theta > sp.theta_R).all()
            assert (theta <= sp.theta_S).all()
        else:
            assert theta > sp.theta_R
            assert theta <= sp.theta_S
    except:
        print('theta <= sp.theta_R or theta > sp.theta_S', theta, sp.theta_R, sp.theta_S)
        raise Exception
    theta = np.minimum(theta, sp.theta_S)  # saturated water conent is the maximum
    return -pow(pow((sp.theta_S - sp.theta_R) / (theta - sp.theta_R), (1. / sp.m)) - 1., 1. / sp.n) / sp.alpha


def specific_moisture_storage(h, sp):  # soil water retention function (Cw(h))
    """ returns the specific moisture storage according to the van genuchten model """
#     h = np.minimum(h, np.zeros(h.shape))
#     h = np.maximum(h, np.ones(h.shape) * -30000.)
    C = -sp.alpha * sp.n * np.sign(h) * (1. / sp.n - 1.) * pow(sp.alpha * abs(h), sp.n - 1.) * (sp.theta_R - sp.theta_S) * pow(pow(sp.alpha * abs(h), sp.n) + 1., 1. / sp.n - 2.)
    return C


def water_diffusivity(TH, theta_i, theta_sur, sp):
    """ returns the water diffusivity (Eqn 11) """
    theta = TH * (theta_i - theta_sur) + theta_sur
    Se = (theta - sp.theta_R) / (sp.theta_S - sp.theta_R)
    m = sp.m
    D = (1 - m) * sp.Ksat / (sp.alpha * m * (sp.theta_S - sp.theta_R)) * pow(Se, 0.5 - 1. / m) * (pow(1 - pow(Se, 1. / m), -m) + pow(1 - pow(Se, 1 / m), m) - 2)
    return D


def water_content(h, sp):
    """ returns the volumetric water content [1] at a given matric potential [cm] according to the VanGenuchten model (Eqn 21) """
    try:
        if isinstance(h, (type(list()), type(np.array([])))):
            assert (h <= 0).all()
            assert (h > -np.inf).all()
        else:
            assert h <= 0
            assert h > -np.inf
    except:
        print('water_content, sp out of range', h)
        raise Exception
    return sp.theta_R + (sp.theta_S - sp.theta_R) / pow(1. + pow(sp.alpha * abs(h), sp.n), sp.m)


def effective_saturation(h, sp):
    """ returns the effective saturation [1] at a given matric potential [cm] according to the VanGenuchten model (dimensionless water content, Eqn 2) """
    # h = min(h, 0)  # pressure head is negative, zero the maximum
    theta = water_content(h, sp)
    se = (theta - sp.theta_R) / (sp.theta_S - sp.theta_R)
    return se


def hydraulic_conductivity(h, sp):
    """ returns the hydraulic conductivity [cm/day] at a given matric potential [cm] according to the van genuchten model (Eqn 8) """
    se = effective_saturation(h, sp)
    K = sp.Ksat * (se ** 0.5) * ((1. - pow(1. - pow(se, 1. / sp.m), sp.m)) ** 2)
    return K


def matric_flux_potential(h, sp):
    """ returns the matric flux potential [cm2/day] for a matric potential [cm]"""
    hmin = -16000
    K = lambda h: hydraulic_conductivity(h, sp)  # integrand
    MFP, _ = integrate.quad(K, hmin, h)
    return MFP


def matric_potential_mfp(mfp, sp):
    """ returns the matric potential [cm] from the matric flux potential [cm2/day]"""
    hmin = -16000  # needed to make hmin larger or else brent did not find a solution
    mfp_ = lambda psi: matric_flux_potential(psi, sp) - mfp
    h = optimize.brentq(mfp_, hmin, 0)
    return h


fast_mfp = {}
""" fast_mfp[sp](h):returns the matric flux potential [cm2/day] for a matric potential [cm], 
    call create_mfp_lookup first, once for each soil parameter @param sp"""

fast_imfp = {}
""" fast_imfp[sp](mfp): returns the matric potential [cm] from the matric flux potential [cm2/day], 
    call create_mfp_lookup first, once for each soil parameter @param sp"""


def create_mfp_lookup(sp, wilting_point = -16000, n = 16001):
    """ initializes the look up tables for soil parameter to use fast_mfp, and fast_imfp """
    print("initializing look up tables")
    global fast_mfp
    global fast_imfp

    h_ = -np.logspace(np.log10(1.), np.log10(np.abs(wilting_point)), n)
    h_ = h_ + np.ones((n,))

    mfp = np.zeros(h_.shape)
    for i, h in enumerate(h_):
        mfp[i] = matric_flux_potential(h, sp)
    fast_mfp[sp] = interpolate.interp1d(h_, mfp, bounds_error = False, fill_value = (mfp[0], mfp[-1]))  #
#     print("Table")
#     print(h_[0], h_[-1])
#     print(mfp[0], mfp[-1])

    imfp = np.zeros(h_.shape)
    for i, _ in enumerate(mfp):
        imfp[i] = h_[i]
    fast_imfp[sp] = interpolate.interp1d(mfp, imfp, bounds_error = False, fill_value = (imfp[0], imfp[-1]))  #

    print("done")

# fast_specific_moisture_storage = {}
# fast_water_content = {}
# fast_hydraulic_conductivity = {}
#
#
# def create_lookups(sp, wilting_point=-15000, n=15001):
#     """ good luck with that..."""
#     global fast_specific_moisture_storage
#     global fast_water_content
#     global fast_hydraulic_conductivity
#
#     h_ = -np.logspace(np.log10(1.), np.log10(np.abs(wilting_point)), n)
#     h_ = h_ + np.ones((n,))
#     print("Creating Van Genuchten Look Up Tables")
#     print("Table")
#     print(h_[0], h_[-1])
#
#     sms = np.zeros(h_.shape)
#     theta = np.zeros(h_.shape)
#     k = np.zeros(h_.shape)
#     for i, h in enumerate(h_):
#         sms[i] = specific_moisture_storage(h, sp)
#         theta[i] = water_content(h, sp)
#         k[i] = water_content(h, sp)
#     fast_specific_moisture_storage[sp] = interpolate.interp1d(h_, sms, bounds_error=False, fill_value=(sms[0], sms[-1]))
#     fast_water_content[sp] = interpolate.interp1d(h_, theta, bounds_error=False, fill_value=(theta[0], theta[-1]))  #
#     fast_hydraulic_conductivity[sp] = interpolate.interp1d(h_, k, bounds_error=False, fill_value=(k[0], k[-1]))  #
#     input()

