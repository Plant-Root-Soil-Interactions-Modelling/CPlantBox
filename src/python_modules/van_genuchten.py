#
# Mualem - van Genuchten model, equations from van Genuchten, MT 1980
#
import numpy as np
from scipy import optimize
from scipy import integrate
from scipy import interpolate
import matplotlib.pyplot as plt

def pa2head(pa_, pnref=1.e5, g=9.81):
    """ pascal to pressure head, converts a numpy array """
    h = np.zeros(len(pa_))
    for i, p in enumerate(pa_):
        h[i] = (p - pnref) * 100. / 1000. / g
    return h


def head2pa(h_, pnref=1.e5, g=9.81):
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


def pressure_head(theta, sp):
    """ returns pressure head at a given volumetric water content according to the van genuchten model """
    theta = min(theta, sp.theta_S)  # saturated water conent is the maximum
    return -pow(pow((sp.theta_S - sp.theta_R) / (theta - sp.theta_R), (1. / sp.m)) - 1., 1. / sp.n) / sp.alpha


def specific_moisture_storage(h, sp):
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
    K = lambda h: hydraulic_conductivity(h, sp)  # integrand 
    MFP, _ = integrate.quad(K, -15000, h)
    return MFP


def matric_potential_mfp(mfp, sp):
    """ returns the matric potential [cm] from the matric flux potential [cm2/day]"""
    mfp_ = lambda psi: matric_flux_potential(psi, sp) - mfp
    h = optimize.brentq(mfp_, -15000, 0)
    return h

fast_mfp = {}
""" fast_mfp[sp](h):returns the matric flux potential [cm2/day] for a matric potential [cm], 
    call create_mfp_lookup first, once for each soil parameter @param sp"""

fast_imfp = {}
""" fast_imfp[sp](mfp): returns the matric potential [cm] from the matric flux potential [cm2/day], 
    call create_mfp_lookup first, once for each soil parameter @param sp"""
    
def create_mfp_lookup(sp, wilting_point = -15000, n = 15001):
    """ initializes the look up tables for soil parameter to use fast_mfp, and fast_imfp """
    print("initializing look up tables")
    global fast_mfp 
    global fast_imfp

    h_ = -np.logspace(np.log10(1.), np.log10(np.abs(wilting_point)), n) 
    h_ = h_ + np.ones((n,))
    
    mfp = np.zeros(h_.shape)
    for i, h in enumerate(h_):
        mfp[i] = matric_flux_potential(h,sp)
    fast_mfp[sp] = interpolate.interp1d(h_, mfp, bounds_error=False, fill_value = (mfp[0], mfp[-1])) # 
#     print("Table")
#     print(h_[0], h_[-1])
#     print(mfp[0], mfp[-1])
    
    imfp = np.zeros(h_.shape)
    for i, _ in enumerate(mfp):
        imfp[i] = h_[i]   
    fast_imfp[sp] = interpolate.interp1d(mfp, imfp, bounds_error=False, fill_value = (imfp[0], imfp[-1])) # 
    
    print("done")
