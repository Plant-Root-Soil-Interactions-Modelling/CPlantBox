import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")

import numpy as np
from scipy.optimize import minimize


def negexp_length(t, r, k):
    return k * (1 - np.exp(-(r / max(k, 1.e-9)) * t))


def negexp_age(l, r, k):
    try:
        return -k / r * np.log(1 - l / k)
    except:
        return 0.


def negexp_rate(l, k, t):
    try:
        return -k / max(t, 1.e-9) * np.log(max(1 - l / max(k, 1.e-9), 1.e-9))
    except:
        return 1.e-9  # if k = 0 growth is probably slow


def target_rate(rate, lengths:np.array, r:float, lmax:float, times:float):
    """ target function for estimating the linear base root production rate [day-1],
    @param rate: linear base root production rate (delay between basal root emergence) [day-1] 
    @param lengths: basal root lengths as numpy array sorted ascending [cm], list (per measurement) of list of sorted root lengths
    @param r: initial growth rate [cm/day]
    @param lmax: maximal root length [cm] 
    @param times: maximal measurement, time list (per measurement)
    @return error 
     """
    rate = max(rate, 0.)
    n = len(lengths)  # number of measurements
    ages, lengths2 = [], []  # flattened ages and lenghts
    for ii in range(0, n):
        nn = len(lengths[ii])
        for i in range(0, nn):
            ages.append(min(max(times[ii] - (i + 1) * rate, 0.), times[ii]))
            lengths2.append(lengths[ii][i])
    # print('target rate', rate, np.array(lengths2), np.array(ages))
    x = target_length(r, lmax, np.array(lengths2), np.array(ages))
    return x


def target_length(r:float, k:float, lengths:np.array, ages:np.array):
    """ target function for optimization root target length [cm],
    @param r initial root length [cm], @param k maximal root length [cm]
    @param lengths root lengths [cm], @param ages root ages [day] corresponding to root lengths"""
    assert lengths.shape == ages.shape, "target_length: number of root lengths must equal number of root ages"
    sum = 0.
    for i in range(0, lengths.shape[0]):
        l = negexp_length(ages[i], r, k)
        sum += (l - lengths[i]) ** 2
    return np.sqrt(sum / (np.max([lengths.shape[0] - 2, 1])))  # -(k+1)


def fit_taproot_r(length, times, k):
    """ fits initial growth rate r, assumes maximal root lenght k as fixed (e.g. literature value) """
    assert(len(length) == len(times))
    f = lambda x0: target_length(x0, k, np.array(length), np.array(times))
    x0 = [1.]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    return res.x[0], f(res.x[0])


def fit_taproot_rk(length, times):
    """ fits initial growth rate r, and maximal root lenght k """
    assert(len(length) == len(times))
    f = lambda x: target_length(x[0], x[1], np.array(length), np.array(times))
    x0 = [5., 200]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    return res.x[0], res.x[1], f(res.x)

# def estiamte_emergance_order0(lengths:np.array, ages:np.array, r:float, k:float):
#     """ fits the emergance time of 2nd basal roots """
#     f = lambda x: target_length(r, k, lengths, ages - np.ones(ages.shape) * x[0])
#     x0 = [np.max(ages)]
#     res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)
#     return res, f


def estimate_order0_rate(lengths:np.array, r:float, k:float, times:float):
    """ fits basal prodcution rate [day-1] for given initial growth rate and maximal root length
    @param lengths list of root lengths [cm] list (per measurement) of list of root sorted lengths
    @param r initial predefined growth rate [cm/day] 
    @param k maximal root length [cm] 
    @param times maximal measurement, time list (per measurement) [day]"""
    assert len(lengths) == len(times), "estimate_order0_rate: size of measuered lengths list must equal measuring times"
    f = lambda x: target_rate(x[0], lengths, r, k, times)
    n = len(lengths)
    x0 = np.max(times) / n
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    ages = [None] * n
    for ii in range(0, n):
        ages[ii] = []
        nn = len(lengths[ii])
        for i in range(0, nn):
            ages[ii].append(min(max(times[ii] - (i + 1) * res.x[0], 0.), times[ii]))
    return res, f, ages


def estimate_order0_rrate(lengths:np.array, r0:float, k:float, times:float):
    """ fits basal prodcution rate [day-1] for given initial growth rate and maximal root length
    @param lengths list of root lengths [cm] list (per measurement) of list of root sorted lengths
    @param r0 initial growth rate for fitting [cm/day] 
    @param k maximal root length [cm] 
    @param times maximal measurement, time list (per measurement) [day]"""
    assert len(lengths) == len(times), "estimate_order0_rate: size of measuered lengths list must equal measuring times"
    f = lambda x: target_rate(x[0], lengths, x[1], k, times)
    n = len(lengths)
    x0 = [ np.max(times) / n, r0]  # initial value
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    ages = [None] * n
    for ii in range(0, n):
        ages[ii] = []
        nn = len(lengths[ii])
        for i in range(0, nn):
            ages[ii].append(min(times[ii] - (i + 1) * res.x[0], times[ii]))
    return res, f, ages


def polyline_length(i0, i1, polyline):
    """ length between two indices i0 and i1 """
    l = 0.
    for i in range(int(i0), int(i1)):
        l += np.linalg.norm(np.array(polyline[i]) - np.array(polyline[i + 1]))
    return l

