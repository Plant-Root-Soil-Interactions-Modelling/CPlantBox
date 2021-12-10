import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")
import numpy as np
from scipy.optimize import minimize


def negexp_length(t, r, k):
    return k * (1 - np.exp(-(r / k) * t))


def negexp_age(l, r, k):
    return -k / r * np.log(1 - l / k)


def negexp_rate(l, k, t):
    return -k / t * np.log(1 - l / k)


def target_length(r:float, k:float, lengths:np.array, ages:np.array):
    """ target function for optimization root target length [cm],
    @param r initial root length [cm], @param k maximal root length [cm]
    @param lengths root lengths [cm], @param ages root ages [day] corresponding to root lengths"""
    assert lengths.shape == ages.shape, "target_length: number of root lengths must equal number of root ages"
    sum = 0.
    for i in range(0, lengths.shape[0]):
        l = negexp_length(ages[i], r, k)
        sum += (l - lengths[i]) ** 2
    return np.sqrt(sum / (lengths.shape[0] - 2))  # -(k+1)


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
    f = lambda x0: target_length(x0[0], x0[1], np.array(length), np.array(times))
    x0 = [5., 200]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    return res.x[0], res.x[1], f(res.x)


def estiamte_emergance_order0(lengths:np.array, ages:np.array, r:float, k:float):
    """ fits the emergance time of 2nd basal roots """
    f = lambda x: target_length(r, k, lengths, ages - np.ones(ages.shape) * x[0])
    x0 = [np.max(ages)]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)
    return res, f


def estimate_order0_rate(lengths:np.array, r:float, k:float, time:float):
    """ fits basal prodcution rate [day-1] for given initial growth rate and maximal root length
    @param lengths list of root lengths [cm] 
    @param r initial root length [cm] 
    @param k maximal root length [cm]
    @param time maximal measurement time """
    f = lambda x: target_rate(x[0], lengths, r, k, time)
    x0 = [time / lengths.shape[0]]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    n = lengths.shape[0]
    ages = np.zeros(lengths.shape)
    for i in range(0, n):
        ages[n - i - 1] = max(time - i * res.x[0], 0.)
    return res, f, ages


def polyline_length(i0, i1, polyline):
    """ length between two indices i0 and i1 """
    l = 0.
    for i in range(int(i0), int(i1)):
        l += np.linalg.norm(np.array(polyline[i]) - np.array(polyline[i + 1]))
    return l

