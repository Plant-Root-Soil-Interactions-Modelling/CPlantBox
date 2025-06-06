import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb

import math
import numpy as np
from scipy.optimize import minimize


def find_base_index(poly_indices, i):
    """
    The index of the base, root this root with inde i is connected to
    
    @param poly_indices      the indices of all parent roots, usually properties["parent-poly"] 
    @param i                 the index of this polyline
     
    @return index of the base root this root i is part of 
    """
    if poly_indices[i] == -1:
        return i
    else:
        return find_base_index(poly_indices, poly_indices[i])


def find_order(poly_indices, i, o = 1):
    """
    The root order of the root with index i starting with 1
    
    @param poly_indices      the indices of all parent roots, usually properties["parent-poly"] 
    @param i                 the index of this polyline
    
    @return root order of root with index i 
    """
    if poly_indices[i] == -1:
        return o
    else:
        return find_order(poly_indices, poly_indices[i], o + 1)


def measurement_time(polylines, properties, functions, time):
    """
    Truncates the rsml data to a certain maximal time
    """
    pl, ps, fs = [], {}, {}
    fet = functions["emergence_time"]
    for i, p in enumerate(polylines):
        if fet[i][0] <= time:  # add
            p_ = []
            for j, node in enumerate(p):
                if fet[i][j] <= time:
                    p_.append(node)
            pl.append(p_)
            for k in properties.keys():
                ps.setdefault(k, []).append(properties[k][i])
            for k in functions.keys():
                fs.setdefault(k, []).append(functions[k][i])
    return pl, ps, fs


def base_roots(polylines, properties):
    """
    return only the base roots plus first order lateral roots
    """
    pp = properties["parent-poly"]
    pl, ps = [], {}
    for i, p in enumerate(polylines):
        if pp[i] == -1:  # add if base root
            # print(properties.keys())
            pl.append(p)
            for k in properties.keys():
                try:
                    dummy = properties[k][i]
                except IndexError:
                    dummy = None
                if dummy is None:
                    print('Found: Base roots that miss some properties found')
                    ps.setdefault(k, [])
                else:
                    ps.setdefault(k, []).append(properties[k][i])
    return pl, ps


def lateral_roots(polylines, properties):
    """
    return only the first order lateral roots
    """
    pp = properties["order"]
    pl, ps = [], {}
    for i, p in enumerate(polylines):
        if pp[i] == 2:  # add if base root
            # print(properties.keys())
            pl.append(p)
            for k in properties.keys():
                try:
                    dummy = properties[k][i]
                except IndexError:
                    dummy = None
                if dummy is None:
                    print('Found: Base roots that miss some properties found')
                    ps.setdefault(k, [])
                else:
                    ps.setdefault(k, []).append(properties[k][i])
    return pl, ps


def laterals1(polylines, properties):
    """
    return only the base roots plus first order lateral roots
    """
    pp = properties["order"]
    pl, ps = [], {}
    for i, p in enumerate(polylines):
        if pp[i] == 2 or pp[i] == 1:  # add if first order lateral or base root
            # print(properties.keys())
            pl.append(p)
            for k in properties.keys():
                try:
                    dummy = properties[k][i]
                except IndexError:
                    dummy = None
                if dummy is None:
                    print('Found: First order laterals that miss some properties')
                    ps.setdefault(k, [])
                else:
                    ps.setdefault(k, []).append(properties[k][i])
    return pl, ps


def laterals(polylines, properties, base_polyline):
    """
    return only laterals that are attached to a given base root
    """

    pp = properties["order"]
    pl, ps = [], {}
    for i, p in enumerate(polylines):
        if pp[i] == 2:  # add if first order lateral
            # print(properties.keys())
            pl.append(p)
            for k in properties.keys():
                try:
                    dummy = properties[k][i]
                except IndexError:
                    dummy = None
                if dummy is None:
                    print('Found: First order laterals that miss some properties')
                    ps.setdefault(k, [])
                else:
                    ps.setdefault(k, []).append(properties[k][i])
    return pl, ps


def create_length(polylines, properties):
    """
    recalculate root lengths, and adds it to properties
    """
    properties["length"] = []
    for i, p in enumerate(polylines):
        l = 0.
        for j, n1 in enumerate(p[:-1]):
            n2 = p[j + 1]
            l += np.linalg.norm(np.array(n2) - np.array(n1))
        properties["length"].append(l)


def create_order(polylines, properties):
    """
    recalculate root orders, and adds it to properties
    """
    properties["order"] = []
    poly_indices = properties["parent-poly"]  # generated by read_rsml
    for i, p in enumerate(polylines):
        o = find_order(poly_indices, i)
        properties["order"].append(o)


def negexp_length(t, r, k):
    return k * (1 - np.exp(-(r / k) * t))


def negexp_age(l, r, k):
    return -k / r * np.log(1 - l / k)


def target(r, k, length, times):
    """ target function for optimization for fit_taproot_r, fit_taproot_rk """
    sum = 0.
    for i, t in enumerate(times):
        sim_l = negexp_length(t, r, k)
        sum_t = 0.
        if isinstance(length[i], float):
            sum_t += (sim_l - length[i]) ** 2
        else:
            for l in length[i]:
                sum_t += (sim_l - l) ** 2

        sum += np.sqrt(sum_t)
    return sum


def fit_taproot_r(length, times, k):
    """ fits initial growth rate r, assumes maximal root lenght k as fixed (e.g. literature value) """
    assert(len(length) == len(times))
    f = lambda x0: target(x0, k, length, times)
    x0 = [1.]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    return res.x[0], f(res.x[0])


def fit_taproot_rk(length, times):
    """ fits initial growth rate r, and maximal root lenght k """
    assert(len(length) == len(times))
    f = lambda x0: target(x0[0], x0[1], length, times)
    x0 = [5., 200]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    return res.x[0], res.x[1], f(res.x)


def target2(delay, times, numbers, i_n):
    """ target function for optimization for fit_number_of_roots"""
    sum = 0.
    for i, t in enumerate(times):
        sim_n = i_n + t / delay  # missing round (continious seems easier to optimize)
        s = 0.
        for j in range(0, len(numbers[i])):
            s += (numbers[i][j] - sim_n) ** 2
        sum += np.sqrt(s)
    return sum


def fit_number_of_roots(times, numbers, initial_number):
    """ we fit the delay between emergence of seminals with, 
        number_of_roots = inital_number + round(t/delay)
    """
    assert(len(numbers) == len(times))
    f = lambda x0: target2(x0, times, numbers, initial_number)
    x0 = [1.]  # days
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # linear regression would be enough in this case
    return res.x[0]

# def fit_seminal_rk(length, times):
#     """ fits initial growth rate r, and maximal root lenght k """
#     assert(len(length.shape[0] == len(times))
#     f = lambda x0: target(x0[0], x0[1], length, times)
#     x0 = [5., 200]
#     res = minimize(f, x0, method='Nelder-Mead', tol=1e-6)  # bounds and constraints are possible, but method dependent
#     # print(x)
#     return res.x[0], res.x[1]


def connect(node, base_polyline):
    """
    connects the node to the closest point in base_polyline
    """
    min_dist = 1.e8  # much
    for i, n in enumerate(base_polyline):
        dist = np.linalg.norm(n - node)
        if dist < min_dist:
            min_dist = dist
            min_i = i
    return min_i, min_dist


def reconstruct_laterals(polylines, properties, base_polyline, snap_radius = 0.5):
    """
    connect all lateral roots to the base root (all roots other than tap need will have order 1)
    """
    s = 0
    npl, nprop = [], {}
    pp = properties["parent-poly"]  # generated by read_rsml
    pn = properties["parent-node"]  # generated by read_rsml
    for i, pi in enumerate(pp):
        if pi == -1:  # tap root
            npl.append(polylines[i])
            for k in properties.keys():
                nprop.setdefault(k, []).append(properties[k][i])
        elif pi == 0:  # lateral root
            npl.append(polylines[i])
            for k in properties.keys():
                try:
                    dummy = properties[k][i]
                except IndexError:
                    dummy = None
                if dummy is None:
                    print('Found in reconstruct_laterals: First order laterals miss some properties')
                else:
                    nprop.setdefault(k, []).append(properties[k][i])
        elif pi > 0:  # attach to base root
            n = np.array(polylines[i][0])
            nni, dist = connect(n, base_polyline)
            if dist < snap_radius:  # TODO criteria
                npl.append(polylines[i])
                for k in properties.keys():
                    try:
                        dummy = properties[k][i]
                    except IndexError:
                        dummy = None
                    if dummy is None:
                        print('Found in reconstruct_laterals: First order laterals that miss some properties found')
                    else:
                        nprop.setdefault(k, []).append(properties[k][i])
                        nprop["parent-poly"][-1] = 0
                        nprop["parent-node"][-1] = nni
            else:
                s += 1

    print("removed", s, "roots")
    return npl, nprop


def polyline_length(i0, i1, polyline):
    """ length between two indices i0 and i1 """
    l = 0.
    for i in range(int(i0), int(i1)):
        l += np.linalg.norm(polyline[i] - polyline[i + 1])
    return l


def analyze_zones(polylines, properties, base_polyline):
    """
    determines la, lb, ln on given basal root
    """
    pni = properties["parent-node"]
    ii = pni[1:]  # lateral parent node indices
    ii.sort()
#     print(ii)
    n0 = 0
    nl0 = ii[0]
    lb = polyline_length(n0, nl0, base_polyline)  # length of basal zone
    ne = ii[-1]
    nle = len(base_polyline) - 1
    la = polyline_length(ne, nle, base_polyline)  # length of apical zone
    ln_ = []
    for i in range(0, len(ii) - 1):  # laterals
        i0 = ii[i]
        i1 = ii[i + 1]
        ln_.append(polyline_length(i0, i1, base_polyline))
    return lb, ln_, la


def analyze_theta(polylines, properties, base_polyline):
    """
    assumes polylines[0] is tap root and determines theta
    """
#     print("analyze_theta")
#     angle = np.arccos(np.clip(np.dot(np.array([0, -1]), np.array([0, -1])), -1.0, 1.0))
#     print("down", angle)
#     angle = np.arccos(np.clip(np.dot(np.array([0, -1]), np.array([1, 0])), -1.0, 1.0))
#     print("90", angle)
#     angle = np.arccos(np.clip(np.dot(np.array([0, -1]), np.array([-1, 0])), -1.0, 1.0))
#     print("90", angle)
    theta = []
    pni = properties["parent-node"]
    for i in range(0, len(polylines)):
        j = int(pni[i])
        if j > 0:
            p = base_polyline[j];  # polylines[0][j]
            p0 = base_polyline[j - 1];  # polylines[0][j - 1]
            if len(polylines[i]) > 2:
                p2 = polylines[i][2]
            else:
                p2 = polylines[i][0]
            v1 = p - p0
            v2 = p2 - p
            v1 = v1 / np.linalg.norm(v1)
            # v1 = [0, 0, -1] # also makes sense
            v2 = v2 / np.linalg.norm(v2)
            angle = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
            theta.append(angle)
    return theta


def create_diameters(polylines, properties, functions):
    """ 
    adds diameter property as mean value of the functions
    """
    properties["diameter"] = []
    for i, p in enumerate(polylines):
        d = 0.
        for j in range(0, len(p)):
#             print(len(p))
#             print(len(functions["diameter"]))
            d += functions["diameter"][i][j] / len(p) / 10.  # mm -> cm
        properties["diameter"].append(d)


def create_age_delay(polylines, properties, time, r0, lmax0, lad0):
    """
    adds estimated lateral root age to the properties, based on given tap root paraeters 
    """
    pni = properties["parent-node"]
    properties.setdefault("age", []).append(time)  # tap root et =  0
    for i in range(1, len(polylines)):
        il = pni[i]
        bl = polyline_length(0, il, polylines[0])  # latearal base length
        et = negexp_age(bl, r0, lmax0) + lad0  # time the lateral emerged
        age = max(time - et, 0.)
        age = min(age, time)
        properties["age"].append(age)  # tap root et =  0

