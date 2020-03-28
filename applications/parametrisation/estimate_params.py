import sys
sys.path.append("../..")
import plantbox as pb

import math
import numpy as np


def find_base_index(poly_indices, i):   
    """
    returns the index of the base root the root is connected to
    @param i the index of the polyline 
    @return the index of the base root, this lateral or base root (==i) is connected to 
    """
    if poly_indices[i] == -1:
        return i  
    else: 
        return find_base_index(poly_indices, poly_indices[i])


def find_order(poly_indices, i, o=0):   
    """
    returns the index of the base root the root is connected to
    @param i the index of the polyline 
    @return the index of the base root, this lateral or base root (==i) is connected to 
    """
    if poly_indices[i] == -1:
        return o  
    else: 
        return find_base_index(parent_poly, parent_poly[i], o + 1)


def measurement_time(polylines, properties, functions, time):
    """
    truncates the rsml data to a certain maximal time
    todo functions
    """
    pl, ps = [], {},
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
    return pl, ps


def base_roots(polylines, properties, functions):
    """
    return only the base roots
    """


def create_length(polylines, properties):
    """
    recalculate root lengths, and adds it to properties
    """
    properties["length"] = 0.
    for i, p in enumerate(polylines):
        l = 0.
        for j, n1 in enumerate(p[:-1]):
            n2 = p[j + 1]
            l += np.linalg.norm(n2 - n1)
        p["length"].append(l)


def create_order(polylines, properties, time):
    """
    recalculate root orders, and adds it to properties
    """
    properties["order"] = 0
    for i, p in enumerate(polylines):
        l = 0.
        for j, n1 in enumerate(p[:-1]):
            n2 = p[j + 1]
            l += np.linalg.norm(n2 - n1)
        p["length"].append(l)


def estimate():
    pass
    
