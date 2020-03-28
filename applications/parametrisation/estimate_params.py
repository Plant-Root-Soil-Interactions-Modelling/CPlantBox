import sys
sys.path.append("../..")
import plantbox as pb

import math
import numpy as np


def find_base(parent_poly, i):   
    """
    returns the index of the base root the root is connected to
    @param i the index of the polyline 
    @return the index of the base root, this lateral or base root (==i) is connected to 
    """
    if parent_poly[i] == -1:
        return i  
    else: 
        return find_base(parent_poly, parent_poly[i])


def measurement_time(polylines, properties, functions, time):
    """
    truncates the rsml data to a certain maximal time
    todo functions
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
#             for k in functions.keys():               
#                 f = functions[k]                 
#                 fs[k].setdefault(k, []).append(f) 
    return pl, ps, fs


def estimate():
    pass
    
