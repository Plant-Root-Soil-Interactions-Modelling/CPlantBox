#PSP_ThomasAlgorithm.py 
from __future__ import division
import numpy as np

def Thomas(a, b, c, d):
    n=len(d)
    x =np.zeros(n, float)
    
    for i in range(0, n-1):
        c[i] /=  b[i]
        d[i] /=  b[i]
        b[i + 1] -=  a[i + 1] * c[i]
        d[i + 1] -=  a[i + 1] * d[i]
    
    # back substitution    
    x[n-1] = d[n-1] / b[n-1]
    for i in range(n-2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    
    return(x)

def ThomasBoundaryCondition(a, b, c, d, x, first, last):
    
    for i in range(first, last):
        c[i] /=  b[i]
        d[i] /=  b[i]
        b[i + 1] -=  a[i + 1] * c[i]
        d[i + 1] -=  a[i + 1] * d[i]
    
    # back substitution    
    x[last] = d[last] / b[last]
    for i in range(last-1, first-1, -1):
        x[i] = d[i] - c[i] * x[i + 1]