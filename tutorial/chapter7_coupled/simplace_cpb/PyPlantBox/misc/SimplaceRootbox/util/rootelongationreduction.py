'''
Created on 25.04.2018

@author: Gunther Krauss <guntherkrauss@uni-bonn.de>
'''

import math
import operator

def vangenuchtenContentToPressure (ta, ts, tr, n, alpha):
    m = 1. - 1./n
    if ta >= ts:
        return 0
    tta = ta
    if ta <= tr:
        tta = tr + 0.0001
    h = pow((ts-tr)/(tta - tr), 1./m)
    h = pow(h-1.,1./n)
    h = -h / alpha
    return h

def vangenuchtenPressureToContent (h, ts, tr, n, alpha):
    m = 1. - 1./n
    th = pow(abs(h)*alpha, n) + 1
    th = tr + ((ts-tr)/pow(th,m))
    return th

def ReductionREWater(h, h1, h2, h3, h4):
    (g, g1, g2, g3, g4) = (abs(h), abs(h1), abs(h2), abs(h3), abs(h4))
    if g <= g1:
        return 0
    elif g1 < g and g <= g2:
        return (g1-g)/(g1-g2)
    elif g2 < g and g <= g3:
        return 1
    elif g3 < g and g <= g4:
        return (g4-g) / (g4-g3)
    else:
        return 0
    
def ReductionRESoilStrengthContPores(q):
    return math.exp(-.3 * q)

def ReductionRESoilStrengthNoContPores(q):
    return math.exp(-0.0025 * q)  # set to -.0025 for sbarley and -0.005 for winter wheat,   # default:-0.4325

def ReductionRESoilStrength(q, cont_mp = False):
    if cont_mp:
        return ReductionRESoilStrengthContPores(q) 
    else: 
        return ReductionRESoilStrengthNoContPores(q)


def ReductionRE(q, 
                              ta, ts, tr, 
                              h1, h2, h3, h4,
                              n, alpha, 
                              fun = lambda a, b : a*b,
                              cont_mp = False): 
       
    return [fun(
        ReductionRESoilStrength(q[i], cont_mp),
        ReductionREWater(vangenuchtenContentToPressure(ta[i], ts[i], tr[i], n[i], alpha[i]), h1, h2, h3, h4))
      for i in range(0,len(q))]

   
def ReductionREMultiplicative(q, 
                              ta, ts, tr, 
                              h1, h2, h3, h4,
                              n, alpha, 
                              cont_mp = False): 
       
    return ReductionRE(q, ta, ts, tr, h1, h2, h3, h4, n, alpha, operator.mul, cont_mp)
    

def ReductionREMin(q, 
                              ta, ts, tr, 
                              h1, h2, h3, h4,
                              n, alpha, 
                              cont_mp = False): 
       
    return ReductionRE(q, ta, ts, tr, h1, h2, h3, h4, n, alpha, min, cont_mp)

def SoilStrengthFromBDandWC(ta, bd, a, b ,c):
    return [a*pow(bd[i],b)*pow(ta[i],c) for i in range(0, len(bd))]
    