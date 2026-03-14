
# this is a small test file by Erik Kopp to test the alternative implementation of the perirhizal resistances and the perirhizal diffusion


import sys; sys.path.append("../src/functional")
#from plantbox import Perirhizal
from plantbox.functional.Perirhizal import PerirhizalPython
import plantbox.functional.van_genuchten as vg
import Perirhizal
import pandas as pd
import numpy as np

# generate random inputs for one perirhizal segment

"""
        rx             xylem matric potential [cm]
        sx             bulk soil matric potential [cm]
        inner_kr       root radius times hydraulic conductivity [cm/day] 
        rho            geometry factor [1] (outer_radius / inner_radius)
        sp             soil parameter: van Genuchten parameter set (type vg.Parameters)
        sp_theta_r
        sp_theta_s
        sp_n
        sp_alpha
        sp_Ksat
"""
labels = ["rx","sx","inner_kr","rho","sp_theta_r","sp_theta_s","sp_n","sp_alpha","sp_Ksat","hsr_base","hsr_simp","hsr_alt"]

#kr is assumed to lie between .5 and 5e-6 cm/hPa/d = 0.5 to 5e-6 1/d 
#r is assumed to be 0.03
#multiply this number by 10 in order to not fall under the threshhold at which the soil hydraulic potential is chosen as a default

Intervals = {
        "rx": [-14000,-10],
        "sx": [-1000,-10],
        "inner_kr": [1.5e-7,1.5e-6], 
        "rho": [(0.5/0.03),(10.0)/0.03],
        "sp_theta_r": [0.01,0.025],
        "sp_theta_s": [0.4,0.43],
        "sp_n": [1.25,1.38],
        "sp_alpha": [0.0083,0.0383],
        "sp_Ksat": [2,60]
        }
Intervals = pd.DataFrame(Intervals, index = ["min","max"])

ntests = 2 
tests = pd.DataFrame(index = range(ntests),columns = labels)

for one_label in Intervals.columns:
    min_val = Intervals.loc["min",one_label]
    max_val = Intervals.loc["max",one_label]

    testvalues = (max_val - min_val) * np.random.rand(ntests) + min_val * np.ones(ntests)
    tests[one_label] = testvalues

peri = PerirhizalPython()
#peri = PerirhizalPython(Perirhizal)

for ind in range(ntests):
    rx = tests.loc[ind,"rx"] #* (-1)
    sx = tests.loc[ind,"sx"] #* (-1)
    inner_kr = tests.loc[ind,"inner_kr"]
    #print(inner_kr)
    rho = tests.loc[ind,"rho"]
    #rho = 10
    sp_theta_r = tests.loc[ind,"sp_theta_r"]
    sp_theta_s = tests.loc[ind,"sp_theta_s"]
    sp_n = tests.loc[ind,"sp_n"]
    sp_alpha = tests.loc[ind,"sp_alpha"]
    sp_Ksat = tests.loc[ind,"sp_Ksat"]
    
    hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]
    #sp = vg.Parameters([sp_theta_r,sp_theta_s,sp_n,sp_alpha,sp_Ksat])
    sp = vg.Parameters(hydrus_loam)

    #peri.sp = sp
    peri.set_soil(sp)

    #the standard implementation
    hsr = peri.soil_root_interface_(rx, sx, inner_kr, rho, sp)
    tests.loc[ind,"hsr_base"] = hsr

    #the first simplified implementation
    hsr_simp = peri.soil_root_interface_simp(peri, rx, sx, inner_kr, rho, sp)
    tests.loc[ind,"hsr_simp"] = hsr_simp

    #the alternative implementation
    hsr_alt = peri.soil_root_interface_alt(peri, rx, sx, inner_kr, rho, sp)
    #print(hsr_alt)
    tests.loc[ind,"hsr_alt"] = hsr_alt


print(tests)

