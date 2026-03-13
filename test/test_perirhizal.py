
# this is a small test file by Erik Kopp to test the alternative implementation of the perirhizal resistances and the perirhizal diffusion


import sys; sys.path.append("../src/functional")
from Perirhizal import *
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

Intervals = {
        "rx": [-14000,-10],
        "sx": [-1000,-10],
        "inner_kr": [1.5e-8,1.5e-7], 
        "rho": [(0.5/0.03),(10.0)/0.03],
        "sp_theta_r": [0.01,0.10],
        "sp_theta_s": [0.2,0.5],
        "sp_n": [1.0,1.5],
        "sp_alpha": [0.005,0.05],
        "sp_Ksat": [1,50]
        }
Intervals = pd.DataFrame(Intervals, index = ["min","max"])

ntests = 100 
tests = pd.DataFrame(index = range(ntests),columns = labels)

for one_label in Intervals.columns:
    min_val = Intervals.loc["min",one_label]
    max_val = Intervals.loc["max",one_label]

    testvalues = (max_val - min_val) * np.random.rand(ntests) + min_val * np.ones(ntests)
    tests[one_label] = testvalues

for ind in range(ntests):
    rx = tests.loc[ind,"rx"]
    sx = tests.loc[ind,"sx"]
    inner_kr = tests.loc[ind,"inner_kr"]
    rho = tests.loc[ind,"rho"]
    sp_theta_r = tests.loc[ind,"sp_theta_r"]
    sp_theta_s = tests.loc[ind,"sp_theta_s"]
    sp_n = tests.loc[ind,"sp_n"]
    sp_alpha = tests.loc[ind,"sp_alpha"]
    sp_Ksat = tests.loc[ind,"sp_Ksat"]

    #the standard implementation
    hsr = soil_root_interface_(rx, sx, inner_kr, rho, sp)
    tests.loc[ind,"hsr_base"] = hsr

    #the first simplified implementation
    hsr_simp = soil_root_interface_simp(rx, sx, inner_kr, rho, sp)
    tests.loc[ind,"hsr_simp"] = hsr_simp

    #the alternative implementation
    hsr_alt = soil_root_interface_alt(rx, sx, inner_kr, rho, sp)
    tests.loc[ind,"hsr_alt"] = hsr_alt




