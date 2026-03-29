
# this is a small test file by Erik Kopp to test the alternative implementation of the perirhizal resistances and the perirhizal diffusion


import sys; sys.path.append("../src/functional")
#from plantbox import Perirhizal
from plantbox.functional.Perirhizal import PerirhizalPython
import plantbox.functional.van_genuchten as vg
import Perirhizal
import pandas as pd
import numpy as np
import time

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
labels = ["rx","sx","inner_kr","rho","hsr","hsr_base","hsr_simp","hsr_global"]

#kr is assumed to lie between .5 and 5e-6 cm/hPa/d = 0.5 to 5e-6 1/d 
#r is assumed to be 0.03
#multiply this number by 10 in order to not fall under the threshhold at which the soil hydraulic potential is chosen as a default

Intervals = {
        "rx": [-14000,-10],
        "sx": [-1000,-10],
        "inner_kr": [1.5e-7,1.5e-6], 
        "rho": [10,199.0],
        }
Intervals = pd.DataFrame(Intervals, index = ["min","max"])

ntests = 30 
tests = pd.DataFrame(index = range(ntests),columns = labels)

for one_label in Intervals.columns:
    min_val = Intervals.loc["min",one_label]
    max_val = Intervals.loc["max",one_label]

    testvalues = (max_val - min_val) * np.random.rand(ntests) + min_val * np.ones(ntests)
    tests[one_label] = testvalues

peri = PerirhizalPython()
#peri = PerirhizalPython(Perirhizal)


rx = np.array(tests.loc[:,"rx"]) #* (-1)
sx = np.array(tests.loc[:,"sx"]) #* (-1)
inner_kr = np.array(tests.loc[:,"inner_kr"])
rho = np.array(tests.loc[:,"rho"])

hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]
sp = vg.Parameters(hydrus_loam)
peri.set_soil(sp)

#create the big lookup table
#t = time.time()
#peri.create_lookup_mpi("hydrus_loam", sp)
#elapsed_time = time.time() - t
#print("Creating the big lookup table took about ", elapsed_time) #

peri.open_lookup("hydrus_loam")

#create the small lookup table
#t = time.time()
#peri.create_lookup_simp("hydrus_loam", sp)
#elapsed_time = time.time() - t

peri.open_simp_lookup("hydrus_loam")
#print("Creating the simp lookup table took about ", elapsed_time) #

#create the global lookup table
#t = time.time()
#peri.create_lookup_global(peri.water_filename, sp)
#elapsed_time = time.time() - t

peri.open_global_lookup(peri.water_filename)
#print("Creating the global lookup table took about ", elapsed_time) #

hsr, hsr_base, hsr_simp, hsr_global = peri.soil_root_interface_potentials(rx, sx, inner_kr, rho)

tests.loc[:,"hsr"] = hsr

#the standard implementation
tests.loc[:,"hsr_base"] = hsr_base

#the first simplified implementation
tests.loc[:,"hsr_simp"] = hsr_simp

#the alternative implementation
tests.loc[:,"hsr_global"] = hsr_global



print(tests)

