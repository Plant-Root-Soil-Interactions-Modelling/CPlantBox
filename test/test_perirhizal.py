
# this is a small test file by Erik Kopp to test the alternative implementation of the perirhizal resistances and the perirhizal diffusion


import sys; sys.path.append("../src/functional")
#from plantbox import Perirhizal
from plantbox.functional.Perirhizal import PerirhizalPython
from numpy import linalg as LA
import plantbox.functional.van_genuchten as vg
from rosi.richards_flat import RichardsFlatWrapper as RichardsWrapper  # Python part, macroscopic soil model
#from richards import RichardsWrapper  # Python part
#from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial
from rosi.rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding), macroscopic soil model
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
        c_sol          concentration of solutes in the bulk soil
        Vmax           Michaelis Menten Kinetics maximal solute uptake rate
        Km             Michaelis Menten Kinetics half saturation
        sp             soil parameter: van Genuchten parameter set (type vg.Parameters)
        sp_theta_r
        sp_theta_s
        sp_n
        sp_alpha
        sp_Ksat
"""
labels = ["rx","sx","inner_kr","rho","c_sol","Vmax","Km","hsr","hsr_lookup","hsr_global","waterflow","sol_c","sol_c_ss","sol_c_sr","sol_U","sol_U_ss","sol_U_sr"]

inputs = ["rx","sx","inner_kr","rho","c_sol","Vmax","Km"]
outputs = ["hsr","hsr_lookup","hsr_global","sol_c","sol_c_ss","sol_c_sr","sol_U","sol_U_ss","sol_U_sr"]
outputs = ["rx","sx","waterflow","sol_c","sol_c_ss","sol_c_sr","sol_U","sol_U_ss","sol_U_sr"]

#kr is assumed to lie between .5 and 5e-6 cm/hPa/d = 0.5 to 5e-6 1/d 
#r is assumed to be 0.03
#multiply this number by 10 in order to not fall under the threshhold at which the soil hydraulic potential is chosen as a default

Intervals = {
        "rx": [-14000,-10],
        "sx": [-1000,-10],
        "inner_kr": [1.5e-7,1.5e-6], 
        "rho": [10,199.0],
        "c_sol": [1.0e-9,1.0e-5],
        "Vmax": [1.0e-11,1.0e-10],
        "Km": [1.0e-8,1.0e-6],
        }
Intervals = pd.DataFrame(Intervals, index = ["min","max"])

ntests = 10 
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
c_sol = np.array(tests.loc[:,"c_sol"])
Vmax = np.array(tests.loc[:,"Vmax"])
Km = np.array(tests.loc[:,"Km"])
Ds = 1.0e1
r_root = 0.03

hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]
sp = vg.Parameters(hydrus_loam)
peri.set_soil(sp)

#create the lookup table
#t = time.time()
#peri.create_lookup("results/hydrus_loam", sp)
#elapsed_time = time.time() - t
#print("Creating the lookup table took about ", elapsed_time) #

peri.open_lookup("results/hydrus_loam")


#create the global lookup table
#t = time.time()
#peri.create_lookup_global("results/"+peri.water_filename, sp)
#elapsed_time = time.time() - t
#print("Creating the global lookup table took about ", elapsed_time) #

peri.open_global_lookup("results/"+peri.water_filename)

# create lookup tables for the solute flow
Ds_=np.logspace(np.log10(1.0e0), np.log10(1.0e1), 2)
#peri.create_integralDiffusion_lookup("results/hydrus_loam_ss_solutes", sp)
#peri.create_integralconcentration_lookup("results/hydrus_loam_sr_solutes", Ds_, sp)

# open lookup tables for the solute flow
peri.open_lookup_solutes("results/hydrus_loam_ss_solutes")
peri.open_lookup_sr_solutes("results/hydrus_loam_sr_solutes")


#numerically solve each root segment, this is the reference solution 
hsr = np.array([peri.soil_root_interface_(rx[i], sx[i], inner_kr[i], rho[i], peri.sp) for i in range(0, len(rx))])
tests.loc[:,"hsr"] = hsr
print("norm of the root soil matrix potentials: ", LA.norm(hsr))

#the standard implementation
hsr1= peri.soil_root_interface_potentials_table(rx, sx, inner_kr, rho)
tests.loc[:,"hsr_lookup"] = hsr1
print("Norm of the difference to basic lookup table:", LA.norm(hsr-hsr1))

#the alternative global implementation
hsr3 = peri.soil_root_interface_potentials_table_global(rx, sx, inner_kr, rho)
tests.loc[:,"hsr_global"] = hsr3
print("Norm of the difference to global lookup table:", LA.norm(hsr-hsr3))

#waterflow = 2*3.14*np.multiply((hsr-rx),inner_kr)/r_root
waterflow = 2*3.14*np.multiply((sx-rx),inner_kr)#/r_root
Phi_root = np.array([vg.fast_mfp[sp](hsr[i]) for i in range(len(hsr))])
Phi_soil = np.array([vg.fast_mfp[sp](sx[i]) for i in range(len(sx))])
tests.loc[:,"waterflow"] = waterflow
Ds = 1.0e1

#base solutes
tests.loc[:,"sol_c"] = c_sol
tests.loc[:,"sol_U"] = np.array([Vmax[i]*c_sol[i]/(Km[i]+c_sol[i])*1e6 for i in range(len(c_sol))])

# steady state solutes
solutes_ss = peri.soil_root_solutes_ss_(Phi_root, Phi_soil, c_sol, Vmax, Km, Ds, waterflow, peri.sp)
tests.loc[:,"sol_c_ss"] = solutes_ss
tests.loc[:,"sol_U_ss"] = np.array([Vmax[i]*solutes_ss[i]/(Km[i]+solutes_ss[i])*1e6 for i in range(len(solutes_ss))])

# steady rate solutes
solutes_sr = peri.soil_root_solutes_sr_(Phi_root, Phi_soil, rho, c_sol, Vmax, Km, Ds, waterflow, peri.sp)
tests.loc[:,"sol_c_sr"] = solutes_sr
tests.loc[:,"sol_U_sr"] = np.array([Vmax[i]*solutes_sr[i]/(Km[i]+solutes_sr[i])*1e6 for i in range(len(solutes_sr))])




print("All inputs:")
print(tests[inputs])
print("All outputs:")
print(tests[outputs])




