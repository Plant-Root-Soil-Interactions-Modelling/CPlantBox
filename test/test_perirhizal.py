
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
labels = ["rx","sx","inner_kr","rho","c_sol","Vmax","Km","hsr","hsr_base","hsr_simp","hsr_global","sol_c","sol_c_ss","sol_c_sr","sol_U","sol_U_ss","sol_U_sr"]
#for now only water flow:
#labels = ["rx","sx","inner_kr","rho","c_sol","Vmax","Km","hsr","hsr_1d","hsr_base","hsr_simp","hsr_global"]

inputs = ["rx","sx","inner_kr","rho","c_sol","Vmax","Km"]
outputs = ["hsr","hsr_base","hsr_simp","hsr_global","sol_c","sol_c_ss","sol_c_sr","sol_U","sol_U_ss","sol_U_sr"]

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
c_sol = np.array(tests.loc[:,"c_sol"])
Vmax = np.array(tests.loc[:,"Vmax"])
Km = np.array(tests.loc[:,"Km"])
Ds = 1.0e2

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
peri.create_lookup_simp("hydrus_loam", sp)
#elapsed_time = time.time() - t

peri.open_simp_lookup("hydrus_loam")
#print("Creating the simp lookup table took about ", elapsed_time) #

#create the global lookup table
#t = time.time()
#peri.create_lookup_global(peri.water_filename, sp)
#elapsed_time = time.time() - t

peri.open_global_lookup(peri.water_filename)
#print("Creating the global lookup table took about ", elapsed_time) #

#hsr = peri.soil_root_interface_potentials(rx, sx, inner_kr, rho)

hsr = np.array([peri.soil_root_interface_(rx[i], sx[i], inner_kr[i], rho[i], peri.sp) for i in range(0, len(rx))])
tests.loc[:,"hsr"] = hsr
print("norm of the root soil matrix potentials: ", LA.norm(hsr))

#the standard implementation
hsr1= peri.soil_root_interface_potentials_table(rx, sx, inner_kr, rho)
tests.loc[:,"hsr_base"] = hsr1
print("Norm of the difference to basic lookup table:", LA.norm(hsr-hsr1))

#the first simplified implementation
hsr2 = peri.soil_root_interface_potentials_table_simp(rx, sx, inner_kr, rho)
tests.loc[:,"hsr_simp"] = hsr2
print("Norm of the difference to simp lookup table:", LA.norm(hsr-hsr2))

#the alternative implementation
hsr3 = peri.soil_root_interface_potentials_table_global(rx, sx, inner_kr, rho)
tests.loc[:,"hsr_global"] = hsr3
print("Norm of the difference to global lookup table:", LA.norm(hsr-hsr3))

waterflow = 2*3.14*np.multiply((hsr-rx),inner_kr)
Phi_root = hsr
Phi_soil = np.array([vg.fast_mfp[sp](sx[i]) for i in range(len(sx))])

#base solutes
tests.loc[:,"sol_c"] = c_sol
tests.loc[:,"sol_U"] = np.array([Vmax[i]*c_sol[i]/(Km[i]+c_sol[i])*1e6 for i in range(len(c_sol))])

# steady state solutes
solutes_ss = peri.soil_root_solutes_ss_(Phi_root, Phi_soil, c_sol, Vmax, Km, Ds, waterflow, peri.sp)
tests.loc[:,"sol_c_ss"] = solutes_ss
tests.loc[:,"sol_U_ss"] = np.array([Vmax[i]*solutes_ss[i]/(Km[i]+solutes_ss[i])*1e6 for i in range(len(solutes_ss))])

# steady rate solutes
#solutes_sr = peri.soil_root_solutes_sr_(Phi_root, Phi_soil, rho, c_sol, Vmax, Km, Ds, waterflow, peri.sp)
#tests.loc[:,"sol_c_sr"] = solutes_sr
#tests.loc[:,"sol_U_sr"] = np.array([Vmax[i]*solutes_sr[i]/(Km[i]+solutes_sr[i])*1e6 for i in range(len(solutes_sr))])

# generate the values for the solute (and water?) uptake with dumux rosi:

#print(tests)

# for i in range(ntests):
    # cyl = RichardsWrapper(RichardsNCCylFoam())
    
    # points = [0.03,0.03*tests["rho"].iloc[0]] # radius of discretisation
    # NC = 100 # number of discretisations + 1
    # seg_length = 1
    # cyl.createGrid1d(points)# cm
    # cyl.setParameter( "Soil.Grid.Cells",str( NC-1))
    # cyl.seg_length = seg_length
    # cyl.setTopBC("noFlux")
    # cyl.setBotBC("constantPressure")
    # cyl.setTopBC_solute("noFlux")
    # cyl.setBotBC_solute("michaelisMenten")
    # cyl.setParameter("RootSystem.Uptake.Vmax",tests.loc[str(i),"Vmax"])
    # cyl.setParameter("RootSystem.Uptake.Km",tests.loc[str(i),"Km"])
    # Cells = cyl.getCellCenters_().reshape(-1)
    # CellsStr = cyl.dumux_str(Cells/100)#cm -> m #changed by Erik
    # cyl.setParameter("Soil.IC.Z",CellsStr)# m
    # cyl.setParameter("Soil.IC.C"+str(j)+"Z",CellsStr) # m
    # cyl.setCriticalPressure(-15000)  # cm pressure head

print("All inputs:")
print(tests[inputs])
print("All outputs:")
print(tests[outputs])




































# #print("solute tests")

# #peri.create_integralDiffusion_lookup("hydrus_loam_solutes", peri.sp)
# peri.open_lookup_solutes("hydrus_loam_solutes")

# labels_sol = ["Phi_root", "Phi_soil", "c_bulk", "Vmax", "Km", "Ds", "waterflow", "c_rootsoil"]

# Intervals_sol = {
        # "Phi_root": [0.0,50],
        # "Phi_soil": [60,140],
        # "c_bulk": [1.5e-8,1.5e-7], 
        # "Vmax": [3.5e-6,3.6e-6],
        # "Km": [1.4e-7,1.5e-7],
        # "waterflow": [0.0,0.0],
        # }
# Intervals_sol = pd.DataFrame(Intervals_sol, index = ["min","max"])

# tests_sol = pd.DataFrame(index = range(ntests),columns = labels_sol)

# for one_label in Intervals_sol.columns:
    # min_val = Intervals_sol.loc["min",one_label]
    # max_val = Intervals_sol.loc["max",one_label]

    # testvalues = (max_val - min_val) * np.random.rand(ntests) + min_val * np.ones(ntests)
    # tests_sol[one_label] = testvalues

# Phi_root = np.array(tests_sol.loc[:,"Phi_root"])
# Phi_soil = np.array(tests_sol.loc[:,"Phi_soil"])
# c_bulk = np.array(tests_sol.loc[:,"c_bulk"])
# Vmax = np.array(tests_sol.loc[:,"Vmax"])
# Km = np.array(tests_sol.loc[:,"Km"])
# waterflow = np.array([(2*3.14*inner_kr[i]*(hsr[i]-rx[i])) for i in range(ntests)])
# Ds = 1.9e-5

# #c_rootsoil = peri.soil_root_solutes(Phi_root, Phi_soil, c_bulk, Vmax, Km, Ds, waterflow)

# tests_sol.loc[:,"c_rootsoil"] = c_rootsoil




# #print(tests_sol)
