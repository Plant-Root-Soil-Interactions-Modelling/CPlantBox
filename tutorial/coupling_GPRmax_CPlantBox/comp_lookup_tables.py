""" how to create a look up table for nonlinear perirhzal resistances, 
    run on multiple threads using: mpiexec -n 4 python3 example7_3_lookup.py """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

from functional.Perirhizal import PerirhizalPython
import functional.van_genuchten as vg

peri = PerirhizalPython()  # |\label{l73l:peri}|
hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]  
hydrus_clay = [0.068, 0.38, 0.008, 1.09, 4.8]
hydrus_sand = [0.045, 0.43, 0.145, 2.68, 712.8]
hydrus_sandyloam = [0.065, 0.41, 0.075, 1.89, 106.1]

#Lena
# subsoil = [0.002, 0.25, 0.04, 1.12, 200]
# topsoil = [0.009, 0.39, 0.009, 1.36, 160]

#layered soil, Jana Bauer et al. 2011, works with CPlantBox / dumux
l1 = [0.008, 0.389, 0.012, 1.97, 91.68] # 0−20 cm
l2 = [0.008, 0.389, 0.023, 1.23, 63.36] # 20−33 cm
l3 = [0.008, 0.389, 0.01, 1.1, 10] # 33−150 cm

# soils = [hydrus_clay, hydrus_sand, hydrus_sandyloam]
# soils_ = ['hydrus_clay', 'hydrus_sand', 'hydrus_sandyloam']

soils = [hydrus_sandyloam, hydrus_loam]
soils_ = ['hydrus_sandyloam', 'hydrus_loam']

for i in range(0, len(soils)): 
    filename = soils_[i]
    sp = vg.Parameters(soils[i])  # |\label{l73l:soil_end}|
    vg.create_mfp_lookup(sp)  # |\label{l73l:mfp}|
    peri.create_lookup_mpi("lookup/" + filename, sp)  # |\label{l73l:lookup}|
