""" how to create a look up table for nonlinear perirhzal resistances, 
    run on multiple threads using: mpiexec -n 4 python3 example7_3_lookup.py """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

from functional.Perirhizal import PerirhizalPython
import functional.van_genuchten as vg

peri = PerirhizalPython()  # |\label{l73l:peri}|
hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]  # |\label{l73l:soil}|
filename = "hydrus_loam"
sp = vg.Parameters(hydrus_loam)  # |\label{l73l:soil_end}|
vg.create_mfp_lookup(sp)  # |\label{l73l:mfp}|
peri.create_lookup_mpi("results/" + filename, sp)  # |\label{l73l:lookup}|
