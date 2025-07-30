""" how to create a look up table for nonlinear perirhzal resistances, 
    run on multiple threads using: mpiexec -n 4 python3 example7_3_lookup.py """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

from functional.Perirhizal import PerirhizalPython
import functional.van_genuchten as vg

peri = PerirhizalPython()
hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]
filename = "hydrus_loam"
sp = vg.Parameters(hydrus_loam)
vg.create_mfp_lookup(sp)
peri.create_lookup_mpi(filename, sp)
# peri.open_lookup(filename)
