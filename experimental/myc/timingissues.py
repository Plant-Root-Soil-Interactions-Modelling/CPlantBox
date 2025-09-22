import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import matplotlib.pyplot as plt
import numpy as np
import math
import time as t

# input("Only run this script if secondary infection is disabled.\nPress Enter to continue...")
mycp = pb.MycorrhizalPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Glycine_max"
# name = "Heliantus_Pagès_2013"
mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

hyphae_parameter = pb.HyphaeRandomParameter(mycp)
hyphae_parameter.subType = 1
# hyphae_parameter.a = 0.01
hyphae_parameter.vs = 0.01
mycp.setOrganRandomParameter(hyphae_parameter)
# print(hyphae_parameter)

root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 1;

infbox = pb.SDF_PlantBox(3, 3, 3)
infbox = pb.SDF_RotateTranslate(infbox, 0, 0, pb.Vector3d(0, 0, -10))
local = False
for i in range(0, len(root)):
    if local:
        root[i].f_inf = pb.SoilLookUpSDF(infbox, 1, 0.0, 0.1)
        print("ATTENTION: The infection radius is not 0, the infection will be local")
    root[i].dx = 0.05
    



# --- Simulationseinstellungen ---
simtime = 10            # Gesamtdauer (Tage)
fpd = 25                # Schritte pro Tag
N = simtime * fpd       # Gesamtanzahl Schritte
dt = simtime / N        # Schrittweite in Tagen
timespan = np.linspace(0, simtime, N)
ratio = False            # Set to True for ratio plots, False for absolute values


def runsimulation(seed):
    mycp.setSeed(seed)
    mycp.initialize(True)

    for t in range(0, len(timespan)):
        mycp.simulate(dt, False)
        
    return ()

start = t.time()
for seed in range(1, 51):
    print(f"Running simulation with seed {seed}")
    runsimulation(seed)
end = t.time()
print(f"Total simulation time for {seed} runs: {end - start:.2f} seconds")