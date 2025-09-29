import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import matplotlib.pyplot as plt
import numpy as np
import math
import time

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
    # rp.hyphalEmergenceDensity = 1
    rp.highresolution = 1

local = False


# --- Simulationseinstellungen ---
simtime = 10            # Gesamtdauer (Tage)
fpd = 25                # Schritte pro Tag
N = simtime * fpd       # Gesamtanzahl Schritte
dt = simtime / N        # Schrittweite in Tagen
timespan = np.linspace(0, simtime, N)
ratio = False            # Set to True for ratio plots, False for absolute values


def runsimulation(seed):
    start = time.perf_counter()
    mycp.setSeed(seed)
    mycp.initialize(True)

    for t in range(0, len(timespan)):
        mycp.simulate(dt, False)
    end = time.perf_counter()
    return (end - start)

runtimes = []
seed = 10

for runs in range(1,101):
    runtimes.append(runsimulation(seed))

# print(f"Total simulation time for {seed} runs: {end - start:.2f} seconds")
plt.hist(runtimes, bins=20, edgecolor='black')
plt.title('Runtime Distribution of Simulation')
plt.xlabel('Runtime (Seconds)')
plt.ylabel('Number of Runs')
plt.show()

# Statistik
print(f"Average Run: {sum(runtimes)/len(runtimes):.6f} Sekunden")
print(f"Fastest Run: {min(runtimes):.6f} Sekunden")
print(f"Slowest Run: {max(runtimes):.6f} Sekunden")