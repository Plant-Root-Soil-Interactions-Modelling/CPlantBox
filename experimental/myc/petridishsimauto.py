import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import numpy as np
import time
import AMFAnalysis as amf

# I changed the c++ code, so set might not be valid anymore
# good_seeds = [10, 25, 45, 50, 70, 105, 115, 125]
good_seeds = [10 + i for i in range(100)]
# good_seeds = [25, 10]
path = "tomatoparameters/"
name = "TwoHyphaePlusBAS"
animation = False
show_3d = False      # vtk window after each simulation (blocks until closed)
hours_hyphae = 60   ### HIER VERÄNDERUNG DAUER HYPHEN SIMULATION

## Setting up petri dish
diameter = 9.4
radius = diameter / 2
height = 1.6
# introduce parameters for barrier and opening
barrier_thickness = 0.16
barrier_height = height
opening_length = 5.0
opening_height = 0.2

nRings = 15
petri_dish, small_hyphae_dish, half_dish, rings = amf.makedishes(diameter, height, barrier_thickness, barrier_height, opening_length, opening_height, 0.1, nRings)
location = amf.ringLocations(radius, nRings)
vol = amf.ringVolume(radius, nRings, height, barrier_thickness)[:, None] # one volume per ring (column vector)

dt, nSteps = amf.setUpSimulationTime(5, 24)

simulations = []
allsims = time.perf_counter()
for seed in good_seeds:
    filename = "petridish_nowait_seed_" + str(seed)
    start = time.perf_counter()
    mycp, time_star, times = amf.makeSimulation(seed, path, name, height, petri_dish, small_hyphae_dish, half_dish, dt, nSteps, hours_hyphae, filename, animation=animation)
    print("Time for simulation: ", time.perf_counter() - start)
    if show_3d:
        vp.plot_roots_and_container(mycp, petri_dish)
    if not animation:
        amf.getMycSegmentAnalyser(mycp, stdwrite=True, filename="animation/" + filename, step=nSteps)

    # Ringe -> Zeile // Zeitabschnitt -> Spalte
    mats = amf.getParasPerRing(["nodeTips", "anastomosis", "length"], times, mycp, rings)
    simulations.append({
        "times": times[1:],   # getParasPerRing returns len(times)-1 columns
        "tip_dens": mats["nodeTips"] / vol,
        "ana_densities": mats["anastomosis"],
        "length_densities": mats["length"] / vol,
    })
print("Time for all simulations: ", time.perf_counter() - allsims)

for sim, seed in zip(simulations, good_seeds):
    if (seed % 10 == 0):
        amf.plot_ring_heatmap(sim["times"], location, sim["tip_dens"], "Tip Density [mm$^{-3}$]", f"plots/tip_densities_seed_{seed:03d}.png")
        amf.plot_ring_heatmap(sim["times"], location, sim["ana_densities"], "Anastomoses per ring", f"plots/ana_densities_seed_{seed:03d}.png")

t_common, tip_dens_interp = amf.interpolate_to_common_time(simulations, "tip_dens")
amf.plot_radial_profiles(location, t_common, tip_dens_interp, indices=[0, 30, 60, 90], ylabel="Tip density")
