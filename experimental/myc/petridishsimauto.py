import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import numpy as np
import os
import time
import AMFAnalysis as amf

# Usage: python3 petridishsimauto.py <seed>
# one seed per run, e.g. one per cluster job (see run_petridish.sh). Seed 0 means random seed in CPlantBox, so use seeds > 0.
# Results go to results_dir, the radial profile over all seeds is made by petridishcombine.py
# I changed the c++ code, so set might not be valid anymore
# good_seeds = [10, 25, 45, 50, 70, 105, 115, 125]
if len(sys.argv) != 2:
    sys.exit("Usage: python3 petridishsimauto.py <seed>")
seed = int(sys.argv[1])
results_dir = "petridish_results/densities/"
os.makedirs(results_dir, exist_ok=True)
path = "tomatoparameters/"
name = "TwoHyphaePlusBAS"
animation = False
show_3d = False      # vtk window after each simulation (blocks until closed)
hours_hyphae = 60   ### HIER VERÄNDERUNG DAUER HYPHEN SIMULATION
max_wait_days = 20  # stop if less than 10 cm hyphae are behind the barrier after this many simulated days (None = no limit)

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

filename = "petridish_nowait_seed_" + str(seed)
stats = {}
mycp, time_star, times = amf.makeSimulation(seed, path, name, height, petri_dish, small_hyphae_dish, half_dish, dt, nSteps, hours_hyphae, filename,
                                            animation=animation, max_wait_days=max_wait_days, stats=stats)
stats["max_wait_days"] = max_wait_days if max_wait_days is not None else np.inf
if stats["status"] != "ok": # store the run statistics anyway, petridishcombine.py skips this seed, petridishtiming.py counts it
    amf.save_simulation(results_dir + f"seed_{seed:03d}.npz", {"seed": seed, **stats})
    sys.exit(0)
t_analysis = time.perf_counter()
if show_3d:
    vp.plot_roots_and_container(mycp, petri_dish)
if not animation:
    amf.getMycSegmentAnalyser(mycp, stdwrite=True, filename="animation/" + filename, step=nSteps)

    # Ringe -> Zeile // Zeitabschnitt -> Spalte
mats = amf.getParasPerRing(["nodeTips", "anastomosis", "length"], times, mycp, rings)
sim = {
    "seed": seed,
    "time_star": time_star,
    "location": location,
    "times": times[1:],   # getParasPerRing returns len(times)-1 columns
    "tip_dens": amf.getActiveTipsPerRing(times, mycp, rings) / vol,   # actively growing tips at each time
    "stopped_tip_dens": mats["nodeTips"] / vol,                        # tips that stopped growing in each interval (old "tip_dens")
    "ana_densities": mats["anastomosis"],
    "length_densities": mats["length"] / vol,
    **stats,
    "wall_analysis": time.perf_counter() - t_analysis,
}
amf.save_simulation(results_dir + f"seed_{seed:03d}.npz", sim)
if (seed % 10 == 0):
    amf.plot_ring_heatmap(sim["times"], amf.ringCenters(location), sim["tip_dens"], "Tip Density [cm$^{-3}$]", f"plots/tip_densities_seed_{seed:03d}.png")
    amf.plot_ring_heatmap(sim["times"], amf.ringCenters(location), sim["ana_densities"], "Anastomoses per ring", f"plots/ana_densities_seed_{seed:03d}.png")
