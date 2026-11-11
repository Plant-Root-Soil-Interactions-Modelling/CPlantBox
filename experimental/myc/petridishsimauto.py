import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import numpy as np
import inspect
import os
import resource
import socket
import time
import AMFAnalysis as amf
import AMFDensities as amfd
import AMFIO as amfio
import AMFPlots as amfp
import AMFSimulationSetup as amfs

def peak_memory_mb():
    # peak resident memory of this process so far (ru_maxrss is in kB on Linux)
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

# Usage: python3 petridishsimauto.py <seed> [name=value ...]
# one seed per run, e.g. one per cluster job (see run_petridish.sh). Seed 0 means random seed in CPlantBox, so use seeds > 0.
# name=value: hyphal parameters set for all hyphal sub types instead of the values of the xml, e.g. distTH=0.02 or ana=0
# (no anastomosis, the network grows without fusion), see amfs.setHyphaeParameters; stored in run_info.json
# Results go below $PETRIDISH_OUT (set in run_petridish.sh, else petridish_results) into folders named by the run
# ($PETRIDISH_RUN, set by submit_petridish.sh, else date and time now), with the parameters in run_info.json, see AMFIO.py.
# The radial profile over all seeds is made by petridishcombine.py
# I changed the c++ code, so set might not be valid anymore
# good_seeds = [10, 25, 45, 50, 70, 105, 115, 125]
if len(sys.argv) < 2 or not all(a.count("=") == 1 for a in sys.argv[2:]):
    sys.exit("Usage: python3 petridishsimauto.py <seed> [name=value ...]")
seed = int(sys.argv[1])
hyphae_parameters = {k: float(v) for k, v in (a.split("=") for a in sys.argv[2:])}
tipdist_R_max = 0.1  # [cm] smallest distance of each growing tip to another hyphal tree is stored up to this distance (amfd.getTipDistances)
out = amfio.petridishOutputDirs()
results_dir = out["densities"]
path = "tomatoparameters/"
name = "TwoHyphaePlusBAS"
animation = False
show_3d = False      # vtk window after each simulation (blocks until closed)
hours_hyphae = 80   ### HIER VERÄNDERUNG DAUER HYPHEN SIMULATION
max_wait_days = 30  # stop if less than 10 cm hyphae are behind the barrier after this many simulated days (None = no limit)
# stop after t* if the front stalls (subtrees of dead hyphae are not simulated anymore): fewer than min_growing_tips growing tips
# behind the barrier, or fewer than stall_fraction times the maximum so far (None = no check). Such seeds get status "stalled".
min_growing_tips = 10
stall_fraction = 0.5
progress_every = 12  # one line of progress output every 12 time steps (= 12 h), None = no output

## Setting up petri dish (geometry of the experiment, height from the filled medium volumes, see amfs.defaultPetridish)
hyphae_layer = 0.1  # hyphae grow on the cellophane in the fungal compartment, i.e. in a thin layer at the top
nRings = 15
dish = amfs.defaultPetridish(nRings=nRings, hyphae_layer=hyphae_layer)
radius, height, barrier_thickness = dish["radius"], dish["height"], dish["barrier_thickness"]
petri_dish, small_hyphae_dish, half_dish, rings = dish["petri_dish"], dish["small_hyphae_dish"], dish["half_dish"], dish["rings"]
location = amfs.ringLocations(radius, nRings)
# densities per area as in Oyarte Galvez et al. 2025 (tips mm^-2, length mm mm^-2), the hyphae grow in a thin layer on the cellophane
area = amfs.ringVolume(radius, nRings, 1.0, barrier_thickness)[:, None] * 100 # area of each ring [mm^2] (column vector)

root_days, steps_per_day = 5, 24
dt, nSteps = amfs.setUpSimulationTime(root_days, steps_per_day)

filename = "petridish_nowait_seed_" + str(seed)
# options of makeSimulation, also stored in run_info.json (together with the defaults of all options not given here)
sim_options = {"animation": animation, "max_wait_days": max_wait_days, "min_growing_tips": min_growing_tips, "stall_fraction": stall_fraction,
               "progress_every": progress_every, "turn_off_root_growth": True, "hyphae_parameters": hyphae_parameters}
make_defaults = {k: p.default for k, p in inspect.signature(amfs.makeSimulation).parameters.items() if p.default is not inspect.Parameter.empty}
amfio.writeRunInfo(out, {
    "parameter_file": path + name + ".xml",
    "root_days": root_days, "steps_per_day": steps_per_day, "dt": dt, "nSteps_root": nSteps, "hours_hyphae": hours_hyphae,
    "makeSimulation": {**make_defaults, **sim_options, "stats": None, "animation_dir": None},
    "hyphae_parameters": hyphae_parameters, "tipdist_R_max": tipdist_R_max,
    "geometry": {k: v for k, v in dish.items() if not k in ("petri_dish", "small_hyphae_dish", "half_dish", "rings")},
    "slurm_job": os.environ.get("SLURM_ARRAY_JOB_ID") or os.environ.get("SLURM_JOB_ID"),
}, parameter_file=path + name + ".xml")
seed_info = {"seed": seed, "run": out["run"], "host": socket.gethostname(), "started": time.strftime("%Y-%m-%d %H:%M:%S")}
stats = {}
mycp, time_star, times = amfs.makeSimulation(seed, path, name, height, petri_dish, small_hyphae_dish, half_dish, dt, nSteps, hours_hyphae, filename,
                                             stats=stats, animation_dir=out["animation"], **sim_options)
stats["max_wait_days"] = max_wait_days if max_wait_days is not None else np.inf
if stats["status"] != "ok": # store the run statistics anyway, petridishcombine.py skips this seed, petridishtiming.py counts it
    amfio.save_simulation(results_dir + f"seed_{seed:03d}.npz", {**seed_info, **stats, "peak_memory_mb": peak_memory_mb()})
    sys.exit(0)
t_analysis = time.perf_counter()
if show_3d:
    vp.plot_roots_and_container(mycp, petri_dish)
if not animation:
    amf.getMycSegmentAnalyser(mycp, stdwrite=True, filename=out["animation"] + filename, step=nSteps)

    # Ringe -> Zeile // Zeitabschnitt -> Spalte
subtypes = amfd.HYPHAL_SUBTYPES # {1: "runner", 2: "bas"}
mats = amfd.getParasPerRing(["nodeTips", "anastomosis", "length"], times, mycp, rings, subTypes=list(subtypes), increments=["anastomosis"])
active_tips = amfd.getActiveTipsPerRingPerSubtype(times, mycp, rings)
zeros = np.zeros((len(rings), len(times) - 1))
per_subtype = {}
for st, key in subtypes.items():
    per_subtype["length_" + key] = mats[f"length_st{st}"] * 10 / area          # filament density of the sub type [mm mm^-2]
    per_subtype["active_tips_" + key] = active_tips.get(st, zeros)            # growing tips of the sub type (number per ring)
    per_subtype["ana_new_" + key] = mats[f"anastomosis_st{st}_new"]            # new anastomoses in each interval, sub type of the tip that fused
# smallest distance of each growing tip to a hypha of another tree, for the encounter fraction with any R = distTH afterwards (amfd.tipDistanceC)
t_tipdist = time.perf_counter()
tipdist = amfd.getTipDistances(times, mycp, rings, R_max=tipdist_R_max)
print(f"tip distances: {len(tipdist['dist'])} tips at {len(times) - 1} times in {time.perf_counter() - t_tipdist:.0f} s", flush=True)
sim = {
    **seed_info,             # seed, run, host, start time
    "time_star": time_star,
    "location": location,
    "height": height,
    "hyphae_layer": hyphae_layer,
    "density_per": "mm^2",   # files without this key have densities per cm^3
    "times": times[1:],   # getParasPerRing returns len(times)-1 columns
    "area": area[:, 0],                                                 # ring areas [mm^2]
    "tip_dens": sum(active_tips.values()) / area,                       # actively growing tips at each time [mm^-2]
    "active_tips": sum(active_tips.values()),                           # actively growing tips (number per ring)
    "ana_new": mats["anastomosis_new"],                                 # new anastomoses in each interval (number per ring)
    **per_subtype,
    "stopped_tip_dens": mats["nodeTips"] / area,                        # tips that stopped growing in each interval (old "tip_dens") [mm^-2]
    "ana_densities": mats["anastomosis"],
    "length_densities": mats["length"] * 10 / area,                     # [mm mm^-2]
    "tipdist_col": tipdist["col"].astype(np.int32),                     # one entry per growing tip and time: column (as active_tips),
    "tipdist_ring": tipdist["ring"].astype(np.int16),                   # ring,
    "tipdist_subtype": tipdist["subtype"].astype(np.int8),              # sub type and
    "tipdist_d": tipdist["dist"].astype(np.float32),                    # distance to the nearest hypha of another tree [cm] (inf: > tipdist_R_max)
    "tipdist_R_max": tipdist_R_max,
    **stats,
    "wall_analysis": time.perf_counter() - t_analysis,
    "peak_memory_mb": peak_memory_mb(),
}
amfio.save_simulation(results_dir + f"seed_{seed:03d}.npz", sim)
if (seed % 10 == 0):
    title = f"seed {seed}, run {out['run']}"
    amfp.plot_ring_heatmap(sim["times"], amfs.ringCenters(location), sim["tip_dens"], "Tip density [mm$^{-2}$]",
                           amfio.plotFile(out, f"tip_densities_seed_{seed:03d}.png"), title=title)
    amfp.plot_ring_heatmap(sim["times"], amfs.ringCenters(location), sim["ana_densities"], "Anastomoses per ring",
                           amfio.plotFile(out, f"ana_densities_seed_{seed:03d}.png"), title=title)
