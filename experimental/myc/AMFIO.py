import numpy as np
import glob
import os

# Output folders, saving and loading of the per-seed results (only numpy, no plantbox needed)

def petridishOutputDirs():
    # output folders of the petri dish scripts (petridishsimauto.py, petridishcombine.py, petridishtiming.py, petridishstats.py):
    # below the environment variable PETRIDISH_OUT (set by run_petridish.sh on the cluster), else below petridish_results in the
    # current folder. Returns a dict with "base" (combined results), "densities" (seed_*.npz), "plots" and "animation" (vtp), each ending with "/"
    base = os.environ.get("PETRIDISH_OUT", "petridish_results")
    dirs = {"base": base, "densities": os.path.join(base, "densities"), "plots": os.path.join(base, "plots"), "animation": os.path.join(base, "animation")}
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)
    return {k: os.path.join(d, "") for k, d in dirs.items()}

def save_simulation(filename, sim):
    # stores one simulation dict (only arrays/numbers, e.g. "times", "tip_dens", ...) as .npz
    np.savez(filename, **sim)

def load_simulations(pattern):
    # loads all .npz files matching pattern (e.g. "results/seed_*.npz") into a list of dicts, sorted by filename
    simulations = []
    for f in sorted(glob.glob(pattern)):
        with np.load(f) as data:
            simulations.append({key: data[key] for key in data.files})
    # print("Loaded " + str(len(simulations)) + " simulations from " + pattern)
    return simulations

def load_petridish_results(pattern, verbose=True):
    # loads the per-seed results of petridishsimauto.py (e.g. "petridish_results/densities/seed_*.npz") that can be combined:
    # seeds stopped before t* (status "not_crossed", see max_wait_days) or because the front stalled after t* (status "stalled",
    # see min_growing_tips, stall_fraction) have no densities and are skipped, older files have no status.
    # Older files have densities per cm^3: with hyphae_layer (thin layer geometry) they are converted to mm^-2
    # (cm^-3 * layer [cm] = cm^-2 = mm^-2 * 100), without it (old geometry, densities per volume of the whole height) they are skipped.
    # "times" are shifted to days since t* (the absolute times of different seeds do not overlap)
    simulations = load_simulations(pattern)
    stopped = [int(sim["seed"]) for sim in simulations if str(sim.get("status", "ok")) != "ok"]
    simulations = [sim for sim in simulations if str(sim.get("status", "ok")) == "ok"]
    for sim in simulations:
        if "density_per" not in sim and "hyphae_layer" in sim:
            h = float(sim["hyphae_layer"])
            for key in ["tip_dens", "stopped_tip_dens"]:
                sim[key] = sim[key] * h / 100
            sim["length_densities"] = sim["length_densities"] * h / 10 # cm cm^-3 -> cm cm^-2 -> mm mm^-2
            sim["density_per"] = "mm^2"
    old = [int(sim["seed"]) for sim in simulations if "density_per" not in sim]
    simulations = [sim for sim in simulations if "density_per" in sim]
    for sim in simulations:
        sim["times"] = sim["times"] - sim["time_star"]
    if verbose:
        if old:
            print("skipping " + str(len(old)) + " seeds of the old geometry (rerun them): " + str(old))
        print(str(len(simulations)) + " simulations reached t*, " + str(len(stopped)) + " were stopped: " + str(stopped))
    return simulations
