import numpy as np
import glob
import json
import os
import shutil
import socket
import subprocess
import time

# Output folders, saving and loading of the per-seed results, run metadata (only numpy, no plantbox needed)
#
# Every run (one set of seeds, e.g. one SLURM array job plus its combine job) has a name, the date and time it was started,
# e.g. "2026-11-11_14-30-05". All outputs of the run go into subfolders of that name:
#   petridish_results/densities/<run>/seed_010.npz, ...   per-seed results, run_info.json, copy of the parameter file
#   petridish_results/plots/<run>/<plot>_<run>.png         all plots (per seed and combined), run_info.json, wave fits csv
#   petridish_results/animation/<run>/*.vtp                vtp files, run_info.json
# The run name comes from the environment variable PETRIDISH_RUN (set once by submit_petridish.sh, so that all array tasks
# and the combine job write into the same folders). Without it, a simulation starts a new run (to run several seeds locally
# into one run: export PETRIDISH_RUN=$(date +%Y-%m-%d_%H-%M-%S)) and the combine scripts use the newest run.

RUN_INFO = "run_info.json"

def newRunName():
    # date and time now, sorts chronologically
    return time.strftime("%Y-%m-%d_%H-%M-%S")

def latestRun(base):
    # newest run folder below base/densities (last modified), None if there is none
    runs = [d for d in glob.glob(os.path.join(base, "densities", "*")) if os.path.isdir(d)]
    if not runs:
        return None
    return os.path.basename(max(runs, key=os.path.getmtime))

def petridishOutputDirs(run = None, latest = False, create = ("densities", "plots", "animation")):
    # output folders of the petri dish scripts (petridishsimauto.py, petridishcombine.py, petridishtiming.py, petridishstats.py):
    # below the environment variable PETRIDISH_OUT (set by run_petridish.sh on the cluster), else below petridish_results in the
    # current folder. Returns a dict with "base", "run" (name of the run) and the run folders "densities" (seed_*.npz),
    # "plots" and "animation" (vtp), each ending with "/"
    # run: name of the run, default: $PETRIDISH_RUN; without both, latest = True uses the newest existing run (combine scripts),
    #      latest = False starts a new run named by the current date and time (simulation)
    # create: run folders that are created if missing (the others are only returned)
    base = os.environ.get("PETRIDISH_OUT", "petridish_results")
    if run is None:
        run = os.environ.get("PETRIDISH_RUN") or None
    if run is None:
        run = latestRun(base) if latest else newRunName()
        if run is None:
            raise SystemExit("no run folder in " + os.path.join(base, "densities") + ": set PETRIDISH_RUN, or move older seed_*.npz into a subfolder densities/<name>/")
    dirs = {"densities": os.path.join(base, "densities", run), "plots": os.path.join(base, "plots", run), "animation": os.path.join(base, "animation", run)}
    for k in create:
        os.makedirs(dirs[k], exist_ok=True)
    out = {k: os.path.join(d, "") for k, d in dirs.items()}
    out.update({"base": os.path.join(base, ""), "run": run})
    return out

def plotFile(out, name):
    # file name of a plot in the plots folder of the run, with the run name (date and time) appended: "TipDens.png" -> plots/<run>/TipDens_<run>.png
    stem, ext = os.path.splitext(name)
    return out["plots"] + stem + "_" + out["run"] + (ext or ".png")

def gitInfo():
    # commit, branch and locally changed (tracked) files of the CPlantBox repository, so the C++ state of a run is known
    def git(*args):
        try:
            r = subprocess.run(["git", "-C", os.path.dirname(os.path.abspath(__file__))] + list(args), capture_output=True, text=True, timeout=30)
            return r.stdout.rstrip() if r.returncode == 0 else None # no strip(): status lines start with a space
        except (OSError, subprocess.SubprocessError):
            return None
    changed = git("status", "--porcelain", "--untracked-files=no")
    return {"commit": git("rev-parse", "HEAD"), "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
            "changed_files": [l[3:] for l in changed.splitlines()] if changed else []}

def _jsonable(o):
    # numpy values and other objects for json.dump
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.generic):
        return o.item()
    return str(o)

def _writeJson(filename, data):
    # write to a temporary file first, then rename: several array tasks may write the same file at the same time
    tmp = filename + "." + str(os.getpid()) + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2, default=_jsonable)
    os.replace(tmp, filename)

def readRunInfo(folder):
    # run_info.json of a run folder, {} if there is none
    try:
        with open(os.path.join(folder, RUN_INFO)) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}

def writeRunInfo(out, info, folders = ("densities", "plots", "animation"), parameter_file = None):
    # stores the metadata of the run (parameters, git state, ...) as run_info.json in the run folders, once per run
    # (the first array task writes it, the others find it), parameter_file (xml) is copied into densities and plots
    info = {"run": out["run"], "written": time.strftime("%Y-%m-%d %H:%M:%S"), "git": gitInfo(), **info}
    for key in folders:
        filename = os.path.join(out[key], RUN_INFO)
        if not os.path.exists(filename):
            _writeJson(filename, info)
        if parameter_file is not None and key in ("densities", "plots"):
            target = os.path.join(out[key], os.path.basename(parameter_file))
            if not os.path.exists(target):
                shutil.copyfile(parameter_file, target)

def updatePlotInfo(out, section, data):
    # run_info.json in the plots folder: metadata of the simulations (from densities) plus data of a combine script under section
    filename = os.path.join(out["plots"], RUN_INFO)
    info = readRunInfo(out["plots"]) or readRunInfo(out["densities"]) or {"run": out["run"]}
    info[section] = {"written": time.strftime("%Y-%m-%d %H:%M:%S"), "host": socket.gethostname(), **data}
    _writeJson(filename, info)

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
    # loads the per-seed results of petridishsimauto.py (e.g. "petridish_results/densities/<run>/seed_*.npz") that can be combined:
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
