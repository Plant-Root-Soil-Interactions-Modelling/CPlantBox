import sys; sys.path.append("../.."); sys.path.append("../../src/")

import numpy as np
import AMFDensities as amfd
import AMFIO as amfio
import AMFPlots as amfp

# Compares several runs of petridishsimauto.py, e.g. with different anastomosis distances (submit_petridish.sh distTH=0.02, ...)
# and one without anastomoses (ana=0): C of f = C rho against R = distTH, from the realized anastomoses per hour of each run and from
# the distances of the growing tips to other hyphal trees (any R, see plot_tip_distances), with the random straight hyphae of hyphal_network.py.
# Usage: python3 petridishcompare.py <run> <run> ...   (names of the run folders in petridish_results/densities)
# Writes C_against_R_<date>.png and run_info.json (the runs and their C) into petridish_results/plots/compare_<date>
if len(sys.argv) < 2:
    sys.exit("Usage: python3 petridishcompare.py <run> <run> ...")
runs = sys.argv[1:]
out = amfio.petridishOutputDirs(run="compare_" + amfio.newRunName(), create=("plots",)) # no densities folder, it would count as newest run
min_tips = 5

results, layer = [], None
for run in runs:
    run_out = amfio.petridishOutputDirs(run=run, create=())
    simulations = amfio.load_petridish_results(run_out["densities"] + "seed_*.npz", verbose=False)
    if not simulations:
        print(f"run {run}: no simulations, skipped")
        continue
    info = amfio.readRunInfo(run_out["densities"])
    changed = info.get("hyphae_parameters") or {}
    distTH = float(np.max(simulations[0].get("hyphae_distTH", [np.nan])))
    ana = float(np.max(simulations[0].get("hyphae_ana", [np.nan])))
    layer = float(simulations[0]["hyphae_layer"])
    label = ", ".join(f"{k}={v:g}" for k, v in changed.items()) or "xml"
    res = {"run": run, "label": f"{label} ({len(simulations)} seeds)", "distTH": distTH, "ana": ana, "n_seeds": len(simulations),
           "R": None, "C_geometric": None, "C_realized": np.nan}
    if ana > 0 and all(k in simulations[0] for k in ["ana_new", "active_tips"]):
        res["C_realized"] = amfd.encounterC(simulations, min_tips)["total"]
    with_tipdist = [sim for sim in simulations if "tipdist_d" in sim]
    if with_tipdist:
        R_max = min(float(sim["tipdist_R_max"]) for sim in with_tipdist)
        res["R"] = np.geomspace(0.001, R_max, 30)
        res["C_geometric"], _ = amfd.tipDistanceC(with_tipdist, res["R"], min_tips)
    print(f"run {run}: {res['label']}, distTH = {distTH:g} cm, ana = {ana:g}, realized C = {res['C_realized']:.3e} mm^2 um^-1 h^-1"
          + (f", geometric C at R = distTH: {np.interp(distTH, res['R'], res['C_geometric']):.3e} mm^2 um^-1" if res["R"] is not None and distTH > 0 else ""))
    results.append(res)

if not results:
    sys.exit("no run with simulations found")
amfp.plot_C_against_R(results, layer, filename=amfio.plotFile(out, "C_against_R.png"))
amfio.updatePlotInfo(out, "petridishcompare", {"runs": results, "min_tips": min_tips})
print("written " + amfio.plotFile(out, "C_against_R.png"))
