"""
Interface to the simulation: Petri dish as petridishsimauto.py, with the parameter file of the run.

1. simulate(xml_path, seed, outdir) -> list[Path]
   Root phase, waiting for t* and HOURS_HYPHAE hours of hyphae (AMFSimulationSetup.makeSimulation),
   then the quantities per ring and hour as in petridishsimauto.py. Stored in outdir/raw.npz,
   from which evaluate.py computes all outputs (new outputs without re-simulating).
   Without t* after MAX_WAIT_DAYS an error is raised: the run counts as failed.

2. load_raw(path) -> dict
   Reads raw.npz, "times" in days since t* (as AMFIO.load_petridish_results).

While the environment variable SA_DUMMY=1 is set, a toy model is used
(dummy_model.py), so that the pipeline can be tested without the real simulation.
"""
import os
import sys
from pathlib import Path

import numpy as np

import config

USE_DUMMY = os.environ.get("SA_DUMMY") == "1"


def import_amf():
    """Make plantbox and the AMF* modules from experimental/myc importable."""
    for p in (config.MYC_DIR, config.CPLANTBOX_DIR / "src", config.CPLANTBOX_DIR):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))


def simulate(xml_path: Path, seed: int, outdir: Path, hours: float | None = None) -> list[Path]:
    """hours: simulation time after t* [h], default config.HOURS_HYPHAE (for the 120 h test)."""
    hours = int(hours) if hours else config.HOURS_HYPHAE
    if USE_DUMMY:
        import dummy_model
        return dummy_model.simulate(xml_path, seed, outdir, hours)

    import_amf()
    import AMFAnalysis as amf
    import AMFDensities as amfd
    import AMFSimulationSetup as amfs

    dish = amfs.defaultPetridish(nRings=config.N_RINGS, hyphae_layer=config.HYPHAE_LAYER, verbose=False)
    location = amfs.ringLocations(dish["radius"], config.N_RINGS)
    area = amfs.ringVolume(dish["radius"], config.N_RINGS, 1.0, dish["barrier_thickness"])[:, None] * 100  # [mm^2]
    dt, nSteps = amfs.setUpSimulationTime(config.ROOT_DAYS, config.STEPS_PER_DAY)
    stats = {}
    mycp, time_star, times = amfs.makeSimulation(
        seed, str(xml_path.parent) + "/", xml_path.stem, dish["height"], dish["petri_dish"], dish["small_hyphae_dish"],
        dish["half_dish"], dt, nSteps, hours, "sim", max_wait_days=config.MAX_WAIT_DAYS, stats=stats,
        min_growing_tips=config.MIN_GROWING_TIPS, stall_fraction=config.STALL_FRACTION, turn_off_root_growth=True)
    if time_star is None:
        raise RuntimeError(f"status {stats['status']} (wait_days {stats.get('wait_days', np.nan):.1f})")
    if config.SAVE_VTP:
        amf.getMycSegmentAnalyser(mycp, stdwrite=True, filename=str(outdir / "final"), step="", timestamp=False)

    # as petridishsimauto.py: rings -> rows, time steps -> columns
    rings = dish["rings"]
    mats = amfd.getParasPerRing(["nodeTips", "anastomosis", "length"], times, mycp, rings,
                                subTypes=list(amfd.HYPHAL_SUBTYPES), increments=["anastomosis"])
    active_tips = amfd.getActiveTipsPerRingPerSubtype(times, mycp, rings)
    zeros = np.zeros((len(rings), len(times) - 1))
    raw = {"seed": seed, "time_star": time_star, "location": location, "area": area[:, 0], "times": times[1:],
           "tip_dens": sum(active_tips.values()) / area, "active_tips": sum(active_tips.values()),
           "ana_new": mats["anastomosis_new"], "length_densities": mats["length"] * 10 / area,
           "branch_new": amfd.getBranchingPerRing(times, mycp, rings),
           "branch_new_runner": amfd.getBranchingPerRing(times, mycp, rings, child_subtype=config.RUNNER),
           "wait_days": stats["wait_days"]}
    if config.ROOT_EFFICIENCY:   # network at the end of the simulation, planar like the image graphs of the experiment (crossings = nodes)
        import time
        import traceback
        t0 = time.time()
        try:   # an error here (including missing networkx) should only cost this output, not the run
            import AMFNetwork as amfn
            eff = amfn.relativeRootEfficiency(mycp, dish["small_hyphae_dish"], dish["barrier_thickness"], dim=2)
            raw.update({"root_eff_" + k: v for k, v in eff.items()})
        except Exception:
            raw["root_eff_E_rel"] = np.nan
            raw["root_eff_error"] = traceback.format_exc()[-2000:]
        raw["root_eff_seconds"] = time.time() - t0
    for st, key in amfd.HYPHAL_SUBTYPES.items():
        raw["length_" + key] = mats[f"length_st{st}"] * 10 / area
        raw["active_tips_" + key] = active_tips.get(st, zeros)
        raw["ana_new_" + key] = mats[f"anastomosis_st{st}_new"]
    path = outdir / "raw.npz"
    np.savez_compressed(path, **raw)
    return [path]


def load_raw(path: Path) -> dict:
    with np.load(path) as d:
        raw = {k: d[k] for k in d.files}
    raw["times"] = raw["times"] - raw["time_star"]
    return raw
