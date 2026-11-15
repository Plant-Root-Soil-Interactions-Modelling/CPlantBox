"""
Schnittstelle zur Simulation: Petrischale wie petridishsimauto.py, mit dem Parameterfile des Laufs.

1. simulate(xml_path, seed, outdir) -> list[Path]
   Wurzelphase, Warten auf t* und HOURS_HYPHAE Stunden Hyphen (AMFSimulationSetup.makeSimulation),
   danach die Größen je Ring und Stunde wie in petridishsimauto.py. Gespeichert in outdir/raw.npz,
   daraus berechnet evaluate.py alle Outputs (neue Outputs ohne Neusimulation).
   Ohne t* nach MAX_WAIT_DAYS wird ein Fehler ausgelöst: der Lauf zählt als fehlgeschlagen.

2. load_raw(path) -> dict
   Liest raw.npz, "times" in Tagen seit t* (wie AMFIO.load_petridish_results).

Solange die Umgebungsvariable SA_DUMMY=1 gesetzt ist, wird ein Spielzeugmodell
verwendet (dummy_model.py), damit die Pipeline ohne die echte Simulation testbar ist.
"""
import os
import sys
from pathlib import Path

import numpy as np

import config

USE_DUMMY = os.environ.get("SA_DUMMY") == "1"


def import_amf():
    """plantbox und die AMF*-Module aus experimental/myc importierbar machen."""
    for p in (config.MYC_DIR, config.CPLANTBOX_DIR / "src", config.CPLANTBOX_DIR):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))


def simulate(xml_path: Path, seed: int, outdir: Path) -> list[Path]:
    if USE_DUMMY:
        import dummy_model
        return dummy_model.simulate(xml_path, seed, outdir)

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
        dish["half_dish"], dt, nSteps, config.HOURS_HYPHAE, "sim", max_wait_days=config.MAX_WAIT_DAYS, stats=stats,
        min_growing_tips=config.MIN_GROWING_TIPS, stall_fraction=config.STALL_FRACTION, turn_off_root_growth=True)
    if time_star is None:
        raise RuntimeError(f"status {stats['status']} (wait_days {stats.get('wait_days', np.nan):.1f})")
    if config.SAVE_VTP:
        amf.getMycSegmentAnalyser(mycp, stdwrite=True, filename=str(outdir / "final"), step="", timestamp=False)

    # wie petridishsimauto.py: Ringe -> Zeile, Zeitabschnitt -> Spalte
    rings = dish["rings"]
    mats = amfd.getParasPerRing(["nodeTips", "anastomosis", "length"], times, mycp, rings,
                                subTypes=list(amfd.HYPHAL_SUBTYPES), increments=["anastomosis"])
    active_tips = amfd.getActiveTipsPerRingPerSubtype(times, mycp, rings)
    zeros = np.zeros((len(rings), len(times) - 1))
    raw = {"seed": seed, "time_star": time_star, "location": location, "area": area[:, 0], "times": times[1:],
           "tip_dens": sum(active_tips.values()) / area, "active_tips": sum(active_tips.values()),
           "ana_new": mats["anastomosis_new"], "length_densities": mats["length"] * 10 / area,
           "wait_days": stats["wait_days"]}
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
