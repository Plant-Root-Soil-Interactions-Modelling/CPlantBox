"""
Wertet die Rohausgabe EINES Laufs aus -> runs/<run_id>/outputs.csv

    python evaluate.py results/design_morris.csv 17
    python evaluate.py results/design_morris.csv 0-9

Getrennt von der Simulation, damit neue Outputs, Zeitpunkte oder Teilbereiche
später ohne Neusimulation berechnet werden können (dann mit --force).

Ein Output wird über (time, region, output) identifiziert; region "all" ist der
gesamte Zielbereich. Neue Outputs: Funktion in OUTPUTS ergänzen.
"""
import json
import sys

import numpy as np
import pandas as pd

import config
import model_hook
from run_one import parse_rows


def _box_volume(box):
    return float(np.prod(np.subtract(box["max"], box["min"])))


def _in_box(pts, box):
    return np.all((pts >= box["min"]) & (pts < box["max"]), axis=1)


def region_outputs(snap, box):
    mid = 0.5 * (snap["seg_start"] + snap["seg_end"])
    m = _in_box(mid, box)
    seg_len = np.linalg.norm(snap["seg_end"] - snap["seg_start"], axis=1)[m]
    vol = _box_volume(box)
    L = float(seg_len.sum())
    n_tips = int(_in_box(snap["tips"], box).sum()) if len(snap["tips"]) else 0
    n_ana = int(_in_box(snap["anastomoses"], box).sum()) if len(snap["anastomoses"]) else 0
    return {
        "total_length": L,                       # [cm]
        "length_density": L / vol,               # [cm/cm³]
        "n_tips": n_tips,
        "tip_density": n_tips / vol,             # [1/cm³]
        "n_anastomoses": n_ana,
        "carbon_cost": model_hook.carbon_cost(seg_len, snap["radius"][m], snap["subtype"][m]),
    }


def evaluate_run(run_id, force=False):
    rdir = config.RUNS_DIR / run_id
    meta_path = rdir / "meta.json"
    if not meta_path.exists():
        print(f"{run_id}: noch nicht simuliert"); return
    meta = json.loads(meta_path.read_text())
    out_path = rdir / "outputs.csv"
    if meta.get("sim_status") != "ok":
        print(f"{run_id}: Simulation fehlgeschlagen, keine Auswertung"); return
    if out_path.exists() and not force:
        print(f"{run_id}: bereits ausgewertet"); return
    rows = []
    regions = {"all": config.DOMAIN, **config.REGIONS}
    for t, name in zip(config.SNAPSHOT_TIMES, meta["snapshots"]):
        snap = model_hook.load_snapshot(rdir / name)
        for rname, box in regions.items():
            for k, v in region_outputs(snap, box).items():
                rows.append({"run_id": run_id, "time": t, "region": rname, "output": k, "value": v})
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"{run_id}: {len(rows)} Werte")


if __name__ == "__main__":
    design = pd.read_csv(sys.argv[1])
    force = "--force" in sys.argv
    for i in parse_rows(sys.argv[2]):
        if i < len(design):
            evaluate_run(design.iloc[i]["run_id"], force)
