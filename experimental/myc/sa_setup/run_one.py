"""
Führt EINEN Simulationslauf aus einer Designtabelle aus.

    python run_one.py results/design_morris.csv 17        # Zeile 17
    python run_one.py results/design_morris.csv 0-9       # Zeilen 0..9 (lokal testen)

Unter SLURM wird die Zeile aus $SLURM_ARRAY_TASK_ID übergeben (s. slurm/).
Neustartfähig: Läufe mit status "ok" in meta.json werden übersprungen.
"""
import json
import sys
import time
import traceback

import pandas as pd

import config
import model_hook
from xml_tools import write_xml


def run_row(row):
    rdir = config.RUNS_DIR / row["run_id"]
    meta_path = rdir / "meta.json"
    if meta_path.exists() and json.loads(meta_path.read_text()).get("sim_status") == "ok":
        print(f"{row['run_id']}: bereits fertig, übersprungen")
        return
    rdir.mkdir(parents=True, exist_ok=True)
    names = [f.name for f in config.FACTORS]
    values = {n: float(row[n]) for n in names}
    dx = None if pd.isna(row.get("dx")) else float(row["dx"])
    xml_path = write_xml(values, rdir / "params.xml", dx=dx)

    meta = {"run_id": row["run_id"], "seed": int(row["seed"]), "dx": dx,
            "values": values, "dev_mode": config.DEV_MODE}
    t0 = time.time()
    try:
        snaps = model_hook.simulate(xml_path, int(row["seed"]), rdir)
        meta.update(sim_status="ok", snapshots=[p.name for p in snaps])
    except Exception:
        meta.update(sim_status="failed", error=traceback.format_exc())
    meta["sim_seconds"] = round(time.time() - t0, 1)
    meta_path.write_text(json.dumps(meta, indent=1))
    print(f"{row['run_id']}: {meta['sim_status']} ({meta['sim_seconds']} s)")


def parse_rows(spec):
    if "-" in spec:
        a, b = spec.split("-")
        return range(int(a), int(b) + 1)
    return [int(spec)]


if __name__ == "__main__":
    design = pd.read_csv(sys.argv[1])
    for i in parse_rows(sys.argv[2]):
        if i < len(design):
            run_row(design.iloc[i])
