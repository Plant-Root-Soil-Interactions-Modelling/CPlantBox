"""
Sammelt alle outputs.csv einer Designtabelle in einer Datei und meldet fehlende/
fehlgeschlagene Läufe.

    python collect.py results/design_morris.csv   -> results/outputs_morris.csv
"""
import json
import sys
from pathlib import Path

import pandas as pd

import config


def collect(design_path):
    design = pd.read_csv(design_path)
    tag = Path(design_path).stem.replace("design_", "")
    frames, status = [], []
    for rid in design["run_id"]:
        rdir = config.RUNS_DIR / rid
        meta = json.loads((rdir / "meta.json").read_text()) if (rdir / "meta.json").exists() else {}
        st = meta.get("sim_status", "missing")
        if st == "ok" and not (rdir / "outputs.csv").exists():
            st = "not_evaluated"
        status.append({"run_id": rid, "status": st, "sim_seconds": meta.get("sim_seconds")})
        if st == "ok":
            frames.append(pd.read_csv(rdir / "outputs.csv"))
    status = pd.DataFrame(status)
    status.to_csv(config.RESULTS_DIR / f"status_{tag}.csv", index=False)
    print(status["status"].value_counts().to_string())
    if frames:
        out = pd.concat(frames, ignore_index=True)
        path = config.RESULTS_DIR / f"outputs_{tag}.csv"
        out.to_csv(path, index=False)
        print(f"-> {path} ({len(out)} Zeilen)")
    t = status["sim_seconds"].dropna()
    if len(t):
        print(f"Laufzeit Simulation: Median {t.median():.0f} s, Max {t.max():.0f} s")


if __name__ == "__main__":
    collect(sys.argv[1])
