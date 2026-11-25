"""
Calibration targets pooled per design point over its replicates, as the calibration pools over its seeds (B11).
Per run the wave fits fail too often for Morris; pooled over the replicates of a point they are defined almost always.

    python evaluate_pooled.py results/design_morris.csv     -> results/pooled_morris.csv
    python evaluate_pooled.py results/design_inoculum.csv   -> results/pooled_inoculum.csv

One row per (traj, step, output): value (evaluate.wave_quantities over the finished replicates, window
max(SNAPSHOT_HOURS) after t*), se (bootstrap over the replicates) and n_runs. Reads raw.npz, no re-simulation.
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import config
import evaluate as ev
import model_hook


def point_outputs(preps, rng):
    if len(preps) < config.POOL_MIN_RUNS:
        nan = dict.fromkeys(config.POOLED_OUTPUTS, np.nan)
        return nan, nan
    w = ev.wave_quantities(preps)
    boot = [ev.wave_quantities([preps[i] for i in rng.integers(0, len(preps), len(preps))]) for _ in range(config.POOL_BOOT)]
    se = {o: float(np.nanstd([b[o] for b in boot])) for o in config.POOLED_OUTPUTS}
    return {o: w[o] for o in config.POOLED_OUTPUTS}, se


def main(design_path):
    design = pd.read_csv(design_path)
    tag = Path(design_path).stem.replace("design_", "")
    hour = max(config.SNAPSHOT_HOURS)
    rng = np.random.default_rng(config.BASE_SEED)
    rows = []
    groups = design.groupby(["traj", "step"])
    for i, ((t, s), d) in enumerate(groups):
        preps = []
        for rid in d.run_id:
            meta_path = config.RUNS_DIR / rid / "meta.json"
            if not meta_path.exists():
                continue
            meta = json.loads(meta_path.read_text())
            if meta.get("sim_status") != "ok":
                continue
            raw = model_hook.load_raw(config.RUNS_DIR / rid / meta["snapshots"][0])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                preps.append(ev.prepare(ev.truncate(raw, hour), hour))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            val, se = point_outputs(preps, rng)
        rows += [{"traj": int(t), "step": int(s), "time": hour, "region": "all", "output": o, "value": val[o],
                  "se": se[o], "n_runs": len(preps)} for o in config.POOLED_OUTPUTS]
        if (i + 1) % 50 == 0:
            print(f"{i + 1} of {groups.ngroups} design points", flush=True)
    df = pd.DataFrame(rows)
    path = config.RESULTS_DIR / f"pooled_{tag}.csv"
    df.to_csv(path, index=False)
    print(f"-> {path} ({groups.ngroups} design points)")
    defined = df.groupby("output").value.apply(lambda v: f"{v.notna().sum()}/{len(v)}")
    print("Defined at design points:\n" + defined.to_string())


if __name__ == "__main__":
    main(sys.argv[1])
