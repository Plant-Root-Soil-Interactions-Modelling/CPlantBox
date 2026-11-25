"""
Runs ONE simulation from a design table.

    python run_one.py results/design_morris.csv 17        # row 17
    python run_one.py results/design_morris.csv 0-9       # rows 0..9 (local testing)

Under SLURM the row is passed from $SLURM_ARRAY_TASK_ID (see slurm/).
Restartable: runs with status "ok" in meta.json are skipped.
"""
import json
import sys
import time
import traceback

import numpy as np
import pandas as pd

import config
import model_hook
from xml_tools import write_xml


def run_row(row):
    rdir = config.RUNS_DIR / row["run_id"]
    meta_path = rdir / "meta.json"
    if meta_path.exists() and json.loads(meta_path.read_text()).get("sim_status") in ("ok", "excluded"):
        print(f"{row['run_id']}: already done, skipped")   # excluded (≤ 1 cm without t*) would give the same again
        return
    rdir.mkdir(parents=True, exist_ok=True)
    # factor list of the design: new designs store factor_set = "v2", the first Morris round has none (v1)
    fs = row.get("factor_set", np.nan)
    factors = config.factor_set("v1" if pd.isna(fs) else str(fs))
    values = {f.name: float(row[f.name]) for f in factors}
    # inoculum factors (only in design_inoculum.csv): root scaling into the XML, the others as arguments of the simulation
    inoc = [(f, float(row[f.name])) for f in config.INOCULUM_FACTORS if f.name in row and not pd.isna(row[f.name])]
    values.update({f.name: v for f, v in inoc if f.target == "root_scale"})
    sim_args = {f.param: v for f, v in inoc if f.target == "sim"}
    if not pd.isna(row.get("hyphae_layer", np.nan)):                     # only in the layer test of Phase 0
        sim_args["hyphae_layer"] = float(row["hyphae_layer"])
    dx = None if pd.isna(row.get("dx")) else float(row["dx"])
    als = row.get("ana_L_scale", np.nan)                                  # only in the anaMinDist test of Phase 0
    xml_path = write_xml(values, rdir / "params.xml", dx=dx, ana_L_scale=None if pd.isna(als) else float(als),
                         factors=factors)

    hours = None if pd.isna(row.get("hours", np.nan)) else float(row["hours"])   # only set in the 120 h test
    spd = None if pd.isna(row.get("steps_per_day", np.nan)) else int(row["steps_per_day"])   # only in the dt test
    ro = row.get("roots_on", np.nan)                                       # only in the long test
    roots_on = False if pd.isna(ro) else str(ro).strip().lower() in ("true", "1", "1.0")
    meta = {"run_id": row["run_id"], "seed": int(row["seed"]), "dx": dx, "hours": hours,
            "steps_per_day": spd or config.STEPS_PER_DAY, "roots_on": roots_on, "sim_args": sim_args,
            "values": values, "dev_mode": config.DEV_MODE}
    t0 = time.time()
    try:
        snaps = model_hook.simulate(xml_path, int(row["seed"]), rdir, hours=hours, steps_per_day=spd, roots_on=roots_on,
                                    sim_args=sim_args)
        meta.update(sim_status="ok", snapshots=[p.name for p in snaps])
    except model_hook.NoTstar as e:
        # as the experiment: no more than EXCLUDE_CROSSED_CM hyphae on the fungal side = unhealthy colony, excluded
        excluded = e.crossed_cm <= config.EXCLUDE_CROSSED_CM
        meta.update(sim_status="excluded" if excluded else "failed", crossed_cm=e.crossed_cm, error=traceback.format_exc())
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
