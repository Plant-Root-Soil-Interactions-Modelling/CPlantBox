"""
Ring resolution: do the outputs depend on the 15 rings of equal area? Every GROUP neighbouring rings are merged into one
(15 -> 5 with GROUP = 3; the merged rings again have equal area) and all outputs are computed again from raw.npz.
No re-simulation needed.

    python analyze_rings.py                                  # Phase 0 noise runs (nominal point), window max(SNAPSHOT_HOURS) = 120 h
    python analyze_rings.py results/design_long.csv long 180 # design, group, window [h after t*]

    -> results/rings_compare_<design>.csv   per output: medians and paired relative difference 5 vs. 15 rings
    Per run (as the SA) and pooled over all runs (as the calibration: runner tips in the wave fits).
"""
import json
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

import config
import evaluate as ev
import model_hook

GROUP = 3
DENSITY = ("length_densities", "tip_dens", "length_runner", "length_bas")   # per area: area-weighted mean, others summed


def merge_rings(raw, g=GROUP):
    n = len(raw["location"])
    if n % g:
        raise ValueError(f"{n} rings cannot be merged in groups of {g}")
    area = np.asarray(raw["area"], dtype=float).ravel()
    area_new = area.reshape(-1, g).sum(1)
    out = dict(raw, location=np.asarray(raw["location"])[g - 1::g], area=area_new)
    for k, v in raw.items():
        if isinstance(v, np.ndarray) and v.ndim == 2 and v.shape[0] == n:
            if k in DENSITY:
                out[k] = (v * area[:, None]).reshape(-1, g, v.shape[1]).sum(1) / area_new[:, None]
            else:
                out[k] = v.reshape(-1, g, v.shape[1]).sum(1)
    return out


def per_run(raw, hour):
    """Outputs of the SA that can depend on the rings (others, e.g. total_length, are sums over all rings)."""
    snap = ev.snapshot_outputs(raw, hour)
    wave = ev.wave_outputs(raw, hour)
    keep = ["front_position", "spread_r90", "max_tip_density"]
    return {**{k: snap[k] for k in keep},
            **{k: wave[k] for k in ("front_velocity", "v_tips", "K1", "K2", "wave_width", "alpha", "beta", "encounter_C",
                                    "n_wave_rings", "n_K1_rings", "n_tip_rings")}}


POOLED = ["alpha", "beta", "K1", "K2", "wave_width", "front_velocity", "n_wave_rings", "n_K1_rings", "n_tip_rings"]


def pooled(raws, hour):
    """Wave quantities pooled over all runs, as the calibration and evaluate_pooled.py (evaluate.wave_quantities)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w = ev.wave_quantities([ev.prepare(ev.truncate(r, hour), hour) for r in raws])
    return {k: w[k] for k in POOLED}


def main():
    design_path = sys.argv[1] if len(sys.argv) > 1 else str(config.RESULTS_DIR / "design_phase0.csv")
    group = sys.argv[2] if len(sys.argv) > 2 else "noise"
    hour = float(sys.argv[3]) if len(sys.argv) > 3 else max(config.SNAPSHOT_HOURS)
    design = pd.read_csv(design_path)
    if "group" in design:
        design = design[design.group == group]
    raws = []
    for rid in design.run_id:
        meta_path = config.RUNS_DIR / rid / "meta.json"
        if meta_path.exists() and json.loads(meta_path.read_text()).get("sim_status") == "ok":
            raws.append((rid, model_hook.load_raw(config.RUNS_DIR / rid / json.loads(meta_path.read_text())["snapshots"][0])))
    if not raws:
        raise SystemExit("No finished runs.")
    n15 = len(raws[0][1]["location"])
    n5 = n15 // GROUP
    print(f"{len(raws)} runs ({group}), window {hour:g} h after t*: {n15} vs. {n5} rings")

    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for rid, raw in raws:
            for n_rings, r in ((n15, raw), (n5, merge_rings(raw))):
                rows += [{"run_id": rid, "rings": n_rings, "output": k, "value": v} for k, v in per_run(r, hour).items()]
    df = pd.DataFrame(rows)
    res = []
    for o, s in df.groupby("output", sort=False):
        w = s.pivot(index="run_id", columns="rings", values="value")
        diff = (w[n5] - w[n15]).dropna()
        m = w[n15].mean()
        if len(diff) >= 3 and m:
            rel = diff.mean() / m
            hw = stats.t.ppf(0.975, len(diff) - 1) * diff.std() / np.sqrt(len(diff)) / abs(m)
        else:
            rel = hw = np.nan
        res.append({"kind": "per run", "output": o, f"valid_{n15}": int(w[n15].notna().sum()), f"valid_{n5}": int(w[n5].notna().sum()),
                    f"median_{n15}": w[n15].median(), f"median_{n5}": w[n5].median(),
                    "rel_diff": rel, "rel_ci_lo": rel - hw, "rel_ci_hi": rel + hw})
    p15, p5 = pooled([r for _, r in raws], hour), pooled([merge_rings(r) for _, r in raws], hour)
    for o in p15:
        res.append({"kind": "pooled", "output": o, f"median_{n15}": p15[o], f"median_{n5}": p5[o],
                    "rel_diff": (p5[o] - p15[o]) / p15[o] if p15[o] else np.nan})
    tab = pd.DataFrame(res)
    tag = design_path.rsplit("/", 1)[-1].replace("design_", "").replace(".csv", "")
    path = config.RESULTS_DIR / f"rings_compare_{tag}.csv"
    tab.to_csv(path, index=False)
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        print("\nPer run (SA outputs), paired over runs; pooled over all runs (calibration targets):")
        print(tab.round(3).to_string(index=False))
    print("\nn_* = number of rings behind the fits (fewer rings: MIN_RINGS = 3 is harder to reach).")
    print(f"-> {path}")


if __name__ == "__main__":
    main()
