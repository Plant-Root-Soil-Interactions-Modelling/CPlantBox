"""
Smoke test: a few real runs end to end before Phase 0 goes to the cluster.
Does not check whether the results are right, but whether everything works together.

    python smoke_test.py nominal            # one point: nominal | low (all factors low) | high (all high)
    python smoke_test.py nominal --quick    # only 12 h after t*: checks the wiring quickly, wave fits mostly empty
    python smoke_test.py --report           # summary of all smoke runs -> results/smoke_report.csv
    python smoke_test.py --report --gate    # plus exit code 1 if the nominal run is not clean (submit_all.sh)

On the cluster from the home directory: bash submit_smoke.sh (all three points in parallel,
then the report as a separate job).

Checked per point: the parameter file is written, makeSimulation runs with it and reaches t*,
raw.npz has all quantities in the right shape, evaluate.py delivers all outputs (which are NaN),
runtime and memory. The runs are in runs/smoke_<point> and are recomputed on every call.
"""
import argparse
import json
import resource
import shutil
import time

import numpy as np
import pandas as pd

import config
import evaluate
import model_hook
from run_one import run_row
from xml_tools import nominal_values

POINTS = ["nominal", "low", "high"]
RAW_2D = ["tip_dens", "active_tips", "ana_new", "length_densities", "length_runner", "length_bas",
          "active_tips_runner", "active_tips_bas", "ana_new_runner", "ana_new_bas", "branch_new", "branch_new_runner"]


def point_values(point):
    if point == "nominal":
        return nominal_values()
    return {f.name: (f.low if point == "low" else f.high) for f in config.FACTORS}


def check_raw(raw, problems):
    n_rings, n_cols = np.shape(raw["length_densities"])
    if n_rings != config.N_RINGS:
        problems.append(f"{n_rings} rings instead of {config.N_RINGS}")
    for k in RAW_2D:
        if k not in raw:
            problems.append(f"raw.npz: {k} missing")
        elif np.shape(raw[k]) != (n_rings, n_cols):
            problems.append(f"raw.npz: {k} has shape {np.shape(raw[k])} instead of {(n_rings, n_cols)}")
    if len(raw["times"]) != n_cols or np.any(np.diff(raw["times"]) <= 0):
        problems.append("raw.npz: times does not match the columns or is not increasing")
    hours = raw["times"][-1] * 24
    if hours < max(config.SNAPSHOT_HOURS) - 1:
        problems.append(f"only {hours:.1f} h after t* simulated, last evaluation time {max(config.SNAPSHOT_HOURS)} h")
    if raw["length_densities"][:, -1].sum() <= 0:
        problems.append("no hyphal length behind the barrier")
    return n_cols, hours


def check_extras(raw, warnings):
    """Outputs for the calibration: errors here are warnings, not an abort (the SA does not need them)."""
    if config.ROOT_EFFICIENCY:
        if "root_eff_error" in raw:
            warnings.append("root efficiency failed: " + str(raw["root_eff_error"]).strip().splitlines()[-1])
        elif not np.isfinite(float(raw.get("root_eff_E_rel", np.nan))):
            warnings.append("root efficiency NaN (no tips or root nodes)")
    if "branch_new" in raw and raw["branch_new"].sum() <= 0:
        warnings.append("no branching events behind the barrier")


def run_point(point, quick):
    if quick:
        config.HOURS_HYPHAE = 12
        config.SNAPSHOT_HOURS = [6, 12]
    run_id = f"smoke_{point}" + ("_quick" if quick else "")
    rdir = config.RUNS_DIR / run_id
    shutil.rmtree(rdir, ignore_errors=True)
    row = pd.Series({"run_id": run_id, "seed": config.BASE_SEED, "dx": np.nan, "factor_set": config.FACTOR_SET,
                     **point_values(point)})
    t0 = time.time()
    run_row(row)
    meta = json.loads((rdir / "meta.json").read_text())
    report = {"run_id": run_id, "point": point, "quick": quick, "sim_status": meta["sim_status"],
              "sim_h": meta["sim_seconds"] / 3600, "problems": [], "warnings": []}
    if meta["sim_status"] != "ok":
        report["problems"].append("simulation failed: " + meta.get("error", "").strip().splitlines()[-1])
        print(meta.get("error", ""))
    else:
        raw = model_hook.load_raw(rdir / meta["snapshots"][0])
        report["wait_days"] = float(raw["wait_days"])
        report["n_cols"], report["hours_after_tstar"] = check_raw(raw, report["problems"])
        check_extras(raw, report["warnings"])
        report["root_eff_E_rel"] = float(raw.get("root_eff_E_rel", np.nan))
        evaluate.evaluate_run(run_id, force=True)
        out = pd.read_csv(rdir / "outputs.csv")
        sa = out[out.output.isin(config.SNAPSHOT_OUTPUTS + config.WAVE_OUTPUTS)]
        expected = {(o, t) for o in config.SNAPSHOT_OUTPUTS for t in config.SNAPSHOT_HOURS} \
            | {(o, max(config.SNAPSHOT_HOURS)) for o in config.WAVE_OUTPUTS}
        missing = expected - set(zip(sa.output, sa.time))
        if missing:
            report["problems"].append(f"outputs missing: {sorted(missing)}")
        nan = sa[sa.value.isna()]
        report["nan_outputs"] = ", ".join(f"{o}@{t:g}h" for o, t in zip(nan.output, nan.time))
        print("\nOutputs:")
        print(out.pivot_table(index="output", columns="time", values="value", aggfunc="first").to_string())
    report["total_h"] = (time.time() - t0) / 3600
    report["peak_mem_gb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2   # ru_maxrss in kB (Linux)
    report["problems"] = "; ".join(report["problems"])
    report["warnings"] = "; ".join(report["warnings"])
    (rdir / "smoke.json").write_text(json.dumps(report, indent=1))
    print("\nSmoke test " + run_id + ": " + ("OK" if not report["problems"] else "PROBLEMS: " + report["problems"]))
    if report["warnings"]:
        print("WARNING: " + report["warnings"])
    if report.get("nan_outputs"):
        print("NaN (fit failed or too few rings): " + report["nan_outputs"])
    print(f"Simulation {report['sim_h']:.2f} h, total {report['total_h']:.2f} h, peak memory {report['peak_mem_gb']:.2f} GB")


def summary(gate=False, quick=False):
    rows = [json.loads(p.read_text()) for p in sorted(config.RUNS_DIR.glob("smoke_*/smoke.json"))]
    if not rows:
        print("no smoke runs yet" + ("\nGATE FAILED: no smoke runs (aborted?) -> slurmOut/sa_smoke-*.out" if gate else ""))
        return 1 if gate else 0
    df = pd.DataFrame(rows)
    config.RESULTS_DIR.mkdir(exist_ok=True)
    df.to_csv(config.RESULTS_DIR / "smoke_report.csv", index=False)
    cols = [c for c in ["run_id", "sim_status", "wait_days", "hours_after_tstar", "sim_h", "peak_mem_gb", "root_eff_E_rel",
                        "nan_outputs", "problems", "warnings"] if c in df]
    with pd.option_context("display.max_colwidth", 80, "display.width", 200):
        print(df[cols].round(2).to_string(index=False))
    print(f"\n-> {config.RESULTS_DIR / 'smoke_report.csv'}")
    print("For sim.sbatch: --time above the longest simulation (with margin), --mem above the largest memory.")
    if not gate:
        return 0
    # Gate: the nominal run of this smoke test must be clean. low/high may fail (edge of the factor space,
    # e.g. t* not reached): that is a statement about the ranges, not a wiring error -> note only
    rid = "smoke_nominal" + ("_quick" if quick else "")
    nom = df[df.run_id == rid]
    if nom.empty:
        print(f"\nGATE FAILED: {rid} missing (run aborted, time limit?) -> slurmOut/sa_smoke-*_0.out")
        return 1
    r = nom.iloc[0]
    if r.sim_status != "ok" or r.problems:
        print(f"\nGATE FAILED: {rid}: {r.problems or r.sim_status}. The following jobs are cancelled.")
        return 1
    edge = df[df.run_id.isin([p + ("_quick" if quick else "") for p in ("smoke_low", "smoke_high")]) & (df.problems != "")]
    for _, e in edge.iterrows():
        print(f"Note: {e.run_id}: {e.problems} (edge of the factor space, check; not aborting)")
    print(f"\nGATE OK: {rid} clean, the following jobs start.")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Smoke test: a few real runs end to end before Phase 0 goes to the cluster.")
    ap.add_argument("point", nargs="?", choices=POINTS)
    ap.add_argument("--quick", action="store_true", help="only 12 h after t* (check the wiring)")
    ap.add_argument("--report", action="store_true", help="summary of all smoke runs")
    ap.add_argument("--gate", action="store_true", help="with --report: exit code 1 if the nominal run is not clean")
    args = ap.parse_args()
    if args.report:
        raise SystemExit(summary(args.gate, args.quick))
    elif args.point:
        run_point(args.point, args.quick)
    else:
        ap.error("give a point (nominal, low, high) or --report")
