"""
Smoke-Test: wenige echte Läufe von Anfang bis Ende, bevor Phase 0 auf den Cluster geht.
Prüft nicht, ob die Ergebnisse stimmen, sondern ob alles zusammenspielt.

    python smoke_test.py nominal            # ein Punkt: nominal | low (alle Faktoren unten) | high (alle oben)
    python smoke_test.py nominal --quick    # nur 12 h nach t*: prüft die Verkabelung schnell, Wellen-Fits meist leer
    python smoke_test.py --report           # Zusammenfassung aller Smoke-Läufe -> results/smoke_report.csv

Auf dem Cluster: slurm/smoke.sbatch ins Home-Verzeichnis kopieren, dort sbatch smoke.sbatch
(alle drei Punkte parallel), danach hier --report.

Geprüft wird je Punkt: Parameterfile wird geschrieben, makeSimulation läuft damit und erreicht t*,
raw.npz hat alle Größen in passender Form, evaluate.py liefert alle Outputs (welche NaN sind),
Laufzeit und Speicher. Die Läufe liegen in runs/smoke_<punkt> und werden bei jedem Aufruf neu gerechnet.
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
          "active_tips_runner", "active_tips_bas", "ana_new_runner", "ana_new_bas"]


def point_values(point):
    if point == "nominal":
        return nominal_values()
    return {f.name: (f.low if point == "low" else f.high) for f in config.FACTORS}


def check_raw(raw, problems):
    n_rings, n_cols = np.shape(raw["length_densities"])
    if n_rings != config.N_RINGS:
        problems.append(f"{n_rings} Ringe statt {config.N_RINGS}")
    for k in RAW_2D:
        if k not in raw:
            problems.append(f"raw.npz: {k} fehlt")
        elif np.shape(raw[k]) != (n_rings, n_cols):
            problems.append(f"raw.npz: {k} hat Form {np.shape(raw[k])} statt {(n_rings, n_cols)}")
    if len(raw["times"]) != n_cols or np.any(np.diff(raw["times"]) <= 0):
        problems.append("raw.npz: times passt nicht zu den Spalten oder ist nicht aufsteigend")
    hours = raw["times"][-1] * 24
    if hours < max(config.SNAPSHOT_HOURS) - 1:
        problems.append(f"nur {hours:.1f} h nach t* simuliert, letzter Auswertezeitpunkt {max(config.SNAPSHOT_HOURS)} h")
    if raw["length_densities"][:, -1].sum() <= 0:
        problems.append("keine Hyphenlänge hinter der Barriere")
    return n_cols, hours


def run_point(point, quick):
    if quick:
        config.HOURS_HYPHAE = 12
        config.SNAPSHOT_HOURS = [6, 12]
    run_id = f"smoke_{point}" + ("_quick" if quick else "")
    rdir = config.RUNS_DIR / run_id
    shutil.rmtree(rdir, ignore_errors=True)
    row = pd.Series({"run_id": run_id, "seed": config.BASE_SEED, "dx": np.nan, **point_values(point)})
    t0 = time.time()
    run_row(row)
    meta = json.loads((rdir / "meta.json").read_text())
    report = {"run_id": run_id, "point": point, "quick": quick, "sim_status": meta["sim_status"],
              "sim_h": meta["sim_seconds"] / 3600, "problems": []}
    if meta["sim_status"] != "ok":
        report["problems"].append("Simulation fehlgeschlagen: " + meta.get("error", "").strip().splitlines()[-1])
        print(meta.get("error", ""))
    else:
        raw = model_hook.load_raw(rdir / meta["snapshots"][0])
        report["wait_days"] = float(raw["wait_days"])
        report["n_cols"], report["hours_after_tstar"] = check_raw(raw, report["problems"])
        evaluate.evaluate_run(run_id, force=True)
        out = pd.read_csv(rdir / "outputs.csv")
        sa = out[out.output.isin(config.SNAPSHOT_OUTPUTS + config.WAVE_OUTPUTS)]
        expected = {(o, t) for o in config.SNAPSHOT_OUTPUTS for t in config.SNAPSHOT_HOURS} \
            | {(o, max(config.SNAPSHOT_HOURS)) for o in config.WAVE_OUTPUTS}
        missing = expected - set(zip(sa.output, sa.time))
        if missing:
            report["problems"].append(f"Outputs fehlen: {sorted(missing)}")
        nan = sa[sa.value.isna()]
        report["nan_outputs"] = ", ".join(f"{o}@{t:g}h" for o, t in zip(nan.output, nan.time))
        print("\nOutputs:")
        print(out.pivot_table(index="output", columns="time", values="value", aggfunc="first").to_string())
    report["total_h"] = (time.time() - t0) / 3600
    report["peak_mem_gb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2   # ru_maxrss in kB (Linux)
    report["problems"] = "; ".join(report["problems"])
    (rdir / "smoke.json").write_text(json.dumps(report, indent=1))
    print("\nSmoke-Test " + run_id + ": " + ("OK" if not report["problems"] else "PROBLEME: " + report["problems"]))
    if report.get("nan_outputs"):
        print("NaN (Fit gescheitert oder zu wenige Ringe): " + report["nan_outputs"])
    print(f"Simulation {report['sim_h']:.2f} h, insgesamt {report['total_h']:.2f} h, Speicher max. {report['peak_mem_gb']:.2f} GB")


def summary():
    rows = [json.loads(p.read_text()) for p in sorted(config.RUNS_DIR.glob("smoke_*/smoke.json"))]
    if not rows:
        print("noch keine Smoke-Läufe"); return
    df = pd.DataFrame(rows)
    config.RESULTS_DIR.mkdir(exist_ok=True)
    df.to_csv(config.RESULTS_DIR / "smoke_report.csv", index=False)
    cols = [c for c in ["run_id", "sim_status", "wait_days", "hours_after_tstar", "sim_h", "peak_mem_gb", "nan_outputs", "problems"] if c in df]
    with pd.option_context("display.max_colwidth", 80, "display.width", 200):
        print(df[cols].round(2).to_string(index=False))
    print(f"\n-> {config.RESULTS_DIR / 'smoke_report.csv'}")
    print("Für sim.sbatch: --time über der längsten Simulation (mit Puffer), --mem über dem größten Speicher.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("point", nargs="?", choices=POINTS)
    ap.add_argument("--quick", action="store_true", help="nur 12 h nach t* (Verkabelung prüfen)")
    ap.add_argument("--report", action="store_true", help="Zusammenfassung aller Smoke-Läufe")
    args = ap.parse_args()
    if args.report:
        summary()
    elif args.point:
        run_point(args.point, args.quick)
    else:
        ap.error("Punkt (nominal, low, high) oder --report angeben")
