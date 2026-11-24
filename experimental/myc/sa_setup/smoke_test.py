"""
Smoke-Test: wenige echte Läufe von Anfang bis Ende, bevor Phase 0 auf den Cluster geht.
Prüft nicht, ob die Ergebnisse stimmen, sondern ob alles zusammenspielt.

    python smoke_test.py nominal            # ein Punkt: nominal | low (alle Faktoren unten) | high (alle oben)
    python smoke_test.py nominal --quick    # nur 12 h nach t*: prüft die Verkabelung schnell, Wellen-Fits meist leer
    python smoke_test.py --report           # Zusammenfassung aller Smoke-Läufe -> results/smoke_report.csv
    python smoke_test.py --report --gate    # dazu Exit-Code 1, wenn der Nominal-Lauf nicht sauber ist (submit_all.sh)

Auf dem Cluster aus dem Home-Verzeichnis: bash submit_smoke.sh (alle drei Punkte parallel,
danach der Bericht als eigener Job).

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
          "active_tips_runner", "active_tips_bas", "ana_new_runner", "ana_new_bas", "branch_new"]


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


def check_extras(raw, warnings):
    """Outputs für die Kalibrierung: Fehler hier sind Warnungen, kein Abbruch (die SA braucht sie nicht)."""
    if config.ROOT_EFFICIENCY:
        if "root_eff_error" in raw:
            warnings.append("Root Efficiency fehlgeschlagen: " + str(raw["root_eff_error"]).strip().splitlines()[-1])
        elif not np.isfinite(float(raw.get("root_eff_E_rel", np.nan))):
            warnings.append("Root Efficiency NaN (keine Spitzen oder Wurzelknoten)")
    if "branch_new" in raw and raw["branch_new"].sum() <= 0:
        warnings.append("keine Verzweigungsereignisse hinter der Barriere")


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
              "sim_h": meta["sim_seconds"] / 3600, "problems": [], "warnings": []}
    if meta["sim_status"] != "ok":
        report["problems"].append("Simulation fehlgeschlagen: " + meta.get("error", "").strip().splitlines()[-1])
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
            report["problems"].append(f"Outputs fehlen: {sorted(missing)}")
        nan = sa[sa.value.isna()]
        report["nan_outputs"] = ", ".join(f"{o}@{t:g}h" for o, t in zip(nan.output, nan.time))
        print("\nOutputs:")
        print(out.pivot_table(index="output", columns="time", values="value", aggfunc="first").to_string())
    report["total_h"] = (time.time() - t0) / 3600
    report["peak_mem_gb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2   # ru_maxrss in kB (Linux)
    report["problems"] = "; ".join(report["problems"])
    report["warnings"] = "; ".join(report["warnings"])
    (rdir / "smoke.json").write_text(json.dumps(report, indent=1))
    print("\nSmoke-Test " + run_id + ": " + ("OK" if not report["problems"] else "PROBLEME: " + report["problems"]))
    if report["warnings"]:
        print("WARNUNG: " + report["warnings"])
    if report.get("nan_outputs"):
        print("NaN (Fit gescheitert oder zu wenige Ringe): " + report["nan_outputs"])
    print(f"Simulation {report['sim_h']:.2f} h, insgesamt {report['total_h']:.2f} h, Speicher max. {report['peak_mem_gb']:.2f} GB")


def summary(gate=False, quick=False):
    rows = [json.loads(p.read_text()) for p in sorted(config.RUNS_DIR.glob("smoke_*/smoke.json"))]
    if not rows:
        print("noch keine Smoke-Läufe" + ("\nGATE FEHLGESCHLAGEN: keine Smoke-Läufe (abgebrochen?) -> slurmOut/sa_smoke-*.out" if gate else ""))
        return 1 if gate else 0
    df = pd.DataFrame(rows)
    config.RESULTS_DIR.mkdir(exist_ok=True)
    df.to_csv(config.RESULTS_DIR / "smoke_report.csv", index=False)
    cols = [c for c in ["run_id", "sim_status", "wait_days", "hours_after_tstar", "sim_h", "peak_mem_gb", "root_eff_E_rel",
                        "nan_outputs", "problems", "warnings"] if c in df]
    with pd.option_context("display.max_colwidth", 80, "display.width", 200):
        print(df[cols].round(2).to_string(index=False))
    print(f"\n-> {config.RESULTS_DIR / 'smoke_report.csv'}")
    print("Für sim.sbatch: --time über der längsten Simulation (mit Puffer), --mem über dem größten Speicher.")
    if not gate:
        return 0
    # Gate: der Nominal-Lauf dieses Smoke-Tests muss sauber sein. low/high dürfen scheitern (Rand des Faktorraums,
    # z. B. t* nicht erreicht): das ist eine Aussage über die Bereiche, kein Verkabelungsfehler -> nur Hinweis
    rid = "smoke_nominal" + ("_quick" if quick else "")
    nom = df[df.run_id == rid]
    if nom.empty:
        print(f"\nGATE FEHLGESCHLAGEN: {rid} fehlt (Lauf abgebrochen, Zeitlimit?) -> slurmOut/sa_smoke-*_0.out")
        return 1
    r = nom.iloc[0]
    if r.sim_status != "ok" or r.problems:
        print(f"\nGATE FEHLGESCHLAGEN: {rid}: {r.problems or r.sim_status}. Folgende Jobs werden abgebrochen.")
        return 1
    edge = df[df.run_id.isin([p + ("_quick" if quick else "") for p in ("smoke_low", "smoke_high")]) & (df.problems != "")]
    for _, e in edge.iterrows():
        print(f"Hinweis: {e.run_id}: {e.problems} (Rand des Faktorraums, prüfen; kein Abbruch)")
    print(f"\nGATE OK: {rid} sauber, die folgenden Jobs starten.")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("point", nargs="?", choices=POINTS)
    ap.add_argument("--quick", action="store_true", help="nur 12 h nach t* (Verkabelung prüfen)")
    ap.add_argument("--report", action="store_true", help="Zusammenfassung aller Smoke-Läufe")
    ap.add_argument("--gate", action="store_true", help="mit --report: Exit-Code 1, wenn der Nominal-Lauf nicht sauber ist")
    args = ap.parse_args()
    if args.report:
        raise SystemExit(summary(args.gate, args.quick))
    elif args.point:
        run_point(args.point, args.quick)
    else:
        ap.error("Punkt (nominal, low, high) oder --report angeben")
