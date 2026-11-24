"""
Welche Morris-Läufe fehlen, warum, und ist das Fehlen zufällig? Rechnet nichts nach, liest nur
results/design_morris.csv, results/status_morris.csv und zum Vergleich mit dem alten Modell
(Baumausschluss) results_treeexcl/status_morris.csv.

    sbatch sa_run.sbatch analyze_missing.py                     # aus dem Home-Verzeichnis
    sbatch sa_run.sbatch analyze_missing.py results_treeexcl    # anderer Ordner für den alten Lauf

1. Status neu gegen alt: wie viele Läufe fehlen in beiden, wie viele nur im neuen Modell
2. Laufzeit der Simulation neu gegen alt (gleiche run_id, gleiche Faktoren und Seeds)
3. Faktorstufen der fehlenden gegen die fertigen Läufe (u in [0, 1]): sind die Timeouts auf
   bestimmte Faktorwerte konzentriert, fehlen gerade die teuren Ecken -> Morris verzerrt
4. Welche Faktoren treiben die Laufzeit (Spearman zwischen Stufe und log-Laufzeit der fertigen Läufe)
5. Vollständige Trajektorien je Replikat jetzt und wenn alle Timeouts nachgerechnet würden
"""
import sys

import numpy as np
import pandas as pd

import config

OK = "ok"


def load_status(folder):
    path = config.ROOT / folder / "status_morris.csv"
    if not path.exists():
        print(f"(no {path})")
        return None
    return pd.read_csv(path)


def main():
    old_folder = sys.argv[1] if len(sys.argv) > 1 else "results_treeexcl"
    design = pd.read_csv(config.RESULTS_DIR / "design_morris.csv")
    new = load_status("results")
    old = load_status(old_folder)
    names = [f.name for f in config.FACTORS]
    d = design.merge(new.rename(columns={"status": "st_new", "sim_seconds": "t_new"}), on="run_id")
    if old is not None:
        d = d.merge(old.rename(columns={"status": "st_old", "sim_seconds": "t_old"}), on="run_id", how="left")

    print("\n=== 1. Status new (rows) vs. old (columns) ===")
    if old is not None:
        print(pd.crosstab(d["st_new"], d["st_old"], margins=True).to_string())
    else:
        print(d["st_new"].value_counts().to_string())
    for st in ("failed", "missing", "not_evaluated"):
        ids = d.loc[d["st_new"] == st, "run_id"].tolist()
        if 0 < len(ids) <= 10:
            print(f"{st}: {', '.join(ids)}")

    print("\n=== 2. Simulation runtime [h] ===")
    t_new = d.loc[d["st_new"] == OK, "t_new"] / 3600
    print(f"new: n={len(t_new)}  median {t_new.median():.2f}  90% {t_new.quantile(.9):.2f}  max {t_new.max():.2f}")
    if old is not None:
        t_old = d.loc[d["st_old"] == OK, "t_old"] / 3600
        print(f"old: n={len(t_old)}  median {t_old.median():.2f}  90% {t_old.quantile(.9):.2f}  max {t_old.max():.2f}")
        both = d[(d["st_new"] == OK) & (d["st_old"] == OK)]
        r = both["t_new"] / both["t_old"]
        print(f"Ratio new/old (same runs, n={len(both)}): median {r.median():.2f}  "
              f"10% {r.quantile(.1):.2f}  90% {r.quantile(.9):.2f}")
        print("Note: the runtime includes the root phase up to t*, which is about the same in both models")

    print("\n=== 3. Mean factor level u: missing/failed vs. finished ===")
    bad = d["st_new"] != OK
    rows = []
    for n in names:
        u = d[f"u_{n}"]
        rows.append({"factor": n, "u_ok": u[~bad].mean(), "u_missing": u[bad].mean(),
                     "diff": u[bad].mean() - u[~bad].mean(),
                     "missing_at_u>=0.67": (bad & (u >= 2 / 3)).sum(), "missing_at_u<=0.33": (bad & (u <= 1 / 3)).sum()})
    t3 = pd.DataFrame(rows).set_index("factor").sort_values("diff", key=np.abs, ascending=False)
    print(t3.round(2).to_string())
    print(f"(missing: {bad.sum()} runs; |diff| > 0.15 indicates clustering at one end of the range)")

    print("\n=== 4. Runtime drivers: Spearman(u, log runtime) of the finished runs ===")
    ok = d[~bad]
    lt = np.log(ok["t_new"])
    rho = {n: ok[f"u_{n}"].rank().corr(lt.rank()) for n in names}
    print(pd.Series(rho).sort_values(key=np.abs, ascending=False).round(2).to_string())

    print("\n=== 5. Complete trajectories (all points ok) per replicate ===")
    k1 = d.groupby("traj")["step"].nunique().max()
    def complete(ok_mask):
        g = d.assign(ok=ok_mask).groupby(["traj", "rep"])["ok"].sum()
        return (g == k1).groupby("rep").sum()
    now = complete(~bad)
    if_rerun = complete(d["st_new"] != "failed")   # Timeouts/nicht ausgewertet nachgerechnet, failed bleibt
    print(pd.DataFrame({"now": now, "after rerunning": if_rerun}).to_string())
    lost = d[d["st_new"] == "failed"].groupby(["traj", "rep"]).size()
    if len(lost):
        print(f"permanently incomplete because of failed (not_crossed): {len(lost)} (trajectory, replicate) pairs")
    by_traj = d.assign(bad=bad).groupby("traj")["bad"].sum()
    print("missing runs per trajectory:", by_traj[by_traj > 0].to_dict())


if __name__ == "__main__":
    main()
