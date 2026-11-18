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
        print(f"(kein {path})")
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

    print("\n=== 1. Status neu (Zeilen) gegen alt (Spalten) ===")
    if old is not None:
        print(pd.crosstab(d["st_new"], d["st_old"], margins=True).to_string())
    else:
        print(d["st_new"].value_counts().to_string())
    for st in ("failed", "missing", "not_evaluated"):
        ids = d.loc[d["st_new"] == st, "run_id"].tolist()
        if 0 < len(ids) <= 10:
            print(f"{st}: {', '.join(ids)}")

    print("\n=== 2. Laufzeit Simulation [h] ===")
    t_new = d.loc[d["st_new"] == OK, "t_new"] / 3600
    print(f"neu: n={len(t_new)}  Median {t_new.median():.2f}  90% {t_new.quantile(.9):.2f}  Max {t_new.max():.2f}")
    if old is not None:
        t_old = d.loc[d["st_old"] == OK, "t_old"] / 3600
        print(f"alt: n={len(t_old)}  Median {t_old.median():.2f}  90% {t_old.quantile(.9):.2f}  Max {t_old.max():.2f}")
        both = d[(d["st_new"] == OK) & (d["st_old"] == OK)]
        r = both["t_new"] / both["t_old"]
        print(f"Verhältnis neu/alt (gleiche Läufe, n={len(both)}): Median {r.median():.2f}  "
              f"10% {r.quantile(.1):.2f}  90% {r.quantile(.9):.2f}")
        print("Hinweis: Laufzeit enthält die Wurzelphase bis t*, die in beiden Modellen etwa gleich ist")

    print("\n=== 3. Mittlere Faktorstufe u: fehlend/fehlgeschlagen gegen fertig ===")
    bad = d["st_new"] != OK
    rows = []
    for n in names:
        u = d[f"u_{n}"]
        rows.append({"factor": n, "u_ok": u[~bad].mean(), "u_fehlend": u[bad].mean(),
                     "diff": u[bad].mean() - u[~bad].mean(),
                     "fehlend_bei_u>=0.67": (bad & (u >= 2 / 3)).sum(), "fehlend_bei_u<=0.33": (bad & (u <= 1 / 3)).sum()})
    t3 = pd.DataFrame(rows).set_index("factor").sort_values("diff", key=np.abs, ascending=False)
    print(t3.round(2).to_string())
    print(f"(fehlend: {bad.sum()} Läufe; |diff| > 0.15 deutet auf Häufung an einem Ende des Bereichs)")

    print("\n=== 4. Laufzeittreiber: Spearman(u, log Laufzeit) der fertigen Läufe ===")
    ok = d[~bad]
    lt = np.log(ok["t_new"])
    rho = {n: ok[f"u_{n}"].rank().corr(lt.rank()) for n in names}
    print(pd.Series(rho).sort_values(key=np.abs, ascending=False).round(2).to_string())

    print("\n=== 5. Vollständige Trajektorien (alle Punkte ok) je Replikat ===")
    k1 = d.groupby("traj")["step"].nunique().max()
    def complete(ok_mask):
        g = d.assign(ok=ok_mask).groupby(["traj", "rep"])["ok"].sum()
        return (g == k1).groupby("rep").sum()
    now = complete(~bad)
    if_rerun = complete(d["st_new"] != "failed")   # Timeouts/nicht ausgewertet nachgerechnet, failed bleibt
    print(pd.DataFrame({"jetzt": now, "nach Nachrechnen": if_rerun}).to_string())
    lost = d[d["st_new"] == "failed"].groupby(["traj", "rep"]).size()
    if len(lost):
        print(f"durch failed (not_crossed) dauerhaft unvollständig: {len(lost)} (Trajektorie, Replikat)-Paare")
    by_traj = d.assign(bad=bad).groupby("traj")["bad"].sum()
    print("fehlende Läufe je Trajektorie:", by_traj[by_traj > 0].to_dict())


if __name__ == "__main__":
    main()
