"""
Phase 0: Seed-Rauschen und dx-Konvergenz am Nominalpunkt.

    python analyze_phase0.py
      -> results/phase0_noise.csv   Mittel, SD, CV, Anteil Nullen je Output
                                    + Rauschniveau der Elementary Effects für R = 1..5
      -> results/phase0_dx.csv      relative Abweichung je dx gegenüber feinstem dx
"""
import numpy as np
import pandas as pd
from scipy import stats

import config

KEY = ["time", "region", "output"]


def main():
    design = pd.read_csv(config.RESULTS_DIR / "design_phase0.csv")
    out = pd.read_csv(config.RESULTS_DIR / "outputs_phase0.csv").merge(
        design[["run_id", "group", "seed"]], on="run_id")

    # ---------- Rauschen ----------
    noise = out[out.group == "noise"]
    g = noise.groupby(KEY)["value"]
    tab = pd.DataFrame({"n": g.size(), "mean": g.mean(), "sd": g.std(),
                        "frac_zero": g.apply(lambda v: np.mean(v == 0))})
    tab["cv"] = tab["sd"] / tab["mean"].replace(0, np.nan)
    # Rauschniveau eines Elementary Effects (normierte Einheiten), gemittelt über R Replikate,
    # konservativ ohne Gewinn durch gemeinsame Seeds: sqrt(2)*sd / (Δ*sqrt(R)).
    p = config.MORRIS_LEVELS
    delta = p / (2 * (p - 1))
    for R in range(1, 6):
        tab[f"EE_noise_R{R}"] = np.sqrt(2) * tab["sd"] / (delta * np.sqrt(R))
    tab.reset_index().to_csv(config.RESULTS_DIR / "phase0_noise.csv", index=False)

    # ---------- dx-Konvergenz (gepaart über Seeds) ----------
    dxo = out[out.group.str.startswith("dx=")].copy()
    dxo["dx"] = dxo.group.str[3:].astype(float)
    finest = dxo.dx.min()
    rows = []
    for key, sub in dxo.groupby(KEY):
        ref = sub[sub.dx == finest].set_index("seed")["value"]
        for dx, s in sub.groupby("dx"):
            if dx == finest:
                continue
            v = s.set_index("seed")["value"]
            diff = (v - ref).dropna()
            m = ref.mean()
            rel = diff.mean() / m if m else np.nan
            if len(diff) > 1 and diff.std() > 0:
                hw = stats.t.ppf(0.975, len(diff) - 1) * diff.std() / np.sqrt(len(diff)) / abs(m) if m else np.nan
            else:
                hw = 0.0
            rows.append(dict(zip(KEY, key), dx=dx, ref_dx=finest, rel_diff=rel,
                             rel_ci_lo=rel - hw, rel_ci_hi=rel + hw))
    dxt = pd.DataFrame(rows)
    dxt.to_csv(config.RESULTS_DIR / "phase0_dx.csv", index=False)

    # ---------- Kurzbericht ----------
    t_end = max(config.SNAPSHOT_TIMES)
    print(f"\nSeed-Rauschen am Nominalpunkt (t = {t_end}, gesamter Zielbereich):")
    print(tab.loc[(t_end, "all")][["mean", "sd", "cv", "frac_zero"]].round(3).to_string())
    nom_dx = max(config.PHASE0_DX_VALUES)
    d = dxt[(dxt.time == t_end) & (dxt.region == "all") & (dxt.dx == nom_dx)]
    print(f"\ndx-Konvergenz: dx = {nom_dx:g} gegenüber dx = {finest:g} (rel. Abweichung, 95%-KI):")
    for _, r in d.iterrows():
        flag = "  <-- prüfen" if (r.rel_ci_lo > 0.02 or r.rel_ci_hi < -0.02) else ""
        print(f"  {r.output:16s} {r.rel_diff:+.3f} [{r.rel_ci_lo:+.3f}, {r.rel_ci_hi:+.3f}]{flag}")
    sparse = tab[tab.frac_zero > 0.2].reset_index()
    if len(sparse):
        print(f"\n{len(sparse)} Output/Region/Zeit-Kombinationen mit >20 % Nullen "
              "(spärlich besetzte Teilbereiche, s. phase0_noise.csv)")


if __name__ == "__main__":
    main()
