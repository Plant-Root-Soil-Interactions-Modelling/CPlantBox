"""
Test mit längerer Simulation (design.py long): Bildet sich die Welle mit 120 h statt 60 h so weit aus, dass die
Wellen-Fits in fast jedem Lauf definiert sind? Für Morris muss ein Output in allen k+1 Punkten einer Trajektorie
vorliegen, also in fast allen Läufen.

    python analyze_long.py      -> results/long_summary.csv
"""
import numpy as np
import pandas as pd

import config

SHOW = ["front_speed", "spread_r90", "max_tip_density", "front_velocity", "v_length", "v_length_sat", "K1", "K2", "wave_width", "encounter_C", "beta",
        "n_wave_rings", "n_K1_rings", "n_tip_rings", "front_position", "total_length", "active_tips"]


def main():
    status = pd.read_csv(config.RESULTS_DIR / "status_long.csv")
    out = pd.read_csv(config.RESULTS_DIR / "outputs_long.csv")
    n = len(status)
    ok = int((status.status == "ok").sum())
    t = status.sim_seconds.dropna() / 3600
    print(f"Runs: {ok} of {n} ok, simulation median {t.median():.2f} h, max {t.max():.2f} h")

    k1 = len(config.FACTORS) + 1
    rows = []
    for (o, h), s in out[out.output.isin(SHOW)].groupby(["output", "time"]):
        v = s.value
        p = v.notna().sum() / n
        rows.append({"output": o, "time": h, "valid": f"{v.notna().sum()}/{n}", "median": v.median(),
                     "cv": v.std() / v.mean() if v.mean() else np.nan,
                     # Anteil vollständiger Morris-Replikate, wenn jeder Punkt mit Wahrscheinlichkeit p einen Wert hat
                     "p_traj": p ** k1})
    tab = pd.DataFrame(rows)
    tab.to_csv(config.RESULTS_DIR / "long_summary.csv", index=False)
    wide = tab.pivot(index="output", columns="time", values=["valid", "median", "cv", "p_traj"])
    wide = wide.reindex([o for o in SHOW if o in wide.index])
    with pd.option_context("display.width", 220, "display.max_columns", 50):
        print("\nPer output and time window (h after t*): valid runs, median, CV, "
              f"p_traj = share of complete trajectories ({k1} points):")
        print(wide.round(3).to_string())
    print(f"\nReference: tip speed {0.576 * 10:.1f} mm/d. An output is usable for Morris if p_traj is close to 1.")
    print(f"-> {config.RESULTS_DIR / 'long_summary.csv'}")


if __name__ == "__main__":
    main()
