"""
Test mit längerer Simulation (design.py long): Bildet sich die Welle mit 120 h statt 60 h so weit aus, dass die
Wellen-Fits in fast jedem Lauf definiert sind? Für Morris muss ein Output in allen k+1 Punkten einer Trajektorie
vorliegen, also in fast allen Läufen.

    python analyze_long.py      -> results/long_summary.csv
"""
import numpy as np
import pandas as pd

import config

SHOW = ["front_velocity", "v_length", "v_length_sat", "K1", "K2", "wave_width", "encounter_C",
        "n_wave_rings", "n_K1_rings", "n_tip_rings", "front_position", "total_length", "active_tips"]


def main():
    status = pd.read_csv(config.RESULTS_DIR / "status_long.csv")
    out = pd.read_csv(config.RESULTS_DIR / "outputs_long.csv")
    n = len(status)
    ok = int((status.status == "ok").sum())
    t = status.sim_seconds.dropna() / 3600
    print(f"Läufe: {ok} von {n} ok, Simulation Median {t.median():.2f} h, Max {t.max():.2f} h")

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
        print("\nJe Output und Zeitfenster (h nach t*): gültige Läufe, Median, CV, "
              f"p_traj = Anteil vollständiger Trajektorien ({k1} Punkte):")
        print(wide.round(3).to_string())
    print(f"\nRichtwert: Tipgeschwindigkeit {0.576 * 10:.1f} mm/d. Für Morris taugt ein Output mit p_traj nahe 1.")
    print(f"-> {config.RESULTS_DIR / 'long_summary.csv'}")


if __name__ == "__main__":
    main()
