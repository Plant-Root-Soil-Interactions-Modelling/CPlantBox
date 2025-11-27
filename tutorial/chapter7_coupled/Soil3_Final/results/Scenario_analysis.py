import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = BASE_DIR  # nicht zwingend nötig, aber ok
OUT_DIR = os.path.join(BASE_DIR, "comparison_plots_groups")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------
# Timeseries-NPZs finden
# ---------------------------------------------------------
pattern = os.path.join(BASE_DIR, "*", "*_timeseries_data.npz")
npz_files = sorted(glob.glob(pattern))

if not npz_files:
    print(f"Keine NPZ-Dateien gefunden unter {pattern}")
    # debug help:
    print("BASE_DIR enthält:")
    for x in os.listdir(BASE_DIR):
        print("  ", x)
    raise SystemExit

# ---------------------------------------------------------
# Gruppen-Definition wie in den Simulationsskripten
# ---------------------------------------------------------
groups = [
    (0, 7),  # boundary edge
    (1, 2),  # off-furrow left
    (3, 4),  # on furrow
    (5, 6),  # off-furrow right
]
group_labels = [
    "G1 boundary (0,7)",
    "G2 off-furrow L (1,2)",
    "G3 on furrow (3,4)",
    "G4 off-furrow R (5,6)",
]

# ---------------------------------------------------------
# Daten einlesen
# ---------------------------------------------------------
runs = {}  # run_name -> dict

for fpath in npz_files:
    # Ordnername = Szenario-Name (z.B. Soil3_Kompost_Scenario1)
    run_name = os.path.basename(os.path.dirname(fpath))
    d = np.load(fpath)

    time = d["sim_time_days"]  # shape: (T_soil,)
    y_all = d["transpiration_per_plant_cm3_day"]  # shape: (N_plants, T_soil)
    trl_all = d["trl_per_plant_cm"]              # shape: (N_plants, T_root) oder (T_root,)

    # Robustheit: 1D → [1, T_root]
    if trl_all.ndim == 1:
        trl_all = trl_all[np.newaxis, :]

    # --- TRL-Zeitvektor berechnen ---
    T_soil = time.shape[0]
    T_root = trl_all.shape[1]

    growth_every = T_soil // T_root
    print(f"{run_name}: T_soil={T_soil}, T_root={T_root}, growth_every≈{growth_every}")

    # Zeitvektor für TRL: jeder growth_every-te Zeitschritt
    trl_time = time[::growth_every][:T_root]

    group_transp = {}
    group_trl = {}

    for (g, label) in zip(groups, group_labels):
        g_idx = list(g)

        # Transpiration: Mittel pro Gruppe (2 Pflanzen), volle Auflösung
        g_T = np.mean([y_all[pidx] for pidx in g_idx], axis=0)  # shape: (T_soil,)

        # TRL: Mittel der beiden Pflanzen, nur T_root Schritte
        if trl_all.shape[0] >= max(g_idx) + 1:
            g_trl = np.mean([trl_all[pidx] for pidx in g_idx], axis=0)  # shape: (T_root,)
        else:
            g_trl = trl_all[0]

        group_transp[label] = g_T
        group_trl[label] = g_trl

    runs[run_name] = {
        "time": time,          # für Transpiration
        "time_trl": trl_time,  # Zeitvektor für TRL
        "group_transp": group_transp,
        "group_trl": group_trl,
    }

# ---------------------------------------------------------
# Farben für Szenarien / Baseline
# ---------------------------------------------------------
# Alle Runs, die "Scenario" im Namen tragen
scenario_names = sorted([rn for rn in runs.keys() if "Scenario" in rn])
# Colormap über Anzahl Szenarien
if scenario_names:
    cmap = cm.get_cmap("viridis", len(scenario_names))
    scenario_colors = {run: cmap(i) for i, run in enumerate(scenario_names)}
else:
    scenario_colors = {}

# Baseline-Farbe (und alles, was nicht "Scenario" enthält)
baseline_color = "black"

# ---------------------------------------------------------
# Plots: pro Gruppe ein Vergleich über alle Runs
# ---------------------------------------------------------

# 1) Transpiration-Vergleich (negative Werte nach oben → -gT)
for label in group_labels:
    plt.figure(figsize=(8, 5))
    for run_name, rdata in sorted(runs.items()):
        t = rdata["time"]
        gT = rdata["group_transp"][label]

        # Farbe wählen: Szenario aus Colormap, sonst Baseline-Farbe
        color = scenario_colors.get(run_name, baseline_color)

        plt.plot(t, -gT, label=run_name, color=color)

    plt.xlabel("Time [days]")
    plt.ylabel("Transpiration per plant group [cm³/day]")
    plt.title(f"Transpiration – {label}")
    plt.grid(True, linestyle=":", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    fname = f"comparison_transp_{label.split()[0]}.png"
    plt.savefig(os.path.join(OUT_DIR, fname), dpi=300)
    plt.close()

# 2) TRL-Vergleich
for label in group_labels:
    plt.figure(figsize=(8, 5))
    for run_name, rdata in sorted(runs.items()):
        t_trl = rdata["time_trl"]
        gL = rdata["group_trl"][label]

        color = scenario_colors.get(run_name, baseline_color)

        plt.plot(t_trl, gL, label=run_name, color=color)

    plt.xlabel("Time [days]")
    plt.ylabel("Total root length per group [cm]")
    plt.title(f"Total root length – {label}")
    plt.grid(True, linestyle=":", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    fname = f"comparison_trl_{label.split()[0]}.png"
    plt.savefig(os.path.join(OUT_DIR, fname), dpi=300)
    plt.close()

print(f"Fertig. Plots liegen in: {OUT_DIR}")
