#!/usr/bin/env python3
"""
Standalone script to recreate Soil3_Kompost_Baseline_uptake_depth_groups.png
when placed in the same folder as the *_timeseries_data.npz file.
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt


def main():
    # --------------------------------------------------------
    # Find the NPZ in the same directory
    # --------------------------------------------------------
    npz_files = glob.glob("*_timeseries_data.npz")
    if not npz_files:
        raise FileNotFoundError("No '*_timeseries_data.npz' file found in this folder.")
    if len(npz_files) > 1:
        print("Warning: multiple *_timeseries_data.npz found, using the first one.")
    npz_path = npz_files[0]
    print(f"Using NPZ: {npz_path}")

    # --------------------------------------------------------
    # Load data
    # --------------------------------------------------------
    data = np.load(npz_path)

    cum_uptake = data["cum_uptake_per_plant_per_layer_cm3"]
    depth_centers_cm = data["layer_depth_centers_cm"]
    transpiration = data["transpiration_per_plant_cm3_day"]

    # --------------------------------------------------------
    # Recalculate total transpiration per plant exactly like original script
    # --------------------------------------------------------
    dt = 150.0 / (24.0 * 3600.0)  # [days]
    plant_total_T = np.sum(transpiration * dt, axis=1)

    # --------------------------------------------------------
    # Plant group definitions (must match the main script)
    # --------------------------------------------------------
    group_defs = [
        ("auf dem Meliorationsstreifen",            (3, 4)),
        ("nahe am Streifen",   (2, 5)),
        ("mittlerer Abstand zum Streifen",  (1, 6)),
        ("maximaler Abstand zum Streifen", (0, 7)),
    ]

    # --------------------------------------------------------
    # Generate plot
    # --------------------------------------------------------
    plt.figure()
    plt.plot([], [], color='none', label="Pflanzenposition:")
    for name, plants in group_defs:
        profile = np.zeros_like(cum_uptake[0])
        soil_total = 0.0
        plant_total = 0.0

        for pidx in plants:
            profile += cum_uptake[pidx]
            soil_total += cum_uptake[pidx].sum()
            plant_total += plant_total_T[pidx]

        label = f"{name}"
        plt.plot(profile, depth_centers_cm, marker="o", label=label)

    plt.gca().invert_yaxis()
    plt.xlabel("Kumulierte Wasseraufnahme [cm³]")
    plt.ylabel("Tiefe [cm]")
    #plt.title("Kumulierte Aufnahme pro Bodentiefenschicht nach Position")
    plt.legend(loc="lower right")
    plt.tight_layout()

    # Output name
    out_png = npz_path.replace("_timeseries_data.npz", "_uptake_depth_groups_german.png")
    plt.savefig(out_png, dpi=300)
    print(f"Saved: {out_png}")




if __name__ == "__main__":
    main()
