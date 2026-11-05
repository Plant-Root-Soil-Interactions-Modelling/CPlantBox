"""
Compute first-order lateral growth time and plot lateral length by E and F.

Outputs:
- results_lat1_growth_time.csv
- lateral_length_vs_growth_time_by_EF.png
- n, mean, std of elongation rate for total, E, and F
"""

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ===== SETTINGS =====

R = 2.197

LAT_FILE = Path("results_lat_origin_nearest_main_pos.csv")
MAIN_FILE = Path("results_main_laterals_metrics.csv")
OUT_FILE = Path("results_lat1_growth_time.csv")
PLOT_DIR = Path("plots")
PLOT_DIR.mkdir(exist_ok=True)

LAT1_PATTERN = re.compile(r"^lat_\d+(?:_\d+)?$")

plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 12,
    "figure.titlesize": 20
})


# ===== HELPERS =====

def is_lat_first_order(name):
    if name is None:
        return False
    s = str(name).strip()
    if s.startswith("lat_lat"):
        return False
    return bool(LAT1_PATTERN.match(s))


def classify_ef(image_name):
    s = str(image_name).strip()
    if re.search(r"\bE\.\d+\b", s, flags=re.IGNORECASE):
        return "E"
    if re.search(r"\bF\.\d+\b", s, flags=re.IGNORECASE):
        return "F"
    return "Other"


# ===== MAIN =====

def main():
    if not LAT_FILE.exists() or not MAIN_FILE.exists():
        print("Required CSV files not found.")
        sys.exit(1)

    lat_df = pd.read_csv(LAT_FILE)
    main_df = pd.read_csv(MAIN_FILE)

    required_lat = {"image", "main", "lat_root", "position_from_tip", "length"}
    required_main = {"image", "main_name", "la"}

    if not required_lat.issubset(lat_df.columns):
        raise ValueError(f"Missing columns in {LAT_FILE.name}")
    if not required_main.issubset(main_df.columns):
        raise ValueError(f"Missing columns in {MAIN_FILE.name}")

    lat_df["lat_root"] = lat_df["lat_root"].astype(str).str.strip()
    lat_df["image"] = lat_df["image"].astype(str).str.strip()
    lat_df["main"] = lat_df["main"].astype(str).str.strip()

    main_df["image"] = main_df["image"].astype(str).str.strip()
    main_df["main_name"] = main_df["main_name"].astype(str).str.strip()

    lat1_df = lat_df[lat_df["lat_root"].map(is_lat_first_order)].copy()

    main_lookup = (
        main_df
        .drop_duplicates(subset=["image", "main_name"], keep="first")
        .loc[:, ["image", "main_name", "la"]]
    )

    merged = lat1_df.merge(
        main_lookup,
        how="left",
        left_on=["image", "main"],
        right_on=["image", "main_name"]
    )

    merged["position_from_tip"] = pd.to_numeric(merged["position_from_tip"], errors="coerce")
    merged["la"] = pd.to_numeric(merged["la"], errors="coerce")
    merged["length"] = pd.to_numeric(merged["length"], errors="coerce")

    merged["growth_length"] = merged["position_from_tip"] - merged["la"]
    merged["growth_time_days"] = merged["growth_length"] / R

    out_cols = [
        "image", "main", "lat_root",
        "position_from_tip", "la",
        "growth_length", "growth_time_days",
        "position_from_base", "length"
    ]
    out_cols = [c for c in out_cols if c in merged.columns]

    out_df = merged[out_cols].copy()
    out_df.to_csv(OUT_FILE, index=False)

    plot_df = merged.dropna(subset=["growth_time_days", "length"]).copy()
    plot_df = plot_df[plot_df["growth_time_days"] > 0].copy()
    plot_df["EF_group"] = plot_df["image"].apply(classify_ef)
    plot_df = plot_df[plot_df["EF_group"].isin(["E", "F"])].copy()

    if len(plot_df) == 0:
        print("\n===== OUTPUT =====")
        print(f"Results saved to: {OUT_FILE.resolve()}")
        print("No valid E/F points for plot or statistics.")
        print("Finished.")
        return

    styles = {
        "E": {"color": "blue", "marker": "o"},
        "F": {"color": "green", "marker": "^"}
    }

    plt.figure(figsize=(8, 6))

    for group in ["E", "F"]:
        sub = plot_df[plot_df["EF_group"] == group]
        if len(sub) == 0:
            continue
        plt.scatter(
            sub["growth_time_days"],
            sub["length"],
            color=styles[group]["color"],
            marker=styles[group]["marker"],
            s=70,
            alpha=0.85,
            label=group
        )

    plt.xlabel("Growth time [days]")
    plt.ylabel("Lateral length [cm]")
    plt.legend()
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "lateral_length_vs_growth_time_by_EF.png", dpi=300)
    plt.close()

    rate_df = plot_df.copy()
    rate_df["elong_rate"] = rate_df["length"] / rate_df["growth_time_days"]
    rate_df = rate_df.replace([np.inf, -np.inf], np.nan)
    rate_df = rate_df.dropna(subset=["elong_rate"])

    total_n = len(rate_df)
    total_mean = rate_df["elong_rate"].mean()
    total_std = rate_df["elong_rate"].std()

    rate_E = rate_df[rate_df["EF_group"] == "E"]
    rate_F = rate_df[rate_df["EF_group"] == "F"]

    print("\n===== OUTPUT =====")
    print(f"Results saved to: {OUT_FILE.resolve()}")
    print("Plot saved: lateral_length_vs_growth_time_by_EF.png")
    print(f"\nTotal: n = {total_n}, mean = {total_mean:.3f} cm/day, std = {total_std:.3f} cm/day")

    if len(rate_E) > 0:
        print(
            f"E: n = {len(rate_E)}, "
            f"mean = {rate_E['elong_rate'].mean():.3f} cm/day, "
            f"std = {rate_E['elong_rate'].std():.3f} cm/day"
        )
    else:
        print("E: n = 0")

    if len(rate_F) > 0:
        print(
            f"F: n = {len(rate_F)}, "
            f"mean = {rate_F['elong_rate'].mean():.3f} cm/day, "
            f"std = {rate_F['elong_rate'].std():.3f} cm/day"
        )
    else:
        print("F: n = 0")

    print("Finished.")


if __name__ == "__main__":
    main()