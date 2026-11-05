"""
Fit basal root count over time with a linear model.

Outputs:
- initiation rate
- DelayB
- firstB
- scatter plot with group colors and markers
"""

import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ===== SETTINGS =====

BASAL_CSV = Path("Basal_Root.csv")
GROUPS_CSV = Path("RootIntensity2.csv")

PLOTS_DIR = Path("plots")
PLOTS_DIR.mkdir(exist_ok=True)


# ===== IO =====

def read_csv_robust(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        try:
            sep = csv.Sniffer().sniff(sample, delimiters=";,").delimiter
        except csv.Error:
            sep = ";"
    df = pd.read_csv(path, sep=sep, encoding="utf-8-sig")
    df.columns = [c.strip().replace("\ufeff", "") for c in df.columns]
    return df


def extract_id(x):
    nums = re.findall(r"\d+", str(x))
    return int(nums[-1]) if nums else np.nan


# ===== MAIN =====

def main():
    basal = read_csv_robust(BASAL_CSV)
    groups = read_csv_robust(GROUPS_CSV)

    basal["ID"] = basal["ID"].apply(extract_id)
    groups["ID"] = groups["ID"].apply(extract_id)

    try:
        merged = basal.merge(groups[["ID", "VAR", "IRR"]], on="ID")
        sign_col = "IRR"
    except KeyError:
        merged = basal.merge(groups[["ID", "VAR", "TRR"]], on="ID")
        sign_col = "TRR"

    merged["Age"] = pd.to_numeric(merged["Age"], errors="coerce")
    merged["BasalrootCount"] = pd.to_numeric(merged["BasalrootCount"], errors="coerce") 
    merged["VAR"] = merged["VAR"].astype(str).str.strip()
    merged[sign_col] = merged[sign_col].astype(str).str.strip()

    merged = merged.dropna(subset=["Age", "BasalrootCount"])

    x = merged["Age"].to_numpy(float)
    y = merged["BasalrootCount"].to_numpy(float)

    m, b = np.polyfit(x, y, 1)

    delay_b = 1.0 / m
    first_b = (1.0 - b) / m
    x_zero = -b / m

    x_fit = np.linspace(min(x_zero, np.min(x)) - 1, np.max(x) + 1, 400)
    y_fit = m * x_fit + b

    colors = {"S": "orange", "D": "green", "M": "blue"}
    markers = {"-": "^", "+": "o"}

    plt.figure()

    seen = set()
    for (var, sign), sub in merged.groupby(["VAR", sign_col]):
        label = f"{var}{sign}"
        plt.scatter(
            sub["Age"],
            sub["BasalrootCount"],
            color=colors.get(var, "gray"),
            marker=markers.get(sign, "o"),
            s=70,
            label=label if label not in seen else None
        )
        seen.add(label)

    plt.plot(x_fit, y_fit, color="black", lw=2.5, label="Linear fit")
    plt.axvline(first_b, color="gray", linestyle="--", lw=1.8, label="First basal root")
    plt.axhline(1, color="gray", linestyle="--", lw=1.8, label="One basal root")

    handles, labels = plt.gca().get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    plt.legend(unique.values(), unique.keys())

    plt.xlabel("Age [days]")
    plt.ylabel("Basal Root Count")

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "Basalroot_vs_Age.png", dpi=300)
    plt.close()

    print("\n===== RESULTS =====")
    print(f"Initiation rate = {m:.3f} roots/day")
    print(f"DelayB = {delay_b:.2f} days/root")
    print(f"firstB = {first_b:.2f} days")


if __name__ == "__main__":
    main()