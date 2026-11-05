"""
Fit deepest root depth with:
- capped linear growth
- sigmoid (logistic) growth

Outputs:
- elongation rates
- R² values
- scatter plot with group colors and markers
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re


# ===== SETTINGS =====

LENGTHS_CSV = Path("deepest_tip_per_id.csv")
AGES_CSV = Path("RootIntensity2.csv")

PLOTS_DIR = Path("plots_growth")
PLOTS_DIR.mkdir(exist_ok=True)

USE_IDS = set(range(81, 87)) | set(range(41, 47)) | set(range(91, 97))
USE_IDS.discard(44)

ID_COL = "ID"
T_COL = "DaS"
L_COL = "length_path_cm"
LMAX = 83.0


# ===== MODELS =====

def capped_linear(t, a):
    return np.minimum(a * t, LMAX)


def logistic(t, k, t0):
    return LMAX / (1 + np.exp(-k * (t - t0)))


# ===== METRICS =====

def r2(y, yhat):
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1 - ss_res / ss_tot


# ===== FITS =====

def fit_linear(t, l):
    t_pos = t[t > 0]
    a_max = np.max(l) / np.min(t_pos) * 3
    grid = np.linspace(1e-6, a_max, 2000)

    best = None
    for a in grid:
        pred = capped_linear(t, a)
        sse = np.sum((l - pred) ** 2)
        if best is None or sse < best["sse"]:
            best = {"a": a, "sse": sse}

    return best["a"]


def fit_sigmoid(t, l):

    k_grid = np.linspace(0.01, 0.5, 150)
    t0_grid = np.linspace(np.min(t), np.max(t), 150)

    best = None

    for k in k_grid:
        for t0 in t0_grid:
            pred = logistic(t, k, t0)
            sse = np.sum((l - pred) ** 2)
            if best is None or sse < best["sse"]:
                best = {"k": k, "t0": t0, "sse": sse}

    return best["k"], best["t0"]


# ===== MAIN =====

def main():

    lengths = pd.read_csv(LENGTHS_CSV)
    ages = pd.read_csv(AGES_CSV, sep=";")

    def extract_id(x):
        nums = re.findall(r"\d+", str(x))
        return int(nums[-1]) if nums else np.nan

    lengths[ID_COL] = lengths[ID_COL].apply(extract_id)
    ages[ID_COL] = ages[ID_COL].apply(extract_id)

    lengths = lengths[lengths[ID_COL].isin(USE_IDS)]
    ages = ages[ages[ID_COL].isin(USE_IDS)]

    merged = lengths.merge(
        ages[[ID_COL, T_COL, "VAR", "IRR"]],
        on=ID_COL
    )

    merged = merged.dropna(subset=[T_COL, L_COL])

    t = merged[T_COL].to_numpy(float)
    l = merged[L_COL].to_numpy(float)

    # linear
    a = fit_linear(t, l)
    l_lin = capped_linear(t, a)
    r2_lin = r2(l, l_lin)

    # sigmoid
    k, t0 = fit_sigmoid(t, l)
    l_sig = logistic(t, k, t0)
    r2_sig = r2(l, l_sig)

    elong_sig = k * LMAX / 4.0

    # plot
    t_grid = np.linspace(0, 40, 600)

    colors = {"S": "orange", "D": "green", "M": "blue"}
    markers = {"-": "^", "+": "o"}

    plt.figure()

    for (var, irr), df in merged.groupby(["VAR", "IRR"]):

        plt.scatter(
            df[T_COL],
            df[L_COL],
            color=colors.get(var, "gray"),
            marker=markers.get(irr, "o"),
            s=70
        )

    plt.plot(t_grid, capped_linear(t_grid, a), color="black", lw=2.5)
    plt.plot(t_grid, logistic(t_grid, k, t0), color="purple", lw=2.5)

    plt.xlabel("Age [days]")
    plt.ylabel("Deepest Root Tip [cm]")

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "fit_ALL.png", dpi=200)
    plt.close()

    # output
    print("\n===== RESULTS =====")
    print(f"Linear elongation rate = {a:.3f} cm/day")
    print(f"Linear R² = {r2_lin:.4f}")
    print(f"\nSigmoid max elongation rate = {elong_sig:.3f} cm/day")
    print(f"Sigmoid R² = {r2_sig:.4f}")


if __name__ == "__main__":
    main()