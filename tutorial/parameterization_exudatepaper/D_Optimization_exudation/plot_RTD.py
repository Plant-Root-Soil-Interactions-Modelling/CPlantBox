import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import sys

fontsize = 24

mpl.rcParams.update({
    "font.size": fontsize,
    "axes.titlesize": fontsize,
    "axes.labelsize": fontsize,
    "xtick.labelsize": fontsize,
    "ytick.labelsize": fontsize,
    "legend.fontsize": fontsize,
    "figure.titlesize": fontsize,
    "mathtext.default": "regular"
})

# ---------------------------------------------------------
# Read data
# ---------------------------------------------------------

df = pd.read_csv("data/column_experiment_mean.csv")

# Select WT, DAS 42 and 63
df = df[
    (df["genotype"] == "WT") &
    (df["DAS"].isin([42, 63])) &
    (df["substrate"].isin(["L", "S"]))
].copy()

df1 = pd.read_csv("data/root_results.csv")

# ---------------------------------------------------------
# Calculate root tissue density
# ---------------------------------------------------------

# Root diameter: mm -> cm
RD_cm = df["RD_mean"] * 1e-1

# Root volume [cm3]
df["Root volume"] = (
    df1["comp_vol"].values
)

# Root tissue density [g/cm3]
df["Root tissue density"] = (
    df["Root dry weight"] / df["Root volume"]
)

df["Root tissue density SE"] = (
    df["Root dry weight SE"] / df["Root volume"]
)

print(df[[
    "DAS",
    "substrate",
    "RL_mean",
    "RD_mean",
    "Root dry weight",
    "Root volume",
    "Root tissue density"
]])


# ---------------------------------------------------------
# Prepare data for plotting
# ---------------------------------------------------------

DAS = [42, 63]

loam = [
    df[(df["DAS"] == das) & (df["substrate"] == "L")]
    ["Root tissue density"].iloc[0]
    for das in DAS
]

loam_SE = [
    df[(df["DAS"] == das) & (df["substrate"] == "L")]
    ["Root tissue density SE"].iloc[0]
    for das in DAS
]


sand = [
    df[(df["DAS"] == das) & (df["substrate"] == "S")]
    ["Root tissue density"].iloc[0]
    for das in DAS
]

sand_SE = [
    df[(df["DAS"] == das) & (df["substrate"] == "S")]
    ["Root tissue density SE"].iloc[0]
    for das in DAS
]

# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------

x = np.arange(len(DAS))
bar_width = 0.1

fig, ax = plt.subplots(figsize=(8, 6))

ax.bar(
    x - bar_width / 2,
    loam,
    width=bar_width,
    yerr=loam_SE,
    color="b"
)

ax.bar(
    x + bar_width / 2,
    sand,
    width=bar_width,
    yerr=sand_SE,
    color="cornflowerblue"
)

ax.set_xticks(x)
ax.set_xticklabels([f"{das}" for das in DAS])

ax.set_xlabel("Days after sowing")
ax.set_ylabel("Root tissue density \n(g dry root cm$^{-3}$ root volume)")

ax.grid(axis="y", alpha=0.3)
ax.set_axisbelow(True)

# ---------------------------------------------------------
# Legend
# ---------------------------------------------------------

legend_handles = [
    Patch(
        facecolor="b",
        edgecolor="b",
        label="Loam"
    ),
    Patch(
        facecolor="cornflowerblue",
        edgecolor="cornflowerblue",
        label="Sand"
    )
]

fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.28, 0.92),
    ncol=2,
    frameon=False
)

plt.subplots_adjust(
    left=0.15,
    right=0.4,
    bottom=0.3,
    top=0.82
)

plt.show()