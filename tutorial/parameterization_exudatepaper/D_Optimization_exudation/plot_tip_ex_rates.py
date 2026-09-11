import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

fontsize = 24

mpl.rcParams.update({
    "font.size": fontsize,
    "axes.titlesize": fontsize,
    "axes.labelsize": fontsize,
    "xtick.labelsize": fontsize,
    "ytick.labelsize": fontsize,
    "legend.fontsize": fontsize,
    "mathtext.default": "regular"
})

# Read data
df = pd.read_csv("data/exudation_rates.csv")

DAS = df["DAS"].values
loam = df["Loam"].values
sand = df["Sand"].values

# Bar positions
x = np.arange(len(DAS))
bar_width = 0.1

# Plot
fig, ax = plt.subplots(figsize=(8, 6))

ax.bar(
    x - bar_width / 2,
    loam,
    width=bar_width,
    color="b"
)

ax.bar(
    x + bar_width / 2,
    sand,
    width=bar_width,
    color="cornflowerblue"
)

# Axes
ax.set_xticks(x)
ax.set_xticklabels([f"{d}" for d in DAS])

ax.set_xlabel("Days after sowing")
ax.set_ylabel("Root tip exudation rate \n(mol cm$^{-2}$ d$^{-1}$)")

ax.ticklabel_format(
    axis="y",
    style="sci",
    scilimits=(-2, 2)
)

ax.grid(axis="y", alpha=0.3)
ax.set_axisbelow(True)

# Legend
legend_handles = [
    Patch(facecolor="b", edgecolor="b", label="Loam"),
    Patch(facecolor="cornflowerblue", edgecolor="cornflowerblue", label="Sand")
]

fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.3, 0.92),
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