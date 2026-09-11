import re
import sys
sys.path.append("../../../")
sys.path.append("../../../src")
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plantbox as pb
import matplotlib as mpl
from matplotlib.patches import Patch

# ============================================================
# Settings
# ============================================================

fontsize = 20

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

left  = 0.175  # the left side of the subplots of the figure
right = 0.4    # the right side of the subplots of the figure
bottom = 0.15   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.18  # the amount of width reserved for blank space between subplots
hspace = 0.5   # the amount of height reserved for white space between subplots

PATH = Path("rootsystem")

NAMES = ["RS_optimized_L_WT", "RS_optimized_S_WT"]
SUBSTRATES = ["Loam", "Sand"]
COLORS = ["darkblue", "cornflowerblue"]

times = [42, 63]
depth_labels = ["0–20 cm", "20–40 cm", "40–60 cm"]

# ============================================================
# Read experimental data
# ============================================================

df = pd.read_csv("data/RLD_RD_field_experiment.csv")

measured = {}
measured_se = {}

for sub, substrate in zip(SUBSTRATES, ["L", "S"]):

    data = df[
        (df["substrate"] == substrate) &
        (df["genotype"] == "WT")
    ]

    measured[sub] = (
        data.pivot(index="DAS", columns="depth", values="RD_mean")
        .reindex(index=times, columns=["D1", "D2", "D3"])
        .to_numpy()
        .T * 1e-3
    )

    measured_se[sub] = (
        data.pivot(index="DAS", columns="depth", values="RD_SE")
        .reindex(index=times, columns=["D1", "D2", "D3"])
        .to_numpy()
        .T * 1e-3
    )


# ============================================================
# Simulate root systems
# ============================================================

simulated = {}

r = 2.5
soil_depth = 20

cores = [
    pb.SDF_RotateTranslate(
        pb.SDF_PlantContainer(r, r, soil_depth, False),
        pb.Vector3d(10, 0, z)
    )
    for z in [0, -20, -40]
]

for sub, name in enumerate(NAMES):

    rs = pb.RootSystem()
    rs.readParameters(str(PATH / f"{name}.xml"))
    rs.setSeed(0)
    rs.initializeLB(5, 4)

    simulated[sub] = np.zeros((3, len(times)))

    for day in range(1, max(times) + 1):

        rs.simulate(1, True)

        if day not in times:
            continue

        i = times.index(day)

        for d, core in enumerate(cores):

            ana = pb.SegmentAnalyser(rs)
            ana.mapPeriodic(45, 20)
            ana.crop(core)

            length = np.asarray(ana.getParameter("length"))
            radius = np.asarray(ana.getParameter("radius"))

            simulated[sub][d, i] = (
                np.sum(length * radius) / np.sum(length) * 20
            )


# ============================================================
# Plot
# ============================================================

fig, axes = plt.subplots(
    3, 1,
    figsize=(9, 8),
    sharex=True,
    sharey=True
)

# Slim bars
bar_width = 0.04
group_spacing = 0.5

# Four positions:
# Loam measured, Loam simulated,
# Sand measured, Sand simulated
offsets = [-1.5, -0.5, 0.5, 1.5]

x = np.arange(len(times))/2*group_spacing

for d, ax in enumerate(axes):

    for sub, substrate in enumerate(SUBSTRATES):

        color = COLORS[sub]

        # ----------------------------------------------------
        # Measured
        # ----------------------------------------------------

        # Only plot measured data if available
        if not np.all(np.isnan(measured[substrate][d])):

            ax.bar(
                x + offsets[2 * sub] * bar_width,
                measured[substrate][d],
                yerr=measured_se[substrate][d],
                width=bar_width,
                color=color,
                edgecolor=color,
                capsize=4
            )

        # ----------------------------------------------------
        # Simulated
        # ----------------------------------------------------

        ax.bar(
            x + offsets[2 * sub + 1] * bar_width,
            simulated[sub][d],
            width=bar_width,
            color="white",
            edgecolor=color,
            hatch="///"
        )

    ax.set_xticks(x)
    ax.set_xticklabels([f"DAS {t}" for t in times])

    ax.grid(axis="y", alpha=0.3)
    ax.set_axisbelow(True)

    # Depth label on right
    ax.text(
        1.02, 0.5,
        depth_labels[d],
        transform=ax.transAxes,
        va="center",
        ha="left"
    )


# ============================================================
# Labels
# ============================================================

fig.supylabel("Mean root diameter (mm)", x = 0.12)
axes[-1].set_xlabel("Days after sowing")


# ============================================================
# Legend
# ============================================================

legend_handles = [
    Patch(
        facecolor="darkblue",
        edgecolor="darkblue",
        label="Loam"
    ),
    Patch(
        facecolor="cornflowerblue",
        edgecolor="cornflowerblue",
        label="Sand"
    ),
    Patch(
        facecolor="white",
        edgecolor="black",
        label="Measured"
    ),
    Patch(
        facecolor="white",
        edgecolor="black",
        hatch="///",
        label="Simulated"
    )
]

fig.legend(
    handles=legend_handles,
    loc="upper center",
    ncol=1,
    bbox_to_anchor=(0.55, 0.9),
    frameon=False
)
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show()