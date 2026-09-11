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

PATH = Path("rootsystem")
NAMES = ["RS_optimized_L_WT", "RS_optimized_S_WT"]
# NAMES = ["RS_base_trial", "RS_base_trial"]
SUBSTRATES = ["Loam", "Sand"]
COLORS = ["darkblue", "cornflowerblue"]
NAMES_ABB = ["L", "S"]

depths = [-20, -40, -60]
depth_labels = ["0–20", "20–40", "40–60"]
times = [42, 63]

bar_height = 1
group_spacing = 5


# ============================================================
# Figure
# ============================================================

fig, axes = plt.subplots(
    2, 1,
    figsize=(12, 9),
    sharex=True,
    sharey=True
)


# ============================================================
# Store data
# ============================================================

measured = {}
measured_SE = {}
simulated = {}
TRL = np.zeros((len(times), len(SUBSTRATES)))


# ============================================================
# Loop over substrates
# ============================================================

for sub, name in enumerate(NAMES):

    # --------------------------------------------------------
    # Experimental data
    # --------------------------------------------------------

    df = pd.read_csv("data/RLD_RD_field_experiment.csv")

    data = df[
        (df["substrate"] == NAMES_ABB[sub]) &
        (df["genotype"] == "WT")
    ]

    measured[sub] = (
        data.pivot(index="DAS", columns="depth", values="RLD_mean")
        .reindex(index=times, columns=["D1", "D2", "D3"])
        .to_numpy()
        .T
    )

    measured_SE[sub] = (
        data.pivot(index="DAS", columns="depth", values="RLD_SE")
        .reindex(index=times, columns=["D1", "D2", "D3"])
        .to_numpy()
        .T
    )

    # --------------------------------------------------------
    # Simulate root system
    # --------------------------------------------------------

    rs = pb.RootSystem()
    rs.readParameters(str(PATH / f"{name}.xml"))
    rs.setSeed(0)
    rs.initializeLB(5, 4)

    r = 2.5
    soil_depth = 20
    sc_vol = np.pi * r**2 * soil_depth

    cores = [
        pb.SDF_RotateTranslate(
            pb.SDF_PlantContainer(r, r, soil_depth, False),
            pb.Vector3d(10, 0, z)
        )
        for z in [0, -20, -40]
    ]

    simulated[sub] = np.zeros((3, len(times)))

    for day in range(1, max(times) + 1):

        rs.simulate(1, True)
        

        # rs.write('test_plots/'+name+'_day'+str(day)+'.vtp')
        if day not in times:
            continue

        i = times.index(day)
        
        #make total length table 
        TRL[i,sub] = np.sum(np.asarray(rs.getParameter("length")))
        
        for d, core in enumerate(cores):
            

            ana = pb.SegmentAnalyser(rs)
            if d == 0: 
                ana.write('test_plots/'+name+'_day'+str(day)+'.vtp')
            ana.mapPeriodic(45, 20)
            if d == 0: 
                ana.write('test_plots/mapped_'+name+'_day'+str(day)+'.vtp')
            ana.crop(core)
            ana.write('test_plots/core'+str(d)+'_'+name+'_day'+str(day)+'.vtp')


            simulated[sub][d, i] = (
                np.sum(ana.getParameter("length")) / sc_vol
            )

times = np.asarray(times).reshape(-1, 1)
TRL = np.hstack((times, TRL))
headers = ['DAS', 'Loam', 'Sand']
np.savetxt(
    "data/simulated_TRL.csv",
    TRL,
    delimiter=",",
    header=",".join(headers),
    comments=""
)

# ============================================================
# Plot
# ============================================================

for i, ax in enumerate(axes):

    # Position of the four bars within each depth group:
    # Loam measured
    # Loam simulated
    # Sand measured
    # Sand simulated
    offsets = [1.5, 0.5, -0.5, -1.5]

    y_positions = []

    # Create positions for the three depth groups
    for d in range(3):

        center = -d * group_spacing

        # ----------------------------------------------------
        # Loam
        # ----------------------------------------------------

        ax.barh(
            center + offsets[0] * bar_height,
            measured[0][d, i],
            xerr=measured_SE[0][d, i],
            height=bar_height,
            color=COLORS[0],
            edgecolor=COLORS[0]
        )

        ax.barh(
            center + offsets[1] * bar_height,
            simulated[0][d, i],
            height=bar_height,
            facecolor="white",
            edgecolor=COLORS[0],
            hatch="///"
        )

        # ----------------------------------------------------
        # Sand
        # ----------------------------------------------------

        ax.barh(
            center + offsets[2] * bar_height,
            measured[1][d, i],
            xerr=measured_SE[1][d, i],
            height=bar_height,
            color=COLORS[1],
            edgecolor=COLORS[1]
        )

        ax.barh(
            center + offsets[3] * bar_height,
            simulated[1][d, i],
            height=bar_height,
            facecolor="white",
            edgecolor=COLORS[1],
            hatch="///"
        )

        y_positions.append(center)

    # --------------------------------------------------------
    # Depth labels
    # --------------------------------------------------------

    ax.set_yticks(y_positions)
    ax.set_yticklabels(depth_labels)

    # DAS label
    ax.text(
        0.95, 0.08,
        f"DAS {times[i]}",
        transform=ax.transAxes,
        ha="right"
    )

    ax.grid(axis="x", linestyle=":", alpha=0.6)
    ax.set_axisbelow(True)


# ============================================================
# Labels
# ============================================================

fig.supxlabel("Root length density (cm cm$^{-3}$)")
fig.supylabel("Soil depth (cm)")


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
    ncol=4,
    bbox_to_anchor=(0.5, 0.98),
    frameon=False
)

plt.show()