import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl


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

# Read CSV
df = pd.read_csv("data/simulated_upscaled_TRL.csv")

# Extract data
DAS = df["DAS"].values
loam_sim = df["Loam simulation approach"].values/100
sand_sim = df["Sand simulation approach"].values/100
loam_up = df["Loam upscaled"].values/100
sand_up = df["Sand upscaled"].values/100

# Bar positions
x = range(len(DAS))
width = 0.18

fig, ax = plt.subplots(figsize=(6, 5))

# Simulated = hatched
ax.bar(
    [i - 1.5 * width for i in x],
    loam_sim,
    width,
    color="w",
    hatch="//",
    edgecolor="blue"
)

ax.bar(
    [i + 0.5 * width for i in x],
    sand_sim,
    width,
    color="w",
    hatch="//",
    edgecolor="cornflowerblue"
)

# Upscaled = filled
ax.bar(
    [i - 0.5 * width for i in x],
    loam_up,
    width,
    color="blue",
    edgecolor="black"
)

ax.bar(
    [i + 1.5 * width for i in x],
    sand_up,
    width,
    color="cornflowerblue",
    edgecolor="black"
)

# X axis
ax.set_xticks(list(x))
ax.set_xticklabels([f"DAS {int(das)}" for das in DAS])

ax.set_ylabel("Total root length (m)")
ax.set_xlabel("")

# Legend
legend_elements = [
    Patch(facecolor="blue", edgecolor="black", label="Loam"),
    Patch(facecolor="cornflowerblue", edgecolor="black", label="Sand"),
    Patch(facecolor="white", edgecolor="black", hatch="//", label="Simulated"),
    Patch(facecolor="white", edgecolor="black", label="Upscaled"),
]

ax.legend(handles=legend_elements)

plt.tight_layout()
plt.show()