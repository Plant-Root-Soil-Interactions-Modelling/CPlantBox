"""Optimize root length / diameter"""

import sys
sys.path.append("../../../../CPlantBox")
sys.path.append("../../../../CPlantBox/src")

import plantbox as pb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import interpolate


# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------

WRITE_ROOTSYS = False

fontsize = 20
mpl.rcParams.update({
    'font.size': fontsize,
    'axes.titlesize': fontsize,
    'axes.labelsize': fontsize,
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,
    'legend.fontsize': fontsize,
    'figure.titlesize': fontsize,
    'mathtext.default': 'regular'
})

left  = 0.1  # the left side of the subplots of the figure
right = 0.95   # the right side of the subplots of the figure
bottom = 0.35   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.5  # the amount of width reserved for blank space between subplots
hspace = 0.5   # the amount of height reserved for white space between subplots

soiltypes = ['L', 'S']
titles = ['Loam', 'Sand']

# Colors
colors = {
    'Loam': 'blue',
    'Sand': 'cornflowerblue'
}

CSV_OUTPUT = "data/root_results.csv"


# ---------------------------------------------------------------------
# Function to update CSV
# ---------------------------------------------------------------------

def update_csv(results, filename):
    """Create or update the results CSV."""

    new_data = pd.DataFrame(results)

    try:
        old_data = pd.read_csv(filename)

        # Remove existing rows for the same DAS + soil type
        keys = new_data[['DAS', 'soil type']]

        old_data = old_data.merge(
            keys,
            on=['DAS', 'soil type'],
            how='left',
            indicator=True
        )

        old_data = old_data[old_data['_merge'] == 'left_only']
        old_data = old_data.drop(columns='_merge')

        data = pd.concat(
            [old_data, new_data],
            ignore_index=True
        )

    except FileNotFoundError:
        data = new_data

    data = data.sort_values(['DAS', 'soil type'])

    data.to_csv(filename, index=False)


# ---------------------------------------------------------------------
# Plot setup
# ---------------------------------------------------------------------

all_results = []


# ---------------------------------------------------------------------
# Loop over soil types
# ---------------------------------------------------------------------

for i, soiltype in enumerate(soiltypes):

    # ---------------------------------------------------------------
    # Experimental data
    # ---------------------------------------------------------------

    df = pd.read_csv("data/column_experiment_mean.csv")

    data = df[
        (df['substrate'] == soiltype) &
        (df['genotype'] == "WT")
    ]

    DAS = data["DAS"].values

    real_length_all = data['RL_mean'].values
    real_length_SE = data['RL_SE'].values

    real_diam_all = data['RD_mean'].values
    real_diam_SE = data['RD_SE'].values

    real_exu_all = data["Total exudation"].values * 24 / 1000
    real_exu_SE = data["Total exudation SE"].values * 24 / 1000

    # Target times
    times_ = DAS[:2]

    real_length = real_length_all[:2]
    real_diam = real_diam_all[:2]
    real_exu = real_exu_all[:2]

    real_length_SE_plot = real_length_SE[:2]
    real_diam_SE_plot = real_diam_SE[:2]
    real_exu_SE_plot = real_exu_SE[:2]


    # ---------------------------------------------------------------
    # Load root system
    # ---------------------------------------------------------------

    path = "rootsystem/"
    name = "RS_optimized_" + soiltype + "_WT"

    rs = pb.RootSystem()
    rs.readParameters(path + name + ".xml")
    rs.setSeed(0)


    # ---------------------------------------------------------------
    # Simulation without tube
    # ---------------------------------------------------------------

    rs.initializeLB(5, 4)

    for day in range(1, times_[-1] + 1):
        rs.simulate(1, True)

        if WRITE_ROOTSYS:
            rs.write(
                f"RS_visualisation/notube_{soiltype}_WT_day_{day}.vtp"
            )


    # ---------------------------------------------------------------
    # Simulation with tube
    # ---------------------------------------------------------------

    tube_ = pb.SDF_PlantContainer(10, 10, 80, False)

    tube = pb.SDF_RotateTranslate(
        tube_,
        0,
        pb.SDF_Axis.zaxis,
        pb.Vector3d(0, 0, 10)
    )

    rs.setGeometry(tube)
    rs.initializeLB(5, 4)

    times = np.arange(0, times_[-1] + 1)

    comp_length = np.zeros(len(times))
    comp_diam = np.zeros(len(times))
    comp_vol = np.zeros(len(times))
    comp_exu = np.zeros(len(times))
    


    # ---------------------------------------------------------------
    # Maximum root lengths
    # ---------------------------------------------------------------

    lmax = [
        pp.lmax
        for pp in rs.getRootRandomParameter()
    ]


    # ---------------------------------------------------------------
    # Exudation rates
    # ---------------------------------------------------------------

    exu_df = pd.read_csv("data/exudation_rates.csv")

    exu_DAS = np.insert(
        exu_df["DAS"].values,
        0,
        0
    )

    exu_rates = np.insert(
        exu_df[titles[i]].values,
        0,
        exu_df[titles[i]].values[0]
    )

    exu_function = interpolate.interp1d(
        exu_DAS,
        exu_rates
    )

    tip = 3.5  # cm


    # ---------------------------------------------------------------
    # Simulation
    # ---------------------------------------------------------------

    for day in range(1, times_[-1] + 1):

        rs.simulate(1, True)

        # -----------------------------------------------------------
        # Root length and diameter
        # -----------------------------------------------------------

        lengths = np.asarray(
            rs.getParameter("length")
        )

        radii = np.asarray(
            rs.getParameter("radius")
        )

        comp_length[day] = np.sum(lengths)

        mean_radius = (
            np.sum(lengths * radii) /
            np.sum(lengths)
        )

        comp_diam[day] = mean_radius * 20  # mm

        comp_vol[day] = np.sum(lengths * radii**2 * np.pi)  # cm^3
        
        # -----------------------------------------------------------
        # Exudation
        # -----------------------------------------------------------

        kex_tip = exu_function(day)
        kex_base = kex_tip / 2

        types = np.asarray(
            rs.getParameter("type")
        )

        polylines = rs.getPolylines()

        exudation = []

        for j, polyline in enumerate(polylines):

            radius = radii[j]
            root_type = int(types[j])

            accumulated_length = 0

            for k in range(len(polyline) - 1):

                p0 = np.array([
                    polyline[-1-k].x,
                    polyline[-1-k].y,
                    polyline[-1-k].z
                ])

                p1 = np.array([
                    polyline[-2-k].x,
                    polyline[-2-k].y,
                    polyline[-2-k].z
                ])

                segment_length = np.linalg.norm(
                    p0 - p1
                )

                accumulated_length += segment_length

                # Tip or base
                if (lengths[j] - accumulated_length) < tip:
                    kexu = kex_tip
                else:
                    kexu = kex_base

                # Root has stopped growing
                if lengths[j] >= lmax[root_type] * 0.99:
                    kexu = kex_base

                # Artificial shoot
                if root_type == 0:
                    kexu = 0

                exudation.append(
                    2 * np.pi *
                    radius *
                    segment_length *
                    kexu
                )

        comp_exu[day] = np.sum(exudation)


    # ---------------------------------------------------------------
    # Save values for target times
    # ---------------------------------------------------------------

    for j, day in enumerate(times_):

        all_results.append({
            'DAS': day,
            'soil type': titles[i],

            'comp_length': comp_length[day] / 100,
            'real_length': real_length[j] / 100,
            'real_length_SE': real_length_SE_plot[j] / 100,

            'comp_diam': comp_diam[day],
            'real_diam': real_diam[j],
            'real_diam_SE': real_diam_SE_plot[j],
            
            'comp_vol': comp_vol[day],

            'comp_exu': comp_exu[day],
            'real_exu': real_exu[j],
            'real_exu_SE': real_exu_SE_plot[j]
        })


# ---------------------------------------------------------------------
# Update CSV
# ---------------------------------------------------------------------

update_csv(
    all_results,
    CSV_OUTPUT
)


# ---------------------------------------------------------------------
# Prepare data for plotting
# ---------------------------------------------------------------------

results = pd.DataFrame(all_results)

DAS_values = sorted(results['DAS'].unique())

x = np.arange(len(DAS_values))

# Width of individual bars
bar_width = 0.1


# ---------------------------------------------------------------------
# Plot setup
# ---------------------------------------------------------------------

fig, ax = plt.subplots(
    1, 3,
    figsize=(15, 12),
    sharex=True
)


# ---------------------------------------------------------------------
# Function for plotting one variable
# ---------------------------------------------------------------------

def plot_bars(
    axis,
    measured_col,
    simulated_col,
    ylabel,
    measured_se
):

    # Positions:
    # Loam measurement
    # Loam simulation
    # Sand measurement
    # Sand simulation
    positions = [-1.5, -0.5, 0.5, 1.5]

    for i, soil in enumerate(['Loam', 'Sand']):

        soil_data = results[
            results['soil type'] == soil
        ].set_index('DAS')

        color = colors[soil]

        # -------------------------------------------------------------
        # Measurement
        # -------------------------------------------------------------

        measured_values = [
            soil_data.loc[d, measured_col]
            for d in DAS_values
        ]

        measured_errors = [
            soil_data.loc[d, measured_se]
            for d in DAS_values
        ]

        axis.bar(
            x + positions[i * 2] * bar_width,
            measured_values,
            width=bar_width,
            color=color,
            edgecolor=color,
            yerr=measured_errors,
            capsize=5,
            error_kw={
                'elinewidth': 1.5,
                'capthick': 1.5
            }
        )

        # -------------------------------------------------------------
        # Simulation
        # -------------------------------------------------------------

        simulated_values = [
            soil_data.loc[d, simulated_col]
            for d in DAS_values
        ]

        axis.bar(
            x + positions[i * 2 + 1] * bar_width,
            simulated_values,
            width=bar_width,
            facecolor='white',
            edgecolor=color,
            hatch='///'
        )

    axis.set_ylabel(ylabel)


# ---------------------------------------------------------------------
# Root length
# ---------------------------------------------------------------------

plot_bars(
    ax[0],
    measured_col='real_length',
    simulated_col='comp_length',
    ylabel='Total root\nlength (m)',
    measured_se='real_length_SE'
)

# ax[0].text(
    # 0.02, 0.95,
    # 'Root length',
    # transform=ax[0].transAxes,
    # ha='left',
    # va='top'
# )


# ---------------------------------------------------------------------
# Root diameter
# ---------------------------------------------------------------------

plot_bars(
    ax[1],
    measured_col='real_diam',
    simulated_col='comp_diam',
    ylabel='Mean root\ndiameter (mm)',
    measured_se='real_diam_SE'
)

ax[1].set_ylim([0, 0.8])

# ax[1].text(
    # 0.02, 0.95,
    # 'Root diameter',
    # transform=ax[1].transAxes,
    # ha='left',
    # va='top'
# )


# ---------------------------------------------------------------------
# Total plant exudation
# ---------------------------------------------------------------------

plot_bars(
    ax[2],
    measured_col='real_exu',
    simulated_col='comp_exu',
    ylabel='Total plant exudation\n(mol/day/plant)',
    measured_se='real_exu_SE'
)

ax[2].set_ylim([0, 0.01])

# ax[2].text(
    # 0.02, 0.95,
    # 'Total plant exudation',
    # transform=ax[2].transAxes,
    # ha='left',
    # va='top'
# )


# ---------------------------------------------------------------------
# X-axis
# ---------------------------------------------------------------------

ax[2].set_xticks(x)
ax[2].set_xticklabels(
    [f'{d}' for d in DAS_values]
)

# ax[2].set_xlabel('Time (days after sowing)')


# ---------------------------------------------------------------------
# Legend
# ---------------------------------------------------------------------

from matplotlib.patches import Patch

legend_handles = [

    # Soil types
    Patch(
        facecolor='blue',
        edgecolor='blue',
        label='Loam'
    ),

    Patch(
        facecolor='cornflowerblue',
        edgecolor='cornflowerblue',
        label='Sand'
    ),

    # Measurement / simulation
    Patch(
        facecolor='white',
        edgecolor='black',
        label='Measurement'
    ),

    Patch(
        facecolor='white',
        edgecolor='black',
        hatch='///',
        label='Simulation'
    )
]

fig.legend(
    handles=legend_handles,
    loc='upper center',
    bbox_to_anchor=(0.5, 0.98),
    ncol=4,
    frameon=False
)


# ---------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------

# plt.tight_layout(rect=[0, 0, 1, 0.95])
fig.supxlabel('Days after sowing (d)',y = 0.25)
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show()

