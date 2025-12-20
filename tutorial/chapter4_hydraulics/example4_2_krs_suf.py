"""Whole root system conductance (Krs) for different root architectures and SUF depth profiles"""

import csv  # |\label{l42:imports}|

import matplotlib.pyplot as plt
import numpy as np

import plantbox as pb
from plantbox.functional.PlantHydraulicModel import HydraulicModel_Meunier
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters  # |\label{l42:imports_end}|
from plantbox.visualisation import figure_style 

# Simulation parameters   # |\label{l42:parameters}|
simtime = 70  # simulate from day 1 to 70
dt = 1 

architectures = [  # |\label{l42:architecture}|
    "Heliantus_Pages_2013",
    "Glycine_max_Moraes2020_opt2",
    "Brassica_oleracea_Vansteenkiste_2014",
    "Zea_mays_1_Leitner_2010",
]

path = "../../modelparameter/structural/rootsystem/"  # |\label{l42:architecture_end}|

# Root hydraulic properties  # |\label{l42:roothydraulics}|
param = PlantHydraulicParameters()

kr0 = np.array([[0.0, 2.2e-4], [12.5, 2.2e-4], [20.9, 8.0e-5], [44.6, 8.0e-5], [62.7, 1.9e-5], [100, 1.9e-5]])
kr1 = np.array([[0.0, 1.8e-4], [10, 1.8e-4], [15, 1.7e-5], [25, 1.7e-5]])
param.set_kr_age_dependent(kr0[:, 0], kr0[:, 1], subType=[1, 4])
param.set_kr_age_dependent(kr1[:, 0], kr1[:, 1], subType=[2, 3])

kx0 = np.array([[0, 2.7e-2], [18.3, 2.7e-2], [21, 3.3e-1], [47, 3.3e-1], [61, 4.2], [100, 4.2]])
kx1 = np.array([[0, 1.0e-4], [9, 2.0e-4], [13, 6.0e-4], [20, 1.73e-3], [25, 1.73e-3]])
param.set_kx_age_dependent(kx0[:, 0], kx0[:, 1], subType=[1, 4])
param.set_kx_age_dependent(kx1[:, 0], kx1[:, 1], subType=[2, 3])  # |\label{l42:roothydraulics_end}|

# Simulation loop # |\label{l42:sim_start}|
krs_all = []
lengths = []
surfaces = []
csv_data = []
suf_profiles = []  # To store suf depth profiles at final simulation time

for name in architectures:
    print(f"\nSimulating: {name}")

    plant = pb.MappedPlant()  # |\label{l42:mappedplant}|
    plant.readParameters(path + name + ".xml")
    plant.initialize()

    hm = HydraulicModel_Meunier(plant, param)  # |\label{l42:model}|

    krs_values = []
    arch_lengths = []
    arch_surfaces = []

    for t in range(0, simtime):
        plant.simulate(dt)
        ns = plant.getNumberOfMappedSegments()  # segments extracted
        krs, _ = hm.get_krs(t)  # krs calculated
        krs_values.append(krs)

        total_length = np.sum(np.array(plant.getParameter("length")))
        total_surface = np.sum(np.array(plant.getParameter("surface")))
        arch_lengths.append(total_length)
        arch_surfaces.append(total_surface)

        csv_data.append([name, t, krs, total_length, total_surface])

    krs_all.append(krs_values)
    lengths.append(arch_lengths[-1])  # final length
    surfaces.append(arch_surfaces[-1])  # final surface |\label{l42:sim_end}|

    suf = hm.get_suf(simtime)  # |\label{l42:suf}|
    ana = pb.SegmentAnalyser(plant)
    ana.addData("SUF", suf)  # suf at each segment in the 3d space

    bin_size = 5  # soil layer thickness is defined, for plotting purposes
    z_max = 0
    z_min = -150
    n_bins = int((z_max - z_min) / bin_size)

    suf_dist = ana.distribution("SUF", z_max, z_min, n_bins, True)  # |\label{l42:suf_end}|

    depths = np.linspace(z_max - bin_size / 2, z_min + bin_size / 2, n_bins)  # |\label{l42:suf_filter}|
    suf_array = np.array(suf_dist)
    depth_array = np.array(depths)

    nonzero_indices = np.where(suf_array > 0)[0]
    if len(nonzero_indices) > 0:
        suf_profiles.append((suf_array[nonzero_indices], depth_array[nonzero_indices]))
    else:  # fallback if no roots (shouldn't normally happen)
        suf_profiles.append((suf_array, depth_array))  # |\label{l42:suf_filter_end}|

# Plotting
n_arch = len(architectures)  # |\label{l42:plotting_krs}|
fig_krs, axes_krs = figure_style.subplots12(1, n_arch, sharey=True)

if n_arch == 1:
    axes_krs = [axes_krs]

for i, ax in enumerate(axes_krs):
    ax.plot(range(0, simtime), krs_all[i])
    ax.set_title(architectures[i])
    ax.set_xlabel("Root system age (day)")
    ax.set_yscale("log")
    if i == 0:
        ax.set_ylabel("Krs (cm$^2$ day$^{-1}$")
    ax.grid(True)

plt.tight_layout()
plt.show()  # |\label{l42:plotting_krs_end}|

fig_suf, axes_suf = figure_style.subplots12(1, n_arch, sharey=True)  # |\label{l42:plotting_suf}|

if n_arch == 1:
    axes_suf = [axes_suf]

for i, ax in enumerate(axes_suf):
    suf, depth = suf_profiles[i]
    ax.plot(suf, depth, "k")
    ax.set_xlabel("SUF [-]")
    if i == 0:
        ax.set_ylabel("Depth [cm]")
    ax.grid(True)
    ax.set_title(architectures[i])
    if len(depth) > 0:
        ax.set_ylim(depth.min() - 5, depth.max() + 5)  # optional padding

plt.tight_layout()
plt.show()  # |\label{l42:plotting_suf_krs}|

# Printing summary and saving outputs
print("\nSummary:")  # |\label{l42:summary}|
for i, name in enumerate(architectures):
    print(f"{name:20s}, Total root length: {lengths[i]:8.2f} cm, Surface area: {surfaces[i]:8.2f} cm2")

csv_file = "results/krs_length_surface.csv"  # write CSV
with open(csv_file, mode="w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["architecture", "day", "krs", "length", "surface"])
    writer.writerows(csv_data)

print(f"\nSaved results to: {csv_file}")  # |\label{l42:summary_end}|
