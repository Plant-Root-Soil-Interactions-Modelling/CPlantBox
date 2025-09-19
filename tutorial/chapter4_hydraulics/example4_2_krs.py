""" Whole root system conductance (Krs) for different root architectures """

import sys; sys.path.append("../"); sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.PlantHydraulicParameters import PlantHydraulicParameters  # |\label{l42:imports}|
from functional.PlantHydraulicModel import HydraulicModel_Meunier  # |\label{l42:imports_end}|

import numpy as np
import matplotlib.pyplot as plt
import figure_style
import csv

""" Simulation parameters """  # |\label{l42:parameters}|

simtime = 70  # simulate from day 1 to 70
dt = 1

architectures = [  # |\label{l42:architecture}|
    "Heliantus_Pages_2013",
    "Glycine_max_Moraes2020_opt2",
    "Brassica_oleracea_Vansteenkiste_2014",
    "Zea_mays_1_Leitner_2010"
]

path = "../../modelparameter/structural/rootsystem/"  # |\label{l42:architecture_end}|

""" Root hydraulic properties """  # |\label{l42:roothydraulics}|

param = PlantHydraulicParameters()

kr0 = np.array([[0., 2.2e-4], [12.5, 2.2e-4], [20.9, 8.0e-5], [44.6, 8.0e-5], [62.7, 1.9e-5], [100, 1.9e-5]])
kr1 = np.array([[0., 1.8e-4], [10, 1.8e-4], [15, 1.7e-5], [25, 1.7e-5]])
param.set_kr_age_dependent(kr0[:, 0], kr0[:, 1], subType = [1, 4])  #
param.set_kr_age_dependent(kr1[:, 0], kr1[:, 1], subType = [2, 3])

kx0 = np.array([[0, 2.7e-2], [18.3, 2.7e-2], [21, 3.3e-1], [47, 3.3e-1], [61, 4.2], [100, 4.2]])
kx1 = np.array([[0, 1.e-4], [9, 2.e-4], [13, 6.e-4], [20, 1.73e-3], [25, 1.73e-3]])
param.set_kx_age_dependent(kx0[:, 0], kx0[:, 1], subType = [1, 4])
param.set_kx_age_dependent(kx1[:, 0], kx1[:, 1], subType = [2, 3])  # |\label{l42:roothydraulics_end}|

""" Simulation loop"""  # |\label{l42:sim_start}|

krs_all = []
lengths = []
surfaces = []
csv_data = []  # to store output for CSV

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
        ns = plant.getNumberOfMappedSegments()  # |\label{l42:segments}|
        krs, _ = hm.get_krs(t)
        krs_values.append(krs)  # |\label{l42:krs}|

        total_length = np.sum(np.array(plant.getParameter("length")))
        total_surface = np.sum(np.array(plant.getParameter("surface")))
        arch_lengths.append(total_length)
        arch_surfaces.append(total_surface)

        # Store all data per day
        csv_data.append([name, t, krs, total_length, total_surface])

    krs_all.append(krs_values)
    lengths.append(arch_lengths[-1])  # final length
    surfaces.append(arch_surfaces[-1])  # final surface # |\label{l42:sim_end}|

""" Plotting """  # |\label{l42:plotting}|

n_arch = len(architectures)
fig, axes = plt.subplots(1, n_arch, figsize = (5 * n_arch, 4), sharey = True)

if n_arch == 1:
    axes = [axes]

for i, ax in enumerate(axes):
    ax.plot(range(0, simtime), krs_all[i])
    ax.set_title(architectures[i])
    ax.set_xlabel("Root system age (days)")
    ax.set_yscale("log")
    if i == 0:
        ax.set_ylabel("Krs (cm$^2$/day)")
    ax.grid(True)

plt.tight_layout()
plt.show()  # |\label{l42:plotting_end}|

""" Printing summary and saving outputs """  # |\label{l42:summary}|

print("\nSummary:")
for i, name in enumerate(architectures):
    print(f"{name:20s} | Total root length: {lengths[i]:8.2f} cm | Surface area: {surfaces[i]:8.2f} cm2")

# Write CSV
csv_file = "results/krs_length_surface.csv"
with open(csv_file, mode = "w", newline = "") as f:
    writer = csv.writer(f)
    writer.writerow(["architecture", "day", "krs", "length", "surface"])
    writer.writerows(csv_data)

print(f"\nSaved results to: {csv_file}")  # |\label{l42:summary_end}|

