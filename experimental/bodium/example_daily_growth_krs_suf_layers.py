"""Daily root growth analysis: Krs, SUF, length per layer, and mean radius per layer.

This example follows the tutorial style and simulates a root system with a daily
(time step = 1 day) loop. For each day it computes:
- root system conductance Krs [cm2 day-1]
- SUF profile per soil layer [-]
- root length per soil layer [cm]
- mean radius per soil layer [cm] (computed from surface/length)

Outputs are written to CSV files in experimental/bodium/results.
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np

import plantbox as pb
from plantbox.structural.MappedOrganism import MappedPlantPython
from plantbox.functional.PlantHydraulicModel import (
    HydraulicModel_Doussan,
    HydraulicModel_Meunier,
)
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters


def make_hydraulic_parameters():
    """Use age-dependent hydraulic parameters (same style as tutorial chapter 4)."""
    params = PlantHydraulicParameters()

    kr0 = np.array([[0.0, 2.2e-4], [12.5, 2.2e-4], [20.9, 8.0e-5], [44.6, 8.0e-5], [62.7, 1.9e-5], [100.0, 1.9e-5]])
    kr1 = np.array([[0.0, 1.8e-4], [10.0, 1.8e-4], [15.0, 1.7e-5], [25.0, 1.7e-5]])
    params.set_kr_age_dependent(kr0[:, 0], kr0[:, 1], subType=[1, 4])
    params.set_kr_age_dependent(kr1[:, 0], kr1[:, 1], subType=[2, 3])

    kx0 = np.array([[0.0, 2.7e-2], [18.3, 2.7e-2], [21.0, 3.3e-1], [47.0, 3.3e-1], [61.0, 4.2], [100.0, 4.2]])
    kx1 = np.array([[0.0, 1.0e-4], [9.0, 2.0e-4], [13.0, 6.0e-4], [20.0, 1.73e-3], [25.0, 1.73e-3]])
    params.set_kx_age_dependent(kx0[:, 0], kx0[:, 1], subType=[1, 4])
    params.set_kx_age_dependent(kx1[:, 0], kx1[:, 1], subType=[2, 3])

    return params


def make_constant_hydraulic_parameters(kr_const=1.728e-4, kx_const=4.32e-2):
    """Use constant hydraulic parameters for all root segments.

    kr_const [day-1]
    kx_const [cm3 day-1]
    """
    params = PlantHydraulicParameters()
    params.set_kr_const(kr_const)
    params.set_kx_const(kx_const)
    return params


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Hydraulic parameter mode: "age_dependent" or "constant"
    hydraulic_mode = "age_dependent"

    # Simulation settings
    sim_days = 62
    dt = 1.0  # daily time step

    # Soil layering for depth profiles
    z_top = 0.0
    z_bottom = -120.0
    n_layers = 24
    layer_depths = np.linspace(z_top - 0.5 * abs((z_top - z_bottom) / n_layers), z_bottom + 0.5 * abs((z_top - z_bottom) / n_layers), n_layers)

    # Root architecture file
    plant_name = "Zea_mays_1_Leitner_2010"
    # plant_name = "Anagallis_femina_Leitner_2010"
    xml_path = os.path.join(script_dir, "..", "..", "modelparameter", "structural", "rootsystem", f"{plant_name}.xml")
    out_dir = os.path.join(script_dir, "results")
    os.makedirs(out_dir, exist_ok=True)

    # Root system and hydraulic model
    plant = pb.MappedPlant()
    plant.readParameters(xml_path)
    # plant.setGeometry(pb.SDF_PlantContainer(1.0e6, 1.0e6, -z_bottom, True))  # large container to avoid boundary effects
    plant.initialize()
    plant_python = MappedPlantPython(plant)

    if hydraulic_mode == "constant":
        hydraulic_params = make_constant_hydraulic_parameters(kr_const=1.728e-4, kx_const=4.32e-2)
    else:
        hydraulic_params = make_hydraulic_parameters()
    # hm = HydraulicModel_Doussan(plant, hydraulic_params)
    hm = HydraulicModel_Meunier(plant, hydraulic_params)

    # Storage
    days = []
    krs_series = []
    length_total_series = []

    length_layers_daily = []
    mean_radius_layers_daily = []
    suf_layers_daily = []
    tip_count_layers_daily = []
    tip_positions_daily = []
    tip_layer_edges = np.linspace(z_bottom, z_top, n_layers + 1)

    for day in range(1, sim_days + 1):

        # Simulate
        plant.simulate(dt, True)

        root_tips = plant_python.get_root_tips()
        tip_positions_daily.append(root_tips.copy())

        # Count root tips per soil layer using tip z-coordinates.
        if root_tips.size > 0:
            tip_z = root_tips[:, 2]
            tip_counts = np.histogram(tip_z, bins=tip_layer_edges)[0][::-1]  # match top->bottom layer order
        else:
            tip_counts = np.zeros(n_layers, dtype=np.int64)
        tip_count_layers_daily.append(tip_counts)

        # Krs at current age
        # hm.update(day)  # in case of doussan
        krs, _ = hm.get_krs(day)
        krs_series.append(krs)
        days.append(day)

        # Segment-based analysis at this day
        ana = pb.SegmentAnalyser(plant)

        length_profile = np.array(ana.distribution("length", z_top, z_bottom, n_layers, True))
        surface_profile = np.array(ana.distribution("surface", z_top, z_bottom, n_layers, True))

        # Mean radius from surface = 2*pi*r*length, guard division by zero for empty layers
        mean_radius_profile = np.zeros_like(length_profile)
        nonzero = length_profile > 0.0
        mean_radius_profile[nonzero] = surface_profile[nonzero] / (2.0 * np.pi * length_profile[nonzero])

        # SUF profile by layer
        # hm.update(day)  # in case of Doussan
        suf = hm.get_suf(day)
        ana.addData("SUF", suf)
        suf_profile = np.array(ana.distribution("SUF", z_top, z_bottom, n_layers, False))  # don't cut SUF

        length_layers_daily.append(length_profile)
        mean_radius_layers_daily.append(mean_radius_profile)
        suf_layers_daily.append(suf_profile)
        length_total_series.append(np.sum(length_profile))
        print("Day {}: Krs = {:.4e} cm2 day-1, suf_sum = {:.2f}".format(day, krs, np.sum(suf)))

    length_layers_daily = np.array(length_layers_daily)  # shape: [time, layer]
    mean_radius_layers_daily = np.array(mean_radius_layers_daily)
    suf_layers_daily = np.array(suf_layers_daily)
    tip_count_layers_daily = np.array(tip_count_layers_daily)

    # 1) Daily scalar outputs
    with open(os.path.join(out_dir, "daily_krs_length.csv"), "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["day", "krs_cm2_per_day", "total_length_cm"])
        for day, krs, length in zip(days, krs_series, length_total_series):
            writer.writerow([day, krs, length])

    # 2) Layer outputs per day
    with open(os.path.join(out_dir, "daily_profiles_by_layer.csv"), "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "day",
                "layer_index",
                "layer_center_depth_cm",
                "length_cm",
                "mean_radius_cm",
                "suf_layer_sum",
                "root_tip_count",
            ]
        )
        for ti, day in enumerate(days):
            for li in range(n_layers):
                writer.writerow(
                    [
                        day,
                        li,
                        layer_depths[li],
                        length_layers_daily[ti, li],
                        mean_radius_layers_daily[ti, li],
                        suf_layers_daily[ti, li],
                        int(tip_count_layers_daily[ti, li]),
                    ]
                )

    # 3) Root tip positions (all tips, all days)
    with open(os.path.join(out_dir, "daily_root_tip_positions.csv"), "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["day", "tip_index", "x_cm", "y_cm", "z_cm"])
        for ti, day in enumerate(days):
            root_tips = tip_positions_daily[ti]
            if root_tips.size == 0:
                continue
            for tip_i, tip in enumerate(root_tips):
                writer.writerow([day, tip_i, tip[0], tip[1], tip[2]])

    # Plot summary curves
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharex=True)

    axes[0].plot(days, krs_series)
    axes[0].set_title("Krs over time")
    axes[0].set_xlabel("Day")
    axes[0].set_ylabel("Krs (cm2 day-1)")
    axes[0].grid(True)

    axes[1].plot(days, length_total_series)
    axes[1].set_title("Total root length")
    axes[1].set_xlabel("Day")
    axes[1].set_ylabel("Length (cm)")
    axes[1].grid(True)

    axes[2].plot(days, np.sum(tip_count_layers_daily, axis=1))
    axes[2].set_title("Total root tips over time")
    axes[2].set_xlabel("Day")
    axes[2].set_ylabel("Root tip count (-)")
    axes[2].grid(True)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "daily_summary.png"), dpi=200)
    plt.show()

    # Plot profile evolution over time
    fig2, axes2 = plt.subplots(1, 4, figsize=(19, 5), sharey=True)
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(vmin=days[0], vmax=days[-1])

    for ti, day in enumerate(days):
        color = cmap(norm(day))
        axes2[0].plot(length_layers_daily[ti, :], layer_depths, color=color, alpha=0.9)
        axes2[1].plot(mean_radius_layers_daily[ti, :], layer_depths, color=color, alpha=0.9)
        axes2[2].plot(suf_layers_daily[ti, :], layer_depths, color=color, alpha=0.9)
        axes2[3].plot(tip_count_layers_daily[ti, :], layer_depths, color=color, alpha=0.9)

    axes2[0].set_title("Length profiles (all days)")
    axes2[0].set_xlabel("Length (cm)")
    axes2[0].set_ylabel("Depth (cm)")
    axes2[0].grid(True)

    axes2[1].set_title("Mean radius profiles (all days)")
    axes2[1].set_xlabel("Mean radius (cm)")
    axes2[1].grid(True)

    axes2[2].set_title("SUF profiles (all days)")
    axes2[2].set_xlabel("Layer SUF sum (-)")
    axes2[2].grid(True)

    axes2[3].set_title("Root tip count profiles (all days)")
    axes2[3].set_xlabel("Tip count (-)")
    axes2[3].grid(True)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig2.subplots_adjust(right=0.88, wspace=0.35)
    cax = fig2.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig2.colorbar(sm, cax=cax)
    cbar.set_label("Simulation day")

    plt.savefig(os.path.join(out_dir, "final_day_profiles.png"), dpi=200)
    plt.show()

    # # Plot root tip z-positions as a day-vs-depth scatter.
    # fig3, ax3 = plt.subplots(1, 1, figsize=(8, 4.5))
    # for ti, day in enumerate(days):
    #     root_tips = tip_positions_daily[ti]
    #     if root_tips.size == 0:
    #         continue
    #     tip_z = root_tips[:, 2]
    #     day_vals = np.full(tip_z.shape, day, dtype=float)
    #     ax3.scatter(day_vals, tip_z, s=12, alpha=0.8, color=cmap(norm(day)))

    # ax3.set_title("Root tip vertical positions over time")
    # ax3.set_xlabel("Day")
    # ax3.set_ylabel("Tip depth z (cm)")
    # ax3.grid(True)
    # plt.tight_layout()
    # plt.savefig(os.path.join(out_dir, "daily_root_tip_positions.png"), dpi=200)
    # plt.show()

    print("Finished simulation.")
    print(f"Wrote: {os.path.join(out_dir, 'daily_krs_length.csv')}")
    print(f"Wrote: {os.path.join(out_dir, 'daily_profiles_by_layer.csv')}")
    print(f"Wrote: {os.path.join(out_dir, 'daily_root_tip_positions.csv')}")
    print(f"Wrote: {os.path.join(out_dir, 'daily_summary.png')}")
    print(f"Wrote: {os.path.join(out_dir, 'final_day_profiles.png')}")
    print(f"Wrote: {os.path.join(out_dir, 'daily_root_tip_positions.png')}")


if __name__ == "__main__":
    main()
