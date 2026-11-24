import os

import numpy as np
import pandas as pd
from analysis_and_plotting_utils import (
    _resample_profile,
    _soil_core_layer_centers,
    plot_soil_core,
)
from ruthe_global_variables import (
    RUTHE_1994_95,
    RUTHE_1995_96,
    RUTHE_1996_97,
    RutheConfig,
)
from soil_core_simulation import run_soil_core_simulation


def run_ruthe_simulations(year: int) -> None:

    if year == 1995:
        dataset = RUTHE_1994_95
    elif year == 1996:
        dataset = RUTHE_1995_96
    elif year == 1997:
        dataset = RUTHE_1996_97
    else:
        raise ValueError("Simulation for the specified year is not implemented.")

    print(year)

    # 1. load soilcore and minirhizotron data
    soil_core_stats = load_soil_core_data(year)

    soil_core_simulated_mean, soil_core_simulated_std = run_soil_core_simulation(
        dataset
    )

    soil_core_simulated = {
        "mean": soil_core_simulated_mean,
        "std": soil_core_simulated_std,
    }

    # 2. plot and save results
    plot_soil_core(dataset, soil_core_stats, soil_core_simulated, year)

    save_results_as_csv(
        dataset.OUTPUT_DIRECTORY,
        soil_core_stats,
        soil_core_simulated_mean,
        soil_core_simulated_std,
        year,
        dataset,
    )


def load_soil_core_data(year: int) -> dict:
    """
    Load soil core data for the specified year.

    Args:
        year (int): The year of the experiment.

    Returns:
        dict: dictionary with soil core statistics (mean and std).
    """

    if year == 1995:
        dataset = RUTHE_1994_95

        # Soil core
        soil_cores = pd.read_csv(dataset.PATH_TO_SOIL_CORE_DATA, sep=";")
        soil_cores["Wurzellängendichte"] = pd.to_numeric(
            soil_cores["Wurzellängendichte"], errors="coerce"
        )
        soil_cores_optimal_N = soil_cores[soil_cores["N-Due"] == "Offizialempfehlung"]
        soil_cores_mean = (
            soil_cores_optimal_N.groupby("Tiefe")["Wurzellängendichte"]
            .mean()
            .to_frame()
        )
        soil_cores_standard_derivation = (
            soil_cores_optimal_N.groupby("Tiefe")["Wurzellängendichte"].std().to_frame()
        )

        soil_core_stats = {
            "mean": soil_cores_mean,
            "std": soil_cores_standard_derivation,
        }

    elif year == 1996:
        dataset = RUTHE_1995_96

        # Soil core
        soil_cores = pd.read_csv(dataset.PATH_TO_SOIL_CORE_DATA, sep=";")
        soil_cores["WLD"] = pd.to_numeric(soil_cores["WLD"], errors="coerce")
        soil_cores_optimal_N = soil_cores[soil_cores["N-Due"] == "Offizialempf."]
        soil_cores_mean = soil_cores_optimal_N.groupby("Tiefe")["WLD"].mean().to_frame()
        soil_cores_standard_derivation = (
            soil_cores_optimal_N.groupby("Tiefe")["WLD"].std().to_frame()
        )

        soil_core_stats = {
            "mean": soil_cores_mean,
            "std": soil_cores_standard_derivation,
        }

    elif year == 1997:
        dataset = RUTHE_1996_97

        # Soil core
        soil_cores = pd.read_csv(dataset.PATH_TO_SOIL_CORE_DATA, sep=";")
        soil_cores["WLD"] = pd.to_numeric(soil_cores["WLD"], errors="coerce")
        soil_cores_optimal_N = soil_cores[soil_cores["N-DÜNG."] == "Offizialempf."]
        soil_cores_mean = soil_cores_optimal_N.groupby("TIEFE")["WLD"].mean().to_frame()
        soil_cores_standard_derivation = (
            soil_cores_optimal_N.groupby("TIEFE")["WLD"].std().to_frame()
        )

        soil_core_stats = {
            "mean": soil_cores_mean,
            "std": soil_cores_standard_derivation,
        }

    else:
        raise ValueError("Data for the specified year is not available.")

    return soil_core_stats


def save_results_as_csv(
    output_directory: str,
    soil_core_stats: dict,
    soil_core_simulated_mean: np.ndarray,
    soil_core_simulated_std: np.ndarray,
    year: int,
    dataset: RutheConfig,
) -> None:

    rows = []

    # Soil core
    layers = len(soil_core_simulated_mean)
    depth = _soil_core_layer_centers(dataset, layers)

    measured_mean = _resample_profile(soil_core_stats["mean"], depth)
    measured_std = _resample_profile(soil_core_stats["std"], depth)

    for d in range(layers):
        rows.append(
            {
                "system": "soil_core",
                "year": year,
                "day": None,
                "depth": float(depth[d]),
                "measured_mean": float(measured_mean[d]),
                "measured_std": float(measured_std[d]),
                "simulated_mean": float(soil_core_simulated_mean[d]),
                "simulated_std": float(soil_core_simulated_std[d]),
            }
        )

    # write csv
    os.makedirs(output_directory, exist_ok=True)

    out_path = os.path.join(
        output_directory,
        f"full_data_csv_{year}_cores_{dataset.NUMBER_OF_SOIL_CORE_SAMPLES}_layers_{dataset.NUMBER_OF_LAYERS}.csv",
    )

    pd.DataFrame(rows).to_csv(out_path, index=False)

    print(f"Saved simulation results to {out_path}")


if __name__ == "__main__":
    for year in [1995, 1996, 1997]:
        run_ruthe_simulations(year)
