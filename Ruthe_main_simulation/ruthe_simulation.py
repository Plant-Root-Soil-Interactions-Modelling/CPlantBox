import json
import os

import numpy as np
import pandas as pd
from analysis_and_plotting_utils import (
    model_evaluation,
    plot_minirhizotron,
    plot_soil_core,
)
from mini_rhizotron_simulation import run_minirhizotron_simulation
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

    # 2. load soilcore and minirhizotron data
    soil_core_stats, mini_rhizo_stats = load_soil_core_and_minirhizotron_data(year)

    # 3. run simulations
    minirhizo_simulated_mean, minirhizo_simulated_std = run_minirhizotron_simulation(
    dataset, only_touching_segments=False
    )

    soil_core_simulated_mean, soil_core_simulated_std = run_soil_core_simulation(
        dataset
    )

    soil_core_simulated = {
        "mean": soil_core_simulated_mean,
        "std": soil_core_simulated_std,
    }

    # 4. compute model evaluation statistics & save them
    evaluation_metrics = {}

    evaluation_metrics["soil_core"] = model_evaluation(
        measured_data=soil_core_stats,
        simulated_data=soil_core_simulated,
        number_of_observations=dataset.NUMBER_OF_SOIL_CORE_OBSERVATIONS,
        type="soil_core",
    )

    for i in range(len(dataset.DAYS_TO_ANALYZE_MINIRHIZOTRON)):
        day = dataset.DAYS_TO_ANALYZE_MINIRHIZOTRON[i]
        stats = mini_rhizo_stats[day]
        minirhizo_simulated = {
            "mean": minirhizo_simulated_mean[i],
            "std": minirhizo_simulated_std[i],
        }

        if year == 1995 and day != 191:
            n_obs = 20
        else:
            n_obs = dataset.NUMBER_OF_MINIRHIZO_OBSERVATIONS

        evaluation_metrics[day] = model_evaluation(stats, minirhizo_simulated, n_obs, type="minirhizotron")

    out_file = dataset.OUTPUT_DIRECTORY / f"evaluation_metrics_{year}.json"

    with out_file.open("w") as f:
        json.dump(evaluation_metrics, f, indent=4)

    # 5. plot and save results
    plot_soil_core(dataset, soil_core_stats, soil_core_simulated, year)

    for i in range(len(dataset.DAYS_TO_ANALYZE_MINIRHIZOTRON)):
        day = dataset.DAYS_TO_ANALYZE_MINIRHIZOTRON[i]
        stats = mini_rhizo_stats[day]
        plot_minirhizotron(
            dataset, stats, minirhizo_simulated_mean[i], minirhizo_simulated_std[i], i
        )

    save_results_as_csv(dataset.OUTPUT_DIRECTORY, soil_core_stats, mini_rhizo_stats,
                        soil_core_simulated_mean, soil_core_simulated_std,
                        minirhizo_simulated_mean, minirhizo_simulated_std,
                        year, dataset)


def load_soil_core_and_minirhizotron_data(year: int) -> tuple[dict, dict]:
    """
    Load soil core and minirhizotron data for the specified year.

    Args:
        year (int): The year of the experiment.

    Returns:
        tuple: A tuple containing two dictionaries:
            - soil_core_stats: Dictionary with soil core statistics (mean and std).
            - mini_rhizo_stats: Dictionary with minirhizotron statistics (mean and std).
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

        # Minirhizotron
        mini_rhizo_2604 = pd.read_csv(
            dataset.PATH_TO_MINIRHIZOTRON_DATA_1, sep=";", encoding="ISO-8859-1"
        )
        mini_rhizo_2604_optimal_N = mini_rhizo_2604[
            mini_rhizo_2604["N-DÜNGUNG"] == "N normal"
        ]
        mini_rhizo_2604_mean = (
            mini_rhizo_2604_optimal_N.groupby("Tiefe")["Bonitur"].mean().to_frame()
        )
        mini_rhizo_2604_standard_derivation = (
            mini_rhizo_2604_optimal_N.groupby("Tiefe")["Bonitur"].std().to_frame()
        )

        mini_rhizo_1505 = pd.read_csv(dataset.PATH_TO_MINIRHIZOTRON_DATA_2, sep=";")
        mini_rhizo_1505_optimal_N = mini_rhizo_1505[
            mini_rhizo_1505["N-DÜNGUNG"] == "N normal"
        ]
        mini_rhizo_1505_mean = (
            mini_rhizo_1505_optimal_N.groupby("Tiefe")["Bonitur"].mean().to_frame()
        )
        mini_rhizo_1505_standard_derivation = (
            mini_rhizo_1505_optimal_N.groupby("Tiefe")["Bonitur"].std().to_frame()
        )

        mini_rhizo_1306 = pd.read_csv(dataset.PATH_TO_MINIRHIZOTRON_DATA_3, sep=";")
        mini_rhizo_1306_optimal_N = mini_rhizo_1306[
            mini_rhizo_1306["N-Due"] == "Offizialempfehlung"
        ]
        mini_rhizo_1306_mean = (
            mini_rhizo_1306_optimal_N.groupby("Tiefe")["Bonitur"].mean().to_frame()
        )
        mini_rhizo_1306_standard_derivation = (
            mini_rhizo_1306_optimal_N.groupby("Tiefe")["Bonitur"].std().to_frame()
        )

        mini_rhizo_stats = {
            172: {
                "mean": mini_rhizo_2604_mean,
                "std": mini_rhizo_2604_standard_derivation,
            },
            191: {
                "mean": mini_rhizo_1505_mean,
                "std": mini_rhizo_1505_standard_derivation,
            },
            220: {
                "mean": mini_rhizo_1306_mean,
                "std": mini_rhizo_1306_standard_derivation,
            },
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

        # Minirhizotron
        mini_rhizo_1996 = pd.read_csv(
            dataset.PATH_TO_MINIRHIZOTRON_DATA_1, sep=";", encoding="ISO-8859-1"
        )

        for col in ["N-Due", "Datum"]:
            mini_rhizo_1996[col] = (
                mini_rhizo_1996[col]
                .astype(str)
                .str.strip()
            )

        mini_rhizo_1996["Bonitur"] = pd.to_numeric(
            mini_rhizo_1996["Bonitur"], errors="coerce"
        )

        mini_rhizo_1996_optimal_N = mini_rhizo_1996[
            mini_rhizo_1996["N-Due"] == "Offizialempf."
        ]

        mini_rhizo_3004 = mini_rhizo_1996_optimal_N[
            mini_rhizo_1996_optimal_N["Datum"] == "30.04.1996"
        ]
        mini_rhizo_3004_mean = (
            mini_rhizo_3004.groupby("Tiefe")["Bonitur"].mean().to_frame()
        )
        mini_rhizo_3004_standard_derivation = (
            mini_rhizo_3004.groupby("Tiefe")["Bonitur"].std().to_frame()
        )

        mini_rhizo_1405 = mini_rhizo_1996_optimal_N[
            mini_rhizo_1996_optimal_N["Datum"] == "14.05.1996"
        ]
        mini_rhizo_1405_mean = (
            mini_rhizo_1405.groupby("Tiefe")["Bonitur"].mean().to_frame()
        )
        mini_rhizo_1405_standard_derivation = (
            mini_rhizo_1405.groupby("Tiefe")["Bonitur"].std().to_frame()
        )

        mini_rhizo_2805 = mini_rhizo_1996_optimal_N[
            mini_rhizo_1996_optimal_N["Datum"] == "28.05.1996"
        ]
        mini_rhizo_2805_mean = (
            mini_rhizo_2805.groupby("Tiefe")["Bonitur"].mean().to_frame()
        )
        mini_rhizo_2805_standard_derivation = (
            mini_rhizo_2805.groupby("Tiefe")["Bonitur"].std().to_frame()
        )

        mini_rhizo_1806 = mini_rhizo_1996_optimal_N[
            mini_rhizo_1996_optimal_N["Datum"] == "18.06.1996"
        ]
        mini_rhizo_1806_mean = (
            mini_rhizo_1806.groupby("Tiefe")["Bonitur"].mean().to_frame()
        )
        mini_rhizo_1806_standard_derivation = (
            mini_rhizo_1806.groupby("Tiefe")["Bonitur"].std().to_frame()
        )

        mini_rhizo_stats = {
            181: {
                "mean": mini_rhizo_3004_mean,
                "std": mini_rhizo_3004_standard_derivation,
            },
            195: {
                "mean": mini_rhizo_1405_mean,
                "std": mini_rhizo_1405_standard_derivation,
            },
            209: {
                "mean": mini_rhizo_2805_mean,
                "std": mini_rhizo_2805_standard_derivation,
            },
            230: {
                "mean": mini_rhizo_1806_mean,
                "std": mini_rhizo_1806_standard_derivation,
            },
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

        # Minirhizotron
        mini_rhizo_1997 = pd.read_csv(
            dataset.PATH_TO_MINIRHIZOTRON_DATA_1, sep=";", encoding="utf-8-sig"
        )

        mini_rhizo_1997["Bonitur"] = pd.to_numeric(
            mini_rhizo_1997["Bonitur"], errors="coerce"
        )

        mini_rhizo_1997_optimal_N = mini_rhizo_1997[
            mini_rhizo_1997["N-Due"] == "Offizialempf."
        ]

        mini_rhizo_1504 = mini_rhizo_1997_optimal_N[
            mini_rhizo_1997_optimal_N["Datum"] == "15.04.1997"
        ]
        mini_rhizo_1504_mean = (
            mini_rhizo_1504.groupby("Tiefe")["Bonitur"].mean().to_frame()
        )
        mini_rhizo_1504_standard_derivation = (
            mini_rhizo_1504.groupby("Tiefe")["Bonitur"].std().to_frame()
        )

        mini_rhizo_1605 = mini_rhizo_1997_optimal_N[
            mini_rhizo_1997_optimal_N["Datum"] == "16.05.1997"
        ]
        mini_rhizo_1605_mean = (
            mini_rhizo_1605.groupby("Tiefe")["Bonitur"].mean().to_frame()
        )
        mini_rhizo_1605_standard_derivation = (
            mini_rhizo_1605.groupby("Tiefe")["Bonitur"].std().to_frame()
        )

        mini_rhizo_3005 = mini_rhizo_1997_optimal_N[
            mini_rhizo_1997_optimal_N["Datum"] == "30.05.1997"
        ]
        mini_rhizo_3005_mean = (
            mini_rhizo_3005.groupby("Tiefe")["Bonitur"].mean().to_frame()
        )
        mini_rhizo_3005_standard_derivation = (
            mini_rhizo_3005.groupby("Tiefe")["Bonitur"].std().to_frame()
        )

        mini_rhizo_1906 = mini_rhizo_1997_optimal_N[
            mini_rhizo_1997_optimal_N["Datum"] == "19.06.1997"
        ]
        mini_rhizo_1906_mean = (
            mini_rhizo_1906.groupby("Tiefe")["Bonitur"].mean().to_frame()
        )
        mini_rhizo_1906_standard_derivation = (
            mini_rhizo_1906.groupby("Tiefe")["Bonitur"].std().to_frame()
        )

        mini_rhizo_0207 = mini_rhizo_1997_optimal_N[
            mini_rhizo_1997_optimal_N["Datum"] == "02.07.1997"
        ]
        mini_rhizo_0207_mean = (
            mini_rhizo_0207.groupby("Tiefe")["Bonitur"].mean().to_frame()
        )
        mini_rhizo_0207_standard_derivation = (
            mini_rhizo_0207.groupby("Tiefe")["Bonitur"].std().to_frame()
        )

        mini_rhizo_stats = {
            162: {
                "mean": mini_rhizo_1504_mean,
                "std": mini_rhizo_1504_standard_derivation,
            },
            193: {
                "mean": mini_rhizo_1605_mean,
                "std": mini_rhizo_1605_standard_derivation,
            },
            207: {
                "mean": mini_rhizo_3005_mean,
                "std": mini_rhizo_3005_standard_derivation,
            },
            227: {
                "mean": mini_rhizo_1906_mean,
                "std": mini_rhizo_1906_standard_derivation,
            },
            240: {
                "mean": mini_rhizo_0207_mean,
                "std": mini_rhizo_0207_standard_derivation,
            },
        }

    else:
        raise ValueError("Data for the specified year is not available.")

    return soil_core_stats, mini_rhizo_stats


def save_results_as_csv(
    output_directory: str,
    soil_core_stats: dict,
    mini_rhizo_stats: dict,
    soil_core_simulated_mean: np.ndarray,
    soil_core_simulated_std: np.ndarray,
    minirhizo_simulated_mean: dict,
    minirhizo_simulated_std: dict,
    year: int,
    dataset: RutheConfig,
) -> None:
    
    rows = []

    # Soil core
    layers = dataset.NUMBER_OF_LAYERS

    depth = soil_core_stats["mean"].index.to_numpy()

    sc_mean = soil_core_stats["mean"]
    sc_std = soil_core_stats["std"]

    measured_mean = sc_mean.iloc[:, 0].to_numpy()
    measured_std = sc_std.iloc[:, 0].to_numpy()

    for d in range(layers):
        rows.append(
            {
                "system": "soil_core",
                "year": year,
                "day": None,
                "depth": depth[d],
                "measured_mean": float(measured_mean[d]),
                "measured_std": float(measured_std[d]),
                "simulated_mean": float(soil_core_simulated_mean[d]),
                "simulated_std": float(soil_core_simulated_std[d]),
            }
        )

    # Minirhizotron
    dates = list(mini_rhizo_stats.keys())

    for day_index, day in enumerate(dates):
        stats = mini_rhizo_stats[day]

        measured_mean = stats["mean"][6:] # measured data still has values for 0 to 30 cm
        measured_std = stats["std"][6:] # these have to be cut out

        depths = measured_mean.index.to_numpy()

        measured_mean = measured_mean.iloc[:, 0].to_numpy()
        measured_std = measured_std.iloc[:, 0].to_numpy()

        simulated_mean = minirhizo_simulated_mean[day_index, :]
        simulated_std = minirhizo_simulated_std[day_index, :]

        layers = dataset.NUMBER_OF_LAYERS_MINIRHIZOTRON

        for d in range(layers):
            rows.append(
                {
                    "system": "minirhizo",
                    "year": year,
                    "day": day,
                    "depth": depths[d],
                    "measured_mean": float(measured_mean[d]),
                    "measured_std": float(measured_std[d]),
                    "simulated_mean": float(simulated_mean[d]),
                    "simulated_std": float(simulated_std[d]),
                }
            )
    
    # write csv
    os.makedirs(output_directory, exist_ok=True)

    out_path = os.path.join(output_directory, f"full_data_csv_{year}.csv")
    
    pd.DataFrame(rows).to_csv(out_path, index=False)

    print(f"Saved simulation results to {out_path}")


if __name__ == "__main__":
    for year in [1995, 1996, 1997]:
        run_ruthe_simulations(year)
