import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ruthe_global_variables import RutheConfig


def _soil_core_layer_centers(dataset: RutheConfig, n_layers: int) -> np.ndarray:
    layer_thickness = dataset.SOIL_CORE_DEPTH / n_layers
    return (np.arange(n_layers, dtype=float) + 0.5) * layer_thickness


def _resample_profile(
    profile: pd.DataFrame, target_depths_cm: np.ndarray
) -> np.ndarray:
    series = pd.Series(
        profile.iloc[:, 0].to_numpy(dtype=float),
        index=pd.Index(profile.index.to_numpy(dtype=float), name="depth_cm"),
    ).sort_index()

    target_index = pd.Index(target_depths_cm, name="depth_cm")
    interpolated = (
        series.reindex(series.index.union(target_index))
        .sort_index()
        .interpolate(method="index", limit_direction="both")
        .reindex(target_index)
    )

    return interpolated.to_numpy()

def plot_soil_core(
    dataset: RutheConfig, measured_data: dict, simulated_data: dict, year: int
) -> None:
    simulated_mean = np.asarray(simulated_data["mean"]).squeeze()
    simulated_std = np.asarray(simulated_data["std"]).squeeze()
    depths = -_soil_core_layer_centers(dataset, simulated_mean.size)

    measured_mean = _resample_profile(measured_data["mean"], -depths)
    measured_std = _resample_profile(measured_data["std"], -depths)

    # set axis label and text font sizes
    plt.rc("font", size=18)
    plt.rc("axes", labelsize=18)

    # plot observed data with error bars
    lower_error = np.minimum(measured_std, measured_mean)
    plt.errorbar(
        measured_mean,
        depths,
        xerr=[lower_error, measured_std],
        fmt="o",
        label="observed",
    )

    # plot simulated data
    plt.plot(simulated_mean, depths, label="simulated")
    # upper and lower uncertainty bounds
    plt.plot(simulated_mean + simulated_std, depths, linestyle="--", color="black")
    plt.plot(
        np.maximum(simulated_mean - simulated_std, 0),
        depths,
        linestyle="--",
        color="black",
    )

    # axis labels
    plt.xlabel(r"RLD (cm / cm$^3$)")
    plt.ylabel("Depth (cm)")

    # add title to give date of soil core measurement
    plt.title(f"{dataset.MEASUREMENT_DAY_SOIL_CORE}")

    plt.savefig(
        f"{dataset.OUTPUT_DIRECTORY}/depth_profiles/soil_core_comparison_{year}_cores_{dataset.NUMBER_OF_SOIL_CORE_SAMPLES}_layers_{dataset.NUMBER_OF_LAYERS}.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()
