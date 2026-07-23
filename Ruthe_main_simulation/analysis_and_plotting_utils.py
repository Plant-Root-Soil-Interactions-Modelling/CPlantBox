import matplotlib.pyplot as plt
import numpy as np
from ruthe_global_variables import RutheConfig


def model_evaluation(measured_data: dict, simulated_data: dict, number_of_observations: int, type: str) -> dict:
    # Log-likelihood
    if type == "minirhizotron":
        measured_mean = np.asarray(measured_data["mean"]).squeeze()[6:]
        measured_std = np.asarray(measured_data["std"]).squeeze()[6:]
    else:
        measured_mean = np.asarray(measured_data["mean"]).squeeze()
        measured_std = np.asarray(measured_data["std"]).squeeze()

    simulated_mean = np.asarray(simulated_data["mean"]).squeeze()
    simulated_std = np.asarray(simulated_data["std"]).squeeze()

    n_obs = number_of_observations

    N1 = len(measured_mean)
    N2 = len(measured_std)

    sigma1 = measured_std / np.sqrt(n_obs)
    sigma2 = measured_std / np.sqrt(2.0 * (n_obs - 1))

    e1 = measured_mean - simulated_mean
    e2 = measured_std - simulated_std

    if np.any(sigma1 <= 0) or np.any(sigma2 <= 0):
        log_p = np.nan

    else:
        log_p1 = (
            -(N1 / 2.0) * np.log(2.0 * np.pi)
            - np.sum(np.log(sigma1))
            - 0.5 * np.sum((e1 / sigma1) ** 2.0)
        )

        log_p2 = (
            -(N2 / 2.0) * np.log(2.0 * np.pi)
            - np.sum(np.log(sigma2))
            - 0.5 * np.sum((e2 / sigma2) ** 2.0)
        )

        log_p = log_p1 + log_p2

    # root mean suare error (RSME)
    rmse = np.sqrt((np.mean(e1**2.0)))

    # normalized RSME
    nrmse = (rmse / np.mean(measured_mean)) * 100

    # nash-sutcliffe efficiency
    denom = np.sum((measured_mean - np.mean(measured_mean)) ** 2.0)
    nse = 1.0 - np.sum(e1**2.0) / denom if denom > 0 else np.nan

    # cosine similarity
    cos_sim = (
        np.dot(measured_mean, simulated_mean)
        / (np.linalg.norm(measured_mean) * np.linalg.norm(simulated_mean))
        if np.linalg.norm(measured_mean) and np.linalg.norm(simulated_mean)
        else np.nan
    )

    result = {
        "log_p": log_p,
        "rmse": rmse,
        "nrmse": nrmse,
        "nse": nse,
        "cos_sim": cos_sim,
    }

    return result


def plot_soil_core(
    dataset: RutheConfig, measured_data: dict, simulated_data: dict, year: int
) -> None:
    depths = (
        np.linspace(-15, -dataset.SOIL_CORE_DEPTH, dataset.NUMBER_OF_LAYERS) + 7.5
    )  # mid-layer depths

    measured_mean = np.asarray(measured_data["mean"]).squeeze()
    measured_std = np.asarray(measured_data["std"]).squeeze()

    simulated_mean = np.asarray(simulated_data["mean"]).squeeze()
    simulated_std = np.asarray(simulated_data["std"]).squeeze()

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
        f"{dataset.OUTPUT_DIRECTORY}/depth_profiles/soil_core_comparison_{year}.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()


def plot_minirhizotron(
    dataset: RutheConfig,
    measured_data: dict,
    simulated_mean: np.ndarray,
    simulated_std: np.ndarray,
    day_index: int,
) -> None:
    depths = (
        np.linspace(
            -35, -dataset.SOIL_CORE_DEPTH, dataset.NUMBER_OF_LAYERS_MINIRHIZOTRON
        )
        + 2.5
    )  # mid-layer depths

    measured_mean = np.asarray(measured_data["mean"]).squeeze()[
        6:
    ]  # cut out first 30 cm like in simulation
    measured_std = np.asarray(measured_data["std"]).squeeze()[
        6:
    ]  # cut out first 30 cm like in simulation

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
    plt.xlabel("Root Score (-)")
    plt.ylabel("Depth (cm)")
    plt.xlim(-0.5, 5)

    # add title to give date of minirhizotron measurement
    plt.title(f"{dataset.MEASUREMENT_DAYS_MINIRHIZOTRON[day_index]}")

    plt.savefig(
        f"{dataset.OUTPUT_DIRECTORY}/depth_profiles/minirhizotron_comparison_day_{dataset.DAYS_TO_ANALYZE_MINIRHIZOTRON[day_index]}.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()
