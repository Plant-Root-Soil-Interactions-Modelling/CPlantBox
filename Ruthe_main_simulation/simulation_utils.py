"""Helper functions called by both the soil core and minirhizotron simulations."""

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import plantbox as pb


def make_seed(
    base_seed: int,
    *,
    namespace: str,
    run: int,
    plant_index: int,
    row_index: int,
) -> int:
    """Create a deterministic, collision-resistant seed.

    The previous seed scheme used simple addition, which caused many collisions
    (different plants ending up with the same seed) and high overlap between
    simulation runs. This function derives a stable, well-mixed integer seed
    from the identifying tuple.

    Returns a non-negative 31-bit integer, compatible with common RNG seed APIs.
    """

    payload = f"{int(base_seed)}|{namespace}|{int(run)}|{int(row_index)}|{int(plant_index)}".encode(
        "utf-8"
    )
    digest = hashlib.blake2b(payload, digest_size=8).digest()
    return int.from_bytes(digest, byteorder="little", signed=False) & 0x7FFFFFFF


def load_soil_temperatures(
    filepath: Path, start_date: str, end_date: str
) -> pd.DataFrame:
    """
    Load soil temperature data from CSV or pickle file.

    Args:
        filepath (Path): Path to the CSV or pickle file.
        start_date (str): Start date for filtering (inclusive).
        end_date (str): End date for filtering (inclusive).

    Returns:
        pd.DataFrame: DataFrame containing soil temperature data within the specified date range.
    """

    if filepath.suffix == ".csv":
        # implement version to use csv instead
        pass
    elif filepath.suffix == ".pkl":
        df_soil_temp = pd.read_pickle(filepath)
        return df_soil_temp.loc[start_date:end_date]
    else:
        raise ValueError(
            "Unsupported file format for soil temperatures. Use csv or pickle file instead."
        )


def initialize_soil_temperature(
    df_soil_temp: pd.DataFrame, time_step: int
) -> np.ndarray:
    """
    Initialize soil temperature profile as numpy array from DataFrame.

    Args:
        df_soil_temp (pd.DataFrame): DataFrame containing soil temperature data.

    Returns:
        np.ndarray: Array of soil temperatures at different depths.
    """

    soil_temp = np.zeros((14,))

    soil_temp[0] = df_soil_temp["10"].iloc[time_step]
    soil_temp[1] = df_soil_temp["20"].iloc[time_step]
    soil_temp[2] = df_soil_temp["30"].iloc[time_step]
    soil_temp[3] = df_soil_temp["40"].iloc[time_step]
    soil_temp[4] = df_soil_temp["50"].iloc[time_step]
    # soil_temp[5] = df_soil_temp['60'].iloc[time_step] #only 1995 has values for 60cm (others have a second 10cm column)

    soil_temp[5:] = soil_temp[
        4
    ]  # anything below measured data depth set to value at 50cm

    return soil_temp


def compute_temperature_scaling(soil_temp: np.array, dataset) -> np.ndarray:
    """
    Compute temperature scaling factors based on soil temperature profile.

    Args:
        soil_temp (np.array): Array of soil temperatures at different depths.
        dataset: Module containing global variables for a certain experiment.

    Returns:
        np.ndarray: Array of temperature scaling factors.
    """

    scales = np.ones((14,))

    for depth in range(soil_temp.size):
        if soil_temp[depth] > dataset.TEMP_MAX or soil_temp[depth] < dataset.TEMP_MIN:
            scales[depth] = 0.0
        elif dataset.TEMP_OPT < 0.5 * (dataset.TEMP_MIN + dataset.TEMP_MAX):
            T = (soil_temp[depth] - dataset.TEMP_MIN) / (
                dataset.TEMP_MAX - dataset.TEMP_MIN
            )
            sigma = np.log(0.5) / np.log(
                (dataset.TEMP_OPT - dataset.TEMP_MIN)
                / (dataset.TEMP_MAX - dataset.TEMP_MIN)
            )
            scales[depth] = np.sin((np.pi) * (T**sigma))
        elif dataset.TEMP_OPT >= 0.5 * (dataset.TEMP_MIN + dataset.TEMP_MAX):
            T = (soil_temp[depth] - dataset.TEMP_MAX) / (
                dataset.TEMP_MIN - dataset.TEMP_MAX
            )
            sigma = np.log(0.5) / np.log(
                (dataset.TEMP_OPT - dataset.TEMP_MAX)
                / (dataset.TEMP_MIN - dataset.TEMP_MAX)
            )
            scales[depth] = np.sin((np.pi) * (T**sigma))

    return scales


def create_scale_elongation_grid(scales: np.array) -> pb.EquidistantGrid1D:
    """
    Create a scale elongation grid based on temperature scaling factors.

    Args:
        scales (np.array): Array of temperature scaling factors.

    Returns:
        pb.EquidistantGrid1D: Equidistant grid with scale elongation data.
    """

    scale_elongation = pb.EquidistantGrid1D(0, -150, 15)
    scale_elongation.data = scales
    return scale_elongation
