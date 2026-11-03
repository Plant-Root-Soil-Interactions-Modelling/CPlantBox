"""Helper functions called by both the soil core and minirhizotron simulations."""

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import plantbox as pb

# depth (cm) down to which the temperature profile / scale elongation grid is defined,
# matches the depth of the soil domain (SOIL_SPACE_PARAMS)
PROFILE_DEPTH = 500
PROFILE_RESOLUTION = 10  # cm, spacing of the soil temperature depths


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

    payload = f"{int(base_seed)}|{namespace}|{int(run)}|{int(row_index)}|{int(plant_index)}".encode()
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
        df_soil_temp = pd.read_csv(filepath, index_col="date", parse_dates=True)
        df_soil_temp = df_soil_temp.drop(
            columns=df_soil_temp.columns[df_soil_temp.columns.str.startswith("source_")]
        )
        return df_soil_temp.loc[start_date:end_date]
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

    measured = df_soil_temp.iloc[time_step].to_numpy(dtype=float)

    # one value per 10cm layer (+1 for the lower grid boundary), anything below
    # the deepest data depth set to the value at that depth
    n_values = PROFILE_DEPTH // PROFILE_RESOLUTION + 1
    soil_temp = np.pad(measured, (0, n_values - measured.size), mode="edge")

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

    scales = np.ones_like(soil_temp)

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

    scale_elongation = pb.EquidistantGrid1D(0, -PROFILE_DEPTH, scales.size)
    scale_elongation.data = scales
    return scale_elongation
