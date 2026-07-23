"""Run simulation and take virtual soil core samples to compare to real life data later"""

import random
from pprint import pprint
from pathlib import Path

import numpy as np
import pandas as pd
import plantbox as pb
from ruthe_global_variables import RUTHE_1994_95, RutheConfig
from simulation_utils import (
    compute_temperature_scaling,
    create_scale_elongation_grid,
    initialize_soil_temperature,
    load_soil_temperatures,
    make_seed,
)


def run_soil_core_simulation(
    dataset: RutheConfig = RUTHE_1994_95,
    n_soil_cores: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Run soil core simulation and return mean/std root length density profile.

    Args:
        dataset (RutheConfig): Configuration for the experiment.
        n_soil_cores (int | None): Number of virtual soil cores to sample per run.

    Returns:
        np.ndarray: Mean root length density profile from soil core samples.
    """

    # 1. load soil temperature data
    soil_temp_df = load_soil_temperatures(
        Path(dataset.PATH_TO_SOIL_TEMPERATURE_DATA),
        dataset.START_DATE,
        dataset.END_DATE,
    )

    # 2. initialize soil temperature profile and scaling function
    soil_temp = initialize_soil_temperature(soil_temp_df, 0)
    temp_scales = compute_temperature_scaling(soil_temp, dataset)

    # 3. set up plantbox elongation scaling function
    scale_elongation = create_scale_elongation_grid(temp_scales)

    # 4. define soil geometry (soil cores and simulation space)
    soil_cores, layer_volume, n_layers, depth = define_soil_core_geometry(
        dataset, n_soil_cores
    )

    # 5. run root system simulation
    root_length_density_profiles = simulate_root_growth(
        soil_temp_df,
        scale_elongation,
        soil_cores,
        layer_volume,
        n_layers,
        dataset,
    )

    # 6. aggregate results (across all individual cores)
    root_length_density_arrays = np.asarray(root_length_density_profiles)
    mean_root_length_density = root_length_density_arrays.mean(axis=0)
    std_root_length_density = root_length_density_arrays.std(axis=0)

    return mean_root_length_density, std_root_length_density

# necessary for 10+ cores to not clip out of the domain
def _select_centered_positions(candidates: np.ndarray, n_positions: int) -> np.ndarray:
    """Select safe positions closest to the field centre.

    If more positions are requested than unique safe candidates are available,
    candidates are reused. This allows large virtual sample sizes, e.g. 100 soil
    cores, without enlarging the simulated plant stand. Reused y-positions still
    become distinct soil cores because each core gets its own deterministic
    along-row x-position later in ``simulate_root_growth``.
    """

    if n_positions == 0:
        return np.empty(0, dtype=float)

    candidates = np.sort(np.asarray(candidates, dtype=float))
    if candidates.size == 0:
        raise ValueError("No safe soil-core candidate positions are available.")

    order = np.argsort(np.abs(candidates), kind="stable")
    centered_candidates = candidates[order]

    if n_positions <= centered_candidates.size:
        selected = centered_candidates[:n_positions]
    else:
        repeats, remainder = divmod(n_positions, centered_candidates.size)
        selected = np.concatenate(
            [
                np.tile(centered_candidates, repeats),
                centered_candidates[:remainder],
            ]
        )

    return np.sort(selected)


def _soil_core_y_offsets(dataset: RutheConfig, n_soil_cores: int) -> np.ndarray:
    """Choose safe across-row core centres, mixing on-row and inter-row positions.

    The root segments are later mapped into a periodic y-domain of width
    NUMBER_OF_ROWS_SOIL_CORE * INTER_ROW_DISTANCE. Therefore the full circular
    core has to fit inside that domain; otherwise boundary clipping creates
    artificially empty / low-RLD cores.
    """

    radius = dataset.SOIL_CORE_RADIUS
    inter_row_distance = dataset.INTER_ROW_DISTANCE
    n_rows = dataset.NUMBER_OF_ROWS_SOIL_CORE

    y_half_period = 0.5 * n_rows * inter_row_distance
    y_min = -y_half_period + radius
    y_max = y_half_period - radius
    if y_min >= y_max:
        raise ValueError(
            "Soil-core radius is too large for the periodic row domain. "
            "Increase NUMBER_OF_ROWS_SOIL_CORE or reduce SOIL_CORE_RADIUS."
        )

    row_positions = (
        inter_row_distance * np.arange(n_rows, dtype=float)
        - 0.5 * (n_rows - 1) * inter_row_distance
    )
    between_row_positions = row_positions[:-1] + 0.5 * inter_row_distance

    safe_row_positions = row_positions[
        (row_positions >= y_min) & (row_positions <= y_max)
    ]
    safe_between_row_positions = between_row_positions[
        (between_row_positions >= y_min) & (between_row_positions <= y_max)
    ]

    # For n > 1, keep both sampling situations represented. On-row gets the
    # extra core for odd n because real soil cores are often deliberately taken
    # on rows and because there is no between-row-only design here.
    n_on_row = 1 if n_soil_cores == 1 else (n_soil_cores + 1) // 2
    n_between_row = n_soil_cores - n_on_row

    on_row_offsets = _select_centered_positions(safe_row_positions, n_on_row)
    between_row_offsets = _select_centered_positions(
        safe_between_row_positions, n_between_row
    )

    return np.sort(np.concatenate([on_row_offsets, between_row_offsets]))


def define_soil_core_geometry(
    dataset: RutheConfig = RUTHE_1994_95,
    n_soil_cores: int | None = None,
):
    """
    Define the geometry of the soil cores.

    Args:
        dataset: Module containing global variables for a certain experiment.
        n_soil_cores: Number of virtual soil cores to generate. If omitted,
            dataset.NUMBER_OF_SOIL_CORE_SAMPLES is used.

    Returns:
        Tuple containing soil core geometries and related parameters.
            - soil_cores: list of soil core geometries.
            - layer_volume: Volume of each soil layer.
            - n_layers: Number of soil layers.
            - depth: Depth of the soil core.
    """

    n_soil_cores = (
        dataset.NUMBER_OF_SOIL_CORE_SAMPLES if n_soil_cores is None else n_soil_cores
    )
    if n_soil_cores < 1:
        raise ValueError("n_soil_cores must be at least 1")

    radius = dataset.SOIL_CORE_RADIUS
    depth = dataset.SOIL_CORE_DEPTH
    n_layers = dataset.NUMBER_OF_LAYERS

    soil_core = pb.SDF_PlantContainer(radius, radius, depth, False)
    y_offsets = _soil_core_y_offsets(dataset, n_soil_cores)
    soil_cores = [
        pb.SDF_RotateTranslate(soil_core, pb.Vector3d(0, float(y_offset), 0))
        for y_offset in y_offsets
    ]

    layer_volume = depth / n_layers * radius * radius * np.pi  # cm3

    return soil_cores, layer_volume, n_layers, depth


def simulate_root_growth(
    df: pd.DataFrame,
    scale_elongation: pb.EquidistantGrid1D,
    soil_cores: list,
    soil_layer_volume: float,
    n_layers: int,
    dataset: RutheConfig = RUTHE_1994_95,
) -> list[np.array]:
    """
    Simulate root growth for multiple plants and return per-core samples.

    Args:
        df (pd.DataFrame): DataFrame containing soil temperature data.
        scale_elongation (pb.EquidistantGrid1D): Scale elongation grid.
        soil_cores (list): Soil core geometries to sample.
        soil_layer_volume (float): Volume of each soil layer.
        n_layers (int): Number of soil layers.
        dataset: Module containing global variables for a certain experiment.

    Returns:
        list[np.array]: List of root length density arrays for each virtual soil core
            in each simulation run.
    """

    n_plants_per_row = dataset.NUMBER_OF_PLANTS_PER_ROW_SOIL_CORE
    n_rows = dataset.NUMBER_OF_ROWS_SOIL_CORE
    interrow_spacing = (
        dataset.NUMBER_OF_PLANTS_PER_ROW_SOIL_CORE * dataset.INTER_PLANT_DISTANCE
    )
    row_spacing = dataset.NUMBER_OF_ROWS_SOIL_CORE * dataset.INTER_ROW_DISTANCE
    n_runs = dataset.N_SIMULATION_RUNS

    root_length_densities = []

    x_period = float(interrow_spacing)
    x_half_period = x_period / 2.0
    x_min = -x_half_period + dataset.SOIL_CORE_RADIUS
    x_max = x_half_period - dataset.SOIL_CORE_RADIUS
    if x_min >= x_max:
        raise ValueError(
            "Soil-core radius is too large for the periodic in-row domain. "
            "Increase NUMBER_OF_PLANTS_PER_ROW_SOIL_CORE or reduce SOIL_CORE_RADIUS."
        )

    for run in range(n_runs):
        all_root_systems, all_analyzed_segments = [], None

        for plant_index in range(n_plants_per_row):
            for row_index in range(n_rows):
                root_system = initilize_root_system(
                    scale_elongation, plant_index, row_index, dataset, run
                )
                run_daily_growth_simulation(root_system, df, scale_elongation, dataset)
                all_root_systems.append(root_system)

                if plant_index + row_index == 0:
                    all_analyzed_segments = pb.SegmentAnalyser(
                        root_system
                    )  # initilizes segment analyser object
                else:
                    all_analyzed_segments.addSegments(root_system)

        # Sample deterministic x-positions for each core. Across-row positions are
        # fixed in define_soil_core_geometry() to include on-row and between-row cores.
        # The x-range is restricted so the full core remains inside the periodic box.
        shifted_soil_cores = []
        for core_index, soil_core in enumerate(soil_cores):
            seed_core = make_seed(
                dataset.BASE_SEED_SOIL_CORE,
                namespace=f"soil_core_sampling_core_{core_index + 1}",
                run=run,
                plant_index=0,
                row_index=0,
            )
            x_core = float(np.random.default_rng(seed_core).uniform(x_min, x_max))
            shifted_soil_cores.append(
                pb.SDF_RotateTranslate(soil_core, pb.Vector3d(x_core, 0, 0))
            )

        core_profiles = analyze_root_length_density(
            all_analyzed_segments,
            interrow_spacing,
            row_spacing,
            shifted_soil_cores,
            soil_layer_volume,
            n_layers,
            dataset,
        )

        # Keep all cores as individual samples to quantify core-to-core variability.
        root_length_densities.extend(core_profiles)

    return root_length_densities


def initilize_root_system(
    scale_elongation: pb.EquidistantGrid1D,
    plant_index: int,
    row_index: int,
    dataset: RutheConfig,
    run: int,
) -> pb.RootSystem:
    """
    Initialize a root system with specified parameters.

    Args:
        scale_elongation (pb.EquidistantGrid1D): Scale elongation grid.
        plant_index (int): Index of the plant in the row.
        row_index (int): Index of the row.
        dataset (RutheConfig): Configuration for the experiment.

    Returns:
        pb.RootSystem: Initialized root system.
    """

    seed = make_seed(
        dataset.BASE_SEED_PLANT,
        namespace="soil_core",
        run=run,
        plant_index=plant_index,
        row_index=row_index,
    )
    np.random.seed(seed)
    random.seed(seed)

    root_system = pb.RootSystem()
    root_system.readParameters(str(dataset.PATH_TO_PLANT_PARAMETERS))

    root_system.setSeed(seed)

    # assign elongation scaling
    for root_type in root_system.getRootRandomParameter():
        root_type.f_se = scale_elongation

    # position of the plant in the field
    root_system.getRootSystemParameter().seedPos = pb.Vector3d(
        dataset.INTER_PLANT_DISTANCE * plant_index
        - (
            dataset.INTER_PLANT_DISTANCE
            * dataset.NUMBER_OF_PLANTS_PER_ROW_SOIL_CORE
            / 2
        ),
        dataset.INTER_ROW_DISTANCE * row_index
        - ((dataset.NUMBER_OF_ROWS_SOIL_CORE - 1) * dataset.INTER_ROW_DISTANCE / 2),
        dataset.SOWING_DEPTH,
    )

    soil_space = pb.SDF_PlantContainer(
        dataset.SOIL_SPACE_PARAMS[0],
        dataset.SOIL_SPACE_PARAMS[1],
        dataset.SOIL_SPACE_PARAMS[2],
        dataset.SOIL_SPACE_PARAMS[3],
    )

    root_system.setGeometry(soil_space)

    root_system.initializeDB(4, 5, False)

    return root_system


def run_daily_growth_simulation(
    root_system: pb.RootSystem,
    df: pd.DataFrame,
    scale_elongation: pb.EquidistantGrid1D,
    dataset: RutheConfig,
) -> None:
    """
    Run daily growth simulation for the root system.

    Args:
        root_system (pb.RootSystem): The root system to simulate.
        df (pd.DataFrame): DataFrame containing soil temperature data.
        scale_elongation (pb.EquidistantGrid1D): Scale elongation grid.
        dataset: Module containing global variables for a certain experiment.
    """

    sim_time = dataset.SIMULATION_TIME
    dt = dataset.TIME_STEP

    soil_temp = np.zeros((14,))
    scales = np.ones_like(soil_temp)

    for time_step in range(round(sim_time / dt)):
        soil_temp = initialize_soil_temperature(df, time_step)

        scales = compute_temperature_scaling(soil_temp, dataset)
        scale_elongation.data = scales

        root_system.simulate(dt, False)


def analyze_root_length_density(
    all_analyzed_segments: pb.SegmentAnalyser,
    interrow: float,
    row: float,
    soil_cores: list,
    soil_layer_volume: float,
    n_layers: int,
    dataset: RutheConfig,
) -> list[np.array]:
    """
    Get root length density from the analyzed segments.

    Args:
        all_analyzed_segments (pb.SegmentAnalyser): Analyzed segments from the root systems.
        interrow (float): Distance between rows.
        row (float): Distance between plants in a row.
        soil_cores (list[float]): Geometries of the soil cores to sample.
        soil_layer_volume (float): Volume of each soil layer.
        n_layers (int): Number of soil layers.
        dataset (RutheConfig): Configuration for the experiment.

    Returns:
        list[np.array]: List containing root length density arrays for all soil cores.
    """

    root_length_densities = []

    for soil_core in soil_cores:
        analyzed_segments = pb.SegmentAnalyser(all_analyzed_segments)
        analyzed_segments.mapPeriodic(interrow, row)
        analyzed_segments.crop(soil_core)
        analyzed_segments.pack()

        distribution = analyzed_segments.distribution(
            "length", 0, -dataset.SOIL_CORE_DEPTH, n_layers, True
        )
        root_length_density = np.array(distribution) / soil_layer_volume  # cm/cm3
        root_length_densities.append(root_length_density)

    pprint(root_length_densities)

    return root_length_densities
