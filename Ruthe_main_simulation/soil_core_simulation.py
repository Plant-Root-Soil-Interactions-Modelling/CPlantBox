"""Run simulation and take virtual soil core samples to compare to real life data later"""

import random
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
    dataset: RutheConfig = RUTHE_1994_95
) -> tuple[np.ndarray, np.ndarray]:
    """Run soil core simulation and return mean/std root length density profile.

    Args:
        dataset (RutheConfig): Configuration for the experiment.

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
    soil_core_1, soil_core_2, layer_volume, n_layers, depth = define_soil_core_geometry(
        dataset
    )

    # 5. run root system simulation
    root_length_density_profiles = simulate_root_growth(
        soil_temp_df,
        scale_elongation,
        soil_core_1,
        soil_core_2,
        layer_volume,
        n_layers,
        dataset,
    )

    # 6. aggregate results (across all individual cores)
    root_length_density_arrays = np.asarray(root_length_density_profiles)
    mean_root_length_density = root_length_density_arrays.mean(axis=0)
    std_root_length_density = root_length_density_arrays.std(axis=0)

    return mean_root_length_density, std_root_length_density


def define_soil_core_geometry(dataset: RutheConfig = RUTHE_1994_95):
    """
    Define the geometry of the soil core.

    Args:
        dataset: Module containing global variables for a certain experiment.

    Returns:
        Tuple containing soil core geometries and related parameters.
            - soil_core_1: SDF_PlantContainer for the main soil core.
            - soil_core_2: SDF_RotateTranslate for the adjacent soil core.
            - layer_volume: Volume of each soil layer.
            - n_layers: Number of soil layers.
            - depth: Depth of the soil core.
    """

    radius = dataset.SOIL_CORE_RADIUS
    depth = dataset.SOIL_CORE_DEPTH
    n_layers = dataset.NUMBER_OF_LAYERS
    inter_row_distance = dataset.INTER_ROW_DISTANCE

    soil_core_1 = pb.SDF_PlantContainer(radius, radius, depth, False)
    soil_core_2 = pb.SDF_RotateTranslate(
        soil_core_1, pb.Vector3d(0, inter_row_distance + (inter_row_distance / 2), 0)
    )

    layer_volume = depth / n_layers * radius * radius * np.pi  # cm3

    return soil_core_1, soil_core_2, layer_volume, n_layers, depth


def simulate_root_growth(
    df: pd.DataFrame,
    scale_elongation: pb.EquidistantGrid1D,
    soil_core_1: pb.SDF_PlantContainer,
    soil_core_2: pb.SDF_RotateTranslate,
    soil_layer_volume: float,
    n_layers: int,
    dataset: RutheConfig = RUTHE_1994_95,
) -> list[np.array]:
    """
    Simulate root growth for multiple plants and return per-core samples.

    Args:
        df (pd.DataFrame): DataFrame containing soil temperature data.
        scale_elongation (pb.EquidistantGrid1D): Scale elongation grid.
        soil_core_1 (pb.SDF_PlantContainer): Main soil core geometry.
        soil_core_2 (pb.SDF_RotateTranslate): Adjacent soil core geometry.
        soil_layer_volume (float): Volume of each soil layer.
        n_layers (int): Number of soil layers.
        dataset: Module containing global variables for a certain experiment.

    Returns:
        list[np.array]: List of root length density arrays for each virtual soil core
            (2 cores per run).
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

        # Sample deterministic x-positions for each core (on-row and between-row).
        # This mimics cores being taken at different along-row locations while keeping
        # their across-row positions fixed.
        seed_core_1 = make_seed(
            dataset.BASE_SEED_SOIL_CORE,
            namespace="soil_core_sampling_core_1",
            run=run,
            plant_index=0,
            row_index=0,
        )
        seed_core_2 = make_seed(
            dataset.BASE_SEED_SOIL_CORE,
            namespace="soil_core_sampling_core_2",
            run=run,
            plant_index=0,
            row_index=0,
        )

        x_core_1 = float(np.random.default_rng(seed_core_1).uniform(-x_half_period, x_half_period))
        x_core_2 = float(np.random.default_rng(seed_core_2).uniform(-x_half_period, x_half_period))

        soil_core_1_shifted = pb.SDF_RotateTranslate(soil_core_1, pb.Vector3d(x_core_1, 0, 0))
        soil_core_2_shifted = pb.SDF_RotateTranslate(soil_core_2, pb.Vector3d(x_core_2, 0, 0))

        core_1, core_2 = analyze_root_length_density(
            all_analyzed_segments,
            interrow_spacing,
            row_spacing,
            soil_core_1_shifted,
            soil_core_2_shifted,
            soil_layer_volume,
            n_layers,
            dataset,
        )

        # Keep both cores as individual samples to quantify core-to-core variability.
        root_length_densities.append(core_1)
        root_length_densities.append(core_2)

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
        - (dataset.INTER_PLANT_DISTANCE * dataset.NUMBER_OF_PLANTS_PER_ROW_SOIL_CORE / 2),
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
    soil_core_1: float,
    soil_core_2: float,
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
        soil_core_1 (float): Geometry of the first soil core.
        soil_core_2 (float): Geometry of the second soil core.
        soil_layer_volume (float): Volume of each soil layer.
        n_layers (int): Number of soil layers.
        dataset (RutheConfig): Configuration for the experiment.

    Returns:
        list[np.array]: List containing root length density arrays for both soil cores.
    """

    analyzed_segments = pb.SegmentAnalyser(all_analyzed_segments)
    analyzed_segments.mapPeriodic(interrow, row)

    analyzed_segments.crop(soil_core_1)

    analyzed_segments.pack()

    distribution_core_1 = analyzed_segments.distribution(
        "length", 0, -dataset.SOIL_CORE_DEPTH, n_layers, True
    )
    root_length_density_core_1 = (
        np.array(distribution_core_1) / soil_layer_volume
    )  # cm/cm3

    analyzed_segments = pb.SegmentAnalyser(all_analyzed_segments)
    analyzed_segments.mapPeriodic(interrow, row)
    analyzed_segments.crop(soil_core_2)
    analyzed_segments.pack()
    distribution_core_2 = analyzed_segments.distribution(
        "length", 0, -dataset.SOIL_CORE_DEPTH, n_layers, True
    )
    root_length_density_core_2 = (
        np.array(distribution_core_2) / soil_layer_volume
    )  # cm/cm3

    return [root_length_density_core_1, root_length_density_core_2]
