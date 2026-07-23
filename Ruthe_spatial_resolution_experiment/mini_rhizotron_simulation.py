"""Minirhizotron root system simulation based on Ruthe field data."""

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


def run_minirhizotron_simulation(
    dataset: RutheConfig = RUTHE_1994_95,
    only_touching_segments: bool = False,
) -> np.ndarray:
    """
    Simulate minirhizotron root measurements based on given root system parameters and dataset configuration

    Args:
        dataset (RutheConfig): Configuration parameters for the experiment.
        only_touching_segments (bool): If True, only consider root segments touching the tube. If False, segments within a viewing depth of 0.5 mm are considered.

    Returns:
        np.ndarray: Simulated minirhizotron root length density measurements.
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

    # 4. define soil domain and minirhizotron tube geometry
    (
        inner_tube_geom,
        outer_tube_geom,
        opening_left,
        opening_right,
        opening_left_offset,
        opening_right_offset,
        rhizotube_domain,
        inner_tubes_row_segment,
    ) = define_rhizotron_tube_geometry(dataset)

    # 5. run simulation and extract root scores
    root_score_data = []

    interrow_spacing = (
        dataset.NUMBER_OF_PLANTS_PER_ROW_MINIRHIZOTRON * dataset.INTER_PLANT_DISTANCE
    )
    row_spacing = dataset.NUMBER_OF_ROWS_MINIRHIZOTRON * dataset.INTER_ROW_DISTANCE

    for run in range(dataset.N_SIMULATION_RUNS):
        all_segments = simulate_root_growth(
            soil_temp_df, scale_elongation, rhizotube_domain, dataset, run
        )

        all_segments.mapPeriodic(interrow_spacing, row_spacing)

        if only_touching_segments:
            root_scores = only_touching_tube(
                dataset,
                all_segments,
                inner_tubes_row_segment,
                [
                    opening_left,
                    opening_right,
                    opening_left_offset,
                    opening_right_offset,
                ],
            )

        else:
            root_scores = extract_and_score_valid_segments(
                dataset,
                all_segments,
                outer_tube_geom,
                inner_tubes_row_segment,
                [
                    opening_left,
                    opening_right,
                    opening_left_offset,
                    opening_right_offset,
                ],
            )

        root_score_data.append(root_scores)

    # 6. process simulation results
    root_score_array = np.asarray(root_score_data)  # (run, day, depth, window)

    # filter: all runs, all measurement days, depths from 6th (depth 30cm) onward, all windows
    filtered_root_scores = root_score_array[:, :, 6:, :]

    root_scores_by_tube = []
    for run in range(filtered_root_scores.shape[0]):
        # average over windows to get two tubes
        mr_r = np.zeros(
            (2, len(dataset.DAYS_TO_ANALYZE_MINIRHIZOTRON), 18)
        )  # 2 tubes, measurement days,  18 depth intervals
        mr_r[0, :, :] = filtered_root_scores[run][:, :, :2].mean(
            axis=2
        )  # collapse 4 windows into two tubes
        mr_r[1, :, :] = filtered_root_scores[run][:, :, 2:].mean(
            axis=2
        )  # meaned over the windows
        root_scores_by_tube.append(mr_r)

    root_scores_by_tube = np.asarray(
        root_scores_by_tube
    )  # convert to array (runs, tubes, measurement days, depth intervals)

    mean_root_scores = root_scores_by_tube.mean(axis=(0, 1))  # (days, depth)
    root_score_deviation = root_scores_by_tube.std(axis=(0, 1))  # (days, depth)

    return mean_root_scores, root_score_deviation


def define_rhizotron_tube_geometry(
    dataset: RutheConfig = RUTHE_1994_95,
) -> tuple[
    pb.SDF_PlantContainer,
    pb.SDF_PlantContainer,
    pb.SDF_PlantContainer,
    pb.SDF_PlantContainer,
    pb.SDF_PlantContainer,
    pb.SDF_PlantContainer,
    pb.SDF_PlantContainer,
]:
    """
    Define the geometry of the minirhizotron tubes and openings using plantbox SDFs.

    Args:
        dataset (RutheConfig): Configuration parameters for the experiment.

    Returns:
        tuple: A tuple containing the SDF geometries for inner tube, outer tube, left
               opening, right opening, left offset opening, right offset opening, and rhizotube domain."""
    inner_tube = pb.SDF_PlantContainer(
        dataset.TUBE_INTERNAL_DIAMETER / 2,
        dataset.TUBE_INTERNAL_DIAMETER / 2,
        dataset.TUBE_LENGTH,
        False,
    )  # single tube
    tilted_inner_tube = pb.SDF_RotateTranslate(
        inner_tube, dataset.TUBE_ANGLE, pb.SDF_Axis.yaxis, pb.Vector3d(45, 0, 0)
    )  # rotated and translated
    inner_tubes_row_segment = pb.SDF_Union(
        tilted_inner_tube,
        pb.SDF_RotateTranslate(
            tilted_inner_tube, pb.Vector3d(0, dataset.INTER_ROW_DISTANCE / 2, 0)
        ),
    )  # second tube, half interrow distance away
    inner_tubes_upper_rows = pb.SDF_RotateTranslate(
        inner_tubes_row_segment,
        pb.Vector3d(
            0, dataset.INTER_ROW_DISTANCE * dataset.NUMBER_OF_ROWS_MINIRHIZOTRON, 0
        ),
    )
    inner_tubes_lower_rows = pb.SDF_RotateTranslate(
        inner_tubes_row_segment,
        pb.Vector3d(
            0, -dataset.INTER_ROW_DISTANCE * dataset.NUMBER_OF_ROWS_MINIRHIZOTRON, 0
        ),
    )

    full_inner_tube_array = pb.SDF_Union(
        pb.SDF_Union(inner_tubes_row_segment, inner_tubes_upper_rows),
        inner_tubes_lower_rows,
    )

    outer_tube = pb.SDF_PlantContainer(
        dataset.TUBE_EXTERNAL_DIAMETER / 2,
        dataset.TUBE_EXTERNAL_DIAMETER / 2,
        dataset.TUBE_LENGTH,
        False,
    )  # outer tube, (internal_diameter - outer_diameter)/2 = viewing depth
    tilted_outer_tube = pb.SDF_RotateTranslate(
        outer_tube, dataset.TUBE_ANGLE, pb.SDF_Axis.yaxis, pb.Vector3d(45, 0, 0)
    )  # rotated and translated
    outer_tubes_row_segment = pb.SDF_Union(
        tilted_outer_tube,
        pb.SDF_RotateTranslate(
            tilted_outer_tube, pb.Vector3d(0, dataset.INTER_ROW_DISTANCE / 2, 0)
        ),
    )  # two outer tubes

    opening_cylinder = pb.SDF_PlantContainer(
        dataset.TUBE_OPENING_DIAMETER / 2,
        dataset.TUBE_OPENING_DIAMETER / 2,
        dataset.TUBE_OPENING_HEIGHT,
        False,
    )
    opening_tilted = pb.SDF_RotateTranslate(
        opening_cylinder,
        90,
        pb.SDF_Axis.xaxis,
        pb.Vector3d(0, 0, dataset.INTER_IMAGE_DISTANCE),
    )
    opening_left = pb.SDF_RotateTranslate(
        opening_tilted, 45, pb.SDF_Axis.zaxis, pb.Vector3d(0, 0, 0)
    )
    opening_right = pb.SDF_RotateTranslate(
        opening_tilted, 135, pb.SDF_Axis.zaxis, pb.Vector3d(0, 0, 0)
    )
    opening_left_offset = pb.SDF_RotateTranslate(
        opening_tilted,
        45,
        pb.SDF_Axis.zaxis,
        pb.Vector3d(0, dataset.INTER_ROW_DISTANCE / 2, 0),
    )
    opening_right_offset = pb.SDF_RotateTranslate(
        opening_tilted,
        135,
        pb.SDF_Axis.zaxis,
        pb.Vector3d(0, dataset.INTER_ROW_DISTANCE / 2, 0),
    )

    simulation_box = pb.SDF_PlantBox(5000, 5000, 5000)
    simulation_box_shifted = pb.SDF_RotateTranslate(
        simulation_box, pb.Vector3d(0, 0, 2500)
    )

    rhizotube_domain = pb.SDF_Difference(simulation_box_shifted, full_inner_tube_array)

    return (
        full_inner_tube_array,
        outer_tubes_row_segment,
        opening_left,
        opening_right,
        opening_left_offset,
        opening_right_offset,
        rhizotube_domain,
        inner_tubes_row_segment,
    )


def simulate_root_growth(
    soil_temp_df: pd.DataFrame,
    scale_elongation: pb.EquidistantGrid1D,
    rhizotube_domain: pb.SDF_PlantContainer,
    dataset: RutheConfig,
    run: int,
) -> pb.SegmentAnalyser:
    n_plants_per_row = dataset.NUMBER_OF_PLANTS_PER_ROW_MINIRHIZOTRON
    n_rows = dataset.NUMBER_OF_ROWS_MINIRHIZOTRON

    all_root_systems, all_analyzed_segments = [], None

    for plant_index in range(n_plants_per_row):
        for row_index in range(n_rows):
            root_system = initialize_root_system(
                scale_elongation, plant_index, row_index, run, dataset, rhizotube_domain
            )

            run_daily_growth_simulation(
                root_system, soil_temp_df, scale_elongation, dataset
            )

            nodes = root_system.getNodes()

            if nodes:
                min_z = min([n.z for n in nodes])
                print(f"Plant {plant_index}-{row_index}, Depth reached: {min_z:.2f} cm")

            all_root_systems.append(root_system)
            if plant_index == 0 and row_index == 0:
                all_analyzed_segments = pb.SegmentAnalyser(root_system)
            else:
                all_analyzed_segments.addSegments(root_system)

    return all_analyzed_segments


def initialize_root_system(
    scale_elongation: pb.EquidistantGrid1D,
    plant_index: int,
    row_index: int,
    run: int,
    dataset: RutheConfig,
    rhizotube_domain: pb.SDF_PlantContainer,
) -> pb.RootSystem:
    seed = make_seed(
        dataset.BASE_SEED,
        namespace="minirhizotron",
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
            * dataset.NUMBER_OF_PLANTS_PER_ROW_MINIRHIZOTRON
            / 2
        ),
        dataset.INTER_ROW_DISTANCE * row_index
        - ((dataset.NUMBER_OF_ROWS_MINIRHIZOTRON - 1) * dataset.INTER_ROW_DISTANCE / 2),
        dataset.SOWING_DEPTH,
    )

    soil_space = pb.SDF_PlantContainer(
        dataset.SOIL_SPACE_PARAMS[0],
        dataset.SOIL_SPACE_PARAMS[1],
        dataset.SOIL_SPACE_PARAMS[2],
        dataset.SOIL_SPACE_PARAMS[3],
    )

    root_system.setGeometry(soil_space)
    root_system.setGeometry(rhizotube_domain)

    root_system.initializeDB(4, 5, False)

    return root_system


def run_daily_growth_simulation(
    root_system: pb.RootSystem,
    soil_temp_df: pd.DataFrame,
    scale_elongation: pb.EquidistantGrid1D,
    dataset: RutheConfig,
) -> None:
    sim_time = dataset.SIMULATION_TIME
    dt = dataset.TIME_STEP

    soil_temp = np.zeros((14,))
    scales = np.ones_like(soil_temp)

    for time_step in range(round(sim_time / dt)):
        soil_temp = initialize_soil_temperature(soil_temp_df, time_step)

        scales = compute_temperature_scaling(soil_temp, dataset)
        scale_elongation.data = scales

        root_system.simulate(dt, False)


def extract_and_score_valid_segments(
    dataset: RutheConfig,
    all_segments: pb.SegmentAnalyser,
    outer_tube_geom: pb.SDF_PlantContainer,
    inner_tubes_row_segment: pb.SDF_PlantContainer,
    openings,
) -> pd.DataFrame:
    measurements = []

    for measurement in range(len(dataset.DAYS_TO_ANALYZE_MINIRHIZOTRON)):
        segments_to_analyze = pb.SegmentAnalyser(all_segments)
        segments_to_analyze.filter(
            "creationTime", 0, dataset.DAYS_TO_ANALYZE_MINIRHIZOTRON[measurement]
        )
        segments_to_analyze.crop(
            pb.SDF_Difference(outer_tube_geom, inner_tubes_row_segment)
        )

        scores = np.zeros((24, 4))

        for window_index in range(4):
            for depth_index in range(24):
                v = pb.Vector3d(0, 0, depth_index * (dataset.INTER_IMAGE_DISTANCE))
                window_ = pb.SDF_RotateTranslate(openings[window_index], v)
                window = pb.SDF_RotateTranslate(
                    window_,
                    dataset.TUBE_ANGLE,
                    pb.SDF_Axis.yaxis,
                    pb.Vector3d(45, 0, 0),
                )
                analyzed_segments = pb.SegmentAnalyser(segments_to_analyze)
                analyzed_segments.crop(window)
                analyzed_segments.pack()

                df_seg = pd.DataFrame(
                    {
                        "ids": analyzed_segments.getParameter("id"),
                        "l": analyzed_segments.getParameter("length"),
                    }
                )

                df = df_seg.groupby("ids").sum()
                df = df.drop(df[df.l < 0.3].index)
                df["score"] = (dataset.BONITUR_SCORING_BASE) * round(
                    (df["l"] / 3) / dataset.BONITUR_SCORING_BASE
                )

                if (df["score"] > 0).any() and (df["score"] < 1).all():
                    score = 0.5
                elif df["score"].sum() > 5:
                    score = 5
                else:
                    score = df["score"].sum()

                scores[depth_index, window_index] = score

        measurements.append(scores)

    return measurements


def only_touching_tube(
    dataset: RutheConfig,
    all_segments: pb.SegmentAnalyser,
    inner_tubes_row_segment: pb.SDF_PlantContainer,
    openings: list,
) -> pd.DataFrame:
    measurements = []

    for measurement in range(len(dataset.DAYS_TO_ANALYZE_MINIRHIZOTRON)):
        segments_to_analyze = pb.SegmentAnalyser(all_segments)
        segments_to_analyze.filter(
            "creationTime", 0, dataset.DAYS_TO_ANALYZE_MINIRHIZOTRON[measurement]
        )

        scores = np.zeros((24, 4))

        for window_index in range(4):
            for depth_index in range(24):
                v = pb.Vector3d(0, 0, depth_index * (dataset.INTER_IMAGE_DISTANCE))
                window_ = pb.SDF_RotateTranslate(openings[window_index], v)
                window = pb.SDF_RotateTranslate(
                    window_,
                    dataset.TUBE_ANGLE,
                    pb.SDF_Axis.yaxis,
                    pb.Vector3d(45, 0, 0),
                )
                analyzed_segments = pb.SegmentAnalyser(segments_to_analyze)
                analyzed_segments.crop(window)
                analyzed_segments.pack()

                segments = analyzed_segments.segments
                nodes = analyzed_segments.nodes
                radii = analyzed_segments.getParameter("radius")
                length = analyzed_segments.getParameter("length")
                ids = analyzed_segments.getParameter("id")

                valid_ids = []
                valid_lengths = []

                for i, seg in enumerate(segments):
                    n1 = nodes[seg.x]
                    n2 = nodes[seg.y]

                    d1 = inner_tubes_row_segment.getDist(pb.Vector3d(n1.x, n1.y, n1.z))
                    d2 = inner_tubes_row_segment.getDist(pb.Vector3d(n2.x, n2.y, n2.z))

                    touching_nodes = 0

                    if d1 <= radii[i]:
                        touching_nodes += 1
                    if d2 <= radii[i]:
                        touching_nodes += 1

                    if touching_nodes == 2:
                        valid_ids.append(ids[i])
                        valid_lengths.append(length[i])

                    if touching_nodes == 1:
                        valid_ids.append(ids[i])
                        valid_lengths.append(length[i] * 0.5)

                if measurement == 0 and window_index == 0 and depth_index == 10:
                    print(
                        f"[Sanity Check] Time {measurement}, Win {window_index}: Total segments {len(segments)} -> Touching segments {len(valid_ids)}"
                    )

                df_seg = pd.DataFrame(
                    {
                        "ids": valid_ids,
                        "l": valid_lengths,
                    }
                )

                df = df_seg.groupby("ids").sum()
                df = df.drop(df[df.l < 0.3].index)
                df["score"] = (dataset.BONITUR_SCORING_BASE) * round(
                    (df["l"] / 3) / dataset.BONITUR_SCORING_BASE
                )

                if (df["score"] > 0).any() and (df["score"] < 1).all():
                    score = 0.5
                elif df["score"].sum() > 5:
                    score = 5
                else:
                    score = df["score"].sum()

                scores[depth_index, window_index] = score

        measurements.append(scores)

    return measurements
