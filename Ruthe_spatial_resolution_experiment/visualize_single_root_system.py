import random
from pathlib import Path

import numpy as np
import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import vtk
from plantbox.visualisation.vtk_plot import (
    create_lookup_table,
    create_scalar_bar,
    segs_to_polydata,
)
from ruthe_global_variables import RUTHE_1994_95, RutheConfig
from simulation_utils import (
    compute_temperature_scaling,
    create_scale_elongation_grid,
    initialize_soil_temperature,
    load_soil_temperatures,
    make_seed,
)


def initialize_visualization_root_system(
    dataset: RutheConfig,
    scale_elongation: pb.EquidistantGrid1D,
    seed_pos: pb.Vector3d,
    domain: pb.SDF_PlantContainer,
    *,
    run: int,
    plant_index: int,
    row_index: int,
) -> pb.RootSystem:
    """Initialize one root system for the single-plant VTK visualization."""

    seed = make_seed(
        dataset.BASE_SEED_PLANT,
        namespace="virtual_field",
        run=run,
        plant_index=plant_index,
        row_index=row_index,
    )
    np.random.seed(seed)
    random.seed(seed)

    root_system = pb.RootSystem()
    root_system.readParameters(str(dataset.PATH_TO_PLANT_PARAMETERS))
    root_system.setSeed(seed)

    for root_type in root_system.getRootRandomParameter():
        root_type.f_se = scale_elongation

    root_system.getRootSystemParameter().seedPos = seed_pos
    root_system.setGeometry(domain)
    root_system.initializeDB(4, 5, False)

    return root_system


def run_visualization_growth_simulation(
    root_system: pb.RootSystem,
    soil_temp_df,
    scale_elongation: pb.EquidistantGrid1D,
    dataset: RutheConfig,
) -> None:
    """Run one root system with the current temperature-scaling workflow."""

    for time_step in range(round(dataset.SIMULATION_TIME / dataset.TIME_STEP)):
        soil_temp = initialize_soil_temperature(soil_temp_df, time_step)
        scale_elongation.data = compute_temperature_scaling(soil_temp, dataset)
        root_system.simulate(dataset.TIME_STEP, False)


def main():
    # 1. Setup configuration
    dataset = RUTHE_1994_95

    # Reduced simulation time for interactive viz
    dataset.SIMULATION_TIME = 200

    print("Setting up single-plant visualization...")

    # 2. Load environmental data and initialize the temperature scaling grid
    print("Loading environmental data...")
    soil_temp_df = load_soil_temperatures(
        Path(dataset.PATH_TO_SOIL_TEMPERATURE_DATA),
        dataset.START_DATE,
        dataset.END_DATE,
    )
    soil_temp = initialize_soil_temperature(soil_temp_df, 0)
    temp_scales = compute_temperature_scaling(soil_temp, dataset)
    scale_elongation = create_scale_elongation_grid(temp_scales)

    # 3. Simulate a single, unobstructed root system
    print(f"Running single plant simulation (SimTime={dataset.SIMULATION_TIME})...")
    soil_space = pb.SDF_PlantContainer(*dataset.SOIL_SPACE_PARAMS)
    seed_pos = pb.Vector3d(0, 0, dataset.SOWING_DEPTH)

    rs = initialize_visualization_root_system(
        dataset,
        scale_elongation,
        seed_pos,
        soil_space,
        run=0,
        plant_index=0,
        row_index=0,
    )
    run_visualization_growth_simulation(
        rs,
        soil_temp_df,
        scale_elongation,
        dataset,
    )

    # 4. Visualization (VTK): roots only
    print("Creating VTK plot...")
    all_roots_ana = pb.SegmentAnalyser(rs)

    print("Converting roots to polydata...")
    try:
        ct = all_roots_ana.getParameter("creationTime")
        print(f"Extracted creationTime with {len(ct)} entries.")
    except Exception as e:
        print(f"Error extracting creationTime: {e}")

    root_pd = segs_to_polydata(all_roots_ana, 1.0, ["creationTime", "radius"])

    # Force the active scalar to be creationTime
    scalar_on_points = False
    if root_pd.GetPointData().HasArray("creationTime"):
        root_pd.GetPointData().SetActiveScalars("creationTime")
        scalar_on_points = True
        print("Set creationTime as active scalar (Point Data).")
    elif root_pd.GetCellData().HasArray("creationTime"):
        root_pd.GetCellData().SetActiveScalars("creationTime")
        print("Set creationTime as active scalar (Cell Data).")
    else:
        print("WARNING: creationTime array NOT found in PolyData!")

    root_mapper = vtk.vtkPolyDataMapper()
    root_mapper.SetInputData(root_pd)
    root_mapper.ScalarVisibilityOn()
    if scalar_on_points:
        root_mapper.SetScalarModeToUsePointFieldData()
    else:
        root_mapper.SetScalarModeToUseCellFieldData()
    root_mapper.SelectColorArray("creationTime")

    lut = create_lookup_table()
    root_mapper.SetLookupTable(lut)
    root_mapper.SetScalarRange(0, dataset.SIMULATION_TIME)

    root_actor = vtk.vtkActor()
    root_actor.SetMapper(root_mapper)
    root_actor.GetProperty().SetOpacity(1.0)

    root_bar = create_scalar_bar(lut, None, "Creation Time (days)")

    print("Displaying window...")
    vp.render_window(
        [root_actor],
        "Single Plant Root System",
        [root_bar],
        root_pd.GetBounds(),
    ).Start()


if __name__ == "__main__":
    main()