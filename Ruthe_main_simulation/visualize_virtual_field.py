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
from soil_core_simulation import define_soil_core_geometry


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
    """Initialize one root system for the virtual-field VTK visualization."""

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


def sdf_to_vtk_mesh_custom(sdf, bounds, resolution=60):
    """
    Converts certain SDF to mesh using explicit bounds.
    bounds: [xmin, xmax, ymin, ymax, zmin, zmax]
    """
    xmin, xmax, ymin, ymax, zmin, zmax = bounds

    # Create vtkImageData
    image_data = vtk.vtkImageData()
    image_data.SetDimensions(resolution, resolution, resolution)
    image_data.SetSpacing(
        (xmax - xmin) / (resolution - 1),
        (ymax - ymin) / (resolution - 1),
        (zmax - zmin) / (resolution - 1),
    )
    image_data.SetOrigin(xmin, ymin, zmin)

    # Evaluate SDF
    scalars = vtk.vtkDoubleArray()
    scalars.SetNumberOfComponents(1)
    scalars.SetNumberOfTuples(resolution**3)

    # Python Loop - slow but functional
    count = 0
    dx = (xmax - xmin) / (resolution - 1)
    dy = (ymax - ymin) / (resolution - 1)
    dz = (zmax - zmin) / (resolution - 1)

    print(f"Meshing SDF with resolution {resolution}...")

    for k in range(resolution):
        z = zmin + k * dz
        for j in range(resolution):
            y = ymin + j * dy
            for i in range(resolution):
                x = xmin + i * dx
                val = sdf.getDist(pb.Vector3d(x, y, z))
                scalars.SetTuple1(count, val)
                count += 1

    image_data.GetPointData().SetScalars(scalars)

    # Contour Filter (Marching Cubes)
    contour = vtk.vtkContourFilter()
    contour.SetInputData(image_data)
    contour.SetValue(0, 0.0)
    contour.Update()

    return contour.GetOutput()


def main():
    # 1. Setup Configuration
    dataset = RUTHE_1994_95
    dataset.NUMBER_OF_ROWS_SOIL_CORE = 8
    dataset.NUMBER_OF_PLANTS_PER_ROW_SOIL_CORE = 100
    dataset.NUMBER_OF_ROWS_MINIRHIZOTRON = 8
    dataset.NUMBER_OF_PLANTS_PER_ROW_MINIRHIZOTRON = 100
    dataset.SOWING_DEPTH = 0.0

    # Reduced simulation time for interactive viz
    dataset.SIMULATION_TIME = 200

    print("Setting up virtual field simulation...")

    # 2. Geometry Definitions

    # Soil Cores (Usually at 0,0 and offsets in the row direction)
    soil_core_geometries, _, _, _ = define_soil_core_geometry(dataset)

    # Minirhizotrons: use the same tube dimensions/orientation as the current
    # minirhizotron simulation. The tubes are obstacles for root growth.
    tube_obstacle = pb.SDF_PlantContainer(
        dataset.TUBE_INTERNAL_DIAMETER / 2,
        dataset.TUBE_INTERNAL_DIAMETER / 2,
        dataset.TUBE_LENGTH,
        False,
    )
    tilted_tube = pb.SDF_RotateTranslate(
        tube_obstacle,
        dataset.TUBE_ANGLE,
        pb.SDF_Axis.yaxis,
        pb.Vector3d(45, 0, 0),
    )
    rhizotube_a_b_o = pb.SDF_Union(
        tilted_tube,
        pb.SDF_RotateTranslate(
            tilted_tube, pb.Vector3d(0, dataset.INTER_ROW_DISTANCE / 2, 0)
        ),
    )

    # Shift minirhizotrons to X=50 (Further into field, no overlap with cores at X=0)
    # Field X range: ~ -70 to +70. X=50 is inside.
    shift_vector = pb.Vector3d(50, 0, 0)
    minirhizo_shifted = pb.SDF_RotateTranslate(rhizotube_a_b_o, shift_vector)

    # Combine obstacles
    # Soil cores are for visualization ONLY, they should not stop growth
    soil_cores = soil_core_geometries[0]
    for soil_core_geometry in soil_core_geometries[1:]:
        soil_cores = pb.SDF_Union(soil_cores, soil_core_geometry)

    # Minirhizotubes ARE obstacles for simulation
    simulation_obstacles = minirhizo_shifted

    # Visualization needs everything
    visualization_obstacles = pb.SDF_Union(soil_cores, simulation_obstacles)

    # 3. Simulation Setup
    print("Loading environmental data...")
    soil_temp_df = load_soil_temperatures(
        Path(dataset.PATH_TO_SOIL_TEMPERATURE_DATA),
        dataset.START_DATE,
        dataset.END_DATE,
    )
    soil_temp = initialize_soil_temperature(soil_temp_df, 0)
    temp_scales = compute_temperature_scaling(soil_temp, dataset)
    scale_elongation = create_scale_elongation_grid(temp_scales)

    # 4. Run Simulation
    n_rows = dataset.NUMBER_OF_ROWS_SOIL_CORE
    n_plants = dataset.NUMBER_OF_PLANTS_PER_ROW_SOIL_CORE
    print(
        f"Running simulations for {n_rows} rows x {n_plants} plants (SimTime={dataset.SIMULATION_TIME})..."
    )

    # Compute field bounds
    plant_x_span = n_plants * dataset.INTER_PLANT_DISTANCE
    row_y_span = n_rows * dataset.INTER_ROW_DISTANCE

    field_bounds = [
        -plant_x_span / 2 - 20,
        plant_x_span / 2 + 20,
        -row_y_span / 2 - 20,
        row_y_span / 2 + 20,
        -150,
        10,
    ]

    # Domain with obstacles - Use ONLY simulation_obstacles (tubes)
    sim_box = pb.SDF_PlantBox(plant_x_span * 2, row_y_span * 2, 300)
    domain_with_obstacles = pb.SDF_Difference(sim_box, simulation_obstacles)

    all_roots_ana = pb.SegmentAnalyser()
    all_root_systems = []  # Keep references alive!

    for j in range(n_rows):
        print(f"Simulating Row {j + 1}/{n_rows}...")
        for i in range(n_plants):
            x_pos = dataset.INTER_PLANT_DISTANCE * i - (
                dataset.INTER_PLANT_DISTANCE * n_plants / 2
            )
            y_pos = dataset.INTER_ROW_DISTANCE * j - (
                (n_rows - 1) * dataset.INTER_ROW_DISTANCE / 2
            )
            z_pos = dataset.SOWING_DEPTH
            seed_pos = pb.Vector3d(x_pos, y_pos, z_pos)

            rs = initialize_visualization_root_system(
                dataset,
                scale_elongation,
                seed_pos,
                domain_with_obstacles,
                run=0,
                plant_index=i,
                row_index=j,
            )
            run_visualization_growth_simulation(
                rs,
                soil_temp_df,
                scale_elongation,
                dataset,
            )

            # Add to combined analyser
            all_root_systems.append(rs)  # Store it
            all_roots_ana.addSegments(rs)

    # 5. Visualization (VTK)
    print("Creating VTK plot...")

    # 5.1 Roots
    # Mapping roots
    print("Converting roots to polydata...")
    # Explicitly check for parameters
    try:
        # Check if creationTime works on the analyser
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

    # Lookup Table
    lut = create_lookup_table()
    root_mapper.SetLookupTable(lut)
    root_mapper.SetScalarRange(0, dataset.SIMULATION_TIME)  # Color by creation time

    root_actor = vtk.vtkActor()
    root_actor.SetMapper(root_mapper)
    root_actor.GetProperty().SetOpacity(0.3)  # Semi-transparent roots

    root_bar = create_scalar_bar(lut, None, "Creation Time (days)")

    # 5.2 Obstacles
    print("Meshing obstacles...")
    # Use the union of all obstacles for visualization
    obstacle_pd = sdf_to_vtk_mesh_custom(
        visualization_obstacles, field_bounds, resolution=80
    )

    obs_mapper = vtk.vtkPolyDataMapper()
    obs_mapper.SetInputData(obstacle_pd)
    obs_actor = vtk.vtkActor()
    obs_actor.SetMapper(obs_mapper)
    obs_actor.GetProperty().SetColor(0.6, 0.4, 0.2)
    obs_actor.GetProperty().SetOpacity(0.5)  # Semi-transparent

    # 5.3 Render
    print("Displaying window...")
    # Passing [] for scalar bar to hide it, or root_bar
    vp.render_window(
        [root_actor, obs_actor],
        "Virtual Field Visualization",
        [root_bar],
        root_pd.GetBounds(),
    ).Start()


if __name__ == "__main__":
    main()