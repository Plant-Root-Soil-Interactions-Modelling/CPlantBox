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
from ruthe_global_variables import RUTHE_1994_95
from simulation_utils import (
    compute_temperature_scaling,
    create_scale_elongation_grid,
    initialize_soil_temperature,
    load_soil_temperatures,
    make_seed,
)
from soil_core_simulation import (
    define_soil_core_geometry,
    initilize_root_system,
    run_daily_growth_simulation,
)


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

    dataset = RUTHE_1994_95

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
    soil_core_geometries, _, _, _ = define_soil_core_geometry(dataset)

    n_plants_per_row = dataset.NUMBER_OF_PLANTS_PER_ROW_SOIL_CORE
    n_rows = dataset.NUMBER_OF_ROWS_SOIL_CORE
    interrow_spacing = (
        dataset.NUMBER_OF_PLANTS_PER_ROW_SOIL_CORE * dataset.INTER_PLANT_DISTANCE
    )
    run = 0

    x_period = float(interrow_spacing)
    x_half_period = x_period / 2.0

    all_root_systems = []

    for plant_index in range(n_plants_per_row):
        for row_index in range(n_rows):
            root_system = initilize_root_system(
                scale_elongation, plant_index, row_index, dataset, run
            )
            run_daily_growth_simulation(
                root_system, soil_temp_df, scale_elongation, dataset
            )
            all_root_systems.append(root_system)

            if plant_index + row_index == 0:
                all_analyzed_segments = pb.SegmentAnalyser(
                    root_system
                )  # initilizes segment analyser object
            else:
                all_analyzed_segments.addSegments(root_system)

    # Helpful sanity check: the 4 plants in a row only span ~6.8 cm with the
    # current INTER_PLANT_DISTANCE, so they can visually collapse into a single
    # "line" when looking at the full field extent.
    seed_positions = [
        (
            float(rs.getRootSystemParameter().seedPos.x),
            float(rs.getRootSystemParameter().seedPos.y),
            float(rs.getRootSystemParameter().seedPos.z),
        )
        for rs in all_root_systems
    ]
    unique_x = sorted({round(p[0], 6) for p in seed_positions})
    unique_y = sorted({round(p[1], 6) for p in seed_positions})
    print(
        "Planting grid seedPos: "
        f"{len(unique_x)} unique x over {max(unique_x) - min(unique_x):.2f} cm, "
        f"{len(unique_y)} unique y over {max(unique_y) - min(unique_y):.2f} cm"
    )

    shifted_soil_cores = []
    for core_index, soil_core_geometry in enumerate(soil_core_geometries):
        seed_core = make_seed(
            dataset.BASE_SEED_SOIL_CORE,
            namespace=f"soil_core_sampling_core_{core_index + 1}",
            run=run,
            plant_index=0,
            row_index=0,
        )
        x_core = float(
            np.random.default_rng(seed_core).uniform(-x_half_period, x_half_period)
        )
        shifted_soil_cores.append(
            pb.SDF_RotateTranslate(soil_core_geometry, pb.Vector3d(x_core, 0, 0))
        )

    # Combine obstacles
    # Soil cores are for visualization ONLY, they should not stop growth
    soil_cores = shifted_soil_cores[0]
    for shifted_soil_core in shifted_soil_cores[1:]:
        soil_cores = pb.SDF_Union(soil_cores, shifted_soil_core)

    # Compute field bounds
    plant_x_span = n_plants_per_row * dataset.INTER_PLANT_DISTANCE
    row_y_span = n_rows * dataset.INTER_ROW_DISTANCE

    field_bounds = [
        -plant_x_span / 2 - 20,
        plant_x_span / 2 + 20,
        -row_y_span / 2 - 20,
        row_y_span / 2 + 20,
        -150,
        10,
    ]

    # 5. Visualization (VTK)
    print("Creating VTK plot...")

    # 5.1 Roots
    # Mapping roots
    print("Converting roots to polydata...")
    # Explicitly check for parameters
    try:
        # Check if creationTime works on the analyser
        ct = all_analyzed_segments.getParameter("creationTime")
        print(f"Extracted creationTime with {len(ct)} entries.")
    except Exception as e:
        print(f"Error extracting creationTime: {e}")

    root_pd = segs_to_polydata(all_analyzed_segments, 1.0, ["creationTime", "radius"])

    # Force the active scalar to be creationTime
    if root_pd.GetPointData().HasArray("creationTime"):
        root_pd.GetPointData().SetActiveScalars("creationTime")
        print("Set creationTime as active scalar (Point Data).")
    elif root_pd.GetCellData().HasArray("creationTime"):
        root_pd.GetCellData().SetActiveScalars("creationTime")
        print("Set creationTime as active scalar (Cell Data).")
    else:
        print("WARNING: creationTime array NOT found in PolyData!")

    root_mapper = vtk.vtkPolyDataMapper()
    root_mapper.SetInputData(root_pd)
    root_mapper.ScalarVisibilityOn()
    root_mapper.SetScalarModeToUsePointFieldData()
    root_mapper.SelectColorArray("creationTime")

    # Lookup Table
    lut = create_lookup_table()
    root_mapper.SetLookupTable(lut)
    root_mapper.SetScalarRange(0, dataset.SIMULATION_TIME)  # Color by creation time

    root_actor = vtk.vtkActor()
    root_actor.SetMapper(root_mapper)
    root_actor.GetProperty().SetOpacity(0.3)  # Semi-transparent roots

    root_bar = create_scalar_bar(lut, None, "Creation Time (days)")

    # 5.1b Seed position markers (so the 4×6 grid is visible even if roots overlap)
    seed_points = vtk.vtkPoints()
    for x, y, z in seed_positions:
        seed_points.InsertNextPoint(x, y, z)

    seed_pd = vtk.vtkPolyData()
    seed_pd.SetPoints(seed_points)

    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(max(0.2, 0.25 * float(dataset.INTER_PLANT_DISTANCE)))
    sphere.SetThetaResolution(16)
    sphere.SetPhiResolution(16)

    glyph = vtk.vtkGlyph3D()
    glyph.SetInputData(seed_pd)
    glyph.SetSourceConnection(sphere.GetOutputPort())
    glyph.ScalingOff()
    glyph.Update()

    seed_mapper = vtk.vtkPolyDataMapper()
    seed_mapper.SetInputConnection(glyph.GetOutputPort())

    seed_actor = vtk.vtkActor()
    seed_actor.SetMapper(seed_mapper)
    seed_actor.GetProperty().SetColor(0.9, 0.1, 0.1)
    seed_actor.GetProperty().SetOpacity(1.0)

    # 5.2 Obstacles
    print("Meshing obstacles...")
    # Use the union of all obstacles for visualization
    obstacle_pd = sdf_to_vtk_mesh_custom(soil_cores, field_bounds, resolution=80)

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
        [root_actor, seed_actor, obs_actor],
        "Virtual Field Visualization",
        [root_bar],
        root_pd.GetBounds(),
    ).Start()


if __name__ == "__main__":
    main()
