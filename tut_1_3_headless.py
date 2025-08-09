"""
Headless renderer for CPlantBox example (self-contained).

This script:
- Builds and simulates a plant from packaged parameters
- Renders the plant off-screen (no window) using VTK
- Saves a PNG image to results/example_plant_headless.png

It embeds a minimal, dependency-free version of plot_plant and helpers
required for off-screen rendering, so it does not rely on the local
vtk_plot.py or vtk_tools.py modules.
"""

import os
from typing import Any, Iterable, List, cast

import numpy as np
import vtk  # type: ignore

import plantbox as pb  # type: ignore

# Make type checkers happy: plantbox is a C++ extension without type hints
pb = cast(Any, pb)


def _vtk_points_from_nodes(nodes: Iterable) -> vtk.vtkPoints:
    """Create vtkPoints array from SegmentAnalyser nodes.

    We mimic vtk_tools.np_convert behavior to avoid relying on field names.
    """
    # Convert nodes to an Nx3 numpy array via array-coercion of pybind objects
    nodes_np = np.array(list(map(np.array, nodes)))
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(nodes_np.shape[0])
    for i in range(nodes_np.shape[0]):
        p = nodes_np[i]
        points.SetPoint(i, float(p[0]), float(p[1]), float(p[2]))
    return points


def _vtk_lines_from_segments(segments: Iterable) -> vtk.vtkCellArray:
    """Create vtkCellArray of lines from SegmentAnalyser segments.

    Important: Do NOT shift indices. The pybind objects already provide 0-based
    indices when coerced to numpy (matching vtk_tools.vtk_cells behavior).
    """
    segs_np = np.array(list(map(np.array, segments)))  # shape (N, 2)
    cells = vtk.vtkCellArray()
    for i in range(segs_np.shape[0]):
        s = segs_np[i]
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, int(s[0]))
        line.GetPointIds().SetId(1, int(s[1]))
        cells.InsertNextCell(line)
    return cells


def _add_cell_data(pd: vtk.vtkPolyData, name: str, values: List[float]) -> None:
    """Attach a list/array as vtk cell data with given name."""
    arr = vtk.vtkDoubleArray()
    arr.SetName(name)
    arr.SetNumberOfValues(len(values))
    for i, v in enumerate(values):
        arr.SetValue(i, float(v))
    pd.GetCellData().AddArray(arr)


def segs_to_polydata(rs, zoom_factor: float, param_names: List[str]) -> vtk.vtkPolyData:
    """Create vtkPolyData from a PlantBox object using line segments.

    - rs: PlantBox Plant/RootSystem or an existing SegmentAnalyser
    - zoom_factor: multiply the radius values (visual thickness)
    - param_names: names of segment parameters to attach as cell data
    """
    # Always use SegmentAnalyser to unify handling. Avoid direct attribute access for linters.
    if all(hasattr(rs, a) for a in ("nodes", "segments", "getParameter")):
        ana = rs
    else:
        ana = getattr(pb, "SegmentAnalyser")(rs)

    points = _vtk_points_from_nodes(ana.nodes)
    cells = _vtk_lines_from_segments(ana.segments)

    pd = vtk.vtkPolyData()
    pd.SetPoints(points)
    pd.SetLines(cells)

    # Attach requested segment parameters as cell data
    for n in param_names:
        vals = list(ana.getParameter(n))
        if n == "radius":
            vals = [zoom_factor * float(v) for v in vals]
        if len(vals) == cells.GetNumberOfCells():
            _add_cell_data(pd, n, vals)

    # Also provide those arrays as point data (tube radius uses point data)
    c2p = vtk.vtkCellDataToPointData()
    c2p.SetPassCellData(True)
    c2p.SetInputData(pd)
    c2p.Update()
    return c2p.GetPolyDataOutput()


def _create_lookup_table(number_of_colors: int = 256) -> vtk.vtkLookupTable:
    """Create a pleasant lookup table using a VTK color series."""
    color_series = vtk.vtkColorSeries()
    # Use a default scheme; can be tuned if needed
    color_series.SetColorScheme(vtk.vtkColorSeries.BREWER_DIVERGING_SPECTRAL_11)
    lut_src = color_series.CreateLookupTable(vtk.vtkColorSeries.ORDINAL)
    n = lut_src.GetNumberOfTableValues()

    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(number_of_colors)
    for i in range(number_of_colors):
        t = (n - 1) * (float(i) / number_of_colors)
        i0 = int(np.floor(t))
        theta = t - i0
        c0 = np.array(lut_src.GetTableValue(min(i0, n - 1)))
        c1 = np.array(lut_src.GetTableValue(min(i0 + 1, n - 1)))
        col = (1.0 - theta) * c0 + theta * c1
        lut.SetTableValue(i, col[0], col[1], col[2], col[3])
    return lut


def _add_quad(
    a,
    b,
    c,
    d,
    leaf_points: vtk.vtkPoints,
    leaf_polys: vtk.vtkCellArray,
    offs: int,
) -> int:
    """Insert a quad into the leaf mesh; returns new offset."""
    q = vtk.vtkPolygon()
    q.GetPointIds().SetNumberOfIds(4)
    for j in range(4):
        q.GetPointIds().SetId(j, offs + j)
    leaf_points.InsertNextPoint(float(a.x), float(a.y), float(a.z))
    leaf_points.InsertNextPoint(float(b.x), float(b.y), float(b.z))
    leaf_points.InsertNextPoint(float(c.x), float(c.y), float(c.z))
    leaf_points.InsertNextPoint(float(d.x), float(d.y), float(d.z))
    leaf_polys.InsertNextCell(q)
    return offs + 4


def _create_leaf_polygons(plant) -> vtk.vtkPolyData:
    """Reconstruct leaf surface polygons from PlantBox leaf organ visuals."""
    leaf_points = vtk.vtkPoints()
    leaf_polys = vtk.vtkCellArray()
    offs = 0

    leaves = plant.getOrgans(ot=getattr(pb, "leaf"))
    for leaf in leaves:
        for i in range(0, leaf.getNumberOfNodes() - 1):
            ln1 = leaf.getLeafVis(i)
            ln2 = leaf.getLeafVis(i + 1)
            if len(ln1) == 0 and len(ln2) == 0:
                continue
            n1 = leaf.getNode(i)
            n2 = leaf.getNode(i + 1)

            if len(ln1) == 2 and len(ln2) == 2:  # normal case
                offs = _add_quad(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = _add_quad(n1, ln1[1], ln2[1], n2, leaf_points, leaf_polys, offs)
            elif len(ln1) == 6 and len(ln2) == 6:  # convex case
                offs = _add_quad(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = _add_quad(ln1[1], ln1[2], ln2[2], ln2[1], leaf_points, leaf_polys, offs)
                offs = _add_quad(n1, ln1[3], ln2[3], n2, leaf_points, leaf_polys, offs)
                offs = _add_quad(ln1[4], ln1[5], ln2[5], ln2[4], leaf_points, leaf_polys, offs)
            elif len(ln1) == 2 and len(ln2) == 6:  # normal to convex
                x1 = leaf.getLeafVisX(i)
                x2 = leaf.getLeafVisX(i + 1)
                offs = _add_quad(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = _add_quad(n1, ln1[1], ln2[3], n2, leaf_points, leaf_polys, offs)
                if x2[1] <= x1[0]:
                    offs = _add_quad(ln1[0], ln1[0], ln2[1], ln2[2], leaf_points, leaf_polys, offs)
                    offs = _add_quad(ln1[1], ln1[1], ln2[4], ln2[5], leaf_points, leaf_polys, offs)
            elif len(ln1) == 6 and len(ln2) == 2:  # convex to normal
                x1 = leaf.getLeafVisX(i)
                x2 = leaf.getLeafVisX(i + 1)
                offs = _add_quad(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = _add_quad(n1, ln1[3], ln2[1], n2, leaf_points, leaf_polys, offs)
                if x1[1] <= x2[0]:
                    offs = _add_quad(ln1[1], ln1[2], ln2[0], ln2[0], leaf_points, leaf_polys, offs)
                    offs = _add_quad(ln1[4], ln1[5], ln2[1], ln2[1], leaf_points, leaf_polys, offs)

    poly = vtk.vtkPolyData()
    poly.SetPoints(leaf_points)
    poly.SetPolys(leaf_polys)
    return poly


def plot_plant_headless(
    plant,
    p_name: str = "age",
    image_path: str = "results/example_plant_headless.png",
    width: int = 1000,
    height: int = 1000,
    background_rgb=(1.0, 1.0, 1.0),
    zoom: float = 1.0,
) -> None:
    """Render the given plant off-screen and save a PNG image."""
    # Build polydata with parameters needed for rendering
    pd = segs_to_polydata(plant, zoom_factor=1.0, param_names=[p_name, "radius", "organType"])  # type: ignore[arg-type]

    # Tube filter uses point-data radius
    pd.GetPointData().SetActiveScalars("radius")
    tube = vtk.vtkTubeFilter()
    tube.SetInputData(pd)
    tube.SetNumberOfSides(12)
    tube.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tube.Update()

    lut = _create_lookup_table()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tube.GetOutputPort())
    mapper.ScalarVisibilityOn()
    mapper.SetScalarModeToUseCellFieldData()
    mapper.SelectColorArray(p_name)
    mapper.UseLookupTableScalarRangeOn()
    mapper.SetLookupTable(lut)
    mapper.Update()

    # Set LUT range from data if available
    cell_arr = tube.GetOutput().GetCellData().GetAbstractArray(p_name)
    if cell_arr is not None:
        rng = cell_arr.GetRange()
        # Keep organ type range stable if chosen
        if p_name == "organType":
            rng = (2.0, 4.0)
        lut.SetTableRange(rng)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Create leaf polygons and a surface actor (solid green)
    leaf_poly = _create_leaf_polygons(plant)
    leaf_mapper = vtk.vtkPolyDataMapper()
    leaf_mapper.SetInputData(leaf_poly)
    leaf_mapper.ScalarVisibilityOff()
    leaf_actor = vtk.vtkActor()
    leaf_actor.SetMapper(leaf_mapper)
    leaf_actor.GetProperty().SetColor(0.0, 0.5, 0.0)
    leaf_actor.GetProperty().SetOpacity(1.0)

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(*background_rgb)
    renderer.AddActor(actor)
    renderer.AddActor(leaf_actor)

    # Off-screen render window (create before camera setup to account for aspect ratio)
    ren_win = vtk.vtkRenderWindow()
    ren_win.SetSize(width, height)
    ren_win.SetWindowName(p_name)
    ren_win.AddRenderer(renderer)
    ren_win.OffScreenRenderingOn()

    # Camera setup (deterministic) with correct aspect
    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    bounds = tube.GetOutput().GetBounds()
    z_mid = 0.5 * (bounds[4] + bounds[5]) if bounds is not None else 0.0
    camera.ParallelProjectionOn()
    camera.SetFocalPoint(0.0, 0.0, z_mid)
    camera.SetPosition(200.0, 0.0, z_mid)
    camera.SetViewUp(0.0, 0.0, 1.0)
    camera.Azimuth(30)
    camera.Elevation(30)
    camera.OrthogonalizeViewUp()
    camera.SetClippingRange(1, 1000)

    # Start from ResetCamera-computed scale (accounts for aspect/orientation), then apply zoom
    base_scale = camera.GetParallelScale()
    eff_zoom = float(zoom) if isinstance(zoom, (int, float)) and zoom > 0 else 1.0
    camera.SetParallelScale(base_scale / eff_zoom)

    # Perform initial render
    ren_win.Render()

    # Capture buffer and save PNG
    w2i = vtk.vtkWindowToImageFilter()
    w2i.SetInput(ren_win)
    w2i.SetScale(1)  # increase for super-sampling
    w2i.SetInputBufferTypeToRGB()
    w2i.ReadFrontBufferOff()
    w2i.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(image_path)
    writer.SetInputConnection(w2i.GetOutputPort())
    writer.Write()


def main() -> None:
    # Use a non-zero deterministic seed (0 would trigger time-based seeding)
    plant = pb.Plant(1)  # pyright: ignore

    # Load model parameters from packaged data
    path = os.path.join(pb.data_path(), "structural", "plant") + "/"
    name = "fspm2023"
    plant.readParameters(path + name + ".xml")

    plant.initialize()
    simtime = 40  # days
    plant.simulate(simtime)

    os.makedirs("results", exist_ok=True)
    # Save golden under test directory for CI/golden testing
    out_dir = os.path.join("test", "golden")
    os.makedirs(out_dir, exist_ok=True)
    out_png = os.path.join(out_dir, "example_plant_headless.png")

    # Render and save PNG without opening any window
    plot_plant_headless(plant, p_name="age", image_path=out_png, width=1000, height=int(1000*2.2), zoom=3.0)
    print(f"Saved headless render to: {out_png}")


if __name__ == "__main__":
    main()
