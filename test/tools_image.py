"""
Headless rendering helpers and image comparison tools for golden tests.
"""

from typing import Any, Iterable, List, Tuple, cast

import numpy as np
import vtk  # type: ignore

import plantbox as _pb  # type: ignore

pb = cast(Any, _pb)


def _vtk_points_from_nodes(nodes: Iterable) -> vtk.vtkPoints:
    nodes_np = np.array(list(map(np.array, nodes)))
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(nodes_np.shape[0])
    for i in range(nodes_np.shape[0]):
        p = nodes_np[i]
        points.SetPoint(i, float(p[0]), float(p[1]), float(p[2]))
    return points


def _vtk_lines_from_segments(segments: Iterable) -> vtk.vtkCellArray:
    segs_np = np.array(list(map(np.array, segments)))
    cells = vtk.vtkCellArray()
    for i in range(segs_np.shape[0]):
        s = segs_np[i]
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, int(s[0]))
        line.GetPointIds().SetId(1, int(s[1]))
        cells.InsertNextCell(line)
    return cells


def _add_cell_data(pd: vtk.vtkPolyData, name: str, values: List[float]) -> None:
    arr = vtk.vtkDoubleArray()
    arr.SetName(name)
    arr.SetNumberOfValues(len(values))
    for i, v in enumerate(values):
        arr.SetValue(i, float(v))
    pd.GetCellData().AddArray(arr)


def segs_to_polydata(rs, zoom_factor: float, param_names: List[str]) -> vtk.vtkPolyData:
    ana = (
        getattr(pb, "SegmentAnalyser")(rs)
        if not all(hasattr(rs, a) for a in ("nodes", "segments", "getParameter"))
        else rs
    )

    points = _vtk_points_from_nodes(ana.nodes)
    cells = _vtk_lines_from_segments(ana.segments)

    pd = vtk.vtkPolyData()
    pd.SetPoints(points)
    pd.SetLines(cells)

    for n in param_names:
        vals = list(ana.getParameter(n))
        if n == "radius":
            vals = [zoom_factor * float(v) for v in vals]
        if len(vals) == cells.GetNumberOfCells():
            _add_cell_data(pd, n, vals)

    c2p = vtk.vtkCellDataToPointData()
    c2p.SetPassCellData(True)
    c2p.SetInputData(pd)
    c2p.Update()
    return c2p.GetPolyDataOutput()


def _create_lookup_table(number_of_colors: int = 256) -> vtk.vtkLookupTable:
    color_series = vtk.vtkColorSeries()
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


def create_leaf_polygons(plant) -> vtk.vtkPolyData:
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

            if len(ln1) == 2 and len(ln2) == 2:
                offs = _add_quad(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = _add_quad(n1, ln1[1], ln2[1], n2, leaf_points, leaf_polys, offs)
            elif len(ln1) == 6 and len(ln2) == 6:
                offs = _add_quad(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = _add_quad(ln1[1], ln1[2], ln2[2], ln2[1], leaf_points, leaf_polys, offs)
                offs = _add_quad(n1, ln1[3], ln2[3], n2, leaf_points, leaf_polys, offs)
                offs = _add_quad(ln1[4], ln1[5], ln2[5], ln2[4], leaf_points, leaf_polys, offs)
            elif len(ln1) == 2 and len(ln2) == 6:
                x1 = leaf.getLeafVisX(i)
                x2 = leaf.getLeafVisX(i + 1)
                offs = _add_quad(n1, ln1[0], ln2[0], n2, leaf_points, leaf_polys, offs)
                offs = _add_quad(n1, ln1[1], ln2[3], n2, leaf_points, leaf_polys, offs)
                if x2[1] <= x1[0]:
                    offs = _add_quad(ln1[0], ln1[0], ln2[1], ln2[2], leaf_points, leaf_polys, offs)
                    offs = _add_quad(ln1[1], ln1[1], ln2[4], ln2[5], leaf_points, leaf_polys, offs)
            elif len(ln1) == 6 and len(ln2) == 2:
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


def render_headless_png(
    plant,
    p_name: str,
    image_path: str,
    width: int = 1200,
    height: int = 1000,
    background_rgb: Tuple[float, float, float] = (1.0, 1.0, 1.0),
) -> None:
    pd = segs_to_polydata(plant, zoom_factor=1.0, param_names=[p_name, "radius", "organType"])  # type: ignore[arg-type]
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

    cell_arr = tube.GetOutput().GetCellData().GetAbstractArray(p_name)
    if cell_arr is not None:
        rng = cell_arr.GetRange()
        if p_name == "organType":
            rng = (2.0, 4.0)
        lut.SetTableRange(rng)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    leaf_poly = create_leaf_polygons(plant)
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

    ren_win = vtk.vtkRenderWindow()
    ren_win.SetSize(width, height)
    ren_win.SetWindowName(p_name)
    ren_win.AddRenderer(renderer)
    ren_win.OffScreenRenderingOn()
    ren_win.Render()

    w2i = vtk.vtkWindowToImageFilter()
    w2i.SetInput(ren_win)
    w2i.SetScale(1)
    w2i.SetInputBufferTypeToRGB()
    w2i.ReadFrontBufferOff()
    w2i.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(image_path)
    writer.SetInputConnection(w2i.GetOutputPort())
    writer.Write()


def compare_images_png(path_a: str, path_b: str) -> float:
    """Return mean absolute difference normalized to [0,1] similarity = 1 - MAD.
    1.0 means identical images.
    """
    reader_a = vtk.vtkPNGReader()
    reader_a.SetFileName(path_a)
    reader_a.Update()
    img_a = reader_a.GetOutput()

    reader_b = vtk.vtkPNGReader()
    reader_b.SetFileName(path_b)
    reader_b.Update()
    img_b = reader_b.GetOutput()

    if img_a.GetDimensions() != img_b.GetDimensions():
        return 0.0

    diff = vtk.vtkImageMathematics()
    diff.SetOperationToSubtract()
    diff.SetInput1Data(img_a)
    diff.SetInput2Data(img_b)
    diff.Update()

    abs_img = vtk.vtkImageMathematics()
    abs_img.SetOperationToAbsoluteValue()
    abs_img.SetInputConnection(diff.GetOutputPort())
    abs_img.Update()

    image_to_float = vtk.vtkImageCast()
    image_to_float.SetOutputScalarTypeToDouble()
    image_to_float.SetInputConnection(abs_img.GetOutputPort())
    image_to_float.Update()

    stat = vtk.vtkImageAccumulate()
    stat.SetInputConnection(image_to_float.GetOutputPort())
    stat.Update()

    total = stat.GetOutput().GetPointData().GetScalars().GetTuple1(0)
    dims = img_a.GetDimensions()
    num_pix = dims[0] * dims[1] * max(1, img_a.GetNumberOfScalarComponents())
    mad = total / float(num_pix) / 255.0
    return float(max(0.0, 1.0 - mad))
