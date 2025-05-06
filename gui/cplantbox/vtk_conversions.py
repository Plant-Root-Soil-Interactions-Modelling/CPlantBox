
import vtk
import numpy as np


def vtk_polyline_to_dict(polydata):
    """ converts polylines to a dict """
    points = polydata.GetPoints()
    lines = polydata.GetLines()

    n_points = points.GetNumberOfPoints()
    n_lines = lines.GetNumberOfCells()

    # Points array
    pts = np.array([points.GetPoint(i) for i in range(n_points)], dtype = np.float32)
    pts = pts.flatten().tolist()

    # Lines connectivity array
    lines.InitTraversal()
    id_list = vtk.vtkIdList()
    conn = []

    for _ in range(n_lines):
        lines.GetNextCell(id_list)
        conn.append(id_list.GetNumberOfIds())
        for j in range(id_list.GetNumberOfIds()):
            conn.append(id_list.GetId(j))

    return {
        "points": pts,
        "lines": conn,
    }


def vtk_polydata_to_dashvtk_dict(polydata):
    """ converts tubes to a dict """
    points = polydata.GetPoints()
    polys = polydata.GetPolys()

    n_points = points.GetNumberOfPoints()
    n_polys = polys.GetNumberOfCells()

    # Debug info
    print(f"Number of points: {n_points}, polys: {n_polys}")

    pts = np.array([points.GetPoint(i) for i in range(n_points)], dtype = np.float32)
    pts = pts.flatten().tolist()

    # Extract polygons
    polys.InitTraversal()
    id_list = vtk.vtkIdList()
    conn = []

    for _ in range(n_polys):
        polys.GetNextCell(id_list)
        conn.append(id_list.GetNumberOfIds())
        for j in range(id_list.GetNumberOfIds()):
            conn.append(id_list.GetId(j))

    return {
        "points": pts,
        "polys": conn,
    }


def create_polyline_with_radii(polydata):
    points = vtk.vtkPoints()
    radii = vtk.vtkFloatArray()
    radii.SetName("TubeRadius")

    num_pts = 20
    for i in range(num_pts):
        x = i * 0.2
        y = np.sin(i * 0.3)
        z = 0
        points.InsertNextPoint(x, y, z)
        radii.InsertNextValue(0.05 + 0.03 * np.cos(i * 0.5))  # varying radius

    lines = vtk.vtkCellArray()
    lines.InsertNextCell(num_pts)
    for i in range(num_pts):
        lines.InsertCellPoint(i)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    polydata.GetPointData().AddArray(radii)
    polydata.GetPointData().SetActiveScalars("TubeRadius")
    return polydata

