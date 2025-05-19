from vtk.util import numpy_support
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
    """Converts vtkPolyData to a dictionary with optimized handling for large arrays."""
    points = polydata.GetPoints()
    polys = polydata.GetPolys()

    n_points = points.GetNumberOfPoints()
    n_polys = polys.GetNumberOfCells()
    print(f"Number of points: {n_points}, polys: {n_polys}")

    # Efficiently extract points to a numpy array
    pts_array = vtk.util.numpy_support.vtk_to_numpy(points.GetData()).astype(np.float32)
    pts = pts_array.flatten().tolist()

    # Efficiently extract polygon connectivity
    polys_data = vtk.util.numpy_support.vtk_to_numpy(polys.GetData())
    conn = polys_data.tolist()

    vtk_data = {
        "points": pts,
        "polys": conn,
    }

    return vtk_data


def apply_tube_filter(polydata):
    """ applies the tube filter """
    polydata.GetPointData().SetActiveScalars("radius")
    tube_filter = vtk.vtkTubeFilter()
    tube_filter.SetInputData(polydata)
    tube_filter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tube_filter.SetNumberOfSides(7)
    tube_filter.SetRadius(1.0)
    tube_filter.SetCapping(True)
    tube_filter.Update()
    triangle_filter = vtk.vtkTriangleFilter()  # Convert triangle strips to regular triangles
    triangle_filter.SetInputConnection(tube_filter.GetOutputPort())
    triangle_filter.Update()

    return triangle_filter.GetOutput()

