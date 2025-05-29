from vtk.util import numpy_support
import vtk
import numpy as np

import plotly.graph_objs as go


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
    tube_filter.SetNumberOfSides(5)
    tube_filter.SetRadius(1.0)
    tube_filter.SetCapping(True)
    tube_filter.Update()
    triangle_filter = vtk.vtkTriangleFilter()  # Convert triangle strips to regular triangles
    triangle_filter.SetInputConnection(tube_filter.GetOutputPort())
    triangle_filter.Update()

    return triangle_filter.GetOutput()


def generate_colorbar_image(vmin, vmax, colormap = "Viridis", height = 500, width = 100, discrete = False):

    if discrete:
        n = max(int(vmax - vmin + 1), 2)
        z = np.linspace(vmin - 0.5, vmax + 0.5, n).reshape(-1, 1)
    else:
        n = 256
        z = np.linspace(vmin, vmax, n).reshape(-1, 1)

    # print("vmin", vmin)
    # print("vmax", vmax)
    # print("n", n)

    if height > width:
        x0 = 0
        dx = 1
        y0 = vmin
        dy = (vmax - vmin) / (n - 1)
    else:
        y0 = 0
        dy = 1
        x0 = vmin
        dx = (vmax - vmin) / (n - 1)
        z = np.transpose(z)

    fig = go.Figure(go.Heatmap(
        z = z,
        colorscale = colormap,
        showscale = False,
        x0 = x0, dx = dx,
        y0 = y0, dy = dy,
        colorbar = None
    ))
    fig.update_layout(
        width = width,
        height = height,
        margin = dict(l = 10, r = 0, t = 0, b = 0),  # for the text
        yaxis = dict(
            showticklabels = False,
            showgrid = False,
            zeroline = False,
            visible = False
        )
    )
    return fig

