""" conversions regarding vtk data, e.g. vtkPolyData to dash store data and back, also colorbar generation, D. Leitner 2026  """ 

import base64 
import json
import zlib

import vtk

import numpy as np
import plotly.graph_objs as go
from vtk.util import numpy_support



def encode_array(arr: np.ndarray) -> str:
    """ numpy -> json string """
    payload = {
        "data": base64.b64encode(zlib.compress(arr.tobytes())).decode("ascii"),
        "dtype": str(arr.dtype),
        "shape": arr.shape,
    }
    return json.dumps(payload)


def decode_array(json_str: str) -> np.ndarray:
    """ json string -> numpy"""
    payload = json.loads(json_str)
    arr = np.frombuffer(
        zlib.decompress(base64.b64decode(payload["data"])),
        dtype=payload["dtype"],
    ).reshape(payload["shape"])
    return arr

def vtk_polydata_to_dashvtk_dict(polydata):
    """Converts vtkPolyData to a dictionary with optimized handling for large arrays."""
    points = polydata.GetPoints()
    polys = polydata.GetPolys()

    n_points = points.GetNumberOfPoints()
    n_polys = polys.GetNumberOfCells()
    print(f"vtk_polydata_to_dashvtk_dict(): Number of points: {n_points}, polys: {n_polys}")

    # Efficiently extract points to a numpy array
    pts_array = vtk.util.numpy_support.vtk_to_numpy(points.GetData()).astype(np.float16)

    # Efficiently extract polygon connectivity
    polys_array = vtk.util.numpy_support.vtk_to_numpy(polys.GetData()).astype(np.int32)

    vtk_data = {
        "points": encode_array(pts_array),
        "polys": encode_array(polys_array)
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

