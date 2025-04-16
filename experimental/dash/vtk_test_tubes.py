import dash
from dash import html
import dash_vtk
import vtk
import numpy as np


def create_polyline_with_radii():
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


def apply_tube_filter(polydata):
    tube_filter = vtk.vtkTubeFilter()
    tube_filter.SetInputData(polydata)
    tube_filter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tube_filter.SetNumberOfSides(20)
    tube_filter.SetRadius(1.0)
    tube_filter.SetCapping(True)
    tube_filter.Update()

    # Convert triangle strips to regular triangles
    triangle_filter = vtk.vtkTriangleFilter()
    triangle_filter.SetInputConnection(tube_filter.GetOutputPort())
    triangle_filter.Update()

    return triangle_filter.GetOutput()


def vtk_polydata_to_dashvtk_dict(polydata):
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


# Build geometry
polyline = create_polyline_with_radii()
tube = apply_tube_filter(polyline)
vtk_data = vtk_polydata_to_dashvtk_dict(tube)

# Dash app
app = dash.Dash(__name__)
app.layout = html.Div([
    html.Div("Polyline with Varying Radius via Tube Filter"),
    dash_vtk.View([
        dash_vtk.GeometryRepresentation([
            dash_vtk.PolyData(
                points = vtk_data["points"],
                polys = vtk_data["polys"]
            )
        ])
    ], style = {"height": "500px", "width": "500px", "border": "1px solid #ccc"})
])

if __name__ == "__main__":
    app.run(debug = True)
