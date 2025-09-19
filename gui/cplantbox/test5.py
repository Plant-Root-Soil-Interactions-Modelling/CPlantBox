import dash
import dash_html_components as html
import dash_vtk
from dash_vtk.utils import to_mesh_state
import numpy as np

# Create a simple point cloud (or any other dataset)
x, y = np.mgrid[-5:5:20j, -5:5:20j]
z = np.sin(np.sqrt(x ** 2 + y ** 2))
points = np.c_[x.flatten(), y.flatten(), z.flatten()]
scalars = z.flatten()

# Convert to VTK mesh
mesh = to_mesh_state(points, scalars = [{
    "name": "elevation",
    "values": scalars,
    "location": "POINTS",
}])

app = dash.Dash(__name__)

app.layout = html.Div([
    dash_vtk.View([
        dash_vtk.GeometryRepresentation([
            dash_vtk.Mesh(state = mesh)
        ],
        colorMapPreset = 'Cool to Warm',
        colorDataRange = [scalars.min(), scalars.max()],
        scalarMode = 'default',
        colorByArrayName = 'elevation',
        lookupTable = 'lut'
        ),
        dash_vtk.ColorLegend(
            title = 'Elevation',
            visibility = True,
            position = [0.8, 0.05],
            size = [0.15, 0.6],
            lookupTable = 'lut'
        )
    ]),
])

if __name__ == '__main__':
    app.run_server(debug = True)
