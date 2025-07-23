""" something working with colours """

import dash
from dash import html
import dash_vtk

app = dash.Dash(__name__)


def create_line(points, color):
    flat_points = [coord for pt in points for coord in pt]
    num_points = len(points)
    connectivity = [num_points] + list(range(num_points))

    return dash_vtk.GeometryRepresentation(
        children = [
            dash_vtk.PolyData(
                points = flat_points,
                lines = connectivity
            )
        ],
        property = {
            "color": color,
            "lineWidth": 3
        }
    )


# Define three lines with different colors
line1 = create_line([[0, 0, 0], [1, 0, 0], [1, 1, 0]], color = [1, 0, 0])  # Red
line2 = create_line([[0, 0, 0], [0, 1, 0], [0, 1, 1]], color = [0, 1, 0])  # Green
line3 = create_line([[0, 0, 0], [0, 0, 1], [1, 0, 1]], color = [0, 0, 1])  # Blue

app.layout = html.Div(
    dash_vtk.View(
        children = [line1, line2, line3],
        background = [1, 1, 1]  # White background
    ),
    style = {"width": "100vw", "height": "100vh"}
)

if __name__ == "__main__":
    app.run_server(debug = True)
