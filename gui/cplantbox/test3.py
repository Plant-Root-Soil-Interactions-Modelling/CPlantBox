import dash
import dash_vtk
from dash import html
import numpy as np

# ─── 1) POINT + CELL SETUP ────────────────────────────────────────────────

# 10 points along X, sin(Y)
pts = np.array([[i, np.sin(i), 0] for i in range(10)], dtype = np.float32).flatten().tolist()
# 9 line‑segments: [2, i, i+1] for i in 0..8
lns = [c for i in range(9) for c in (2, i, i + 1)]

# 9 RGB colors in 0–255 (one per segment)
cols = np.array([
    [0, 0, 0],
    [  10, 255, 0],
    [  20, 0, 255],
    [255, 255, 0],
    [255, 0, 255],
    [  0, 255, 255],
    [128, 128, 128],
    [255, 128, 0],
    [  0, 128, 255],
], dtype = np.uint8).flatten().tolist()

# ─── 2) BUILD THE APP ────────────────────────────────────────────────────

app = dash.Dash(__name__)
app.layout = html.Div([
    dash_vtk.View([
        dash_vtk.GeometryRepresentation(
            # tell VTK what array to color by, and that it's per‑cell
            mapper = {
                "colorByArrayName": "Colors",
                "colorMode": "cell",
            },
            # ensure the lookup table covers 0–255
            colorDataRange = [0, 255],
            children = [
                dash_vtk.PolyData(
                    points = pts,
                    lines = lns,
                    children = [
                        dash_vtk.CellData([
                            dash_vtk.DataArray(
                                # 1) registration makes it active scalars
                                registration = "setScalars",
                                # 2) type must be a JS TypedArray
                                # type = "Uint8Array",
                                name = "Colors",
                                numberOfComponents = 3,
                                values = cols,
                            )
                        ])
                    ]
                )
            ]
        )
    ], style = {"height": "500px", "width": "100%"})
])

if __name__ == "__main__":
    app.run_server(debug = True)
