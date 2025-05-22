""" methods creating the figures of the cplantbox webapp (D. Leitner, 2025)"""

import numpy as np

from dash import html, dcc, Input, Output, State, ctx
import plotly.graph_objs as go
from plotly.colors import qualitative
import dash_vtk


def vtk3D_plot(vtk_data, color_pick):
    """ 3D and 3D age plot """
    color_range = [np.min(color_pick), np.max(color_pick)]  # set range from min to max
    geom_rep = dash_vtk.GeometryRepresentation(
        mapper = {
                "colorByArrayName": "Colors",
                "colorMode": "cell",
                },
        colorDataRange = color_range,
        children = [
            dash_vtk.PolyData(
                points = vtk_data["points"],
                polys = vtk_data["polys"],
                children = [
                    dash_vtk.CellData([
                        dash_vtk.DataArray(
                            # 1) registration makes it active scalars
                            registration = "setScalars",
                            # 2) type must be a JS TypedArray
                            # type = "Uint8Array",
                            name = "Colors",
                            numberOfComponents = 1,
                            values = color_pick,
                        )
                    ])
                ]
                )
            ]
        )
    leaf_rep = dash_vtk.GeometryRepresentation(
        children = [
            dash_vtk.PolyData(
                points = vtk_data["leaf_points"],
                polys = vtk_data["leaf_polys"]
            )
        ],
        property = {
            "color": [0, 1, 0],  # Green
            "opacity": 0.85
        }
    )
    content = dash_vtk.View(children = [ geom_rep, leaf_rep ])
    return html.Div(content, style = {"width": "100%", "height": "600px"})


def profile_plot(vtk_data):

    z_ = vtk_data["z"]
    rld = vtk_data["rld"]
    rld = go.Scatter(x = rld, y = z_, mode = 'lines', name = "length fraction")
    traces = [rld]

    content = dcc.Graph(
            id = 'profile-plot',
            figure = {
                'data': traces,
                'layout': go.Layout(
                    # title='Line Plot of Multiple Arrays vs Time',
                    xaxis = {'title': 'Fraction [-]'},
                    yaxis = {'title': 'Depth [cm]'},
                    hovermode = 'closest',
                    colorway = qualitative.Set2,
                ),
            },
            style = {'width': '100%', 'height': '600px'}
        )

    return  html.Div(content, style = {"width": "100%", "height": "600px"})


def dynamics_plot(vtk_data):

    N = 50  # hard coded, see conversions.py
    number_r = vtk_data["number_r"]
    number_s = vtk_data["number_s"]
    time = vtk_data["time"]

    # fetch results
    root_length = np.zeros((number_r, N))
    stem_length = np.zeros((number_s, N))
    for j in range(number_r):
        root_length[j,:] = vtk_data[f"root_length-{j+1}"]
    for j in range(number_s):
        stem_length[j,:] = vtk_data[f"stem_length-{j+1}"]
    leaf_length = vtk_data["leaf_length"]
    rlength = np.sum(root_length, axis = 0)
    slength = np.sum(stem_length, axis = 0)

    # Create line traces
    go_leaf_length = go.Scatter(x = time, y = leaf_length, mode = 'lines', name = "leaf length")  # , line = { 'color': 'blue' }
    go_stem_length = go.Scatter(x = time, y = slength, mode = 'lines', name = "stem length")
    go_root_length = go.Scatter(x = time, y = rlength, mode = 'lines', name = "root length")
    go_stem_length_ = []
    for j in range(number_s):
        go_stem_length_.append(go.Scatter(x = time, y = stem_length[j,:], mode = 'lines', name = f"stem_length-{j+1}"))
    go_root_length_ = []
    for j in range(number_r):
        go_root_length_.append(go.Scatter(x = time, y = root_length[j,:], mode = 'lines', name = f"root_length-{j+1}"))

    traces = [go_leaf_length, go_stem_length] + go_stem_length_ + [go_root_length] + go_root_length_

    content = dcc.Graph(
            id = 'profile-plot',
            figure = {
                'data': traces,
                'layout': go.Layout(
                    # title='Line Plot of Multiple Arrays vs Time',
                    xaxis = {'title': 'Time [day]'},
                    yaxis = {'title': 'Length [cm]'},
                    hovermode = 'closest',
                    # hovertemplate='Time: %{x:.2f}<br>sin(t): %{y:.2f}<extra></extra>'
                    colorway = qualitative.Set2,
                )
            },
            style = {'width': '100%', 'height': '600px'}
        )

    return  html.Div(content, style = {"width": "100%", "height": "600px"})

