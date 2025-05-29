""" methods creating the figures of the cplantbox webapp (D. Leitner, 2025)"""

import numpy as np

from dash import html, dcc, Input, Output, State, ctx
import dash_bootstrap_components as dbc

import plotly.graph_objs as go
from plotly.colors import qualitative

import dash_vtk

from vtk_conversions import *


def vtk3D_plot(vtk_data, color_pick, type_):

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

    if type_ == "Age":
        label_ = "Age [day]"
        cbar = generate_colorbar_image(vmin = color_range[0], vmax = color_range[1], colormap = "Jet_r", height = 40, width = 200)
    else:
        label_ = "Type [1]"
        cbar = generate_colorbar_image(vmin = color_range[0], vmax = color_range[1], colormap = "Jet", height = 40, width = 200, discrete = True)

    cbarcontent = dcc.Graph(
                # id = 'dummy-id',
                figure = cbar,
                style = {'width': '200px', 'height': '60px','marginTop': '10px'},
                config = {
                    'displayModeBar': False  # Hides the entire toolbar
                }
            )
    style_ = {
        'marginTop': '10px',
        'marginRight': '10px',
        'marginBottom': '10px',
        'marginLeft': '10px'
    }

    colbar = html.Div([html.H6(label_, style = style_), cbarcontent], style = {'display': 'flex', 'justifyContent': 'flex-end'})

    return [html.Div(content, style = {"width": "100%", "height": "600px"}), colbar]


def profile_plot(vtk_data):

    time = vtk_data["time"][-1]  # final simtime

    traces = []
    for i in range(0, 5):
        z_ = vtk_data[f"z{i}"]
        rld = vtk_data[f"rld{i}"]
        print(rld)
        traces.append(go.Scatter(x = rld, y = z_, mode = 'lines', name = "Day {:g}".format(time / 5.*(i + 1))))

    content = dcc.Graph(
            id = 'profile-plot',
            figure = {
                'data': traces,
                'layout': go.Layout(
                    # title='Line Plot of Multiple Arrays vs Time',
                    xaxis = {'title': 'Fraction of plant organ length per cm layer [-]'},
                    yaxis = {'title': 'Vertical position [cm]'},
                    hovermode = 'closest',
                    colorway = qualitative.Set2,
                ),
            },
            style = {'width': '100%', 'height': '600px'}
        )

    return  html.Div(content, style = {"width": "100%", "height": "600px"})


def dynamics_plot(vtk_data, typename_data):

    N = 25  # hard coded, see conversions.py
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
    traces = []
    traces.append(go.Scatter(x = time, y = leaf_length, mode = 'lines', name = "Leaf"))
    if number_s > 1:
        traces.append(go.Scatter(x = time, y = slength, mode = 'lines', name = "Stem (total)"))
    for j in range(number_s):  # sub types
        traces.append(go.Scatter(x = time, y = stem_length[j,:], mode = 'lines', name = typename_data[f"stem tab-{j+1}"]))
    if number_r > 1:
        traces.append(go.Scatter(x = time, y = rlength, mode = 'lines', name = "Root (total)"))
    for j in range(number_r):  # sub type
        traces.append(go.Scatter(x = time, y = root_length[j,:], mode = 'lines', name = typename_data[f"root tab-{j+1}"]))

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

