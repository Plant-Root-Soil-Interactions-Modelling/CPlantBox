"""methods creating the figures of the cplantbox webapp (D. Leitner, 2025)"""

import io

import dash_bootstrap_components as dbc
import dash_vtk
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from dash import Input, Output, State, ctx, dcc, html
from plotly.colors import qualitative
from vtk_conversions import *


def vtk_camera_preset(data, distance_scale=1.5):
    """Compute a deterministic camera setup so rotation can be reset as well."""
    points = decode_array(data["points"]).reshape(-1, 3)
    if points.size == 0:
        center = np.zeros(3)
        radius = 1.0
    else:
        center = points.mean(axis=0)
        radius = np.linalg.norm(points - center, axis=1).max()
        if radius <= 0:
            radius = 1.0
    distance = distance_scale * radius
    return {
        "cameraPosition": (center + np.array([-distance, -distance, distance])).tolist(),
        "cameraViewUp": [0, 0, 1],
    }


def vtk3D_plot(data, color_pick, type_, reset_camera_token, reset_camera, buttons):
    """3D and 3D age plot"""
    color_range = [np.min(color_pick), np.max(color_pick)]  # set range from min to max
    geom_rep = dash_vtk.GeometryRepresentation(
        mapper={
            "colorByArrayName": "Colors",
            "colorMode": "cell",
        },
        colorDataRange=color_range,
        children=[
            dash_vtk.PolyData(
                points=decode_array(data["points"]).flatten().tolist(),
                polys=decode_array(data["polys"]).tolist(),
                children=[
                    dash_vtk.CellData(
                        [
                            dash_vtk.DataArray(
                                # 1) registration makes it active scalars
                                registration="setScalars",
                                # 2) type must be a JS TypedArray
                                # type = "Uint8Array",
                                name="Colors",
                                numberOfComponents=1,
                                values=color_pick,
                            )
                        ]
                    )
                ],
            )
        ],
    )
    leaf_rep = dash_vtk.GeometryRepresentation(
        children=[dash_vtk.PolyData(points=decode_array(data["leaf_points"]).flatten().tolist(), polys=decode_array(data["leaf_polys"]).tolist())],
        property={"color": [0, 1, 0], "opacity": 0.85},  # Green
    )

    if reset_camera:
        camera_props = vtk_camera_preset(data)
        content = dash_vtk.View(children=[geom_rep, leaf_rep], triggerResetCamera=reset_camera_token, **camera_props)
    else:
        content = dash_vtk.View(children=[geom_rep, leaf_rep])

    if type_ == "Age":
        label_ = "Age [day]"
        cbar = generate_colorbar_image(vmin=color_range[0], vmax=color_range[1], colormap="Jet_r", height=40, width=200)
    else:
        label_ = "Type [1]"
        cbar = generate_colorbar_image(vmin=color_range[0], vmax=color_range[1], colormap="Jet", height=40, width=200, discrete=True)

    cbarcontent = dcc.Graph(
        # id = 'dummy-id',
        figure=cbar,
        style={"width": "200px", "height": "60px", "marginTop": "0px"},
        config={"displayModeBar": False},  # Hides the entire toolbar
    )
    style_ = {"marginTop": "0px", "marginRight": "10px", "marginBottom": "10px", "marginLeft": "10px"}

    colbar = html.Div([html.H6(label_, style=style_), cbarcontent], style={"display": "flex", "justifyContent": "flex-end"})

    return html.Div(
        [html.Div(className="spacer"), html.Div(content, style={"width": "100%", "height": "600px"}), buttons, colbar],
        style={"width": "100%", "display": "flex", "flex-direction": "column"},
    )


def profile_plot(data):
    """profile plots per organ"""
    traces = []
    for i in range(0, 5):
        z_ = decode_array(data[f"z{i}"])
        rld = decode_array(data[f"rld{i}"])
        t = data[f"time{i}"]
        traces.append(go.Scatter(x=rld, y=z_, mode="lines", name=f"Day {t}"))

    content = dcc.Graph(
        id="profile-plot",
        figure={
            "data": traces,
            "layout": go.Layout(
                xaxis={"title": "Fraction of plant organ length per cm layer [-]"},
                yaxis={"title": "Vertical position [cm]"},
                hovermode="closest",
                colorway=qualitative.Set2,
            ),
        },
        style={"position": "absolute", "top": 0, "left": 0, "width": "100%", "height": "100%"},
    )

    return html.Div(
        content,
        style={
            "position": "relative",
            "width": "60%",  # of parent width
            "padding-top": "80%",
            "margin-left": "auto",
            "margin-right": "auto",  # centers the container
        },
    )


def profile_to_excel(data):
    """
    Convert VTK data used for plotting into an Excel file and return it as a Dash download.

    Parameters
    ----------
    data : dict
        Dictionary containing z0..z4, rld0..rld4, time0..time4 arrays.

    Returns
    -------
    dash.dcc.Download
        A Dash download response containing the Excel file ``plant_profiles.xlsx``.
    """
    rows = []

    for i in range(5):
        z_ = decode_array(data[f"z{i}"])
        rld = decode_array(data[f"rld{i}"])
        t = data[f"time{i}"]
        for z_val, rld_val in zip(z_, rld):
            rows.append({"Day": t, "Vertical Position [cm]": z_val, "Fraction of Organ Length": rld_val})

    df = pd.DataFrame(rows)

    # --- write Excel to memory ---
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer) as writer:
        df.to_excel(writer, index=False, sheet_name="Dynamics")
    buffer.seek(0)

    # --- send to Dash download ---
    return dcc.send_bytes(buffer.getvalue(), "plant_profiles.xlsx")


def dynamics_plot(data, typename_data):
    """Creates a plot describeing the plant development by it's organs length over time"""
    number_r = data["number_r"]
    number_s = data["number_s"]
    time = decode_array(data["time"])
    N = int(time[-1])

    # fetch results
    root_length = np.zeros((number_r, N))
    stem_length = np.zeros((number_s, N))
    for j in range(number_r):
        root_length[j, :] = decode_array(data[f"root_length-{j+1}"])
    for j in range(number_s):
        stem_length[j, :] = decode_array(data[f"stem_length-{j+1}"])
    leaf_length = decode_array(data["leaf_length"])
    rlength = np.sum(root_length, axis=0)
    slength = np.sum(stem_length, axis=0)

    # Create line traces
    traces = []
    traces.append(go.Scatter(x=time, y=leaf_length, mode="lines", name="Leaf"))
    if number_s > 1:
        traces.append(go.Scatter(x=time, y=slength, mode="lines", name="Stem (total)"))
    for j in range(number_s):  # sub types
        traces.append(go.Scatter(x=time, y=stem_length[j, :], mode="lines", name=typename_data[f"stem tab-{j+1}"]))
    if number_r > 1:
        traces.append(go.Scatter(x=time, y=rlength, mode="lines", name="Root (total)"))
    for j in range(number_r):  # sub type
        traces.append(go.Scatter(x=time, y=root_length[j, :], mode="lines", name=typename_data[f"root tab-{j+1}"]))

    content = dcc.Graph(
        id="profile-plot",
        figure={
            "data": traces,
            "layout": go.Layout(
                # title='Line Plot of Multiple Arrays vs Time',
                xaxis={"title": "Time [day]"},
                yaxis={"title": "Length [cm]"},
                hovermode="closest",
                # hovertemplate='Time: %{x:.2f}<br>sin(t): %{y:.2f}<extra></extra>'
                colorway=qualitative.Set2,
            ),
        },
        style={"width": "100%", "height": "600px"},
    )

    return html.Div(content, style={"width": "100%", "height": "600px"})


import numpy as np
import pandas as pd


def dynamics_to_excel(data, typename_data):
    """
    Converts plant organ dynamics data into an Excel file.

    Parameters
    ----------
    data : dict
        Dictionary containing root_length, stem_length, leaf_length, time arrays, and number of roots/stems
    typename_data : dict
        Dictionary mapping 'root tab-{i}' and 'stem tab-{i}' to organ subtype names
    """
    number_r = data["number_r"]
    number_s = data["number_s"]
    time = decode_array(data["time"])
    N = int(time[-1])  # max time as int, or len(time)?

    # fetch results
    root_length = np.zeros((number_r, len(time)))
    stem_length = np.zeros((number_s, len(time)))
    for j in range(number_r):
        root_length[j, :] = decode_array(data[f"root_length-{j+1}"])
    for j in range(number_s):
        stem_length[j, :] = decode_array(data[f"stem_length-{j+1}"])
    leaf_length = decode_array(data["leaf_length"])

    # total lengths
    rlength = np.sum(root_length, axis=0)
    slength = np.sum(stem_length, axis=0)

    data_dict = {"Time [day]": time, "Leaf": leaf_length}  # create a dictionary to build the DataFrame

    if number_s > 1:
        data_dict["Stem Total"] = slength
    for j in range(number_s):
        col_name = typename_data.get(f"stem tab-{j+1}", f"Stem_{j+1}")
        data_dict[col_name] = stem_length[j, :]

    if number_r > 1:
        data_dict["Root Total"] = rlength
    for j in range(number_r):
        col_name = typename_data.get(f"root tab-{j+1}", f"Root_{j+1}")
        data_dict[col_name] = root_length[j, :]

    df = pd.DataFrame(data_dict)  # convert to DataFrame

    # --- write Excel to memory ---
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer) as writer:
        df.to_excel(writer, index=False, sheet_name="Dynamics")
    buffer.seek(0)

    # --- send to Dash download ---
    return dcc.send_bytes(buffer.getvalue(), "plant_dynamics.xlsx")
