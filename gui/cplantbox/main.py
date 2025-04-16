import sys; sys.path.append("../.."); sys.path.append("../../src/")

import numpy as np
import vtk

import dash
from dash import html, dcc, Input, Output, State, ctx
import dash_bootstrap_components as dbc
import dash_vtk
from dash_vtk.utils import to_mesh_state

from conversions import *

""" INITIALIZE """
app = dash.Dash(__name__, external_stylesheets = [dbc.themes.MINTY])
app.title = "CPlantbox Dash App"

param_names = get_parameter_names()
plants = [{'label': name[0], 'value': str(i)} for i, name in enumerate(param_names)]

tropisms_names = { "Plagiotropism": 0, "Gravitropism":1, "Exotropism": 2 }  # "Hydrotropism": 3

tabs_names = ["Tap", "First order", "Second order", "Seminal"]

root_parameter_sliders = get_root_slider_names()
seed_parameter_sliders = get_seed_slider_names()

ROOT_SLIDER_INITIALS = [100, 3, 45, 1, 1, 7, 0.1, 1, 0.2, "Gravitropism"]
SEED_SLIDER_INITIALS = [180, 1, 7, 7, 20, 7, 7, 15, 5]

""" LAYOUT """
app.layout = dbc.Container([

    dcc.Store(id = 'seed-store', data = {"seed": SEED_SLIDER_INITIALS}),
    dcc.Store(id = 'root-store', data = {tab: ROOT_SLIDER_INITIALS for tab in tabs_names}),
    dcc.Store(id = 'stem-store', data = {}),
    dcc.Store(id = 'leaf-store', data = {}),

    dbc.Row([
        # Panel 1
        dbc.Col([
            html.H5("Simulation"),
            html.H6("Plant"),
            dcc.Dropdown(
                id = 'plant-dropdown',
                options = plants,
                value = plants[0]['value'],
                className = 'customDropdown'
            ),
            html.Div(className = "spacer"),
            html.H6("Simulation time [day]"),
            dcc.Slider(id = 'time-slider', min = 1, max = 180, step = 1, value = 30,
                       marks = {1: "1", 180: "180"},
                       tooltip = { "always_visible": False}),
            html.Div(className = "spacer"),
            dbc.Button("Create", id = "create-button"),
        ], width = 2),

        # Panel 2
        dbc.Col([
            html.H5("Parameters"),
            dcc.Tabs(
                id = 'organtype-tabs',
                value = "Root",
                children = [
                            dcc.Tab(label = 'Seed', value = 'Seed', className = 'tab', selected_className = 'tabSelected'),
                            dcc.Tab(label = 'Root', value = 'Root', className = 'tab', selected_className = 'tabSelected'),
                            dcc.Tab(label = 'Stem', value = 'Stem', className = 'tab', selected_className = 'tabSelected'),
                            dcc.Tab(label = 'Leaf', value = 'Leaf', className = 'tab', selected_className = 'tabSelected')
                        ]
                ),
            html.Div(id = 'organtype-tabs-content')
        ], width = 3),

        # Panel 3
        dbc.Col([
            html.Div(className = "largeSpacer"),
            html.Div(id = 'panel3')
        ], width = 6),
    ]),

    html.Div([
            html.Img(src = '/assets/cplantbox.png', className = 'logo'),
            html.Img(src = '/assets/fzj.png', className = 'logo'),
        ], className = 'logoContainer')

], fluid = True)

""" Simulation Panel """


@app.callback(
    # Output('root-store', 'data'),
    Input('plant-dropdown', 'value'),
    State('root-store', 'data')
)
def update_plant(plant_value, data):
    debug_params(plant_value)
    return None


""" Parameters Panel"""


# Display sliders for the selected tab, using stored values
@app.callback(
    Output('organtype-tabs-content', 'children'),
    Input('organtype-tabs', 'value'),
    State('root-store', 'data'),
    State('seed-store', 'data'),
    State('stem-store', 'data'),
    State('leaf-store', 'data'),
)
def render_organtype_tab(tab, root_data, seed_data, stem_data, leaf_data):
    if tab == 'Seed':
        print("render_organtype_tab() seed:", seed_data)
        return seed_layout(seed_data)
    elif tab == 'Root':
        print("render_organtype_tab() root:", root_data)
        return root_layout(root_data)
    elif tab == 'Stem':
        print("render_organtype_tab() stem:", stem_data)
        return stem_layout(stem_data)
    elif tab == 'Leaf':
        print("render_organtype_tab() leaf:", leaf_data)
        return leaf_layout(leaf_data)


# Generate sliders for seed tab from stored values
def generate_seed_sliders(seed_values):
    sliders = [html.Div(className = "spacer"), html.Div(className = "spacer"), ]
    sliders.append(html.Div(className = "spacer"))
    for i, key in enumerate(seed_parameter_sliders.keys()):
        min_ = seed_parameter_sliders[key][0]
        max_ = seed_parameter_sliders[key][1]
        sliders.append(html.H6(key))
        sliders.append(
            dcc.Slider(
                        id = {'type': 'dynamic-seed-slider', 'index': i},
                        min = min_,
                        max = max_,
                        value = seed_values[i],
                        marks = {
                            min_: str(min_),
                            max_: str(max_)
                            },
                        tooltip = { "always_visible": False},
                    )
            )
    return html.Div(sliders)


# Generate sliders for root tabs from stored values
def generate_root_sliders(root_values):
    sliders = []
    sliders.append(html.Div(className = "spacer"))
    for i, key in enumerate(root_parameter_sliders.keys()):
        min_ = root_parameter_sliders[key][0]
        max_ = root_parameter_sliders[key][1]
        sliders.append(html.H6(key))
        sliders.append(
            dcc.Slider(
                        id = {'type': 'dynamic-slider', 'index': i},
                        min = min_,
                        max = max_,
                        value = root_values[i],
                        marks = {
                            min_: str(min_),
                            max_: str(max_)
                            },
                        tooltip = { "always_visible": False},
                    )
            )
    sliders.append(dcc.Dropdown(
                id = {'type': 'dynamic-slider', 'index': i},  # little white lie
                options = list(tropisms_names.keys()),
                value = root_values[i + 1],
                className = 'customDropdown'
            ))
    return html.Div(sliders)


def root_layout(data):
    # todo: derive number and label names from data
    rootPanelChildren = [dcc.Tab(label = f'{name}', value = f'tab-{i}', className = 'tab', selected_className = 'tabSelected') for i, name in enumerate(tabs_names)]
    tab = dcc.Tabs(
        id = 'root-tabs',
        value = 'tab-1',
        children = rootPanelChildren,
    )
    content = html.Div(id = 'root-tabs-content')
    stored_values = data.get('tab-1', ROOT_SLIDER_INITIALS)
    return [tab, content]


def seed_layout(data):
    return generate_seed_sliders(data["seed"])


def stem_layout(data):
    return html.Div([html.H5("not implemented (yet)")])


def leaf_layout(data):
    return html.Div([html.H5("not implemented (yet)")])


# render_organtype_tab - Display sliders for the selected tab, using stored values
@app.callback(
    Output('root-tabs-content', 'children'),
    Input('root-tabs', 'value'),
    State('root-store', 'data'),
    suppress_callback_exceptions = True
)
def render_organtype_tab(selected_tab, data):
    stored_values = data.get(selected_tab, ROOT_SLIDER_INITIALS)
    print("update_tab_content()", selected_tab, data)
    return generate_root_sliders(stored_values)


@app.callback(
    Output('seed-store', 'data', allow_duplicate = True),
    Input({'type': 'dynamic-seed-slider', 'index': dash.ALL}, 'value'),
    State('seed-store', 'data'),
    prevent_initial_call = True,
    suppress_callback_exceptions = True
)
def update_seed_store(slider_values, store_data):
    print("update_seed_store()", slider_values)
    if len(store_data) > 1:  # called empty
        store_data = slider_values
    return store_data


# Update root-store when any slider changes
@app.callback(
    Output('root-store', 'data', allow_duplicate = True),
    Input({'type': 'dynamic-slider', 'index': dash.ALL}, 'value'),
    State('root-tabs', 'value'),
    State('root-store', 'data'),
    prevent_initial_call = True,
    suppress_callback_exceptions = True
)
def update_root_store(slider_values, current_tab, store_data):
    print("update_root_store()", current_tab, slider_values)
    store_data[current_tab] = slider_values
    return store_data


# Create summary output on button press
@app.callback(
    Output('panel3', 'children'),
    Input('create-button', 'n_clicks'),
    State('root-tabs', 'value'),
    State('plant-dropdown', 'value'),
    State('time-slider', 'value'),
    State('root-store', 'data'),
    prevent_initial_call = True,
    suppress_callback_exceptions = True
)
def click_simulate(n_clicks, selected_tab, plant_value, time_slider_value, root_store_data):

    root_slider_values = root_store_data.get(selected_tab, [])

    p_name = "creationTime"
    plant = simulate_plant(plant_value, time_slider_value)
    pd = vp.segs_to_polydata(plant, 1., ["radius", "organType", "creationTime", p_name])  # poly data
    tube_plot_actor, color_bar, lut = vp.plot_roots(pd, p_name, "", render = False, returnLut = True)

    vtk_data = vtk_polyline_to_dict(pd)

    content = dash_vtk.View([
        dash_vtk.GeometryRepresentation([
            dash_vtk.PolyData(
                points = vtk_data["points"],
                lines = vtk_data["lines"]
                )
        ]),
    ])

    return html.Div(content, style = {"width": "100%", "height": "600px"})


if __name__ == '__main__':
    app.run(debug = True)
