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
app = dash.Dash(__name__, suppress_callback_exceptions = True, external_stylesheets = [dbc.themes.SANDSTONE])  # SANDSTONE, MINTY, MORPH
app.title = "CPlantbox Dash App"

param_names = get_parameter_names()
plants = [{'label': name[0], 'value': str(i)} for i, name in enumerate(param_names)]

tropisms_names = { "Plagiotropism": 0, "Gravitropism":1, "Exotropism": 2 }  # "Hydrotropism": 3

root_parameter_sliders = get_root_slider_names()
seed_parameter_sliders = get_seed_slider_names()

ROOT_SLIDER_INITIALS = [100, 3, 45, 1, 1, 7, 0.1, 1, 0.2, "Gravitropism"]
SEED_SLIDER_INITIALS = [180, 1, 7, 7, 20, 7, 7, 15, 5]

""" LAYOUT """
app.layout = dbc.Container([

    dcc.Store(id = 'seed-store', data = {"seed": SEED_SLIDER_INITIALS, "basal-checkbox": False, "shoot-checkbox": False, "tillers-checkbox": False }),
    dcc.Store(id = 'root-store', data = {f"tab-{i}": ROOT_SLIDER_INITIALS for i in range(1, 5)}),
    dcc.Store(id = 'root-typename-store', data = {f"tab-{i}": f"Order {i} root" for i in range(1, 5)}),
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

""" 1. LEFT - Simulation Panel """


# plant-dropdown
@app.callback(
    Output('seed-store', 'data'),  # , allow_duplicate=True
    Output('root-store', 'data'),
    Output('root-typename-store', 'data'),
    Output('organtype-tabs', 'value'),  # to trigger update
    Input('plant-dropdown', 'value'),
    State('seed-store', 'data'),
    State('root-store', 'data'),
    State('root-typename-store', 'data'),
    State('organtype-tabs', 'value'),
    # prevent_initial_call = True,
)
def update_plant(plant_value, seed_data, root_data, typename_data, tabs_value):
    set_data(plant_value, seed_data, root_data, typename_data)
    return (seed_data, root_data, typename_data, tabs_value)


# create-button
@app.callback(
    Output('panel3', 'children'),
    Input('create-button', 'n_clicks'),
    State('plant-dropdown', 'value'),
    State('time-slider', 'value'),
    State('seed-store', 'data'),
    State('root-store', 'data'),
    State('stem-store', 'data'),
    State('leaf-store', 'data'),
    prevent_initial_call = True,
)
def click_simulate(n_clicks, plant_value, time_slider_value, seed_data, root_data, stem_data, leaf_data):
    print("click_simulate()", plant_value)
    p_name = "creationTime"
    plant = simulate_plant(plant_value, time_slider_value, seed_data, root_data, stem_data, leaf_data)
    pd = vp.segs_to_polydata(plant, 1., ["radius", "organType", "creationTime", p_name])  # poly data
    tube_plot_actor, color_bar, lut = vp.plot_roots(pd, p_name, "", render = False, returnLut = True)

    vtk_data = vtk_polyline_to_dict(pd)

    geom_rep = dash_vtk.GeometryRepresentation([
            dash_vtk.PolyData(
                points = vtk_data["points"],
                lines = vtk_data["lines"]
                )
        ])

    content = dash_vtk.View(children = [ geom_rep ])

    return html.Div(content, style = {"width": "100%", "height": "600px"})


""" Parameters Panel"""


# organtype-tabs
@app.callback(
    Output('organtype-tabs-content', 'children'),
    Input('organtype-tabs', 'value'),
    State('seed-store', 'data'),
    State('root-store', 'data'),
    State('root-typename-store', 'data'),
    State('stem-store', 'data'),
    State('leaf-store', 'data'),
)
def render_organtype_tab(tab, seed_data, root_data, root_type_names, stem_data, leaf_data):
    if tab == 'Seed':
        print("render_organtype_tab() seed:", seed_data)
        return generate_seed_sliders(seed_data)
    elif tab == 'Root':
        print("render_organtype_tab() root:", root_data)
        return root_layout(root_data, root_type_names)
    elif tab == 'Stem':
        print("render_organtype_tab() stem:", stem_data)
        return stem_layout(stem_data)
    elif tab == 'Leaf':
        print("render_organtype_tab() leaf:", leaf_data)
        return leaf_layout(leaf_data)


""" Parameters Panel - Seed"""


# Generate sliders for seed tab from stored values
def generate_seed_sliders(data):
    seed_values = data["seed"]
    print("generate_seed_sliders()", seed_values)
    sliders = [
        html.Div(className = "spacer"),
        html.Div(className = "spacer"),
        ]
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
    if data["shoot-checkbox"]:
        v = ['agree']
    else:
        v = []
    sliders.insert(2, dcc.Checklist(
        id = 'shoot-checkbox',
        options = [{'label': html.Span(' Shoot borne roots', className = 'checkbox-label'), 'value': 'agree'}],
        className = 'checkbox',
        value = v
    ))
    if data["basal-checkbox"]:
        v = ['agree']
    else:
        v = []
    sliders.insert(3 + 2 * 2, dcc.Checklist(
        id = 'basal-checkbox',
        options = [{'label': html.Span(' Basal roots', className = 'checkbox-label'), 'value': 'agree'}],
        className = 'checkbox',
        value = v
    ))
    if data["tillers-checkbox"]:
        v = ['agree']
    else:
        v = []
    sliders.insert(4 + 5 * 2, dcc.Checklist(
        id = 'tillers-checkbox',
        options = [{'label': html.Span(' Tillers', className = 'checkbox-label'), 'value': 'agree'}],
        className = 'checkbox',
        value = v
    ))
    return html.Div(sliders)


# dynamic-seed-slider
@app.callback(
    Output('seed-store', 'data', allow_duplicate = True),
    Input({'type': 'dynamic-seed-slider', 'index': dash.ALL}, 'value'),
    State('seed-store', 'data'),
    prevent_initial_call = True,
)
def update_seed_store(slider_values, data):
    # print("update_seed_store()", store_data, slider_values, len(store_data["seed"]))
    if len(slider_values) > 1:  # called empty
        data["seed"] = slider_values
    return data


# checkbox
@app.callback(
    Output({'type': 'dynamic-seed-slider', 'index': 0}, 'disabled'),
    Output('seed-store', 'data', allow_duplicate = True),
    Input('shoot-checkbox', 'value'),
    State('seed-store', 'data'),
    prevent_initial_call = True
)
def toggle_slider(checkbox_value, data):
    data["shoot-checkbox"] = 'agree' in checkbox_value
    return ('agree' not in checkbox_value, data)  # disable if not checked


@app.callback(
    Output({'type': 'dynamic-seed-slider', 'index': 1}, 'disabled'),
    Input('shoot-checkbox', 'value')
)
def toggle_slider(checkbox_value):
    return 'agree' not in checkbox_value  # disable if not checked


@app.callback(
    Output({'type': 'dynamic-seed-slider', 'index': 2}, 'disabled'),
    Output('seed-store', 'data', allow_duplicate = True),
    Input('basal-checkbox', 'value'),
    State('seed-store', 'data'),
    prevent_initial_call = True
)
def toggle_slider(checkbox_value, data):
    data["basal-checkbox"] = 'agree' in checkbox_value
    return ('agree' not in checkbox_value, data)  # disable if not checked


@app.callback(
    Output({'type': 'dynamic-seed-slider', 'index': 3}, 'disabled'),
    Input('basal-checkbox', 'value')
)
def toggle_slider(checkbox_value):
    return 'agree' not in checkbox_value  # disable if not checked


@app.callback(
    Output({'type': 'dynamic-seed-slider', 'index': 4}, 'disabled'),
    Input('basal-checkbox', 'value')
)
def toggle_slider(checkbox_value):
    return 'agree' not in checkbox_value  # disable if not checked


@app.callback(
    Output({'type': 'dynamic-seed-slider', 'index': 5}, 'disabled'),
    Output('seed-store', 'data', allow_duplicate = True),
    Input('tillers-checkbox', 'value'),
    State('seed-store', 'data'),
    prevent_initial_call = True
)
def toggle_slider(checkbox_value, data):
    data["tillers-checkbox"] = 'agree' in checkbox_value
    return ('agree' not in checkbox_value, data)  # disable if not checked


@app.callback(
    Output({'type': 'dynamic-seed-slider', 'index': 6}, 'disabled'),
    Input('tillers-checkbox', 'value')
)
def toggle_slider(checkbox_value):
    return 'agree' not in checkbox_value  # disable if not checked


@app.callback(
    Output({'type': 'dynamic-seed-slider', 'index': 7}, 'disabled'),
    Input('tillers-checkbox', 'value')
)
def toggle_slider(checkbox_value):
    return 'agree' not in checkbox_value  # disable if not checked


""" Parameters Panel - Root"""


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


def root_layout(data, root_type_names):
    """ root tab layout: with root subTypes as sub tabs """
    rootPanelChildren = []
    for i, name in enumerate(root_type_names):
        ctab = dcc.Tab(label = root_type_names[f'tab-{i+1}'], value = f'tab-{i+1}',
                      className = 'tab', selected_className = 'tabSelected')
        rootPanelChildren.append(ctab)

    tab = dcc.Tabs(
        id = 'root-tabs',
        value = 'tab-1',
        children = rootPanelChildren,
    )
    content = html.Div(id = 'root-tabs-content')
    stored_values = data.get('tab-1', ROOT_SLIDER_INITIALS)
    return [tab, content]


# render_organtype_tab - Display sliders for the selected tab, using stored values
@app.callback(
    Output('root-tabs-content', 'children'),
    Input('root-tabs', 'value'),
    State('root-store', 'data'),
)
def render_organtype_tab(selected_tab, data):
    stored_values = data.get(selected_tab, ROOT_SLIDER_INITIALS)
    print("update_tab_content()", selected_tab, data)
    return generate_root_sliders(stored_values)


def stem_layout(data):
    return html.Div([html.H5("not implemented (yet)")])


def leaf_layout(data):
    return html.Div([html.H5("not implemented (yet)")])


# Update root-store when any slider changes
@app.callback(
    Output('root-store', 'data', allow_duplicate = True),  #
    Input({'type': 'dynamic-slider', 'index': dash.ALL}, 'value'),
    State('root-tabs', 'value'),
    State('root-store', 'data'),
    prevent_initial_call = True,
)
def update_root_store(slider_values, current_tab, store_data):
    # print("update_root_store()", current_tab, slider_values)
    store_data[current_tab] = slider_values
    return store_data


if __name__ == '__main__':
    app.run(debug = True)
