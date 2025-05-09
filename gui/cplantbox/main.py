""" CPlantBox Webapp (using Python dash), D. Leitner 2025 """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import numpy as np
import vtk

import dash
from dash import html, dcc, Input, Output, State, ctx
import dash_bootstrap_components as dbc

from conversions import *  # auxiliary stuff
from plots import *  # figures

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
    dcc.Store(id = 'result-store', data = {}),
    dcc.Store(id = 'stem-store', data = {}),
    dcc.Store(id = 'leaf-store', data = {}),

    dbc.Row([
        # Panel 1
        dbc.Col([
            html.H5("Simulation"),
            html.H6("Plant"),
            dcc.Dropdown(id = 'plant-dropdown',
                        options = plants,
                        value = plants[0]['value'],
                        className = 'customDropdown'),
            html.Div(className = "spacer"),
            html.H6("Simulation time [day]"),
            dcc.Slider(id = 'time-slider', min = 1, max = 180, step = 1, value = 30,
                       marks = {1: "1", 180: "180"},
                       tooltip = { "always_visible": False}),
            html.Div(className = "spacer"),
            dbc.Button("Create", id = "create-button"),
            html.Div(className = "largeSpacer"),
            dcc.Loading(id = "loading-spinner",
                        type = "circle",  # Options: "default", "circle", "dot", "cube"
                        children = html.Div(id = "loading-spinner-output"))
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
            html.H5("Results"),
            dcc.Tabs(
                id = 'result-tabs',
                value = "VTK3D",
                children = [
                            dcc.Tab(label = '3D', value = 'VTK3D', className = 'tab', selected_className = 'tabSelected'),
                            dcc.Tab(label = '3D Age', value = 'VTK3DAge', className = 'tab', selected_className = 'tabSelected'),
                            dcc.Tab(label = '1D Profiles', value = 'Profile1D', className = 'tab', selected_className = 'tabSelected'),
                            dcc.Tab(label = 'Dynamics', value = 'Dynamics', className = 'tab', selected_className = 'tabSelected')
                        ]
                ),
            html.Div(className = "spacer"),
            html.Div(id = 'result-tabs-content')
        ], width = 6),
    ]),

    html.Div([
            html.Img(src = '/assets/cplantbox.png', className = 'logo'),
            html.Img(src = '/assets/fzj.png', className = 'logo'),
        ], className = 'logoContainer')

], fluid = True)

""" 1. LEFT - Simulation Panel """


@app.callback(# plant-dropdown
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


@app.callback(# create-button
    Output('result-tabs-content', 'children', allow_duplicate = True),
    Output('result-store', 'data'),
    Output('loading-spinner-output', 'children'),
    Input('create-button', 'n_clicks'),
    State('plant-dropdown', 'value'),
    State('time-slider', 'value'),
    State('seed-store', 'data'),
    State('root-store', 'data'),
    State('stem-store', 'data'),
    State('leaf-store', 'data'),
    State('result-store', 'data'),
    State('result-tabs', 'value'),
    prevent_initial_call = True,
)
def click_simulate(n_clicks, plant_value, time_slider_value, seed_data, root_data, stem_data, leaf_data, vtk_data, result_value):
    print("click_simulate()", plant_value)
    vtk_data = simulate_plant(plant_value, time_slider_value, seed_data, root_data, stem_data, leaf_data)
    content = render_result_tab(result_value, vtk_data)
    return (content, vtk_data, html.H6(""))


""" Parameters Panel"""


@app.callback(# organtype-tabs
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


#
#  Parameters Panel - Seed
#
def generate_seed_sliders(data):  # Generate sliders for seed tab from stored values
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


@app.callback(# dynamic-seed-slider
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


@app.callback(# shoot-checkbox
    Output({'type': 'dynamic-seed-slider', 'index': 0}, 'disabled'),
    Output('seed-store', 'data', allow_duplicate = True),
    Input('shoot-checkbox', 'value'),
    State('seed-store', 'data'),
    prevent_initial_call = True
)
def toggle_slider(checkbox_value, data):
    data["shoot-checkbox"] = 'agree' in checkbox_value
    return ('agree' not in checkbox_value, data)  # disable if not checked


@app.callback(# shoot-checkbox
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


#
# Parameters Panel - Root
#
def generate_root_sliders(root_values):  # Generate sliders for root tabs from stored values
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


# render_root_tab - Display sliders for the selected tab, using stored values
@app.callback(
    Output('root-tabs-content', 'children'),
    Input('root-tabs', 'value'),
    State('root-store', 'data'),
)
def render_root_tab(selected_tab, data):
    stored_values = data.get(selected_tab, ROOT_SLIDER_INITIALS)
    print("render_organtype_tab()", selected_tab, data)
    return generate_root_sliders(stored_values)


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


#
# Parameters Panel - Stem
#
def stem_layout(data):
    return html.Div([html.H5("not implemented (yet)")])


#
# Parameters Panel - Leaf
#
def leaf_layout(data):
    return html.Div([html.H5("not implemented (yet)")])


""" Results Panel """


@app.callback(
    Output('result-tabs-content', 'children', allow_duplicate = True),
    Input('result-tabs', 'value'),
    State('result-store', 'data'),
    prevent_initial_call = True,
)
def render_result_tab(tab, vtk_data):

    print("render_result_tab()", tab)

    if not vtk_data:
        print("no data")
        return html.Div([html.H6("press the create button")])

    if tab == 'VTK3D':
        color_pick = vtk_data["subType"]
        color_pick = np.repeat(color_pick, 24)  # 24 = 3*(7+1) ???
        print("number of cell colors", len(color_pick), "\n", type(color_pick))
        return vtk3D_plot(vtk_data, color_pick)
    elif tab == 'VTK3DAge':
        color_pick = vtk_data["creationTime"]
        color_pick = np.repeat(color_pick, 24)  # 24 = 3*(7+1) ???
        return vtk3D_plot(vtk_data, color_pick)
    elif tab == 'Profile1D':
        return profile_plot(vtk_data)
    elif tab == 'Dynamics':
        return dynamics_plot(vtk_data)


if __name__ == '__main__':
    app.run(debug = True)
