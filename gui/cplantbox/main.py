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

seed_parameter_sliders = get_seed_slider_names()
root_parameter_sliders = get_root_slider_names()
stem_parameter_sliders = get_stem_slider_names()
leaf_parameter_sliders = get_leaf_slider_names()

SEED_SLIDER_INITIALS = [180, 1, 7, 7, 20, 7, 7, 15, 5]
ROOT_SLIDER_INITIALS = [100, 3, 45, 1, 1, 7, 0.1, 1, 0.2, "Gravitropism"]
STEM_SLIDER_INITIALS = [100, 3, 45, 1, 0.1, 7, 14, 180., 1, 0.2, "Gravitropism"]
LEAF_SLIDER_INITIALS = ["Defined", 30, 1, 45., 2., 90, 1., 0.2, "Gravitropism"]

""" LAYOUT """
app.layout = dbc.Container([

    dcc.Store(id = 'seed-store', data = {"seed": SEED_SLIDER_INITIALS, "basal-checkbox": False, "shoot-checkbox": False, "tillers-checkbox": False }),
    dcc.Store(id = 'root-store', data = {f"tab-{i}": ROOT_SLIDER_INITIALS for i in range(1, 5)}),
    dcc.Store(id = 'stem-store', data = {f"tab-{i}": STEM_SLIDER_INITIALS for i in range(1, 5)}),
    dcc.Store(id = 'leaf-store', data = {"leaf": LEAF_SLIDER_INITIALS}),
    dcc.Store(id = 'typename-store', data = {f"tab-{i}": f"Order {i} root" for i in range(1, 5)}),
    dcc.Store(id = 'result-store', data = {}),

    dbc.Row([
        # Panel 1
        dbc.Col([
            html.H5("Simulation"),
            html.H6("Plant"),
            dcc.Dropdown(id = 'plant-dropdown',
                        options = plants,
                        value = plants[0]['value'],
                        clearable = False,
                        className = 'customDropdown'),
            html.Div(className = "spacer"),
            html.H6("Simulation time [day]"),
            dcc.Slider(id = 'time-slider', min = 1, max = 45, step = 1, value = 20,
                       marks = {1: "1", 45: "45"},
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
            html.Div(id = 'organtype-tabs-content'),
            dcc.Tabs(id = 'root-tabs', children = [], value = ""),
            dcc.Tabs(id = 'stem-tabs', children = [], value = ""),
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
        html.A(
            html.Img(src = '/assets/cplantbox.png', className = 'logo'),
            href = 'https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox',  # Replace with actual URL
            target = '_blank'
        ),
        html.A(
            html.Img(src = '/assets/fzj.png', className = 'logo'),
            href = 'https://www.fz-juelich.de/de',  # Replace with actual URL
            target = '_blank'
        ),
    ], className = 'logoContainer')

], fluid = True)

""" 1. LEFT - Simulation Panel """


@app.callback(# plant-dropdown
    Output('seed-store', 'data'),  # , allow_duplicate=True
    Output('root-store', 'data'),
    Output('stem-store', 'data'),
    Output('leaf-store', 'data'),
    Output('typename-store', 'data'),
    Output('organtype-tabs', 'value'),  # to trigger update
    Output('time-slider', 'value'),
    Input('plant-dropdown', 'value'),
    State('seed-store', 'data'),
    State('root-store', 'data'),
    State('stem-store', 'data'),
    State('leaf-store', 'data'),
    State('typename-store', 'data'),
    State('organtype-tabs', 'value'),
    # prevent_initial_call = True,
)
def update_plant(plant_value, seed_data, root_data, stem_data, leaf_data, typename_data, tabs_value):
    print("update_plant()",plant_value)
    set_data(plant_value, seed_data, root_data, stem_data, leaf_data, typename_data)
    return (seed_data, root_data, stem_data, leaf_data, typename_data, tabs_value, seed_data["simulationTime"])


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
    State('typename-store', 'data'),
    State('result-store', 'data'),
    State('result-tabs', 'value'),
    prevent_initial_call = True,
)

def click_simulate(n_clicks, plant_value, time_slider, seed_data, root_data, stem_data, leaf_data, typename_data, vtk_data, result_value):
    print("click_simulate()", plant_value)
    vtk_data = simulate_plant(plant_value, time_slider, seed_data, root_data, stem_data, leaf_data)
    content = render_result_tab(result_value, vtk_data, typename_data)  # call by hand
    return (content, vtk_data, html.H6(""))


""" Parameters Panel"""


@app.callback(# organtype-tabs
    Output('organtype-tabs-content', 'children'),
    Input('organtype-tabs', 'value'),
    State('seed-store', 'data'),
    State('root-store', 'data'),
    State('typename-store', 'data'),
    State('stem-store', 'data'),
    State('leaf-store', 'data'),
)
def render_organtype_tab(tab, seed_data, root_data, type_names, stem_data, leaf_data):
    print()
    if tab == 'Seed':
        print("render_organtype_tab() seed:", seed_data)
        return generate_seed_sliders(seed_data)
    elif tab == 'Root':
        print("render_organtype_tab() root:", root_data)
        return root_layout(root_data, type_names)
    elif tab == 'Stem':
        print("render_organtype_tab() stem:", stem_data)
        return stem_layout(stem_data, type_names)
    elif tab == 'Leaf':
        print("render_organtype_tab() leaf:", leaf_data)
        return generate_leaf_sliders(leaf_data)


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
        if i in [4, 7]:  # for number of basal, number of tillers
            step_ = 1
        else:
            step_ = 0.1
        style = {}
        if i in [0, 1]:
            style = {'display': 'none'}
        min_ = seed_parameter_sliders[key][0]
        max_ = seed_parameter_sliders[key][1]
        sliders.append(html.H6(key, style = style))
        sliders.append(html.Div([
            dcc.Slider(
                        id = {'type': 'dynamic-seed-slider', 'index': i},
                        min = min_,
                        max = max_,
                        value = seed_values[i],
                        marks = {
                            min_ + 1.e-4: str(min_),
                            max_ - 1.e-4: str(max_)
                            },
                        tooltip = { "always_visible": False},
                        step = step_,
                    )], style = style)
            )
    if data["shoot-checkbox"]:
        v = ['agree']
    else:
        v = []
    sliders.insert(2, dcc.Checklist(
        id = 'shoot-checkbox',
        options = [{'label': html.Span(' Shoot borne roots', className = 'checkbox-label'), 'value': 'agree'}],
        className = 'checkbox',
        value = v,
        style = {'display': 'none'}
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
    print("update_seed_store()", data, slider_values, len(data["seed"]))
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
def generate_root_sliders(root_values, tab):  # Generate sliders for root tabs from stored values
    """ @param root_values is root_data[current_tab] """
    print("generate_root_sliders()", root_values)
    style = {}
    successors = root_values[-1]
    sliders = []
    sliders.append(html.Div(className = "spacer"))
    for i, key in enumerate(root_parameter_sliders.keys()):
        style = {}
        if (not successors) and (i in [3, 4, 5]):  # no lb, ln, la for highest order
            style = {'display': 'none'}
        if (tab == 1) and (i == 2):  # no initial growth rate theta for tap
            style = {'display': 'none'}

        min_ = root_parameter_sliders[key][0]
        max_ = root_parameter_sliders[key][1]
        sliders.append(html.H6(key, style = style))
        # print(key, str(min_), str(max_), min_, max_, "value", root_values[i])
        sliders.append(html.Div([
            dcc.Slider(
                        id = {'type': 'dynamic-slider', 'index': i},
                        min = min_,
                        max = max_,
                        value = root_values[i],
                        marks = {
                            min_ + 1.e-4: str(min_),
                            max_ - 1.e-4: str(max_)
                            },
                        tooltip = { "always_visible": False}
                    )], style = style)
            )
    sliders.append(dcc.Dropdown(
                id = {'type': 'dynamic-slider', 'index': i},  # little white lie
                options = list(tropism_names.keys()),
                value = root_values[i + 1],
                clearable = False,
                className = 'customDropdown'
            ))    
    return html.Div(sliders)


def root_layout(data, type_names):
    """ root tab layout: with root subTypes as sub tabs """
    rootPanelChildren = []
    for i in range(0, type_names["number_roottypes"]):
        ctab = dcc.Tab(label = type_names[f'root tab-{i+1}'], value = f'tab-{i+1}',
                      className = 'tab', selected_className = 'tabSelected')
        rootPanelChildren.append(ctab)
    tab = dcc.Tabs(
        id = 'root-tabs',
        value = 'tab-1',
        children = rootPanelChildren,
    )
    content = html.Div(id = 'root-tabs-content')
    return [tab, content]



@app.callback( # render_root_tab - Display sliders for the selected tab, using stored values
    Output('root-tabs-content', 'children'),
    Input('root-tabs', 'value'),
    State('root-store', 'data'), suppress_callback_exceptions = True
)
def render_root_tab(selected_tab, root_data):
    stored_values = root_data.get(selected_tab, ROOT_SLIDER_INITIALS)
    print("render_root_tab()", selected_tab, root_data[selected_tab], stored_values)
    return generate_root_sliders(stored_values, int(selected_tab[-1]))



@app.callback( # Update root-store when any slider changes
    Output('root-store', 'data', allow_duplicate = True),  #
    Input({'type': 'dynamic-slider', 'index': dash.ALL}, 'value'),
    State('root-tabs', 'value'),
    State('root-store', 'data'),
    prevent_initial_call = True,
)
def update_root_store(slider_values, current_tab, store_data):
    print("update_root_store()", store_data[current_tab])
    successors = store_data[current_tab][-1]
    store_data[current_tab] = slider_values
    store_data[current_tab].append(successors)
    # print("update_root_store() to  ", store_data[current_tab])
    return store_data


#
# Parameters Panel - Stem
#
def generate_stem_sliders(stem_values, tab):  # Generate sliders for stem tabs from stored values
    sliders = []
    successors = stem_values[-1]
    sliders.append(html.Div(className = "spacer"))
    for i, key in enumerate(stem_parameter_sliders.keys()):
        style = {}
        #if i in [7]:  # rotBeta (not working)
        #    style = {'display': 'none'}
        if (not successors) and (i in [3]):
            style = {'display': 'none'}
        # if (tab == 1) and (i == 2):  # no initial theta for main stem; makes sense for tillers
        #     style = {'display': 'none'}
        min_ = stem_parameter_sliders[key][0]
        max_ = stem_parameter_sliders[key][1]
        sliders.append(html.H6(key, style = style))
        sliders.append(html.Div([
            dcc.Slider(
                        id = {'type': 'stem-dynamic-slider', 'index': i},
                        min = min_,
                        max = max_,
                        value = stem_values[i],
                        marks = {
                            min_: str(min_),
                            max_: str(max_)
                            },
                        tooltip = { "always_visible": False},
                    )], style = style)
            )
    sliders.append(dcc.Dropdown(
                id = {'type': 'stem-dynamic-slider', 'index': i},  # little white lie
                options = list(tropism_names.keys()),
                value = stem_values[i + 1],
                clearable = False,
                className = 'customDropdown'
            ))
    return html.Div(sliders)


def stem_layout(data, type_names):
    """ stem tab layout: with stem subTypes as sub tabs """
    panelChildren = []
    for i in range(0, type_names["number_stemtypes"]):
        ctab = dcc.Tab(label = type_names[f'stem tab-{i+1}'], value = f'tab-{i+1}',
                      className = 'tab', selected_className = 'tabSelected')
        panelChildren.append(ctab)
    tab = dcc.Tabs(
        id = 'stem-tabs',
        value = 'tab-1',
        children = panelChildren,
    )
    content = html.Div(id = 'stem-tabs-content')
    stored_values = data.get('tab-1', STEM_SLIDER_INITIALS)
    return [tab, content]


# render_stem_tab - Display sliders for the selected tab, using stored values
@app.callback(
    Output('stem-tabs-content', 'children'),
    Input('stem-tabs', 'value'),
    State('stem-store', 'data'),
)
def render_stem_tab(selected_tab, data):
    stored_values = data.get(selected_tab, STEM_SLIDER_INITIALS)
    print("render_stem_tab()", selected_tab, stored_values)
    return generate_stem_sliders(stored_values, int(selected_tab[-1]))


# Update root-store when any slider changes
@app.callback(
    Output('stem-store', 'data', allow_duplicate = True),  #
    Input({'type': 'stem-dynamic-slider', 'index': dash.ALL}, 'value'),
    State('stem-tabs', 'value'),
    State('stem-store', 'data'),
    prevent_initial_call = True,
)
def update_stem_store(slider_values, current_tab, store_data):
    print("update_stem_store()", current_tab, slider_values)
    successors = store_data[current_tab][-1]
    store_data[current_tab] = slider_values
    store_data[current_tab].append(successors)
    return store_data


#
# Parameters Panel - Leaf
#
def generate_leaf_sliders(data):  # Generate sliders for leaf tabs from stored values

    if data["leaf"] is None:
        return [html.Div(className = "spacer"), html.H6("There is no leaf defined in your plant data set") ]

    leaf_values = data["leaf"]
    # print(len(leaf_values), len(leaf_parameter_sliders.keys()))
    sliders = []
    sliders.append(html.Div(className = "spacer"))
    sliders.append(html.H6("Leaf shape"))
    sliders.append(dcc.Dropdown(
                id = {'type': 'leaf-dynamic-slider', 'index': 0},  # little white lie
                options = ["Defined", "Long", "Round", "Maple", "Flower"],
                value = "Defined",
                clearable = False,
                className = 'customDropdown'
    ))
    sliders.append(html.Div(className = "spacer"))
    for i, key in enumerate(leaf_parameter_sliders.keys()):
        min_ = leaf_parameter_sliders[key][0]
        max_ = leaf_parameter_sliders[key][1]
        sliders.append(html.H6(key))
        sliders.append(
            dcc.Slider(
                        id = {'type': 'leaf-dynamic-slider', 'index': i + 1},
                        min = min_,
                        max = max_,
                        value = leaf_values[i + 1],
                        marks = {
                            min_: str(min_),
                            max_: str(max_)
                            },
                        tooltip = { "always_visible": False},
                    )
            )
    sliders.append(dcc.Dropdown(
                id = {'type': 'leaf-dynamic-slider', 'index': i + 1},  # little white lie
                options = list(tropism_names.keys()),
                value = leaf_values[i + 2],
                clearable = False,
                className = 'customDropdown'
            ))
    return html.Div(sliders)

# @app.callback(# leaf-dropdown
#     Output({'type': 'leaf-dynamic-slider', 'index': dash.ALL}, 'value'),
#     Input({'type': 'leaf-dynamic-slider', 'index': 0}, 'value'),  # leaf-dropdown
#     State({'type': 'leaf-dynamic-slider', 'index': dash.ALL}, 'value'),
# )
# def update_leaf_shape(shape, slider_values):
#     print("update_leaf_shape()", slider_values)
#     set_leaf_sliders(slider_values)
#     return slider_values


@app.callback( # Update leaf-store when any slider changes
    Output('leaf-store', 'data', allow_duplicate = True),  #
    Input({'type': 'leaf-dynamic-slider', 'index': dash.ALL}, 'value'),
    State('leaf-store', 'data'),
    prevent_initial_call = True,
)
def update_leaf_store(slider_values, store_data):
    print("update_leaf_store()", slider_values)
    if len(slider_values) > 1:  # called empty
        store_data["leaf"] = slider_values
    return store_data


""" Results Panel """


@app.callback(
    Output('result-tabs-content', 'children', allow_duplicate = True),
    Input('result-tabs', 'value'),
    State('result-store', 'data'),
    State('typename-store', 'data'),
    prevent_initial_call = True,
)
def render_result_tab(tab, vtk_data, typename_data):
    print("render_result_tab()", tab)
    if not vtk_data:
        print("no data")
        return html.Div([html.H6("press the create button")])
    if tab == 'VTK3D':
        color_pick = vtk_data["subType"]
        color_pick = np.repeat(color_pick, 16)  # 24 = 3*(7+1) (für n=7) ??? 16 (für n=5)
        # print("number of cell colors", len(color_pick), "cells", len(vtk_data["subType"]) , "\n", type(color_pick))
        return vtk3D_plot(vtk_data, color_pick, "Type")
    elif tab == 'VTK3DAge':
        color_pick = vtk_data["creationTime"]
        color_pick = np.repeat(color_pick, 16)  # 24 = 3*(7+1) (für n=7) ??? 16 (für n=5)
        return vtk3D_plot(vtk_data, color_pick, "Age")
    elif tab == 'Profile1D':
        return profile_plot(vtk_data)
    elif tab == 'Dynamics':
        return dynamics_plot(vtk_data, typename_data)


if __name__ == '__main__':
    app.run(debug = True)
