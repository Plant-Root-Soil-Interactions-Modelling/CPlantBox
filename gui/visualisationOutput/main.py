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
app.title = "Model Output Visualisation"

param_names = get_parameter_names()
plants = [{'label': name[0], 'value': str(i)} for i, name in enumerate(param_names)]

seed_parameter_sliders = get_seed_slider_names()
root_parameter_sliders = get_root_slider_names()
stem_parameter_sliders = get_stem_slider_names()
leaf_parameter_sliders = get_leaf_slider_names()

SEED_SLIDER_INITIALS = [180, 1, 7, 7, 20, 7, 7, 15, 5]
ROOT_SLIDER_INITIALS = [100, 3, 45, 1, 1, 7, 0.1, 1, 0.2, "Gravitropism"]
STEM_SLIDER_INITIALS = [100, 3, 45, 1, 0.1, 7, 14, 180., 1, 0.2, "Gravitropism"]
LEAF_SLIDER_INITIALS = ["Defined", 30, 1, 45, 2., 90, 1., 0.2, "Gravitropism"]

""" LAYOUT """
app.layout = dbc.Container([
    dcc.Store(id = 'seed-store', data = {"scenarios": ['baseline'],"variable":'psiXyl', "pSets": [44] }),
    dcc.Store(id = 'root-store', data = {"time": 22,"variable_plant0":'psiXyl', "scenarios0": 'baseline', "pSets0": 44,"soil0":-1,"rsi0":-1,
                                            "variable_plant1":'psiXyl', "scenarios1": 'baseline', "pSets1": 44,"soil1":-1,"rsi1":-1}),
    dcc.Store(id = 'stem-store', data = {"time": 22,"variable0":1, "variableC0":-1, "scenarios0": 'baseline', "pSets0": 44,"doLog0":[],
                                            "variable1":1, "variableC1":-1,  "scenarios1": 'baseline', "pSets1": 44,"doLog1": []}),
    dcc.Store(id = 'leaf-store', data = {"leaf": LEAF_SLIDER_INITIALS}),
    dcc.Store(id = 'typename-store', data = {f"tab-{i}": f"Order {i} root" for i in range(1, 5)}),
    dcc.Store(id = 'result-store', data = {}),

    dbc.Row([
        # Panel 1
        dbc.Col([
    html.Div([
        # Top content
        html.Div([
            dcc.Tabs(
                id='organtype-tabs',
                value="plantDataTab",
                children=[
                    dcc.Tab(label='Plant Data', value='plantDataTab', className='tab', selected_className='tabSelected'),
                    dcc.Tab(label='3D', value='VTK3D', className='tab', selected_className='tabSelected'),
                    dcc.Tab(label='1D Profiles', value='Profile1D', className='tab', selected_className='tabSelected'),
                ]
            ),
            html.Div(id='organtype-tabs-content'),
            dcc.Tabs(id='root-tabs', children=[], value=""),
            dcc.Tabs(id='stem-tabs', children=[], value=""),
        ], style={"flex": "1 1 auto", "overflowY": "auto",
        "padding": "5px"}),

        # Logos at the bottom
        html.Div([
    html.A(
        html.Img(src='/assets/cplantbox.png', className='logo'),
        href='https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox',
        target='_blank'  # Opens in a new tab
    ),
            html.Img(src='/assets/fzj.png', className='logo'),
        ], className='logoContainer')
    ])
], md=5, lg=3, style={"padding": "0"})  , # Ensures full viewport height,

        # Panel 3
        dbc.Col([
            #html.H5("Results"),
            #html.Div(className = "spacer"),
            html.Div(id = 'result-tabs-content',style={"height": "100%", "display": "flex", "flexDirection": "column"})
        ], md=7, lg=9, style={"height": "100vh", "padding": "0"})
    ], style={"height": "100vh"})  # Make row full viewport height

    #html.Div([
    #        html.Img(src = '/assets/cplantbox.png', className = 'logo'),
    #        html.Img(src = '/assets/fzj.png', className = 'logo'),
    #    ], className = 'logoContainer')

], fluid = True)

""" 1. LEFT - Simulation Panel """



@app.callback(# create-button
    Output('result-tabs-content', 'children', allow_duplicate = True),
    Output('result-store', 'data'),
    Output('loading-spinner-output', 'children'),
    Input('create-button', 'n_clicks'),
    #State('plant-dropdown', 'value'),
    State('scenario-dropdown', 'value'),
    State('pSet-dropdown', 'value'),
    State('time-slider', 'value'),
    prevent_initial_call = True,
)
def click_simulate(n_clicks,  time_slider_value, scenario_dropdown_value, pSet_dropdown_value):
    print("click_simulate()", plant_value)
    #vtk_data = simulate_plant(plant_value, time_slider_value, seed_data, root_data, stem_data, leaf_data)
    content = render_result_tab(time_slider_value)  # call by hand
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
    print("render_organtype_tab", tab)
    if tab == 'plantDataTab':
        print("render_organtype_tab() seed:", seed_data)
        return generate_seed_sliders(seed_data)
    elif tab == 'VTK3D':
        print("render_organtype_tab() VTK3D:", root_data)
        return html.Div([                      
                html.Div(className = "spacer"),
                html.H5("Simulation time [day]"),
                dcc.Slider(
                    id='time-slider-vtk3d',
                    min=10, max=25, step=0.1, value=22,
                    marks={10: "10", 25: "25"},
                    tooltip={"always_visible": False}
                ),
                html.Div([
                    generate_VTK3D_sliders(root_data, idPlant="0"),
                    generate_VTK3D_sliders(root_data, idPlant="1")
                ], style={"display": "flex", "flexDirection": "row"})
            ])
    elif tab == 'Profile1D':
        print("render_organtype_tab() stem:", stem_data)
        return html.Div([
                html.Div(className = "spacer"),
                html.H5("Simulation time [day]"),
                dcc.Slider(
                    id='time-slider-1d',
                    min=10.013889, max=25, step=0.013889, value=22,
                    marks={10.013889: "10-", 25: "25"},
                    tooltip={"always_visible": False}
                ),
                html.Div([
                    generate_1D_sliders(stem_data, idPlant="0"),
                    generate_1D_sliders(stem_data, idPlant="1")
                ], style={"display": "flex", "flexDirection": "row"})
            ])

def generate_1D_sliders(data, idPlant):  # Generate sliders for seed tab from stored values

    checkboxes = []        
    if idPlant == "0":
        checkboxes.append(html.H5("Scenario"))
    else:
        checkboxes.append(html.H5("Scenario", style={"color": "white"}))
    checkboxes.append(
        dcc.Dropdown(id = 'Scenario-Dropdown-1d'+idPlant,options=[
                    {"label": " Baseline", "value": "baseline"},
                    {"label": " Early dry spell", "value": "earlyDry"},
                    {"label": " Late dry spell", "value": "lateDry"}
                    ], value="baseline",
                clearable = False))
    if idPlant == "0":
        checkboxes.append(html.H5("Microbial traits"))
    else:
        checkboxes.append(html.H5("Microbial traits", style={"color": "white"}))
    checkboxes.append(dcc.Dropdown(id = 'pSet-Dropdown-1d'+idPlant,options=[
                    {"label": " High respirating", "value": 5},
                    {"label": " High development", "value": 44},
                    {"label": " Low activity", "value": 61}], value=5,
                clearable = False))
    
    if idPlant == "0":
        checkboxes.append(html.H5("Perirhizal zone--values"))
    else:
        checkboxes.append(html.H5("Perirhizal zone--values", style={"color": "white"}))
    checkboxes.append(dcc.Dropdown(id = 'variables-Dropdown-1d'+idPlant,options=[
                    #{"label": " Water", "value": 0},
                    {"label": " Small molecules-C", "value": 1},
                    {"label": " Active copiotrophs-C", "value": 5}], value=1,
                clearable = False))
    if idPlant == "0":
        checkboxes.append(html.H5("Perirhizal zone--colors"))
    else:
        checkboxes.append(html.H5("Perirhizal zone--colors", style={"color": "white"}))
    checkboxes.append(dcc.Dropdown(id = 'variablesC-Dropdown-1d'+idPlant,options=[
                    {"label": " Mesh-cell id", "value": -1},
                    #{"label": " Water", "value": 0},
                    {"label": " Small molecules-C", "value": 1},
                    {"label": " Active copiotrophs-C", "value": 5},
                    {"label": " Active-to-total copiotrophs-C", "value": 6}], value=-1,
                clearable = False))
                    
                    
    checkboxes.append(html.Div(className = "spacer"))
    checkboxes.append(dcc.Checklist(id = 'doLog'+idPlant,options=[
                    {"label": " Logarithmic view", "value": 1}], value=[]))
    return html.Div(checkboxes, style={'width': '50%'})


@app.callback(# dynamic-seed-slider
    Output('stem-store', 'data', allow_duplicate = True),
    Input('time-slider-1d', 'value'),
    Input('variables-Dropdown-1d0', 'value'),
    Input('variablesC-Dropdown-1d0', 'value'),
    Input('Scenario-Dropdown-1d0', 'value'),
    Input('pSet-Dropdown-1d0', 'value'),
    Input('doLog0', 'value'),
    Input('variables-Dropdown-1d1', 'value'),
    Input('variablesC-Dropdown-1d1', 'value'),
    Input('Scenario-Dropdown-1d1', 'value'),
    Input('pSet-Dropdown-1d1', 'value'),
    Input('doLog1', 'value'),
    State('stem-store', 'data'),
    prevent_initial_call = True,
)
def update_stem_store(time, variable0,variableC0,scenarios0, pSets0,doLog0, variable1, variableC1, scenarios1, pSets1, doLog1,data):
    data["time"] = time
    data["variable0"] = variable0
    data["variableC0"] = variableC0
    data["scenarios0"] = scenarios0
    data["pSets0"] = pSets0
    data["doLog0"] = doLog0
    data["variable1"] = variable1
    data["variableC1"] = variableC1
    data["scenarios1"] = scenarios1
    data["pSets1"] = pSets1
    data["doLog1"] = doLog1
    return data
    
def generate_VTK3D_sliders(data, idPlant):  # Generate sliders for seed tab from stored values

    checkboxes = []        
    
    html.Div(className = "spacer")
    if idPlant == "0":
        checkboxes.append(html.H5("Scenario"))
    else:
        checkboxes.append(html.H5("Scenario", style={"color": "white"}))
    checkboxes.append(
        dcc.Dropdown(id = 'Scenario-Dropdown'+idPlant,options=[
                    {"label": " Baseline", "value": "baseline"},
                    {"label": " Early dry spell", "value": "earlyDry"},
                    {"label": " Late dry spell", "value": "lateDry"}
                    ], value="baseline",
                clearable = False))
    if idPlant == "0":
        checkboxes.append(html.H5("Microbial traits"))
    else:
        checkboxes.append(html.H5("Microbial traits", style={"color": "white"}))
    checkboxes.append(dcc.Dropdown(id = 'pSet-Dropdown'+idPlant,options=[
                    {"label": " High respirating", "value": 5},
                    {"label": " High development", "value": 44},
                    {"label": " Low activity", "value": 61}], value=5,
                clearable = False))
    
    if idPlant == "0":
        checkboxes.append(html.H5("Plant")) # Physiological Response
    else:
        checkboxes.append(html.H5("Plant", style={"color": "white"}))
                
    vars_ = {
        'psiXyl': "Mean plant water potential",
        'Q_Gr': "C used for growth",
        "Q_Rm": "C used for maintenance respiration",
        'Q_Exud': "C used for exudation"
    }
    #'psiXyl', 'Q_Gr', 'Q_Exud'
    checkboxes.append(dcc.Dropdown(
                id = 'variables-dropdown-vtk3d'+idPlant, 
                options = [{'label': v, 'value': k} for k, v in vars_.items()],
                value = list(vars_.keys())[0],
                clearable = False
            ))
    if idPlant == "0":
        checkboxes.append(html.H5("Soil"))
    else:
        checkboxes.append(html.H5("Soil", style={"color": "white"}))
    checkboxes.append(dcc.Dropdown(id = 'Soil-Dropdown'+idPlant,options=[
                    {"label": " None", "value": -1},
                    {"label": " Water", "value": 0},
                    {"label": " small molecules-C", "value": 1},
                    {"label": " active copiotrophs-C", "value": 5}], value=-1,
                clearable = False))
    if idPlant == "0":
        checkboxes.append(html.H5("Root-soil interface"))
    else:
        checkboxes.append(html.H5("Root-soil interface", style={"color": "white"}))
    checkboxes.append(dcc.Dropdown(id = 'RSI-Dropdown'+idPlant,options=[
                    {"label": " None", "value": -1},
                    {"label": " Water", "value": 0},
                    {"label": " small molecules-C", "value": 1},
                    {"label": " active copiotrophs-C", "value": 5}], value=-1,
                clearable = False))
    return html.Div(checkboxes, style={'width': '50%'})


@app.callback(# dynamic-seed-slider
    Output('root-store', 'data', allow_duplicate = True),
    Input('time-slider-vtk3d', 'value'),
    Input('variables-dropdown-vtk3d0', 'value'),
    Input('Scenario-Dropdown0', 'value'),
    Input('pSet-Dropdown0', 'value'),
    Input('Soil-Dropdown0', 'value'),
    Input('RSI-Dropdown0', 'value'),
    Input('variables-dropdown-vtk3d1', 'value'),
    Input('Scenario-Dropdown1', 'value'),
    Input('pSet-Dropdown1', 'value'),
    Input('Soil-Dropdown1', 'value'),
    Input('RSI-Dropdown1', 'value'),
    State('root-store', 'data'),
    prevent_initial_call = True,
)
def update_root_store(time, variable0,scenarios0, pSets0, soil0, rsi0,variable1, scenarios1, pSets1,soil1,rsi1, data):
    data["time"] = time
    data["variable_plant0"] = variable0
    data["scenarios0"] = scenarios0
    data["pSets0"] = pSets0
    data["soil0"] = soil0
    data["rsi0"] = rsi0
    data["variable_plant1"] = variable1
    data["scenarios1"] = scenarios1
    data["pSets1"] = pSets1
    data["soil1"] = soil1
    data["rsi1"] = rsi1
    return data

#
#  Parameters Panel - Seed
#
def generate_seed_sliders(data):  # Generate sliders for seed tab from stored values

    checkboxes = [
        html.Div(className = "spacer"),
        ]        
    
    
    #checkboxes.append(html.H5("Simulation time [day]"))
    #checkboxes.append(dcc.Slider(id = 'time-slider', min = 10.1, max = 25, step = 0.1, value = 22,
    #           marks = {10.1: "10.1", 25: "25"},
    #           tooltip = { "always_visible": False}))
    #checkboxes.append(html.Div(className = "spacer"))
    #checkboxes.append(dcc.Checklist(id = 'cst-lims',options=[
    #                {"label": " Constant limits", "value": 1}], value=[]))
    html.Div(className = "spacer")
    checkboxes.append(html.H5("Scenario"))
    checkboxes.append(
        dcc.Checklist(id = 'Scenario-Checklist',options=[
                    {"label": " Baseline", "value": "baseline"},
                    {"label": " Early dry spell", "value": "earlyDry"},
                    {"label": " Late dry spell", "value": "lateDry"}
                    ], value=["baseline"]))
    checkboxes.append(html.H5("Microbial traits"))
    checkboxes.append(dcc.Checklist(id = 'pSet-Checklist',options=[
                    {"label": " High respirating", "value": 5},
                    {"label": " High development", "value": 44},
                    {"label": " Low activity", "value": 61}], value=[5]))
                    
    checkboxes.append(html.H5("Plant Physiological Response"))
    vars_ = {
        'psiXyl': "Mean plant water potential",
        'Q_Gr': "C used for growth",
        "Q_Rm": "C used for maintenance respiration",
        'Q_Exud': "C used for exudation"
    }
    #'psiXyl', 'Q_Gr', 'Q_Exud'
    checkboxes.append(dcc.Dropdown(
                id = 'variables-dropdown', 
                options = [{'label': v, 'value': k} for k, v in vars_.items()],#list(vars_.values()),
                value = list(vars_.keys())[0],
                clearable = False,
            ))
    return html.Div(checkboxes)


@app.callback(# dynamic-seed-slider
    Output('seed-store', 'data', allow_duplicate = True),
    Input('variables-dropdown', 'value'),
    #Input('time-slider', 'value'),
    Input('Scenario-Checklist', 'value'),
    Input('pSet-Checklist', 'value'),
    #Input('cst-lims', 'value'),
    State('seed-store', 'data'),
    prevent_initial_call = True,
)
def update_seed_store(variable, scenarios, pSets,  data):
    data["variable"] = variable
    #data["time"] = time
    data["scenarios"] = scenarios
    data["pSets"] = pSets
    #data["cstLims"] = cstLims
    return data

""" Results Panel """


@app.callback(
    Output('result-tabs-content', 'children', allow_duplicate = True),
    Input('organtype-tabs', 'value'),
    Input('root-store', 'data'),
    Input('seed-store', 'data'),
    Input('stem-store', 'data'),
    prevent_initial_call = True,
)
def render_result_tab(tab,rootdata,seeddata, stemdata):# scenarios, pSet,
    
    if tab == 'VTK3D':
        if (rootdata['rsi0'] >= 0) or (rootdata['rsi1'] >= 0):
            thekey = str(time.time())
        else:
            thekey = 'coucou'
        
        return html.Div(
                    [vtk3D_plot(rootdata, 0), vtk3D_plot(rootdata, 1)],
                    style={"display": "flex", "flexDirection": "row", "height": "100%"},
                    key=thekey
                )# <- this forces remount
    elif tab == 'plantDataTab':
        data = seeddata
        toplot = data['variable']
        scenarios = data['scenarios']
        pSet = data['pSets']
        return plantData_plot_plotly([toplot], scenarios, pSet)
    elif tab == 'Profile1D':
    
        return html.Div([profile_plot(stemdata, 0),
                        profile_plot(stemdata, 1)], 
                        style={"display": "flex", "flexDirection": "row"})



if __name__ == '__main__':
    app.run(debug = True)
