""" methods creating the figures of the cplantbox webapp (D. Leitner, 2025)"""

import sys; sys.path.append("../.."); sys.path.append("../../src/")
from visualisation.vtk_tools import *
from visualisation.vtk_plot import *
import vtk

import numpy as np
import pandas as pd
import string

from dash import html, dcc, Input, Output, State, ctx
import plotly.graph_objs as go
from plotly.colors import qualitative
import dash_vtk
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter
import pickle
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from vtk.util.numpy_support import vtk_to_numpy
from dash_vtk.utils import to_mesh_state, to_volume_state
from vtk.util.numpy_support import vtk_to_numpy

with open("data/4Andrea_coordcs.pkl", "rb") as file:
    dataradial = pickle.load( file)
dfcoord, dfcs, dfccat, dfcca, dfcoa = dataradial[0], dataradial[1], dataradial[2], dataradial[3], dataradial[4]
with open("data/plantData.pkl", "rb") as file:
    plantData = pickle.load(file)
    

def read_vtp(filename):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def process_plant_vtp(file_main, file_leaf, pname):
    # Read input
    polydata = read_vtp(file_main)
    leafdata = read_vtp(file_leaf)

    # Convert cell data to point data
    c2p = vtk.vtkCellDataToPointData()
    c2p.SetInputData(polydata)
    c2p.PassCellDataOn()
    c2p.Update()
    pointdata = c2p.GetOutput()

    # Tube filter
    tube = vtk.vtkTubeFilter()
    tube.SetInputData(pointdata)
    tube.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tube.SetRadius(0.1)
    tube.SetNumberOfSides(12)
    tube.SetInputArrayToProcess(0, 0, 0, 0, "radius")
    tube.Update()
    
    array = pointdata.GetPointData().GetArray(pname)
    #if array:
    #    tube.GetOutput().GetPointData().AddArray(array)
    #    tube.GetOutput().GetPointData().SetActiveScalars(pname)
    # Scalar range
    q_exud = vtk_to_numpy(array) if array else np.array([])
    q_range = [float(q_exud.min()), float(q_exud.max())] if q_exud.size > 0 else [0.0, 1e-6]

    # Create tube state and manually set scalar field + lookup table
    tube_state = to_mesh_state(tube.GetOutput(), pname)
    q_range = tube_state["field"]["dataRange"]
    
    # Create solid-colored leaf mesh
    leaf_state = to_mesh_state(leafdata, pname)
    #leaf_state["color"] = [0.333, 0.667, 0.498]  # RGB
    #leaf_state.pop("activeScalars", None)
    #leaf_state.pop("lookupTable", None)

    return tube_state, leaf_state, q_range

def process_rsi_vtp(file_main):
    # Read input
    polydata = read_vtp(file_main)

    # Convert cell data to point data
    c2p = vtk.vtkCellDataToPointData()
    c2p.SetInputData(polydata)
    c2p.PassCellDataOn()
    c2p.Update()
    pointdata = c2p.GetOutput()
    print('get arrays')
    point_data = pointdata.GetPointData()
    for i in range(point_data.GetNumberOfArrays()):
        print(f"- {point_data.GetArrayName(i)}")
    print('get arraysDone')
    pname = point_data.GetArrayName(0)
    # Tube filter
    tube = vtk.vtkTubeFilter()
    tube.SetInputData(pointdata)
    tube.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tube.SetRadius(0.1)
    tube.SetNumberOfSides(12)
    tube.SetInputArrayToProcess(0, 0, 0, 0, "radius")
    tube.Update()
    
    array = pointdata.GetPointData().GetArray(pname)
    q_exud = vtk_to_numpy(array) if array else np.array([])
    q_range = [float(q_exud.min()), float(q_exud.max())] if q_exud.size > 0 else [0.0, 1e-6]

    print('q_exud',pname, q_range)
    # Create tube state and manually set scalar field + lookup table
    tube_state = to_mesh_state(tube.GetOutput(), pname)
    print('tube_state',tube_state.keys())
    q_range = tube_state["field"]["dataRange"]
    print('q_exud',pname, q_range)
    
    return tube_state, q_range
    

def uniform_grid(min_, max_, res):
    """ Creates an uniform grid
    @param min_    minimum of bounding rectangle
    @param max_    maximum of bounding rectangle
    @param res_    cell resolution
    @return A vtkUniformGrid
    """
    grid = vtk.vtkUniformGrid()
    grid.SetDimensions(int(res[0]) + 1, int(res[1]) + 1, int(res[2]) + 1)  # cells to corner points
    grid.SetOrigin(min_[0], min_[1], min_[2])
    s = (max_ - min_) / res
    grid.SetSpacing(s[0], s[1], s[2])
    return grid
def get_bounds_from_volume_state(volume_state):
    image = volume_state["image"]
    origin = np.array(image["origin"])
    spacing = np.array(image["spacing"])
    dimensions = np.array(image["dimensions"])

    # The bounds go from origin to origin + spacing * (dim - 1)
    max_bound = origin + spacing * (dimensions - 1)
    min_bound = origin  # just the origin

    return min_bound, max_bound
    
def plot_mesh_cuts(grid, array_name, nz = 3):
    """ plots orthogonal nz vertical cuts z[:-1] (xy-planes), with z = linspace(min_z, max_z, nz+1),
    and two additonal sclices at x=0 (yz-plane), y=0 (xz-plane)
    @param grid         some vtk grid (structured or unstructured)
    @param array_name       parameter to visualize
    @param nz           number of vertical slices
    """

    eps = 1.e-2
    planes = []  # create the cut planes
    bounds = grid.GetBounds()
    z = np.linspace(bounds[4] + eps, bounds[5], nz + 1)
    for i in range(0, nz):  # z-slices (implicit functions)
        p = vtk.vtkPlane()
        p.SetOrigin(0, 0, z[i])
        p.SetNormal(0, 0, 1)
        planes.append(p)
    for n in [(1, 0, 0), (0, 1, 0)]:
        p = vtk.vtkPlane()
        p.SetOrigin(bounds[0] + eps, bounds[2] + eps, bounds[4])
        p.SetNormal(n[0], n[1], n[2])
        planes.append(p)
    output_states = []
    for p in planes:
        cutter = vtk.vtkCutter()
        cutter.SetInputData(grid)

        cutter.SetCutFunction(p)
        cutter.Update()
        cut_output = cutter.GetOutput()
        # Ensure scalar array is active
        cut_output.GetCellData().SetActiveScalars(array_name)

        mesh_state = to_mesh_state(cut_output, field_to_keep=array_name)
        output_states.append(mesh_state)

    return output_states

def process_soil_vti(file_main, pname):
    # Read input
    image_data  = read_rect_vtu(file_main)
    cell_data = image_data.GetCellData()

    array = cell_data.GetArray(pname)
    data_array = vtk_to_numpy(array).flatten()
    data_range = [float(data_array.min()), float(data_array.max())] 
    cut_states =plot_mesh_cuts(image_data, pname, nz = 3)
    return cut_states, data_range
    
def generate_colorbar_image(vmin, vmax, colormap="Viridis", height=500, width=100):
    # Create gradient data
    z = np.linspace(vmin, vmax, 256).reshape(-1, 1)

    fig = go.Figure(go.Heatmap(
        z=z,
        colorscale=colormap,
        showscale=False,
        x0=0, dx=1,
        y0=vmin, dy=(vmax - vmin) / 256,
        colorbar=None
    ))

    fig.update_layout(
        width=width,
        height=height,
        margin=dict(l=0, r=100, t=20, b=20),  # add enough right margin for text
        xaxis=dict(
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            visible=False  # hide x-axis entirely
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            ticks="outside",
            side="right",  # display ticks on the right
            tickvals=[vmin, vmax],
            ticktext=[f"{vmin:.2e}", f"{vmax:.2e}"],
            tickfont=dict(size=14),
            range=[vmin, vmax],
            automargin=True
        )
    )

    return fig
    
def vtk3D_plot(data, plantid):#): #vtk_data, color_pick):
    """ 3D and 3D age plot """
    time = data['time']
    scenario = data['scenarios'+str(plantid)]
    pSet = data['pSets'+str(plantid)]
    pname = data['variable_plant'+str(plantid)]
    soildata = data['soil'+str(plantid)]
    rsidata = data['rsi'+str(plantid)]
    
    print('vtk3D_plot', scenario, pSet, pname)
    if scenario == "lateDry":
        sc = "lD"
    elif scenario == "earlyDry":
        sc = "eD"
    elif scenario == "baseline":
        sc = "b"
    path = "data/vtpvti/"+sc+str(pSet)+"/vtpvti/"
    filename = "plantAt"+str(int(time*10))
    file_main = path + filename + ".vtp"
    file_leaf = path + filename + "_leaf.vtp"

    tube_state, leaf_state, q_range = process_plant_vtp(file_main, file_leaf, pname)
    viewChildren = [
            
            dash_vtk.GeometryRepresentation(
            id="leaf-rep",
            colorMapPreset="erdc_rainbow_bright",
            colorDataRange=q_range,
            children=[
                dash_vtk.Mesh(state=leaf_state)
            ])
        ]
    q_range_rsi = None
    q_range_soil = None
    if rsidata >= 0:
        if rsidata == 0:
            filename = path + "soil_rx"+str(int(time*10))+ ".vtp"
            #pname = "pressure head"
        else:
            filename = path + "C"+str(int(rsidata))+"_"+str(int(time*10))+ ".vtp"
            #pname = "[C"+str(int(soildata))+"] (mol/cm3)"
        rsi_state, q_range_rsi = process_rsi_vtp(filename)
        
    if soildata >= 0:
        if soildata == 0:
            filename = path + "soil_rx"+str(int(time*10))+ ".vti"
            pname = "pressure head"
        else:
            filename = path + "C"+str(int(soildata))+"_"+str(int(time*10))+ ".vti"
            pname = "[C"+str(int(soildata))+"] (mol/cm3)"
        cut_states, q_range_soil = process_soil_vti(filename, pname)
    
    if (rsidata == soildata) and (q_range_rsi is not None) and (q_range_soil is not None):
        q_range = [min(q_range_rsi[0], q_range_soil[0]), max(q_range_rsi[1], q_range_soil[1])]
        q_range_rsi = q_range
        q_range_soil = q_range
        
    if rsidata >= 0:
        viewChildren.append(
            dash_vtk.GeometryRepresentation(
            id="tubesRSI-rep",
            colorMapPreset="erdc_rainbow_bright",
            colorDataRange=q_range_rsi,
            children=[
                dash_vtk.Mesh(state=rsi_state)
            ])
        )
    else:
        viewChildren.append(dash_vtk.GeometryRepresentation(
            id="tubes-rep",
            colorMapPreset= "rainbow",#"Viridis (matplotlib)" ,#"erdc_rainbow_bright",
            colorDataRange=q_range,
            #showCubeAxes=True,
            children=[
                dash_vtk.Mesh(state=tube_state)
            ]))
    cbar = generate_colorbar_image(vmin=q_range[0], vmax=q_range[1], 
                            colormap="Rainbow", #"Viridis", 
                            height=200, width=150)
        
    
    if soildata >= 0:
        print('adding soil' )
        for state in cut_states:
            viewChildren.append(
               dash_vtk.GeometryRepresentation(
                id="soil-rep",
            children=[
                dash_vtk.Mesh(state=state)
            ],
            colorMapPreset="erdc_rainbow_bright",
            colorDataRange=q_range_soil,
        ))
            
    
    content = dash_vtk.View(
        id="vtk-view",
        children= viewChildren,
        background=[0.9, 0.9, 0.9],  # Light gray background
        cameraPosition=[81, -49, 17],
        #cameraFocalPoint=[4, -7, -1],
        cameraViewUp=[-0.18, 0.09, 0.979],
        style={'flex': '1 1 auto', 'height': '100%', 'width': '100%'})
    cbarcontent = dcc.Graph(
            id = 'profile-plot',
            figure = cbar,
            style = {'width': '25%', 'height': '600px'}
        )


    return html.Div(
                children=[
                    html.Div(content, style={
                        "flex": "1 1 auto",
                        "height": "100%",
                        "width": "100%",
                        "display": "flex",
                        "flexDirection": "column",
                    })
                ],
                style={
                    "flex": "1 1 auto",
                    "height": "100%",
                    "width": "100%",
                    "display": "flex",
                    "flexDirection": "column",
        "padding": "5px", 
                }
            )
'''
    
    return html.Div([
            # First child: VTK 3D view (full width)
            html.Div(
                content,
                style={'width': '100%', 'height': '100%'}
            ),

            # Second child: colorbar
            #html.Div(
            #    cbarcontent,
            #    style={
            #        "width": "100px",         # Narrow fixed width
            #        "height": "20px",
            #        "margin": "20px auto",     # auto horizontal margins centers it
            #        "alignSelf": "center"
            #    }
            #)
        ],
        style={'flex': '1 1 auto', 'height': '100%', 'display': 'flex', 'flexDirection': 'column'})#return  html.Div(content, style = {'width': '100%', 'height': '80vh', 'minHeight': '600px'})
'''
def dotheplot_plotly(fig, toplot, cumsum, df, namesyaxes=None, ncols=3, maxTime=None, 
                     indexlegend1=0, indexlegend2=1, microbes=[5, 44, 61], 
                     scenarios=["baseline", "earlyDry", "lateDry"],cstLims=False):

    nrows = len(toplot)
    cm3tommol = 1000 * (1 / 18.01528)
    ratioTrans = 250
    
    changeUp = 1.01
    changeDown = 0.9
    if cstLims:
        max_y = []
        min_y = []
        for rowid, tipi in enumerate(toplot):
            if tipi == "psiXyl":
                max_y.append(0.)
                min_y.append(min([dftemp[tipi].min(axis=0) for dftemp in df]
                       )*changeUp)
            else:
                min_y.append(0.)
                max_y.append(max([dftemp[tipi].max(axis=0) for dftemp in df]
                       )*changeUp)
                       
    if len(microbes) < 3:
        df = [dftemp for dftemp in df if int(dftemp['pSet'].unique()[0]) in microbes]
    if len(scenarios) < 3:
        df = [dftemp for dftemp in df if dftemp['scenario'].unique()[0] in scenarios]

    for ddid, dd in enumerate(df):
        pset_ = int(dd['pSet'].unique()[0])
        scenario_ = dd['scenario'].unique()[0]
        colid = scenarios.index(scenario_) if ncols > 1 else 0

        for rowid, tipi in enumerate(toplot):
            factor, rot, unit = 1000, 0, " (mmol C)"
            if tipi in ["trans", 'psiXyl']:
                factor = cm3tommol
                unit = "(MPa)"
            if tipi == 'root:shoot':
                factor = 1
            if tipi == "psiXyl":
                factor = 1.e-4
                rot = 0
                unit = "(MPa)"

            ddd = dd[tipi] * factor
            tts = dd['time']
            ddd[tts> maxTime] = np.nan
            if not cumsum:
                tts = dd['time'][1:]
                ddd = np.diff(ddd) * 4 # because have 20mn time step

            index = rowid * ncols + colid
            subplot_id = index + 1  # plotly subplot index starts from 1
            label = '(' + string.ascii_lowercase[index] + ')'

            linestyle = l_styles[pset_] #if len(microbes) == 3 else 'solid'
            color = c_styles[scenario_]

            fig.add_trace(go.Scatter(
                x=tts,
                y=ddd,
                mode='lines',
                name=scenario_,
                line=dict(color=color, dash=linestyle, width=4),
                showlegend=False
            ), row=rowid + 1, col=colid + 1)
            
            subplot_idx = rowid * ncols + colid + 1  # Plotly subplots index from 1
            xref = f'x{subplot_idx}' if subplot_idx > 1 else 'x'
            yref = f'y{subplot_idx}' if subplot_idx > 1 else 'y'
            if False:
                fig.add_annotation(
                    text=label,
                    x=0,
                    y=0.95,
                    xref=f'{xref} domain',
                    yref=f'{yref} domain',
                    showarrow=False,
                    font=dict(size=14, family='serif')
                )

            # add vertical lines for stress periods
            vlines = []
            if scenario_ == "lateDry":
                vlines = [18, 25]
            elif scenario_ == "earlyDry":
                vlines = [11, 18]
            for xval in vlines:
                fig.add_vline(x=xval, line_dash='dash', line_color='black', row=rowid + 1, col=colid + 1)

            if cstLims:#maxTime is not None:
                print('ig.update_yaxes', min_y[rowid]* factor, max_y[rowid]* factor)
                fig.update_xaxes(range=[10.,25.], row=rowid + 1, col=colid + 1)
                fig.update_yaxes(range=[min_y[rowid]* factor, max_y[rowid]* factor], row=rowid + 1, col=colid + 1)
            else:
                fig.update_xaxes(range=[10., maxTime], row=rowid + 1, col=colid + 1)

            # Y-axis label
            if colid == 0:
                ylabel = namesyaxes.get(tipi, tipi) if namesyaxes else tipi
                if tipi == "psiXyl":
                    ylabel += ' ' +unit
                elif tipi == 'root:shoot':
                    ylabel += ' (-)'
                else:
                    ylabel = f"cumulative\nmmol C\nfor {ylabel}"
                fig.update_yaxes(title_text=ylabel, row=rowid + 1, col=colid + 1) # Simulation time [day]
                fig.update_xaxes(title_text='Simulation time [day]', row=rowid + 1, col=colid + 1) # 

            # X-axis label
            if rowid == nrows - 1:
                fig.update_xaxes(title_text="time (day)", row=rowid + 1, col=colid + 1)

l_styles = {5:'solid',44:'dash',61:'dot'} 
# Define the color dictionary
c_styles = {
    'baseline':'#8da0cb',  # Blue
    'earlyDry':  '#66c2a5' ,  # Green
    'lateDry': '#fc8d62'  # red
}
def plantData_plot_plotly(toplot, scenarios, pSet):
    """ Generate a Plotly figure showing plant physiological data across scenarios and parameter sets """
    #print('toplot', toplot)
    ncols = 1
    #toplot = ['psiXyl', 'Q_Gr', 'Q_Exud']
    nrows = len(toplot)

    fig = make_subplots(
        rows=nrows,
        cols=ncols,
        shared_xaxes=False,
        shared_yaxes=False,
        subplot_titles=[None] * (nrows * ncols),
        vertical_spacing=0.07
    )

    namesyaxes = {
        'psiXyl': "mean plant water potential",
        'Resp': 'Resp',
        'root:shoot': 'root:shoot\nC ratio',
        'Q_Gr': " growth",
        "Q_Rm": "maintenance",
        "Q_Ag": "Q_Ag",
        'Q_Exud': " exudation",
        'Q_Exud_mean': "exudation rate\n at root tip",
        'Q_Mucil': "mucilage release"
    }

    dotheplot_plotly(
        fig,
        toplot=toplot,
        cumsum=True,
        df=plantData,
        namesyaxes=namesyaxes,
        ncols=ncols,
        maxTime=26.,
        microbes=pSet,
        scenarios=scenarios,
        indexlegend1=1,
        indexlegend2=1,
        cstLims=True
    )
    # === Custom legend traces ===

    # Color legend (scenarios)
    color_legend_traces = [
        go.Scatter(
            x=[None], y=[None],
            mode='lines',
            line=dict(color=c_styles['baseline'], width=4),
            name='baseline',
            showlegend=True
        ),
        go.Scatter(
            x=[None], y=[None],
            mode='lines',
            line=dict(color=c_styles['earlyDry'], width=4),
            name='earlyDry',
            showlegend=True
        ),
        go.Scatter(
            x=[None], y=[None],
            mode='lines',
            line=dict(color=c_styles['lateDry'], width=4),
            name='lateDry',
            showlegend=True
        )
    ]

    # Line style legend (microbes)
    linestyle_legend_traces = [
        go.Scatter(
            x=[None], y=[None],
            mode='lines',
            line=dict(color='black', dash=l_styles[5], width=4),
            name='High respirating',
            showlegend=True
        ),
        go.Scatter(
            x=[None], y=[None],
            mode='lines',
            line=dict(color='black', dash=l_styles[44], width=4),
            name='High development',
            showlegend=True
        ),
        go.Scatter(
            x=[None], y=[None],
            mode='lines',
            line=dict(color='black', dash=l_styles[61], width=4),
            name='Low activity',
            showlegend=True
        )
    ]

    # Add to the figure
    for trace in color_legend_traces + linestyle_legend_traces:
        fig.add_trace(trace)
    fig.update_layout(
        #height=int(1800/3),
        #width=1000,
        #title_text="",
        autosize=True,
        showlegend=True,
    legend=dict(
        orientation="h",
        yanchor="bottom",
        y=-0.25,
        xanchor="center",
        x=0.5
    ),
        font=dict(size=14)
    )

    content = dcc.Graph(
            id = 'profile-plot',
            figure = fig,
        config={'responsive': True},
        style={'width': '100%', 'height': '100%'}
        )

    return  html.Div(content, style = {'width': '100%', 'height': '80vh', 'minHeight': '600px'})
    
def get_val_along_r_figure(dfx_, dfy_, scenarios, pSets, ylab=None, xlab=None,
                            unitChangex=10, unitChangey=1e6, timeeval=-1.,
                            colordf_=None, doLog=False):

    color_palette = ['#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef',
                     '#e6f5d0','#b8e186','#7fbc41','#4d9221']
    color_mapping = dict(zip(range(9), color_palette))

    if timeeval < 0.:
        timeeval = max(dfx_['time'])
        
    closest_time = dfx_['time'].iloc[(dfx_['time'] - timeeval).abs().argsort().values[0]]
    dfx = dfx_[dfx_['time'] == closest_time].copy()
    dfy = dfy_[dfy_['time'] == closest_time].copy()
    if colordf_ is not None:
        colordf = colordf_[colordf_['time'] == closest_time].copy()#colordf_.loc[[closest_time]]

    # Convert pSet to string once if needed
    pSet_str = str(pSets[0])
    scenario = scenarios[0]

    mask = (dfx['pSet'] == pSet_str) & (dfx['scenario'] == scenario)
    dfx__ = dfx.loc[mask]
    dfy__ = dfy.loc[mask]
    if colordf_ is not None:
        colordf__ = colordf.loc[mask]

    # Preselect column names and values
    xcols = [f'cellS{i}' for i in range(8, -1, -1)]
    ycols = [f'cell{i}' for i in range(8, -1, -1)]

    xdata = dfx__[xcols].values * unitChangex
    ydata = dfy__[ycols].values * unitChangey
    if colordf_ is not None:
        cdata = colordf__[ycols].values

    # Generate Plotly traces
    traces = []
    for i, (x, y) in enumerate(zip(xdata.T, ydata.T)):
        cell_id = 8 - i  # because range is from 8 to 0
        if colordf_ is not None:
            trace = go.Scatter(
                x=x,
                y=y,
                mode='markers',
                marker=dict(
                    color=cdata[:, i],
                    colorscale='Rainbow',
                    colorbar=dict(title='Value'),
                    showscale=False
                ),
                name=f'Cell-id {cell_id}'
            )
        else:
            trace = go.Scatter(
                x=x,
                y=y,
                mode='markers',
                marker=dict(color=color_mapping[cell_id]),
                name=f'Cell {cell_id}'
            )
        traces.append(trace)

    layout = go.Layout(
        xaxis=dict(title=xlab, showgrid=False),
        yaxis=dict(title=ylab, type='log' if doLog else 'linear', showgrid=False),
        font=dict(size=16),
        hovermode='closest',
        showlegend=False
    )

    return go.Figure(data=traces, layout=layout)


def profile_plot(data, index):
    variable = data["variable"+str(index)] 
    variableC = data["variableC"+str(index)] 
    scenarios = data["scenarios"+str(index)] 
    pSets = data["pSets"+str(index)]
    doLog = (1 in data["doLog"+str(index)] )
    time = data["time"] 
    
    if variable == 1:
        df_y = dfcs
    elif variable == 5:
        df_y = dfcca
        
    if variableC == -1:
        colordf = None
    if variableC == 1:
        colordf = dfcs
    elif variableC == 5:
        colordf = dfcca
    elif variableC == 6:
        colordf = dfccat
        
        
    gofig = get_val_along_r_figure( dfcoord, df_y,scenarios=[scenarios],pSets=[pSets],unitChangex=10,unitChangey=1e6,
                       ylab="low weight organic molecules-C (mumol/cm3)",
                       xlab="distance from root surface (mm)",
             doLog=doLog, timeeval=time,
             colordf_=colordf)

    content = dcc.Graph(
            id = 'profile-plot',
            figure = gofig,
            style = {'width': '100%',  'height': '100%'}
        )
    return  html.Div(content, style = {"width": "100%",  'height': '80vh', 'minHeight': '600px'})


    
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

