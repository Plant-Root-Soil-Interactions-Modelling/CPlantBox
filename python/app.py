#!/usr/bin/ python3
# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import py_plantbox as pb
from rb_tools import *
import timeit
import numpy as np
from scipy.integrate import odeint
from werkzeug.contrib.fixers import ProxyFix
#import pandas as pd



app = dash.Dash()
app.server.wsgi_app = ProxyFix(app.server.wsgi_app)
app.config.update({
    # as the proxy server will remove the prefix
    'routes_pathname_prefix': '/',

    # the front-end will prefix this string to the requests
    # that are made to the proxy server
    'requests_pathname_prefix': '/'
})


server = app.server

input_value = 'PMA2018'


BACKGROUND = 'rgb(230, 230, 230)'

COLORSCALE = [ [0, "rgb(244,236,21)"], [0.3, "rgb(249,210,41)"], [0.4, "rgb(134,191,118)"],
                [0.5, "rgb(37,180,167)"], [0.65, "rgb(17,123,215)"], [1, "rgb(54,50,153)"] ]
plant = pb.Plant()
plant.openXML(input_value)
plant.initialize()
plant.simulate(160,True)
nodes = vv2a(plant.getNodes())/100
def scatter_plot_3d(
        x = nodes[:,0],
        y = nodes[:,1],
        z = nodes[:,2],
        size = '1',
        color = 'red',
        xlabel = 'x',
        ylabel = 'y',
        zlabel = 'z',
        plot_type = 'scatter3d',
        markers = [] ,
        ):

    def axis_template_3d( title, type='linear' ):
        return dict(
            showbackground = True,
            backgroundcolor = BACKGROUND,
            gridcolor = 'rgb(255, 255, 255)',
            title = title,
            type = type,
            zerolinecolor = 'rgb(255, 255, 255)'
        )
    
    def axis_template_2d(title):
        return dict(
            xgap = 10, ygap = 10,
            backgroundcolor = BACKGROUND,
            gridcolor = 'rgb(255, 255, 255)',
            title = title,
            zerolinecolor = 'rgb(255, 255, 255)',
            color = '#444'
        )

    def blackout_axis( axis ):
        axis['showgrid'] = False
        axis['zeroline'] = False
        axis['color']  = 'white'
        return axis

    data = [ dict(
        x = x,
        y = y,
        z = z,
        mode = 'markers',
        marker = dict(
                colorscale = COLORSCALE,
                colorbar = dict( title = "one colorbar" ),
                line = dict( color = '#444' ),
                reversescale = True,
                sizeref = 45,
                sizemode = 'diameter',
                opacity = 0.7,
                size = size,
                color = color,
            ),
        text = '...',
        type = plot_type,
    ) ]

    layout = dict(
        font = dict( family = 'Raleway' ),
        hovermode = 'closest',
        margin = dict( r=20, t=0, l=0, b=0 ),
        showlegend = False,
        scene = dict(
            xaxis = axis_template_3d( xlabel ),
            yaxis = axis_template_3d( ylabel ),
            zaxis = axis_template_3d( zlabel ),
            camera = dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=0.08, y=2.2, z=0.08)
            )
        )
    )

    if plot_type in ['histogram2d', 'scatter']:
        layout['xaxis'] = axis_template_2d(xlabel)
        layout['yaxis'] = axis_template_2d(ylabel)
        layout['plot_bgcolor'] = BACKGROUND
        layout['paper_bgcolor'] = BACKGROUND
        del layout['scene']
        del data[0]['z']

    if plot_type == 'histogram2d':
        # Scatter plot overlay on 2d Histogram
        data[0]['type'] = 'scatter'
        data.append( dict(
                x = x,
                y = y,
                type = 'histogram2d',
                colorscale = 'Greys',
                showscale = False
            ) )
        layout['plot_bgcolor'] = 'black'
        layout['paper_bgcolor'] = 'black'
        layout['xaxis'] = blackout_axis(layout['xaxis'])
        layout['yaxis'] = blackout_axis(layout['yaxis'])
        layout['font']['color'] = 'white'

    if len(markers) > 0:
        data = data + add_markers( data, markers, plot_type = plot_type )

    return dict( data=data, layout=layout )


FIGURE = scatter_plot_3d()



app.layout = html.Div([
    dcc.Input(id='my-id', value='PMA2018-refresh=run model', type="text"),
    html.Button('Load Parameterfile', id='button'),
    html.Div([
	dcc.Graph(id='clickable-graph',
	    style=dict(width='700px'),
	    hoverData=dict( points=[dict(pointNumber=0)] ),
	    figure=FIGURE ),
	], className='nine columns', style=dict(textAlign='center')),
    ])
 






    #return u'new coordinates of CPlantBox is "{}" parameterfile is "{}"'.format(nodes, input_value)
    
    # Initialize
#    plant.initialize()
    # Simulate
#    plant.simulate(60, True)
    # Export final result (as vtp) 2 = root 4 = stem 8 = leaf 15 = all
#    plant.write("example_1a.vtp", 15)
external_css = ["https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
                "//fonts.googleapis.com/css?family=Raleway:400,300,600",
                "//fonts.googleapis.com/css?family=Dosis:Medium",
                "https://cdn.rawgit.com/plotly/dash-app-stylesheets/0e463810ed36927caf20372b6411690692f94819/dash-drug-discovery-demo-stylesheet.css"]


for css in external_css:
    app.css.append_css({"external_url": css})

if __name__ == '__main__':
    app.run_server(debug=True)
