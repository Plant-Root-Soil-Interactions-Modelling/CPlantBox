import dash
from dash import html, dcc, Input, Output, State, ctx
import dash_bootstrap_components as dbc

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.MINTY])
app.title = "CPlantbox Dash App"

plants = [{'label': f'{i}', 'value': f'{i}'} for i in ["Anagallis", "Soybean", "Maize"]]
plant_paramter_files = ["todo"]


tropisms_names = {"Plagiotropism":  0, "Gravitropism":1, "Exotropism": 2, "Hydrotropism": 3}

tabs_names = ["Tap", "First order", "Second order", "Seminal"]

parameter_sliders = {
    "Maximal length [cm]": (1,200),
    "Initisal growh rate [cm/day]": (0.5,10),
    "Basal zone [cm]": (0.1, 20),
    "Interlateral distance [cm]": (0.1,20),
    "Apical zone [cm]": (0.1, 20),
    "Radius [cm]": (1.e-3, 0.25),
    "Tropsim strength [1]": (0, 6),
    "Tropsim tortuosity [1]": (0., 1.),
}

SLIDER_INITIALS = [100, 3, 3, 1, 7, 0.1,1,0.2, "Gravitropism"]


# Layout
app.layout = dbc.Container([
    
    dcc.Store(id='slider-store', data={tab: SLIDER_INITIALS for tab in tabs_names}),  # default values

    dbc.Row([
        # Panel 1
        dbc.Col([
            html.H5("Simulation"),
            html.Br(),
            html.H6("Plant"),
            dcc.Dropdown(
                id='combo-box',
                options=plants, 
                value=plants[0]['value']
            ),
            html.Br(),
            html.H6("Simulation time [day]"),
            dcc.Slider(id='main-slider', min=1, max=180, step=1, value=30,
                       marks={1: "1", 180: "180"}),
            html.Br(),
            dbc.Button("Create", id="create-button"),  # , color="primary", className="mt-2"
        ], width=2),

        # Panel 2
        dbc.Col([
            html.H5("Parameter"),
            html.Br(),
            dcc.Tabs(
                id='tabs-container',
                value='tab-1',
                children=[dcc.Tab(label=f'{name}', value=f'tab-{i}') for i,name in enumerate(tabs_names)]
            ),
            html.Div(id='tab-content')
        ], width=3),

        # Panel 3
        dbc.Col([
            html.H5("Summary"),
            html.Br(),
            html.Div(id='summary-output')
        ], width=4),
    ])
], fluid=True)


# Generate sliders from stored values
def generate_sliders(stored_values):
    sliders = []
    for i, key in enumerate(parameter_sliders.keys()):
        min_ = parameter_sliders[key][0]
        max_ = parameter_sliders[key][1]
        sliders.append(html.H6(key))
        sliders.append(            
            dcc.Slider(
                        id={'type': 'dynamic-slider', 'index': i},
                        min=min_, 
                        max=max_, 
                        value=stored_values[i],
                        marks={
                            min_: str(min_), 
                            max_: str(max_)
                            },
                        tooltip={ "always_visible": False},  
                    )                        
            )
    sliders.append(dcc.Dropdown(
                id={'type': 'dynamic-slider', 'index': i}, # little white lie
                options=list(tropisms_names.keys()), 
                value=stored_values[i+1],
            ))
    return html.Div(sliders)
    


# Display sliders for the selected tab, using stored values
@app.callback(
    Output('tab-content', 'children'),
    Input('tabs-container', 'value'),
    State('slider-store', 'data')
)
def update_tab_content(selected_tab, slider_store):
    values = slider_store.get(selected_tab, SLIDER_INITIALS)
    return generate_sliders(values)


# Update slider-store when any slider changes
@app.callback(
    Output('slider-store', 'data', allow_duplicate=True),
    Input({'type': 'dynamic-slider', 'index': dash.ALL}, 'value'),
    State('tabs-container', 'value'),
    State('slider-store', 'data'),
    prevent_initial_call=True
)
def update_store(slider_values, current_tab, store_data):
    store_data[current_tab] = slider_values
    return store_data


# Create summary output on button press
@app.callback(
    Output('summary-output', 'children'),
    Input('create-button', 'n_clicks'),
    State('tabs-container', 'value'),
    State('combo-box', 'value'),
    State('main-slider', 'value'),
    State('slider-store', 'data'),
    prevent_initial_call=True
)
def summarize(n_clicks, selected_tab, combo_value, main_slider_value, store_data):
    slider_values = store_data.get(selected_tab, [])
    summary = f"""Selected Tab: {selected_tab}
    Combo Box Value: {combo_value}
    Main Slider Value: {main_slider_value}
    Slider Values in {selected_tab}: {slider_values}"""
    return summary


if __name__ == '__main__':
    app.run(debug=True)
