import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State, ALL
import dash_bootstrap_components as dbc

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Layout
app.layout = dbc.Container([
    dcc.Store(id='slider-store', data={}),  # Stores slider values per pick
    
    dbc.Row([
        # First Panel: Dropdown and Create Button
        dbc.Col([
            html.H5("Select a Pick"),
            dcc.Dropdown(
                id='combo-box',
                options=[{'label': f'Pick {i}', 'value': f'pick-{i}'} for i in range(1, 6)],
                value='pick-1'
            ),
            html.Button('Create', id='create-button', n_clicks=0, style={'marginTop': '10px'})
        ], width=3, style={'padding': '20px', 'border': '1px solid #ddd'}),
        
        # Second Panel: Sliders (7 per pick)
        dbc.Col([
            html.H5("Adjust Values"),
            html.Div(id='sliders-container')
        ], width=3, style={'padding': '20px', 'border': '1px solid #ddd'}),
        
        # Third Panel: Output Summary
        dbc.Col([
            html.H5("Summary"),
            html.Div(id='output-summary', style={'whiteSpace': 'pre-wrap'})
        ], width=6, style={'padding': '20px', 'border': '1px solid #ddd'})
    ])
], fluid=True)

# Callback to create sliders based on the selected pick
@app.callback(
    Output('sliders-container', 'children'),
    Input('combo-box', 'value'),
    State('slider-store', 'data')
)
def create_sliders(selected_pick, slider_data):
    # Retrieve stored values or use default values
    stored_values = slider_data.get(selected_pick, [50] * 7)
    return [
        dcc.Slider(
            id={'type': 'slider', 'index': i},
            min=0, max=100, step=10, value=stored_values[i - 1]
        ) for i in range(1, 8)
    ]

# Callback to update the summary and store slider values
@app.callback(
    [Output('output-summary', 'children'),
     Output('slider-store', 'data')],
    Input('create-button', 'n_clicks'),
    State('combo-box', 'value'),
    State({'type': 'slider', 'index': ALL}, 'value'),
    State('slider-store', 'data')
)
def update_summary(n_clicks, selected_pick, slider_values, slider_data):
    if n_clicks > 0 and slider_values:
        slider_data[selected_pick] = slider_values  # Store values per pick
        return f'Selected Pick: {selected_pick}\nSlider Values: {slider_values}', slider_data
    return '', slider_data

if __name__ == '__main__':
    app.run(debug=True)
