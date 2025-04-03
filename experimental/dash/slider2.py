import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State, MATCH, ALL
import dash_bootstrap_components as dbc

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Layout
app.layout = dbc.Container([
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

# Callbacks to update sliders and display results
@app.callback(
    Output('sliders-container', 'children'),
    Input('combo-box', 'value')
)
def create_sliders(value):
    return [dcc.Slider(id={'type': 'slider', 'index': i}, min=0, max=100, step=10, value=50) for i in range(1, 8)]

@app.callback(
    Output('output-summary', 'children'),
    Input('create-button', 'n_clicks'),
    State({'type': 'slider', 'index': ALL}, 'value')
)
def update_summary(n_clicks, slider_values):
    if n_clicks > 0 and slider_values:
        return f'Selected Values: {slider_values}'
    return ''

if __name__ == '__main__':
    app.run(debug=True)
