import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State

app = dash.Dash(__name__)

app.layout = html.Div([
    dcc.Slider(id='slider-1', min=0, max=100, step=1, value=50,
               marks={i: str(i) for i in range(0, 101, 20)}),
    dcc.Slider(id='slider-2', min=0, max=100, step=1, value=50,
               marks={i: str(i) for i in range(0, 101, 20)}),
    dcc.Slider(id='slider-3', min=0, max=100, step=1, value=50,
               marks={i: str(i) for i in range(0, 101, 20)}),
    html.Button('Show Values', id='submit-button', n_clicks=0),
    html.Div(id='output-text', style={'marginTop': 20, 'fontSize': 20})
])

@app.callback(
    Output('output-text', 'children'),
    Input('submit-button', 'n_clicks'),
    State('slider-1', 'value'),
    State('slider-2', 'value'),
    State('slider-3', 'value')
)
def update_output(n_clicks, val1, val2, val3):
    if n_clicks > 0:
        return f'Slider Values: {val1}, {val2}, {val3}'
    return ''

if __name__ == '__main__':
    app.run(debug=True)