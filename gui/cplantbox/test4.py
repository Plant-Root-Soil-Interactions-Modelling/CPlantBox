import dash
from dash import html, dcc, Output, Input
import time

app = dash.Dash(__name__)

app.layout = html.Div([
    dcc.Input(id = 'input', value = 'initial value', type = 'text'),
    html.Br(),
    dcc.Loading(
        id = "loading-spinner",
        type = "circle",  # Options: "default", "circle", "dot", "cube"
        children = html.Div(id = "output")
    )
])


@app.callback(
    Output("output", "children"),
    Input("input", "value")
)
def update_output(value):
    time.sleep(2)  # Simulate a delay
    return f"You entered: {value}"


if __name__ == "__main__":
    app.run_server(debug = True)
