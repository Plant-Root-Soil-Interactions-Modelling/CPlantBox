import dash
from dash import dcc, html, Output, Input
import dash_bootstrap_components as dbc

app = dash.Dash(__name__, suppress_callback_exceptions = True, external_stylesheets = [dbc.themes.BOOTSTRAP])

app.layout = dbc.Container([
    html.H2("Nested Tabs Example", className = "my-4"),

    dcc.Tabs(id = 'main-tabs', value = 'A', children = [
        dcc.Tab(label = 'Tab A', value = 'A'),
        dcc.Tab(label = 'Tab B', value = 'B')
    ]),

    html.Div(id = 'main-tabs-content', className = "mt-3")
], fluid = True)


# Top-level main tab callback — generates sub-tab layout dynamically
@app.callback(
    Output('main-tabs-content', 'children'),
    Input('main-tabs', 'value')
)
def render_main_tab(tab):
    if tab == 'A':
        return html.Div([
            dcc.Tabs(id = 'sub-tabs-a', value = 'a', children = [
                dcc.Tab(label = 'Sub-tab a', value = 'a'),
                dcc.Tab(label = 'Sub-tab b', value = 'b'),
                dcc.Tab(label = 'Sub-tab c', value = 'c')
            ]),
            html.Div(id = 'sub-tabs-a-content', className = "mt-3")
        ])
    elif tab == 'B':
        return html.Div([
            dcc.Tabs(id = 'sub-tabs-b', value = 'd', children = [
                dcc.Tab(label = 'Sub-tab d', value = 'd'),
                dcc.Tab(label = 'Sub-tab e', value = 'e')
            ]),
            html.Div(id = 'sub-tabs-b-content', className = "mt-3")
        ])


# Dynamic sub-tab content callbacks must be registered despite conditional layout
@app.callback(
    Output('sub-tabs-a-content', 'children'),
    Input('sub-tabs-a', 'value'),
    prevent_initial_call = True
)
def render_subtab_a(subtab):
    return html.Div([
        html.H5("Summary"),
        html.P(f"Selected Main Tab: A"),
        html.P(f"Selected Sub-tab: {subtab}")
    ])


@app.callback(
    Output('sub-tabs-b-content', 'children'),
    Input('sub-tabs-b', 'value'),
    prevent_initial_call = True
)
def render_subtab_b(subtab):
    return html.Div([
        html.H5("Summary"),
        html.P(f"Selected Main Tab: B"),
        html.P(f"Selected Sub-tab: {subtab}")
    ])


if __name__ == '__main__':
    app.run_server(debug = True)
