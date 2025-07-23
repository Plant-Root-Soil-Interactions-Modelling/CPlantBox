import dash
from dash import dcc, html, Input, Output

app = dash.Dash(__name__)

app.layout = html.Div([
    dcc.Tabs(
        id="tabs",
        value="tab-1",
        children=[
            dcc.Tab(label="Tab 1", value="tab-1"),
            dcc.Tab(label="Tab 2", value="tab-2"),
            dcc.Tab(label="Tab 3", value="tab-3"),
        ],
    ),
    html.Button("Next Tab", id="next-tab-btn", n_clicks=0),
    html.Div(id="tab-content"),
])


@app.callback(
    Output("tabs", "value"),  # Change the active tab
    Input("next-tab-btn", "n_clicks")
)
def switch_tab(n_clicks):
    tab_order = ["tab-1", "tab-2", "tab-3"]
    return tab_order[n_clicks % len(tab_order)]  # Cycle through tabs


if __name__ == "__main__":
    app.run(debug=True)
