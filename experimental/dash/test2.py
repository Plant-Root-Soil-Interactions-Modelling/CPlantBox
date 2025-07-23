import dash
from dash import html, dcc, Input, Output, State
import dash_bootstrap_components as dbc

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SLATE])  # Dark Theme
app.title = "Modern Dash Layout"

# Store for maintaining slider values
TABS = [f"tab-{i}" for i in range(1, 6)]
SLIDER_COUNT = 7

app.layout = html.Div([
    dcc.Store(id='slider-store', data={tab: [5] * SLIDER_COUNT for tab in TABS}),  # Stores slider values

    # 🔹 Sidebar Panel (Panel 1)
    html.Div([
        html.H3("Controls", className="title"),
        dcc.Dropdown(
            id='combo-box',
            options=[{'label': f'Option {i}', 'value': f'Option {i}'} for i in range(1, 6)],
            value='Option 1',
            className="dropdown"
        ),
        dcc.Slider(id='main-slider', min=0, max=100, step=1, value=50,
                   marks={i: str(i) for i in range(0, 101, 10)}, className="slider"),
        dbc.Button("Create", id="create-button", color="primary", className="button"),
    ], className="sidebar"),

    # 🔹 Main Content (Panel 2 & 3)
    html.Div([
        # Tabs and sliders
        html.Div([
            html.H3("Adjust Values", className="subtitle"),
            dcc.Tabs(
                id='tabs-container',
                value='tab-1',
                children=[dcc.Tab(label=f'Tab {i}', value=f'tab-{i}') for i in range(1, 6)],
                className="custom-tabs"
            ),
            html.Div(id='tab-content', className="tab-content")
        ], className="tab-container"),

        # Summary Panel
        html.Div([
            html.H3("Summary", className="subtitle"),
            html.Div(id='summary-output', className="summary-box")
        ], className="summary-container"),
    ], className="main-content")
], className="container")


# 🔹 Generate sliders dynamically
def generate_sliders(tab_index, stored_values):
    return html.Div([
        dcc.Slider(
            id={'type': 'dynamic-slider', 'index': i},
            min=0, max=10, step=1,
            value=stored_values[i],
            marks={i: str(i) for i in range(11)},
            className="slider"
        ) for i in range(SLIDER_COUNT)
    ], className="sliders")


# 🔹 Callback to update tab content dynamically
@app.callback(
    Output('tab-content', 'children'),
    Input('tabs-container', 'value'),
    State('slider-store', 'data')
)
def update_tab_content(selected_tab, slider_store):
    return generate_sliders(selected_tab, slider_store.get(selected_tab, [5] * SLIDER_COUNT))


# 🔹 Store slider values persistently
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


# 🔹 Generate summary on button press
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
    return f"""
    Selected Tab: {selected_tab}
    Combo Box Value: {combo_value}
    Main Slider Value: {main_slider_value}
    Slider Values in {selected_tab}: {slider_values}
    """


if __name__ == '__main__':
    app.run(debug=True)
