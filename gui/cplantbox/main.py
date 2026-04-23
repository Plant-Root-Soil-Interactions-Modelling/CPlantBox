"""CPlantBox Webapp (using Python dash), D. Leitner 2026"""

import base64
import os
import uuid
import webbrowser
from threading import Timer  # to open browser automatically ,see __main__

import bibtexxml_apa  # APA citation helper
import conversions  # auxiliary stuff
import dash
import dash_bootstrap_components as dbc
import numpy as np
import plots  # figures
import simulate_plant  # the simulation loop
import vtk_conversions  # auxiliary stuff
from dash import Input, Output, State, ctx, dcc, html, no_update
from pympler import asizeof

MAX_CACHED_RUNS = 16
RUN_CACHE = {}
RUN_CACHE_ORDER = []


def prepare_vtk_render_data(vtk_data):
    """Decode geometry once so 3D tab switches do not rebuild large Python lists repeatedly."""
    points = vtk_conversions.decode_array(vtk_data["points"])
    point_coords = points.reshape(-1, 3)
    if point_coords.size == 0:
        center = np.zeros(3)
        radius = 1.0
    else:
        center = point_coords.mean(axis=0)
        radius = np.linalg.norm(point_coords - center, axis=1).max()
        if radius <= 0:
            radius = 1.0
    distance = 1.5 * radius

    return {
        "points": points.flatten().tolist(),
        "polys": vtk_conversions.decode_array(vtk_data["polys"]).tolist(),
        "leaf_points": vtk_conversions.decode_array(vtk_data["leaf_points"]).flatten().tolist(),
        "leaf_polys": vtk_conversions.decode_array(vtk_data["leaf_polys"]).tolist(),
        "sub_type_colors": np.repeat(vtk_conversions.decode_array(vtk_data["subType"]), 16).tolist(),
        "age_colors": np.repeat(vtk_conversions.decode_array(vtk_data["creationTime"]), 16).tolist(),
        "camera_props": {
            "cameraPosition": (center + np.array([-distance, -distance, distance])).tolist(),
            "cameraViewUp": [0, 0, 1],
        },
    }


def cache_simulation_run(vtk_data, result_data):
    """Store one simulation result server-side and return its run id."""
    run_id = uuid.uuid4().hex
    RUN_CACHE[run_id] = {
        "vtk_data": vtk_data,
        "result_data": result_data,
        "vtk_render_data": prepare_vtk_render_data(vtk_data),
    }
    RUN_CACHE_ORDER.append(run_id)
    while len(RUN_CACHE_ORDER) > MAX_CACHED_RUNS:
        oldest = RUN_CACHE_ORDER.pop(0)
        RUN_CACHE.pop(oldest, None)
    return run_id


def get_cached_simulation_run(run_id):
    if not run_id:
        return None
    return RUN_CACHE.get(run_id)


def open_browser():
    webbrowser.open_new("http://127.0.0.1:8050")


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MD_PATH = os.path.join(BASE_DIR, "assets", "readme.md")  # Path to your Markdown file in the assets folder
with open(MD_PATH, "r", encoding="utf-8") as f:  # Read the Markdown content
    ABOUT_TEXT = f.read()
XML_INIT_PATH = os.path.join(BASE_DIR, "params", "root_only.xml")
with open(XML_INIT_PATH, "r", encoding="utf-8") as f:
    XML_INIT = f.read()
MAX_XML_SIZE = 200 * 1024  # 200 KB in bytes

#
# INITIALIZE
#
app = dash.Dash(
    __name__,
    suppress_callback_exceptions=True,
    external_stylesheets=[dbc.themes.MINTY],  # SANDSTONE, MINTY, MORPH
    meta_tags=[
        {
            "name": "description",
            "content": (
                "CPlantBox WebApp — an interactive interface to visualize and simulate plant architecture. "
                "CPlantBox creates plant geometry based on parameter sets and supports model coupling for "
                "integrated plant-soil-atmosphere interaction models, including water and carbon flow, "
                "photosynthesis, root-soil interaction, and mycorrhizal associations."
            ),
        },
        {
            "name": "keywords",
            "content": (
                "CPlantBox, plant modelling, root system simulation, plant architecture, "
                "soil-plant-atmosphere, water flow, carbon flow, plant parameters, "
                "functional-structural plant model, dumux-rosi, mycorrhizal, photosynthesis, "
                "CRootBox, root-soil interaction, Python plant simulation"
            ),
        },
        {
            "name": "author",
            "content": "Daniel Leitner, IBG-3, Forschungszentrum Jülich",
        },
        {"property": "og:title", "content": "CPlantBox WebApp"},
        {
            "property": "og:description",
            "content": (
                "Interactive plant architecture simulation — visualize and experiment with " "plant parameters for soil-plant-atmosphere interaction models."
            ),
        },
        {"property": "og:type", "content": "website"},
    ],
)
app.title = "CPlantBox WebApp"
app.index_string = """
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <script data-goatcounter="https://cplantboxwebapp.goatcounter.com/count" async src="//gc.zgo.at/count.js"></script>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
"""

_sep_style = {"fontStyle": "italic", "color": "gray"}
param_names = conversions.get_parameter_names()
plants = []
for i, name in enumerate(param_names):
    if i == 0:
        plants.append({"label": html.Span("- Dummy data -", style=_sep_style), "value": "separator-demo", "disabled": True})
    if i == 4:
        plants.append({"label": html.Span("- Plants -", style=_sep_style), "value": "separator-plants", "disabled": True})
    if i == 6:
        plants.append({"label": html.Span("- Root systems -", style=_sep_style), "value": "separator-roots", "disabled": True})
    plants.append({"label": name[0], "value": str(i)})

seed_parameter_sliders = conversions.get_seed_slider_names()
root_parameter_sliders = conversions.get_root_slider_names()
stem_parameter_sliders = conversions.get_stem_slider_names()
leaf_parameter_sliders = conversions.get_leaf_slider_names()

SEED_SLIDER_INITIALS = [180, 1, 7, 7, 20, 7, 7, 15, 5]
ROOT_SLIDER_INITIALS = [100, 3, 45, 1, 1, 7, 0.1, 1, 0.2, "Gravitropism"]
STEM_SLIDER_INITIALS = [100, 3, 45, 1, 0.1, 7, 14, 180.0, 1, 0.2, "Gravitropism"]
LEAF_SLIDER_INITIALS = ["Defined", 30, 1, 45.0, 2.0, 90, 1.0, 0.2, "Gravitropism"]


def small_button(label, id_, tool_tip="", url="cloud-download.svg"):
    """creates a small button with an icon from the assets folder (from Bootstrap Icons)"""
    return dcc.Button([html.Img(src=app.get_asset_url(url), className="buttonIcon"), label], id=id_, title=tool_tip, className="smallButton")


#
# LAYOUT
#
app.layout = dbc.Container(
    [
        dcc.Store(
            id="seed-store", data={"seed": SEED_SLIDER_INITIALS, "basal-checkbox": False, "shoot-checkbox": False, "tillers-checkbox": False}
        ),  # seed slider values
        dcc.Store(id="root-store", data={f"tab-{i}": ROOT_SLIDER_INITIALS for i in range(1, 5)}),  # root slider values
        dcc.Store(id="stem-store", data={f"tab-{i}": STEM_SLIDER_INITIALS for i in range(1, 5)}),  # stem slider values
        dcc.Store(id="leaf-store", data={"leaf": LEAF_SLIDER_INITIALS}),  # leaf slider values
        dcc.Store(id="typename-store", data={f"tab-{i}": f"Order {i} root" for i in range(1, 5)}),  # root sub type names
        dcc.Store(id="settings-store", data={"token": 0, "reset": True, "random_seed": 0}),  # further settings
        dcc.Store(id="run-id-store", data={"run_id": None}),  # active server-side simulation run id
        dcc.Store(id="xml-store", data={"xml": XML_INIT}),  # for uploaded xml content
        dcc.Download(id="download-xml"),
        dcc.Download(id="download-vtk"),
        dcc.Download(id="download-rsml"),
        dcc.Download(id="download-profiles-xls"),
        dcc.Download(id="download-dynamics-xls"),
        html.Div(
            dbc.Row(
                [
                    #
                    # Panel 1
                    #
                    dbc.Col(
                        [
                            html.H5("Simulation"),
                            conversions.into_panel(
                                [
                                    html.H6("Plant parameters"),
                                    dcc.Dropdown(
                                        id="plant-dropdown",
                                        options=plants,
                                        value=plants[1]["value"],
                                        clearable=False,
                                        searchable=False,
                                        className="dropdown",
                                        style={"fontSize": "12px", "padding-top": "5px"},  # hard coded (should be same as h6) did not work with css allown
                                    ),
                                    html.Div(className="smallSpacer"),
                                    html.Div(
                                        [
                                            html.Div(
                                                children=small_button("xml", "xml-download-button", "Download the XML parameter file")
                                            ),  # ohterwise the layout is a pixel wrong
                                            dcc.Upload(
                                                id="xml-upload-button",
                                                accept=".xml",
                                                multiple=False,
                                                children=small_button(
                                                    "xml", "xml-upload-button-hidden", "Upload your XML parameter file", url="cloud-upload.svg"
                                                ),
                                            ),
                                        ],
                                        style={"display": "flex", "alignItems": "center"},  # vertical alignment
                                    ),
                                ]
                            ),
                            html.Div(className="spacer"),
                            conversions.into_panel(
                                [
                                    html.H6("Simulation time [day]"),
                                    dcc.Slider(id="time-slider", min=1, max=45, step=1, value=20, marks={1: "1", 45: "45"}, tooltip={"always_visible": False}),
                                    html.Div(className="spacer"),
                                    html.Div(
                                        [
                                            dcc.Button("Create", id="create-button", title="Create new plant geometry", className="button"),
                                        ],
                                    ),
                                ]
                            ),
                            html.Div(className="spacer"),
                            html.Div(id="reference-panel"),
                            html.Div(className="spacer"),
                            dcc.Loading(id="loading-spinner", type="circle", children=html.Div(id="loading-spinner-output")),
                            dcc.Loading(id="loading-spinner2", type="circle", children=html.Div(id="loading-spinner-output2")),
                            dcc.Loading(id="loading-spinner3", type="circle", children=html.Div(id="loading-spinner-output3")),
                        ],
                        width=2,
                        style={"minWidth": "160px"},
                    ),
                    #
                    # Panel 2
                    #
                    dbc.Col(
                        [
                            html.H5("Parameters"),
                            dcc.Tabs(
                                id="organtype-tabs",
                                value="Seed",
                                children=[
                                    dcc.Tab(label="Seed", value="Seed", className="tab", selected_className="tabSelected"),
                                    dcc.Tab(label="Root", value="Root", className="tab", selected_className="tabSelected"),
                                    dcc.Tab(label="Stem", value="Stem", className="tab", selected_className="tabSelected"),
                                    dcc.Tab(label="Leaf", value="Leaf", className="tab", selected_className="tabSelected"),
                                ],
                                className="tabs",
                            ),
                            html.Div(id="organtype-tabs-content"),  # , className="tabContentScroll"
                        ],
                        width=3,
                        style={"minWidth": "250px"},
                    ),
                    #
                    # Panel 3
                    #
                    dbc.Col(
                        [
                            html.H5("Results"),
                            dcc.Tabs(
                                id="result-tabs",
                                value="VTK3D",
                                children=[
                                    dcc.Tab(label="3D", value="VTK3D", className="tab", selected_className="tabSelected"),
                                    dcc.Tab(label="3D Age", value="VTK3DAge", className="tab", selected_className="tabSelected"),
                                    dcc.Tab(label="1D Profiles", value="Profile1D", className="tab", selected_className="tabSelected"),
                                    dcc.Tab(label="Dynamics", value="Dynamics", className="tab", selected_className="tabSelected"),
                                    dcc.Tab(label="About", value="About", className="tab", selected_className="tabSelected"),
                                ],
                                className="tabs",
                            ),
                            html.Div(id="result-tabs-content"),
                        ],
                        width=6,
                        style={"minWidth": "320px"},
                    ),
                ],
                className="mainLayoutRow",
            ),
            className="mainLayoutScroll",
        ),
        html.Div(
            [
                html.A(
                    html.Img(src="/assets/cplantbox.png", className="logo"),
                    href="https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox",
                    target="_blank",
                ),
                html.A(html.Img(src="/assets/fzj.png", className="logo"), href="https://www.fz-juelich.de/de", target="_blank"),
            ],
            className="logoContainer",
        ),
    ],
    fluid=True,
)


#
# 1. LEFT - Simulation Panel
#
@app.callback(  # Plant parameters: plant-dropdown
    Output("seed-store", "data"),
    Output("root-store", "data"),
    Output("stem-store", "data"),
    Output("leaf-store", "data"),
    Output("typename-store", "data"),
    Output("organtype-tabs", "value"),  # to trigger update
    Output("create-button", "n_clicks"),  # to simulation
    Output("time-slider", "value"),
    Input("plant-dropdown", "value"),
    State("seed-store", "data"),
    State("root-store", "data"),
    State("stem-store", "data"),
    State("leaf-store", "data"),
    State("typename-store", "data"),
    State("organtype-tabs", "value"),
    State("xml-store", "data"),
)
def plant_dropdown(plant_value, seed_data, root_data, stem_data, leaf_data, typename_data, tabs_value, xml_data):
    triggered = ctx.triggered_id
    print("[plant_dropdown()", plant_value, "triggered:", triggered)  # dynamic-seed-slider
    conversions.set_data(plant_value, seed_data, root_data, stem_data, leaf_data, typename_data, xml_data)
    print("]plant_dropdown()", plant_value)  # seed_data
    return (seed_data, root_data, stem_data, leaf_data, typename_data, tabs_value, None, seed_data["simulationTime"])


@app.callback(  # Create button and time slider
    Output("result-tabs-content", "children"),
    Output("run-id-store", "data"),
    Output("settings-store", "data"),
    Output("loading-spinner-output", "children"),
    Input("create-button", "n_clicks"),
    Input("time-slider", "value"),
    State("plant-dropdown", "value"),
    State("seed-store", "data"),
    State("root-store", "data"),
    State("stem-store", "data"),
    State("leaf-store", "data"),
    State("typename-store", "data"),
    State("result-tabs", "value"),
    State("settings-store", "data"),
    State("xml-store", "data"),
)
def handle_simulation(
    create_clicks, time_slider, plant_value, seed_data, root_data, stem_data, leaf_data, typename_data, result_value, settings_data, xml_data
):
    triggered = ctx.triggered_id
    print("**[handle_simulation()", triggered, create_clicks, time_slider)

    # ---- CREATE BUTTON ----
    if triggered is None or triggered == "create-button":
        settings_data["token"] += 1
        settings_data["reset"] = True
        rng = np.random.default_rng()
        settings_data["random_seed"] = rng.integers(1, 10001)

    # ---- TIME SLIDER ----
    elif triggered == "time-slider":
        settings_data["reset"] = False

    # ---- Run simulation (common part) ----
    vtk_data, result_data = simulate_plant.simulate_plant(
        plant_value, time_slider, seed_data, root_data, stem_data, leaf_data, settings_data["random_seed"], xml_data
    )
    run_id = cache_simulation_run(vtk_data, result_data)

    content = render_result_content(result_value, run_id, typename_data, settings_data)

    print("**]handle_simulation()", "\n\n")

    return (content, {"run_id": run_id}, settings_data, html.H6(""))


@app.callback(  # parameter file download
    Output("download-xml", "data"),
    Input("xml-download-button", "n_clicks"),
    State("plant-dropdown", "value"),
    State("seed-store", "data"),
    State("root-store", "data"),
    State("stem-store", "data"),
    State("leaf-store", "data"),
    State("xml-store", "data"),
    prevent_initial_call=True,
)
def download_xml(n_clicks, plant_value, seed_data, root_data, stem_data, leaf_data, xml_data):
    xml_string = simulate_plant.get_xml(plant_value, seed_data, root_data, stem_data, leaf_data, xml_data)
    return dict(content=xml_string, filename="cplantbox_parameters.xml", type="application/xml")


@app.callback(  # parameter file upload
    Output("xml-store", "data"),
    Output("plant-dropdown", "value"),
    Input("xml-upload-button", "contents"),
    State("xml-store", "data"),
    prevent_initial_call=True,
)
def handle_xml_upload(contents, data):
    print("handle_xml_upload")
    if contents is None:
        return dash.no_update, dash.no_update
    content_type, content_string = contents.split(",")
    decoded = base64.b64decode(content_string)
    if len(decoded) > MAX_XML_SIZE:  # Size check (before decoding to UTF-8 string)
        print("XML file exceeds size limit")
        return dash.no_update, dash.no_update  # or raise PreventUpdate
    xml_string = decoded.decode("utf-8")
    data["xml"] = xml_string  # Store the XML content in the dcc.Store
    user_data_value = next((str(i) for i, name in enumerate(conversions.get_parameter_names()) if name[0] == "User Data"), None)
    return data, user_data_value if user_data_value is not None else dash.no_update


#
# 2. MID - Parameters Panel
#
@app.callback(  # organtype-tabs
    Output("organtype-tabs-content", "children"),
    Input("organtype-tabs", "value"),
    State("seed-store", "data"),
    State("root-store", "data"),
    State("typename-store", "data"),  # root and stem sub-type names
    State("stem-store", "data"),
    State("leaf-store", "data"),
)
def render_organtype_tab(tab, seed_data, root_data, type_names, stem_data, leaf_data):

    triggered = ctx.triggered_id
    # print("render_organtype_tab()", triggered, "tab", tab)
    if triggered is None:
        tab = "Seed"
    if tab == "Seed":
        # print("render_organtype_tab() seed:", seed_data)
        return generate_seed_sliders(seed_data)
    elif tab == "Root":
        # print("render_organtype_tab() root:", root_data)
        return root_layout(root_data, type_names)
    elif tab == "Stem":
        # print("render_organtype_tab() stem:", stem_data)
        return stem_layout(stem_data, type_names)
    elif tab == "Leaf":
        # print("render_organtype_tab() leaf:", leaf_data)
        return generate_leaf_sliders(leaf_data)

    raise ValueError("Unknown tab selected")


#
#  Parameters Panel - Seed
#
def generate_seed_sliders(data):  # Generate sliders for seed tab from stored values
    # print("generate_seed_sliders()", seed_values)
    seed_values = data["seed"]
    sliders = []
    for i, key in enumerate(seed_parameter_sliders.keys()):
        #     if False:  # in case we want to hide certain sliders
        #         style = {"display": "none"}
        #     else:
        #         style = {}  # default
        min_ = seed_parameter_sliders[key][0]
        max_ = seed_parameter_sliders[key][1]
        step_ = seed_parameter_sliders[key][2]
        sliders.append(html.H6(key))  # , style=style
        sliders.append(
            html.Div(
                [
                    dcc.Slider(
                        id={"type": "dynamic-seed-slider", "index": i},
                        min=min_,
                        max=max_,
                        value=seed_values[i],
                        marks={min_ + 1.0e-4: str(min_), max_ - 1.0e-4: str(max_)},
                        tooltip={"always_visible": False},
                        step=step_,
                        className="slider",
                    )
                ],
                # style=style,
            )
        )
    # Insert checkboxes

    v = ["agree"] if data["shoot-checkbox"] else []
    sliders.insert(
        0,
        html.Div(
            [
                dcc.Checklist(
                    id="shoot-checkbox",
                    options=[{"label": "Shoot borne roots", "value": "agree"}],
                    className="checkbox",
                    value=v,
                    style={"fontSize": "12px", "padding-top": "0px"},
                )
            ]
        ),
    )
    v = ["agree"] if data["basal-checkbox"] else []
    sliders.insert(
        1 + 2 * 2,
        dcc.Checklist(
            id="basal-checkbox",
            options=[{"label": "Basal roots", "value": "agree"}],
            className="checkbox",
            value=v,
            style={"fontSize": "12px", "padding-top": "0px"},
        ),
    )
    v = ["agree"] if data["tillers-checkbox"] else []
    sliders.insert(
        2 + 5 * 2,
        dcc.Checklist(
            id="tillers-checkbox",
            options=[{"label": "Tillers", "value": "agree"}],
            className="checkbox",
            value=v,
            style={"fontSize": "12px", "padding-top": "0px"},
        ),
    )
    panel1 = conversions.into_panel(sliders, range(0, 1 + 2 * 2))  # shoot
    panel2 = conversions.into_panel(sliders, range(1 + 2 * 2, 2 + 5 * 2))  # basal
    panel3 = conversions.into_panel(sliders, range(2 + 5 * 2, len(sliders)))  # tillers
    return html.Div([panel1, panel2, panel3])


@app.callback(  # store dynamic-seed-slider in seed-store
    Output("seed-store", "data", allow_duplicate=True),
    Input({"type": "dynamic-seed-slider", "index": dash.ALL}, "value"),
    State("seed-store", "data"),
    prevent_initial_call=True,
)
def update_seed_store(slider_values, data):
    triggered = ctx.triggered_id
    # print("update_seed_store()", triggered, ", ", data, slider_values, len(data["seed"]))
    if len(slider_values) > 1:  # called empty
        data["seed"] = slider_values
    return data


@app.callback(  # shoot-checkbox
    Output({"type": "dynamic-seed-slider", "index": 0}, "disabled"),
    Output({"type": "dynamic-seed-slider", "index": 1}, "disabled"),
    Output("seed-store", "data", allow_duplicate=True),
    Input("shoot-checkbox", "value"),
    State("seed-store", "data"),
    prevent_initial_call=True,
)
def toggle_shoot_checkbox(checkbox_value, data):
    triggered = ctx.triggered_id
    # print("toggle_shoot_checkbox()", triggered, checkbox_value)
    data["shoot-checkbox"] = "agree" in checkbox_value
    disabled = "agree" not in checkbox_value
    return (disabled, disabled, data)  # disable if not checked


@app.callback(  # basal-checkbox
    Output({"type": "dynamic-seed-slider", "index": 2}, "disabled"),
    Output({"type": "dynamic-seed-slider", "index": 3}, "disabled"),
    Output({"type": "dynamic-seed-slider", "index": 4}, "disabled"),
    Output("seed-store", "data", allow_duplicate=True),
    Input("basal-checkbox", "value"),
    State("seed-store", "data"),
    prevent_initial_call=True,
)
def toggle_basal_checkbox(checkbox_value, data):
    # triggered = ctx.triggered_id
    # print("toggle_basal_checkbox()", triggered, checkbox_value)
    data["basal-checkbox"] = "agree" in checkbox_value
    disabled = "agree" not in checkbox_value
    return (disabled, disabled, disabled, data)  # disable if not checked


@app.callback(  # tillers-checkbox
    Output({"type": "dynamic-seed-slider", "index": 7}, "disabled"),
    Output({"type": "dynamic-seed-slider", "index": 5}, "disabled"),
    Output({"type": "dynamic-seed-slider", "index": 6}, "disabled"),
    Output("seed-store", "data", allow_duplicate=True),
    Input("tillers-checkbox", "value"),
    State("seed-store", "data"),
    prevent_initial_call=True,
)
def toggle_tillers_checkbox(checkbox_value, data):
    triggered = ctx.triggered_id
    # print("toggle_tillers_checkbox()", triggered, checkbox_value)
    data["tillers-checkbox"] = "agree" in checkbox_value
    disabled = "agree" not in checkbox_value
    return (disabled, disabled, disabled, data)


#
# Parameters Panel - Root
#
def generate_root_sliders(root_values, tab):  # Generate sliders for root tabs from stored values
    """@param root_values is root_data[current_tab]"""
    # print("generate_root_sliders()", root_values)
    style = {}
    successors = root_values[-1]
    sliders = []
    for i, key in enumerate(root_parameter_sliders.keys()):
        style = {}
        if (not successors) and (i in [3, 4, 5]):  # no lb, ln, la for highest order
            style = {"display": "none"}
        # if (tab == 1) and (i == 2):  # no initial angle for tap
        #     style = {"display": "none"}
        min_ = root_parameter_sliders[key][0]
        max_ = root_parameter_sliders[key][1]
        step_ = root_parameter_sliders[key][2]
        sliders.append(html.H6(key, style=style))
        # print(key, str(min_), str(max_), min_, max_, "value", root_values[i])
        sliders.append(
            html.Div(
                [
                    dcc.Slider(
                        id={"type": "root-dynamic-slider", "index": i},
                        min=min_,
                        max=max_,
                        value=root_values[i],
                        marks={min_ + 1.0e-4: str(min_), max_ - 1.0e-4: str(max_)},
                        tooltip={"always_visible": False},
                        step=step_,
                    )
                ],
                style=style,
            )
        )
    sliders.append(
        dcc.Dropdown(
            id={"type": "root-dynamic-slider", "index": i},  # little white lie
            options=list(conversions.tropism_names.keys()),
            value=root_values[i + 1],
            clearable=False,
            className="dropdown",
            style={"fontSize": "12px", "padding-top": "5px"},
        )
    )
    panel2 = conversions.into_panel(sliders, range(0, len(sliders) - (2 * 2 + 1)))  # all but tropism
    panel3 = conversions.into_panel(sliders, range(len(sliders) - (2 * 2 + 1), len(sliders)))  # tropism
    return html.Div([panel2, panel3])


def root_layout(data, type_names):
    """root tab layout: with root subTypes as sub tabs"""
    rootPanelChildren = []
    for i in range(0, type_names["number_roottypes"]):
        ctab = dcc.Tab(
            label=type_names[f"root tab-{i+1}"],
            value=f"tab-{i+1}",
            className="tab",
            selected_className="tabSelected",
        )
        rootPanelChildren.append(ctab)
    tab = dcc.Tabs(
        id="root-tabs",
        value="tab-1",
        children=rootPanelChildren,
    )
    content = html.Div(id="root-tabs-content")
    return [tab, content]


@app.callback(  # render_root_tab - Display sliders for the selected tab, using stored values
    Output("root-tabs-content", "children"),
    Input("root-tabs", "value"),
    State("root-store", "data"),
    suppress_callback_exceptions=True,
)
def render_root_tab(selected_tab, root_data):
    stored_values = root_data.get(selected_tab, ROOT_SLIDER_INITIALS)
    print("render_root_tab()", selected_tab, root_data[selected_tab], stored_values)
    return generate_root_sliders(stored_values, int(selected_tab[-1]))


@app.callback(  # Update root-store when any slider changes
    Output("root-store", "data", allow_duplicate=True),  #
    Input({"type": "root-dynamic-slider", "index": dash.ALL}, "value"),
    State("root-tabs", "value"),
    State("root-store", "data"),
    prevent_initial_call=True,
)
def update_root_store(slider_values, current_tab, store_data):
    # print("update_root_store()", store_data[current_tab])
    successors = store_data[current_tab][-1]
    store_data[current_tab] = slider_values
    store_data[current_tab].append(successors)
    # print("update_root_store() to  ", store_data[current_tab])
    return store_data


#
# Parameters Panel - Stem
#
def generate_stem_sliders(stem_values, tab):  # Generate sliders for stem tabs from stored values
    sliders = []
    successors = stem_values[-1]
    hidden_without_successor = [3, 5, 6, 7]  # phytomer distance, nodal growth start/time, fixed rotation
    sliders.append(html.Div(className="spacer"))
    for i, key in enumerate(stem_parameter_sliders.keys()):
        style = {}
        if i in [5, 6]:  # nodal growth start/time
            style = {"display": "none"}
        if (not successors) and (i in hidden_without_successor):
            style = {"display": "none"}
        # if (tab == 1) and (i == 2):  # no initial theta for main stem; makes sense for tillers
        #     style = {'display': 'none'}
        min_ = stem_parameter_sliders[key][0]
        max_ = stem_parameter_sliders[key][1]
        step_ = stem_parameter_sliders[key][2]
        sliders.append(html.H6(key, style=style))
        sliders.append(
            html.Div(
                [
                    dcc.Slider(
                        id={"type": "stem-dynamic-slider", "index": i},
                        min=min_,
                        max=max_,
                        value=stem_values[i],
                        marks={min_: str(min_), max_: str(max_)},
                        tooltip={"always_visible": False},
                        step=step_,
                    )
                ],
                style=style,
            )
        )
    sliders.append(
        dcc.Dropdown(
            id={
                "type": "stem-dynamic-slider",
                "index": i,
            },  # little white lie - before we move it top, we need to fix that id
            options=list(conversions.tropism_names.keys()),
            value=stem_values[i + 1],
            clearable=False,
            className="customDropdown",
            style={"fontSize": "12px", "padding-top": "5px"},
        )
    )
    panel2 = conversions.into_panel(sliders, range(0, len(sliders) - (2 * 2 + 1)))  # all but tropism
    panel3 = conversions.into_panel(sliders, range(len(sliders) - (2 * 2 + 1), len(sliders)))  # tropism
    return html.Div([panel2, panel3])


def stem_layout(data, type_names):
    """stem tab layout: with stem subTypes as sub tabs"""
    panelChildren = []
    for i in range(0, type_names["number_stemtypes"]):
        ctab = dcc.Tab(
            label=type_names[f"stem tab-{i+1}"],
            value=f"tab-{i+1}",
            className="tab",
            selected_className="tabSelected",
        )
        panelChildren.append(ctab)
    tab = dcc.Tabs(
        id="stem-tabs",
        value="tab-1",
        children=panelChildren,
    )
    content = html.Div(id="stem-tabs-content")
    stored_values = data.get("tab-1", STEM_SLIDER_INITIALS)
    return [tab, content]


# render_stem_tab - Display sliders for the selected tab, using stored values
@app.callback(
    Output("stem-tabs-content", "children"),
    Input("stem-tabs", "value"),
    State("stem-store", "data"),
)
def render_stem_tab(selected_tab, data):
    stored_values = data.get(selected_tab, STEM_SLIDER_INITIALS)
    # print("render_stem_tab()", selected_tab, stored_values)
    return generate_stem_sliders(stored_values, int(selected_tab[-1]))


# Update root-store when any slider changes
@app.callback(
    Output("stem-store", "data", allow_duplicate=True),  #
    Input({"type": "stem-dynamic-slider", "index": dash.ALL}, "value"),
    State("stem-tabs", "value"),
    State("stem-store", "data"),
    prevent_initial_call=True,
)
def update_stem_store(slider_values, current_tab, store_data):
    print("update_stem_store()", current_tab, slider_values)
    successors = store_data[current_tab][-1]
    store_data[current_tab] = slider_values
    store_data[current_tab].append(successors)
    return store_data


#
# Parameters Panel - Leaf
#
def generate_leaf_sliders(data):  # Generate sliders for leaf tabs from stored values

    if data["leaf"] is None:
        return [
            html.Div(className="spacer"),
            html.H6("There is no leaf defined in your plant data set"),
        ]

    leaf_values = data["leaf"]
    # print(len(leaf_values), len(leaf_parameter_sliders.keys()))
    sliders = []
    sliders.append(html.H6("Leaf shape"))
    sliders.append(
        dcc.Dropdown(
            id={"type": "leaf-dynamic-slider", "index": 0},  # little white lie
            options=["Defined", "Long", "Round", "Maple", "Flower"],
            value="Defined",  # leaf_values[0] ???
            clearable=False,
            className="dropdown",
            style={"fontSize": "12px", "padding-top": "5px"},
        )
    )
    for i, key in enumerate(leaf_parameter_sliders.keys()):
        min_ = leaf_parameter_sliders[key][0]
        max_ = leaf_parameter_sliders[key][1]
        step_ = leaf_parameter_sliders[key][2]
        sliders.append(html.H6(key))
        sliders.append(
            dcc.Slider(
                id={"type": "leaf-dynamic-slider", "index": i + 1},
                min=min_,
                max=max_,
                value=leaf_values[i + 1],
                marks={min_: str(min_), max_: str(max_)},
                tooltip={"always_visible": False},
                step=step_,
            )
        )
    sliders.append(
        dcc.Dropdown(
            id={"type": "leaf-dynamic-slider", "index": i + 1},  # little white lie
            options=list(conversions.tropism_names.keys()),
            value=leaf_values[i + 2],
            clearable=False,
            className="customDropdown",
            style={"fontSize": "12px", "padding-top": "5px"},
        )
    )
    panel1 = conversions.into_panel(sliders, range(0, 2 + 2 * 2))
    panel2 = conversions.into_panel(sliders, range(2 + 2 * 2, len(sliders) - (2 * 2 + 1)))  # all but tropism
    panel3 = conversions.into_panel(sliders, range(len(sliders) - (2 * 2 + 1), len(sliders)))  # tropism
    return html.Div([panel1, panel2, panel3])


# @app.callback(# leaf-dropdown
#     Output({'type': 'leaf-dynamic-slider', 'index': dash.ALL}, 'value'),
#     Input({'type': 'leaf-dynamic-slider', 'index': 0}, 'value'),  # leaf-dropdown
#     State({'type': 'leaf-dynamic-slider', 'index': dash.ALL}, 'value'),
# )
# def update_leaf_shape(shape, slider_values):
#     print("update_leaf_shape()", slider_values)
#     set_leaf_sliders(slider_values)
#     return slider_values


@app.callback(  # Update leaf-store when any slider changes
    Output("leaf-store", "data", allow_duplicate=True),  #
    Input({"type": "leaf-dynamic-slider", "index": dash.ALL}, "value"),
    State("leaf-store", "data"),
    prevent_initial_call=True,
)
def update_leaf_store(slider_values, store_data):
    # print("update_leaf_store()", slider_values)
    if len(slider_values) > 1:  # called empty
        store_data["leaf"] = slider_values
    return store_data


#
# 3. RIGHT - Results Panel
#
def create_geometry_buttons():
    vtk_btn = small_button("vtp", "vkt-download-button", "Download plant architecture as Paraview VTP file")
    rsml_btn = small_button("rsml", "rsml-download-button", "Download plant architecture root system markup language file (RSML)")
    return html.Div([vtk_btn, rsml_btn])


def attach_buttons(graph, buttons):
    return html.Div(
        [graph, buttons],
        style={
            "display": "flex",
            "flex-direction": "column",  # stack vertically
            "align-items": "flex-start",
        },
    )


def render_result_content(tab, run_id, typename_data, settings_data):
    """Build tab content from cached simulation data identified by run id."""

    if tab == "About":
        return html.Div(
            html.Div(dcc.Markdown(ABOUT_TEXT, style={"whiteSpace": "pre-line"}), className="aboutScroll"),
            className="aboutContainer",
        )

    cached = get_cached_simulation_run(run_id)
    if not cached:
        print("no cached data for run", run_id)
        return html.Div([html.H6("press the create button")])

    vtk_data = cached["vtk_data"]
    result_data = cached["result_data"]
    vtk_render_data = cached["vtk_render_data"]

    print("render_result_tab()", tab, "vtk data size:", asizeof.asizeof(vtk_data) / 1e6, "MB", "result data size:", asizeof.asizeof(result_data) / 1e6, "MB")

    if tab == "VTK3D":
        buttons = create_geometry_buttons()
        graph = plots.vtk3D_plot(vtk_render_data, "Type", settings_data["token"], settings_data["reset"], buttons)
        return graph
    elif tab == "VTK3DAge":
        buttons = create_geometry_buttons()
        graph = plots.vtk3D_plot(vtk_render_data, "Age", settings_data["token"], settings_data["reset"], buttons)
        return graph
    elif tab == "Profile1D":
        button = small_button("xls", "xls-profile-download-button", "Download 1D profile data as XLS")
        graph = plots.profile_plot(result_data)
        return attach_buttons(graph, button)
    elif tab == "Dynamics":
        button = small_button("xls", "xls-dynamics-download-button", "Download 1D dynamic data as XLS")
        graph = plots.dynamics_plot(result_data, typename_data)
        return attach_buttons(graph, button)


@app.callback(
    Output("result-tabs-content", "children", allow_duplicate=True),
    Input("result-tabs", "value"),
    State("run-id-store", "data"),
    State("typename-store", "data"),
    State("settings-store", "data"),
    prevent_initial_call=True,
)
def render_result_tab(tab, run_data, typename_data, settings_data):
    run_id = (run_data or {}).get("run_id")
    return render_result_content(tab, run_id, typename_data, settings_data)


@app.callback(
    Output("download-vtk", "data"),
    Output("loading-spinner-output2", "children"),
    Input("vkt-download-button", "n_clicks"),
    State("time-slider", "value"),
    State("plant-dropdown", "value"),
    State("seed-store", "data"),
    State("root-store", "data"),
    State("stem-store", "data"),
    State("leaf-store", "data"),
    State("settings-store", "data"),
    State("xml-store", "data"),
    prevent_initial_call=True,
)
def download_vtp(n_clicks, time_slider, plant_value, seed_data, root_data, stem_data, leaf_data, settings_data, xml_data):
    # print("download_vtp")
    if n_clicks is None:
        triggered = ctx.triggered_id
        return dash.no_update, dash.no_update
    plant, _, _ = simulate_plant.get_plant(plant_value, seed_data, root_data, stem_data, leaf_data, xml_data)
    N = time_slider  # makes dt = 1
    t_ = np.linspace(0.0, time_slider, N + 1)
    plant.setSeed(settings_data["random_seed"])
    plant.initialize()
    for _, dt in enumerate(np.diff(t_)):
        plant.simulate(dt)
    filename = "cplantbox_" + param_names[int(plant_value)][0] + ".vtp"
    vtp_string = plant.write(filename, False)
    return dict(content=vtp_string, filename=filename, type="application/vtp"), html.H6("")


@app.callback(
    Output("download-rsml", "data"),
    Output("loading-spinner-output3", "children"),
    Input("rsml-download-button", "n_clicks"),
    State("time-slider", "value"),
    State("plant-dropdown", "value"),
    State("seed-store", "data"),
    State("root-store", "data"),
    State("stem-store", "data"),
    State("leaf-store", "data"),
    State("settings-store", "data"),
    State("xml-store", "data"),
    prevent_initial_call=True,
)
def download_rsml(n_clicks, time_slider, plant_value, seed_data, root_data, stem_data, leaf_data, settings_data, xml_data):
    # print("download_rsml()")
    if n_clicks is None:
        triggered = ctx.triggered_id
        return dash.no_update, dash.no_update
    plant, _, _ = simulate_plant.get_plant(plant_value, seed_data, root_data, stem_data, leaf_data, xml_data)
    N = time_slider  # makes dt = 1
    t_ = np.linspace(0.0, time_slider, N + 1)
    plant.setSeed(settings_data["random_seed"])
    plant.initialize()
    for _, dt in enumerate(np.diff(t_)):
        plant.simulate(dt)
    filename = "cplantbox_" + param_names[int(plant_value)][0] + ".rsml"
    vtp_string = plant.write(filename, False)
    return dict(content=vtp_string, filename=filename, type="application/rsml"), html.H6("")


@app.callback(
    Output("download-profiles-xls", "data"),
    Input("xls-profile-download-button", "n_clicks"),
    State("run-id-store", "data"),
)
def download_profiles_xls(n_clicks, run_data):
    print("download_profiles_xls()")
    if n_clicks is None:
        return dash.no_update
    run_id = (run_data or {}).get("run_id")
    cached = get_cached_simulation_run(run_id)
    if not cached:
        return dash.no_update
    data = cached["result_data"]
    data_xls = plots.profile_to_excel(data)
    return data_xls


@app.callback(
    Output("download-dynamics-xls", "data"),
    Input("xls-dynamics-download-button", "n_clicks"),
    State("run-id-store", "data"),
    State("typename-store", "data"),
)
def download_dynamics_xls(n_clicks, run_data, typename_data):
    print("download_dynamics_xls()")
    if n_clicks is None:
        return dash.no_update
    run_id = (run_data or {}).get("run_id")
    cached = get_cached_simulation_run(run_id)
    if not cached:
        return dash.no_update
    data = cached["result_data"]
    data_xls = plots.dynamics_to_excel(data, typename_data)
    return data_xls


@app.callback(
    Output("reference-panel", "children"),
    Input("plant-dropdown", "value"),
    State("xml-store", "data"),
)
def update_reference_panel(plant_value, xml_data):
    """Show APA reference from the selected dataset's XML file, or hide if absent."""
    fname = conversions.get_parameter_names()[int(plant_value)][1]

    if fname == "xml-store":
        xml_str = (xml_data or {}).get("xml", "")
    else:
        xml_path = os.path.join(BASE_DIR, "params", fname)
        try:
            with open(xml_path, "r", encoding="utf-8") as f:
                xml_str = f.read()
        except OSError:
            return []

    citations = bibtexxml_apa.bibtexxml_apa_from_string(xml_str)
    if not citations:
        return []

    md_text = "\n\n".join(citations)
    items = [html.H6("Reference"), dcc.Markdown(md_text, style={"fontSize": "11px", "marginBottom": "4px"})]
    return conversions.into_panel(items)


if __name__ == "__main__":

    # Timer(1, open_browser).start()
    app.run(debug=True)
