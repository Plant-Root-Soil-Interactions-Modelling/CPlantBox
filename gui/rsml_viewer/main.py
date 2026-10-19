"""RSML Viewer Webapp (using Python Dash), D. Leitner 2026"""

import base64
import hashlib
import os
import sys
import tempfile
from collections import OrderedDict
from threading import Lock

import dash
import dash_bootstrap_components as dbc
import dash_vtk
import numpy as np
import plotly.graph_objects as go
import vtk
from dash import Input, Output, State, ctx, dcc, html, no_update

# Add the viewer folder so ViewerDataModel and viewer_conductivities can be imported
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "viewer"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cplantbox"))

import viewer_conductivities  # noqa: E402
from viewer_data import ViewerDataModel  # noqa: E402
from vtk_conversions import (  # noqa: E402
    apply_tube_filter,
    decode_array,
    generate_colorbar_image,
    vtk_polydata_to_dashvtk_dict,
)

import plantbox as pb  # noqa: E402
import plantbox.visualisation.vtk_plot as vp  # noqa: E402

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
md_path = os.path.join(BASE_DIR, "assets", "readme.md")
with open(md_path, "r", encoding="utf-8") as f:
    ABOUT_TEXT = f.read()

MAX_RSML_SIZE = 10 * 1024 * 1024  # 10 MB

HYDRAULIC_SCENARIOS = [
    "Constant scenario 1",
    "Constant scenario 2",
    "Dynamic scenario 1",
    "Dynamic scenario 2",
    "Wine",
]

VTK_VIEW_OPTIONS = [
    ("Type", "subType"),
    ("Segment length", "length"),
    ("Creation time", "creationTime"),
    ("SUF", "SUF"),
]

COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]

MAX_GEOMETRY_CACHE_ITEMS = 8
GEOMETRY_CACHE = OrderedDict()
GEOMETRY_CACHE_LOCK = Lock()
PAGE_BACKGROUND_RGB = [1.0, 1.0, 1.0]

#
# INITIALISE
#
app = dash.Dash(
    __name__,
    suppress_callback_exceptions=True,
    external_stylesheets=[dbc.themes.MINTY],
)
app.title = "RSML Viewer"


def into_panel(content):
    """Wraps html content in a styled panel div."""
    return html.Div(content, className="panel")


def small_button(label, id_, tool_tip="", url="cloud-download.svg"):
    """Creates a small button with an icon from the assets folder"""
    return dcc.Button(
        [html.Img(src=app.get_asset_url(url), className="buttonIcon"), label],
        id=id_,
        title=tool_tip,
        className="smallButton",
    )


#
# LAYOUT
#
app.layout = dbc.Container(
    [
        dcc.Store(id="rsml-store", data={}),
        dcc.Store(
            id="state-store",
            data={"shoot_added": False, "ct_added": False, "max_ct": None, "filename": None},
        ),
        dcc.Download(id="download-rsml"),
        dcc.Download(id="download-vtp"),
        dbc.Row(
            [
                #
                # Panel 1 – Controls
                #
                dbc.Col(
                    [
                        html.H5("RSML Viewer"),
                        into_panel(
                            [
                                html.H6("Upload RSML file"),
                                dcc.Upload(
                                    id="rsml-upload",
                                    accept=".rsml,.RSML,.xml",
                                    multiple=False,
                                    children=small_button(
                                        "rsml",
                                        "rsml-upload-btn-hidden",
                                        "Upload an RSML or XML root system file",
                                        url="cloud-upload.svg",
                                    ),
                                ),
                                html.Div(
                                    id="upload-feedback",
                                    style={"fontSize": "10px", "color": "#555", "marginTop": "5px"},
                                ),
                            ]
                        ),
                        html.Div(className="spacer"),
                        into_panel(
                            [
                                html.H6("Edit"),
                                dcc.Checklist(
                                    id="shoot-checkbox",
                                    options=[{"label": "Add artificial shoot", "value": "agree"}],
                                    value=[],
                                    className="checkbox",
                                    style={"fontSize": "12px"},
                                ),
                                html.Div(className="smallSpacer"),
                                html.H6("Add creation times [day]"),
                                dbc.InputGroup(
                                    [
                                        dbc.Input(
                                            id="max-ct-input",
                                            type="number",
                                            placeholder="Final age",
                                            size="sm",
                                            min=0,
                                        ),
                                        dbc.Button(
                                            "Apply",
                                            id="apply-ct-button",
                                            size="sm",
                                            className="button",
                                        ),
                                    ],
                                    size="sm",
                                ),
                                html.Div(
                                    id="ct-feedback",
                                    style={"fontSize": "10px", "color": "#555", "marginTop": "4px"},
                                ),
                            ]
                        ),
                        html.Div(className="spacer"),
                        into_panel(
                            [
                                html.H6("Download"),
                                html.Div(
                                    [
                                        small_button(
                                            "rsml",
                                            "download-rsml-button",
                                            "Download as RSML file",
                                        ),
                                        small_button(
                                            "vtp",
                                            "download-vtp-button",
                                            "Download as Paraview VTP file",
                                        ),
                                    ],
                                    style={"display": "flex", "flexWrap": "wrap", "gap": "4px"},
                                ),
                            ]
                        ),
                        html.Div(className="largeSpacer"),
                        dcc.Loading(
                            id="loading-spinner",
                            type="circle",
                            children=html.Div(id="loading-output"),
                        ),
                    ],
                    width=2,
                    style={"minWidth": "200px"},
                ),
                #
                # Panel 2 – Results
                #
                dbc.Col(
                    [
                        html.H5("Results"),
                        dcc.Tabs(
                            id="result-tabs",
                            value="info",
                            children=[
                                dcc.Tab(label="Information", value="info", className="tab", selected_className="tabSelected"),
                                dcc.Tab(label="3D View", value="vtk", className="tab", selected_className="tabSelected"),
                                dcc.Tab(label="Depth Profile", value="profile", className="tab", selected_className="tabSelected"),
                                dcc.Tab(label="Development", value="development", className="tab", selected_className="tabSelected"),
                                dcc.Tab(label="Hydraulic Properties", value="suf", className="tab", selected_className="tabSelected"),
                                dcc.Tab(label="Hydraulic Development", value="krs", className="tab", selected_className="tabSelected"),
                                dcc.Tab(label="About", value="about", className="tab", selected_className="tabSelected"),
                            ],
                            className="tabs",
                        ),
                        dcc.Loading(
                            type="circle",
                            children=html.Div(id="tabs-content"),
                        ),
                    ],
                    width=9,
                    style={"minWidth": "320px"},
                ),
            ]
        ),
        html.Div(
            [
                html.A(
                    html.Img(src="/assets/cplantbox.png", className="logo"),
                    href="https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox",
                    target="_blank",
                ),
                html.A(
                    html.Img(src="/assets/fzj.png", className="logo"),
                    href="https://www.fz-juelich.de/de",
                    target="_blank",
                ),
            ],
            className="logoContainer",
        ),
    ],
    fluid=True,
)


# ---------------------------------------------------------------------------
# Helper: recreate ViewerDataModel from stored XML + apply modifications
# ---------------------------------------------------------------------------


def load_data(rsml_store, state_store):
    """Load a ViewerDataModel from stored XML content and apply any edits."""
    rsml_xml = rsml_store.get("xml", None)
    if rsml_xml is None:
        return None
    with tempfile.NamedTemporaryFile(suffix=".rsml", mode="w", delete=False, encoding="utf-8") as f:
        f.write(rsml_xml)
        fname = f.name
    try:
        data = ViewerDataModel()
        data.open_rsml(fname)
        if state_store.get("shoot_added", False) and len(data.base_nodes) > 1:
            data.add_artificial_shoot()
        if state_store.get("ct_added", False) and state_store.get("max_ct") is not None:
            if not data.tagnames[1]:  # only interpolate when file has no creation times
                data.max_ct = float(state_store["max_ct"])
                data.add_creation_times()
    finally:
        os.unlink(fname)
    return data


def _make_geometry_cache_key(rsml_store, state_store, view_type):
    """Create a stable cache key from RSML content + edit state + view type."""
    rsml_xml = rsml_store.get("xml")
    if not rsml_xml:
        return None
    xml_hash = hashlib.md5(rsml_xml.encode("utf-8")).hexdigest()
    shoot_added = int(bool(state_store.get("shoot_added", False)))
    ct_added = int(bool(state_store.get("ct_added", False)))
    max_ct = state_store.get("max_ct")
    return f"{xml_hash}|s:{shoot_added}|ct:{ct_added}|m:{max_ct}|v:{view_type}"


def _get_cached_geometry_payload(cache_key):
    """Fetch a cached VTK payload and refresh recency ordering."""
    if cache_key is None:
        return None
    with GEOMETRY_CACHE_LOCK:
        payload = GEOMETRY_CACHE.get(cache_key)
        if payload is not None:
            GEOMETRY_CACHE.move_to_end(cache_key)
        return payload


def _set_cached_geometry_payload(cache_key, payload):
    """Store payload in cache and keep only the newest N entries."""
    if cache_key is None:
        return
    with GEOMETRY_CACHE_LOCK:
        GEOMETRY_CACHE[cache_key] = payload
        GEOMETRY_CACHE.move_to_end(cache_key)
        while len(GEOMETRY_CACHE) > MAX_GEOMETRY_CACHE_ITEMS:
            GEOMETRY_CACHE.popitem(last=False)


# ---------------------------------------------------------------------------
# 1. Upload
# ---------------------------------------------------------------------------


@app.callback(
    Output("rsml-store", "data"),
    Output("state-store", "data"),
    Output("upload-feedback", "children"),
    Input("rsml-upload", "contents"),
    State("rsml-upload", "filename"),
    State("state-store", "data"),
    prevent_initial_call=True,
)
def handle_upload(contents, filename, state_data):
    if contents is None:
        return no_update, no_update, ""
    _content_type, content_string = contents.split(",")
    decoded = base64.b64decode(content_string)
    if len(decoded) > MAX_RSML_SIZE:
        return no_update, no_update, f"File too large (max {MAX_RSML_SIZE // (1024 * 1024)} MB)"
    xml_string = decoded.decode("utf-8")
    new_state = {"shoot_added": False, "ct_added": False, "max_ct": None, "filename": filename}
    return {"xml": xml_string}, new_state, f"Loaded: {filename}"


# ---------------------------------------------------------------------------
# 2. Edit – add artificial shoot
# ---------------------------------------------------------------------------


@app.callback(
    Output("state-store", "data", allow_duplicate=True),
    Input("shoot-checkbox", "value"),
    State("state-store", "data"),
    prevent_initial_call=True,
)
def toggle_shoot(checkbox_value, state_data):
    state_data["shoot_added"] = "agree" in checkbox_value
    return state_data


# ---------------------------------------------------------------------------
# 3. Edit – add creation times
# ---------------------------------------------------------------------------


@app.callback(
    Output("state-store", "data", allow_duplicate=True),
    Output("ct-feedback", "children"),
    Input("apply-ct-button", "n_clicks"),
    State("max-ct-input", "value"),
    State("state-store", "data"),
    State("rsml-store", "data"),
    prevent_initial_call=True,
)
def apply_creation_times(n_clicks, max_ct, state_data, rsml_store):
    if max_ct is None:
        return state_data, "Enter a final age first"
    if rsml_store.get("xml") is None:
        return state_data, "Upload a file first"
    # Quick check: does the file already have creation times?
    data = load_data(rsml_store, {"shoot_added": state_data.get("shoot_added", False), "ct_added": False, "max_ct": None})
    if data is not None and data.tagnames[1]:
        return state_data, f"File already has creation times (tag '{data.tagnames[1]}')"
    state_data["ct_added"] = True
    state_data["max_ct"] = max_ct
    return state_data, f"Creation times set (max {max_ct} days)"


# ---------------------------------------------------------------------------
# 4. Render tab layout
# ---------------------------------------------------------------------------


@app.callback(
    Output("tabs-content", "children"),
    Input("result-tabs", "value"),
    Input("state-store", "data"),
    State("rsml-store", "data"),
)
def render_tab(tab, state_data, rsml_store):
    if tab == "about":
        return html.Div(dcc.Markdown(ABOUT_TEXT), className="aboutContainer")

    data = load_data(rsml_store, state_data)

    if tab == "info":
        return _render_info_tab(data, state_data)
    elif tab == "vtk":
        return _render_vtk_tab(data)
    elif tab == "profile":
        return _render_plot_tab(data, "profile")
    elif tab == "development":
        return _render_plot_tab(data, "development")
    elif tab == "suf":
        return _render_plot_tab(data, "suf")
    elif tab == "krs":
        return _render_plot_tab(data, "krs")


def _no_data_msg(msg="Upload an RSML file to view results."):
    return html.Div(msg, style={"padding": "20px", "color": "#666"})


def _render_info_tab(data, state_data):
    if data is None:
        return _no_data_msg("Upload an RSML file to view information.")

    # Node count
    c = sum(1 + len(pl) for pl in data.polylines)
    metadata = data.metadata
    min_str = str(data.analyser.getMinBounds())
    max_str = str(data.analyser.getMaxBounds())
    filename = state_data.get("filename", "N/A")

    general_rows = [
        ("Software", metadata.software or "—"),
        ("Filename", filename),
        ("Number of plants (base nodes)", str(len(data.base_nodes))),
        ("Number of base roots (segments)", str(len(data.base_segs))),
        ("Number of roots", str(len(data.polylines))),
        ("Number of nodes", str(c)),
        ("Bounding box", f"{min_str} – {max_str} (cm)"),
        ("Unit (length scale)", metadata.unit or "—"),
        ("Resolution", f"{metadata.resolution} (dots per {metadata.unit or '?'})"),
    ]

    prop_rows = []
    for k, v_ in data.properties.items():
        v = np.array(v_, dtype=float)
        unit = metadata.properties[k].unit if (metadata.properties and k in metadata.properties) else ""
        prop_rows.append((k, f"[{np.nanmin(v):.4g}, {np.nanmax(v):.4g}] {unit}".strip()))

    fun_rows = []
    for k, v_ in data.functions.items():
        flat = [x for sub in v_ for x in sub]
        v = np.array(flat, dtype=float)
        unit = metadata.properties[k].unit if (metadata.properties and k in metadata.properties) else ""
        fun_rows.append((k, f"[{np.nanmin(v):.4g}, {np.nanmax(v):.4g}] {unit}".strip()))

    tagnames = data.tagnames
    using_rows = [
        (
            "Radius",
            f"from tag '{tagnames[0]}' within [{np.min(data.radii):.4g}, {np.max(data.radii):.4g}] cm" if tagnames[0] else "not found (set to 0.1 cm)",
        ),
        (
            "Creation time",
            f"from tag '{tagnames[1]}' within [{np.min(data.cts):.4g}, {np.max(data.cts):.4g}] days" if tagnames[1] else "not found",
        ),
        (
            "Types",
            f"from tag '{tagnames[2]}' within [{np.min(data.types):.4g}, {np.max(data.types):.4g}]" if tagnames[2] else "not found",
        ),
    ]

    def _table(rows, title):
        return html.Div(
            [
                html.H6(title),
                dbc.Table(
                    html.Tbody(
                        [
                            html.Tr(
                                [
                                    html.Td(r[0], style={"fontWeight": "bold", "width": "45%", "fontSize": "12px"}),
                                    html.Td(r[1], style={"fontSize": "12px"}),
                                ]
                            )
                            for r in rows
                        ]
                    ),
                    bordered=False,
                    size="sm",
                    style={"marginBottom": "0"},
                ),
            ],
            className="panel",
        )

    sections = [_table(general_rows, "General")]
    if prop_rows:
        sections.append(_table(prop_rows, "Properties (values per root)"))
    if fun_rows:
        sections.append(_table(fun_rows, "Functions (values per node)"))
    sections.append(_table(using_rows, "Using"))

    return html.Div(sections, style={"padding": "10px"})


def _render_plot_tab(data, tab_id):
    """Return the layout for a plot tab (profile / development / suf / krs)."""
    CONFIGS = {
        "profile": {
            "label": "Plot type",
            "options": ["Length", "Surface", "Volume"],
            "dropdown_id": "profile-dropdown",
            "graph_id": "profile-graph",
            "empty_msg": "Upload an RSML file to view depth profiles.",
            "graph_style": {"height": "70vh", "width": "50%", "margin": "0 auto"},
        },
        "development": {
            "label": "Plot type",
            "options": ["Length", "Surface", "Volume"],
            "dropdown_id": "development-dropdown",
            "graph_id": "development-graph",
            "empty_msg": "Upload an RSML file to view development.",
            "graph_style": {"height": "70vh"},
        },
        "suf": {
            "label": "Hydraulic scenario",
            "options": HYDRAULIC_SCENARIOS,
            "dropdown_id": "suf-dropdown",
            "graph_id": "suf-graph",
            "empty_msg": "Upload an RSML file to view hydraulic properties.",
            "graph_style": {"height": "70vh", "width": "50%", "margin": "0 auto"},
        },
        "krs": {
            "label": "Hydraulic scenario",
            "options": HYDRAULIC_SCENARIOS,
            "dropdown_id": "krs-dropdown",
            "graph_id": "krs-graph",
            "empty_msg": "Upload an RSML file to view hydraulic development.",
            "graph_style": {"height": "70vh"},
        },
    }
    cfg = CONFIGS[tab_id]
    if data is None:
        return _no_data_msg(cfg["empty_msg"])

    return html.Div(
        [
            html.Div(
                [
                    html.H6(cfg["label"]),
                    dcc.Dropdown(
                        id=cfg["dropdown_id"],
                        options=[{"label": v, "value": str(i)} for i, v in enumerate(cfg["options"])],
                        value="0",
                        clearable=False,
                        className="dropdown",
                        style={"fontSize": "12px"},
                    ),
                ],
                style={"maxWidth": "320px", "padding": "10px 10px 0 10px"},
            ),
            dcc.Graph(id=cfg["graph_id"], style=cfg["graph_style"]),
        ]
    )


def _render_vtk_tab(data):
    if data is None:
        return _no_data_msg("Upload an RSML file to view 3D geometry.")

    return html.Div(
        [
            html.Div(
                [
                    html.H6("View type"),
                    dcc.Dropdown(
                        id="vtk-view-dropdown",
                        options=[{"label": label, "value": value} for label, value in VTK_VIEW_OPTIONS],
                        value="subType",
                        clearable=False,
                        className="dropdown",
                        style={"fontSize": "12px"},
                    ),
                ],
                style={"maxWidth": "320px", "padding": "10px 10px 0 10px"},
            ),
            html.Div(
                dcc.Loading(type="circle", children=html.Div(id="vtk-content")),
                style={"marginTop": "5px"},
            ),
        ]
    )


def _vtk_scalar_data(polydata, scalar_name):
    """Fetch scalar values from vtkPolyData and report whether values are point- or cell-based."""
    arr = polydata.GetCellData().GetArray(scalar_name)
    if arr is not None:
        values = np.array([arr.GetTuple1(i) for i in range(arr.GetNumberOfTuples())], dtype=np.float32)
        return values, "cell"

    arr = polydata.GetPointData().GetArray(scalar_name)
    if arr is not None:
        values = np.array([arr.GetTuple1(i) for i in range(arr.GetNumberOfTuples())], dtype=np.float32)
        return values, "point"

    # Fallback: render with constant color if scalar array is missing.
    n_cells = polydata.GetNumberOfCells()
    values = np.zeros((n_cells,), dtype=np.float32)
    return values, "cell"


def _build_vtk_payload(data, view_type):
    """Build mesh/scalar payload for the VTK tab (expensive path)."""
    analyser = pb.SegmentAnalyser(data.analyser)

    if view_type == "length":
        seg_lengths = np.array([analyser.getSegmentLength(i) for i in range(len(analyser.segments))], dtype=np.float64)
        analyser.addData("length", seg_lengths)
    elif view_type == "SUF":
        _init_hydraulic_scenario(data.hydraulic_params, 0)
        suf = data.hydraulic_model.get_suf(data.max_ct)
        analyser.addData("SUF", suf)

    polydata = vp.segs_to_polydata(analyser, zoom_factor=1.0, param_names=["radius", view_type])
    tube_polydata = apply_tube_filter(polydata)
    vtk_data = vtk_polydata_to_dashvtk_dict(tube_polydata)
    scalar_values, color_mode = _vtk_scalar_data(tube_polydata, view_type)

    vmin = float(np.min(scalar_values)) if scalar_values.size > 0 else 0.0
    vmax = float(np.max(scalar_values)) if scalar_values.size > 0 else 1.0
    if np.isclose(vmin, vmax):
        vmax = vmin + 1.0

    title = {
        "subType": "Root Type",
        "length": "Segment Length [cm]",
        "creationTime": "Creation Time [day]",
        "SUF": "Surface Uptake Fraction [1]",
    }.get(view_type, view_type)

    return {
        "points": decode_array(vtk_data["points"]).flatten().tolist(),
        "polys": decode_array(vtk_data["polys"]).tolist(),
        "scalar_values": scalar_values.tolist(),
        "color_mode": color_mode,
        "vmin": vmin,
        "vmax": vmax,
        "title": title,
        "discrete": (view_type == "subType"),
    }


def _render_vtk_content(payload):
    """Render cached VTK payload as dash_vtk components."""
    scalar_values = payload["scalar_values"]
    color_mode = payload["color_mode"]
    vmin = payload["vmin"]
    vmax = payload["vmax"]

    scalar_container = dash_vtk.CellData if color_mode == "cell" else dash_vtk.PointData
    scalar_bar = dcc.Graph(
        figure=generate_colorbar_image(vmin=vmin, vmax=vmax, colormap="Jet", height=40, width=220, discrete=payload["discrete"]),
        style={"width": "220px", "height": "70px", "marginTop": "8px"},
        config={"displayModeBar": False},
    )

    view = dash_vtk.View(
        children=[
            dash_vtk.GeometryRepresentation(
                mapper={"colorByArrayName": "Colors", "colorMode": color_mode},
                colorDataRange=[vmin, vmax],
                children=[
                    dash_vtk.PolyData(
                        points=payload["points"],
                        polys=payload["polys"],
                        children=[
                            scalar_container(
                                [
                                    dash_vtk.DataArray(
                                        registration="setScalars",
                                        name="Colors",
                                        numberOfComponents=1,
                                        values=scalar_values,
                                    )
                                ]
                            )
                        ],
                    )
                ],
            )
        ],
        background=PAGE_BACKGROUND_RGB,
    )

    return html.Div(
        [
            html.Div(
                view,
                style={"width": "100%", "height": "620px", "backgroundColor": "var(--bs-body-bg, #ffffff)"},
            ),
            html.Div([html.H6(payload["title"], style={"marginTop": "10px"}), scalar_bar], style={"display": "flex", "justifyContent": "flex-end"}),
        ],
        style={"width": "100%", "display": "flex", "flexDirection": "column", "backgroundColor": "var(--bs-body-bg, #ffffff)"},
    )


# ---------------------------------------------------------------------------
# 5. Plot callbacks
# ---------------------------------------------------------------------------


@app.callback(
    Output("vtk-content", "children"),
    Input("vtk-view-dropdown", "value"),
    State("rsml-store", "data"),
    State("state-store", "data"),
)
def update_vtk_view(view_type, rsml_store, state_data):
    if rsml_store.get("xml") is None:
        return _no_data_msg("Upload an RSML file to view 3D geometry.")

    cache_key = _make_geometry_cache_key(rsml_store, state_data, view_type)
    payload = _get_cached_geometry_payload(cache_key)
    if payload is None:
        data = load_data(rsml_store, state_data)
        if data is None:
            return _no_data_msg("Upload an RSML file to view 3D geometry.")
        payload = _build_vtk_payload(data, view_type)
        _set_cached_geometry_payload(cache_key, payload)

    return _render_vtk_content(payload)


@app.callback(
    Output("profile-graph", "figure"),
    Input("profile-dropdown", "value"),
    State("rsml-store", "data"),
    State("state-store", "data"),
)
def update_profile_graph(j_str, rsml_store, state_data):
    data = load_data(rsml_store, state_data)
    if data is None:
        return go.Figure()
    return _plot_depth_profile(data.analyser, int(j_str))


@app.callback(
    Output("development-graph", "figure"),
    Input("development-dropdown", "value"),
    State("rsml-store", "data"),
    State("state-store", "data"),
)
def update_development_graph(j_str, rsml_store, state_data):
    data = load_data(rsml_store, state_data)
    if data is None:
        return go.Figure()
    return _plot_rootsystem_development(data.analyser, int(j_str))


@app.callback(
    Output("suf-graph", "figure"),
    Input("suf-dropdown", "value"),
    State("rsml-store", "data"),
    State("state-store", "data"),
)
def update_suf_graph(j_str, rsml_store, state_data):
    data = load_data(rsml_store, state_data)
    if data is None:
        return go.Figure()
    return _plot_suf(data, int(j_str))


@app.callback(
    Output("krs-graph", "figure"),
    Input("krs-dropdown", "value"),
    State("rsml-store", "data"),
    State("state-store", "data"),
)
def update_krs_graph(j_str, rsml_store, state_data):
    data = load_data(rsml_store, state_data)
    if data is None:
        return go.Figure()
    return _plot_krs(data, int(j_str))


# ---------------------------------------------------------------------------
# 6. Download callbacks
# ---------------------------------------------------------------------------


@app.callback(
    Output("download-rsml", "data"),
    Output("loading-output", "children"),
    Input("download-rsml-button", "n_clicks"),
    State("rsml-store", "data"),
    State("state-store", "data"),
    prevent_initial_call=True,
)
def download_rsml_file(n_clicks, rsml_store, state_data):
    import plantbox.visualisation.vtk_plot as vp
    import plantbox.visualisation.vtk_tools as vt

    data = load_data(rsml_store, state_data)
    if data is None:
        return no_update, html.H6("")
    with tempfile.NamedTemporaryFile(suffix=".rsml", delete=False) as f:
        fname = f.name
    try:
        pd = vp.segs_to_polydata(data.analyser, zoom_factor=1.0, param_names=["subType", "radius", "creationTime"])
        vt.write_rsml(fname, pd, 0, None, data.base_nodes)
        with open(fname, "r", encoding="utf-8") as f:
            content = f.read()
    finally:
        os.unlink(fname)
    return dict(content=content, filename="rsml_viewer_output.rsml", type="application/xml"), html.H6("")


@app.callback(
    Output("download-vtp", "data"),
    Output("loading-output", "children", allow_duplicate=True),
    Input("download-vtp-button", "n_clicks"),
    State("rsml-store", "data"),
    State("state-store", "data"),
    prevent_initial_call=True,
)
def download_vtp_file(n_clicks, rsml_store, state_data):
    data = load_data(rsml_store, state_data)
    if data is None:
        return no_update, html.H6("")
    with tempfile.NamedTemporaryFile(suffix=".vtp", delete=False) as f:
        fname = f.name
    try:
        data.analyser.write(fname)
        with open(fname, "rb") as f:
            content = base64.b64encode(f.read()).decode("ascii")
    finally:
        os.unlink(fname)
    return (
        dict(base64=True, content=content, filename="rsml_viewer_output.vtp", type="application/octet-stream"),
        html.H6(""),
    )


# ---------------------------------------------------------------------------
# Plotly plot functions
# ---------------------------------------------------------------------------


def _init_hydraulic_scenario(params, j: int):
    if j == 0:
        viewer_conductivities.init_constant_scenario1(params)
    elif j == 1:
        viewer_conductivities.init_constant_scenario2(params)
    elif j == 2:
        viewer_conductivities.init_dynamic_scenario1(params)
    elif j == 3:
        viewer_conductivities.init_dynamic_scenario2(params)
    elif j == 4:
        viewer_conductivities.init_constant_scenario_wine(params)


def _plot_depth_profile(analyser, j: int) -> go.Figure:
    type_str = ["length", "surface", "volume"]
    unit_str = ["(cm)", "(cm²)", "(cm³)"]
    n = int(np.ceil(-analyser.getMinBounds().z))
    z_ = np.linspace(-0.5, -n + 0.5, n)
    d = analyser.distribution(type_str[j], 0.0, float(-n), int(n), True)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=d, y=z_, mode="lines+markers", name="total", line=dict(color=COLORS[0])))
    max_type = int(np.max(analyser.data["subType"]))
    for i in range(max_type + 1):
        ana = pb.SegmentAnalyser(analyser)
        ana.filter("subType", i)
        if len(ana.segments) > 0:
            d = ana.distribution(type_str[j], 0.0, float(-n), int(n), True)
            fig.add_trace(go.Scatter(x=d, y=z_, mode="lines+markers", name=f"type {i}", line=dict(color=COLORS[(i + 1) % len(COLORS)])))
    fig.update_layout(
        xaxis_title=f"Root system {type_str[j]} per 1 cm layer {unit_str[j]}",
        yaxis_title="Depth (cm)",
        legend=dict(x=0.7, y=0.02),
        margin=dict(l=60, r=20, t=30, b=60),
        template="plotly_white",
    )
    return fig


def _plot_rootsystem_development(analyser, j: int) -> go.Figure:
    type_str = ["length", "surface", "volume"]
    unit_str = ["(cm)", "(cm²)", "(cm³)"]
    radii = analyser.data["radius"]
    if j == 0:
        weights = [analyser.getSegmentLength(i) for i in range(len(analyser.segments))]
    elif j == 1:
        weights = [2 * np.pi * radii[i] * analyser.getSegmentLength(i) for i in range(len(analyser.segments))]
    else:
        weights = [np.pi * radii[i] ** 2 * analyser.getSegmentLength(i) for i in range(len(analyser.segments))]
    cts = np.array(analyser.data["creationTime"])
    fig = go.Figure()
    try:
        l_, t_ = np.histogram(cts, 100, weights=weights)
        fig.add_trace(go.Scatter(x=0.5 * (t_[1:] + t_[:-1]), y=np.cumsum(l_), mode="lines", name="total", line=dict(color=COLORS[0])))
        max_type = int(np.max(analyser.data["subType"]))
        for i in range(max_type + 1):
            ana = pb.SegmentAnalyser(analyser)
            ana.filter("subType", i)
            if len(ana.segments) > 0:
                radii_i = ana.data["radius"]
                if j == 0:
                    w = [ana.getSegmentLength(k) for k in range(len(ana.segments))]
                elif j == 1:
                    w = [2 * np.pi * radii_i[k] * ana.getSegmentLength(k) for k in range(len(ana.segments))]
                else:
                    w = [np.pi * radii_i[k] ** 2 * ana.getSegmentLength(k) for k in range(len(ana.segments))]
                cts_i = np.array(ana.data["creationTime"])
                l_, t_ = np.histogram(cts_i, 100, weights=w)
                fig.add_trace(
                    go.Scatter(x=0.5 * (t_[1:] + t_[:-1]), y=np.cumsum(l_), mode="lines", name=f"type {i}", line=dict(color=COLORS[(i + 1) % len(COLORS)]))
                )
    except Exception as e:
        print(f"_plot_rootsystem_development: {e}")
    fig.update_layout(
        xaxis_title="Time (days)",
        yaxis_title=f"Root system {type_str[j]} {unit_str[j]}",
        legend=dict(x=0.02, y=0.98),
        margin=dict(l=60, r=20, t=30, b=60),
        template="plotly_white",
    )
    return fig


def _plot_suf(data, j: int) -> go.Figure:
    _init_hydraulic_scenario(data.hydraulic_params, j)
    krs, _ = data.hydraulic_model.get_krs(data.max_ct)
    suf = data.hydraulic_model.get_suf(data.max_ct)
    data.analyser.addData("SUF", suf)
    n = int(np.ceil(-data.analyser.getMinBounds().z))
    z_ = np.linspace(-0.5, -n + 0.5, n)
    d = data.analyser.distribution("SUF", 0.0, float(-n), int(n), False)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=d, y=z_, mode="lines+markers", name="total", line=dict(color=COLORS[0])))
    max_type = int(np.max(data.analyser.data["subType"]))
    for i in range(max_type + 1):
        ana = pb.SegmentAnalyser(data.analyser)
        ana.filter("subType", i)
        if len(ana.segments) > 0:
            d = ana.distribution("SUF", 0.0, float(-n), int(n), False)
            fig.add_trace(go.Scatter(x=d, y=z_, mode="lines+markers", name=f"type {i}", line=dict(color=COLORS[(i + 1) % len(COLORS)])))
    fig.update_layout(
        title=f"Root system krs {krs:.4g} cm²/day",
        xaxis_title="Surface uptake fraction (SUF) per 1 cm layer (1)",
        yaxis_title="Depth (cm)",
        legend=dict(x=0.7, y=0.02),
        margin=dict(l=60, r=20, t=50, b=60),
        template="plotly_white",
    )
    return fig


def _plot_krs(data, j: int) -> go.Figure:
    _init_hydraulic_scenario(data.hydraulic_params, j)
    t = max(1, int(np.ceil(data.max_ct)))
    t_ = np.linspace(1, t, t)
    krs_ = [data.hydraulic_model.get_krs(float(t_val))[0] for t_val in t_]
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=t_, y=krs_, mode="lines", name="Krs", line=dict(color=COLORS[0])))
    fig.update_layout(
        xaxis_title="Time (days)",
        yaxis_title="Root system conductivity Krs [cm²/day]",
        margin=dict(l=60, r=20, t=30, b=60),
        template="plotly_white",
    )
    return fig


if __name__ == "__main__":
    app.run(debug=True)
