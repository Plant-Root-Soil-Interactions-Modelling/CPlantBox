import os, sys
from pathlib import Path

from PyQt5 import QtWidgets, QtCore, QtGui
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import vtk

# Add CPlantBox build and src
root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / "build" / "Release"))
sys.path.insert(0, str(root / "src"))

# General
from PyQt5 import QtWidgets, QtCore, QtGui
import vtk
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import colorsys

# CPB
import functional.xylem_flux as xylem_flux
import visualisation.vtk_plot_rhizomancer as vp
import visualisation.vtk_tools as vt
from visualisation.vtk_plot_rhizomancer import keypress_callback_
from viewer_data import ViewerDataModel
import viewer_plots, viewer_conductivities


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.bounds = (0,0,0,0,0,0) 
        self.data_properties = {
        

            'subType':      {'title': 'Type (-)', 
                             'format': ticker.StrMethodFormatter('{x:.0f}'), 
                             'type': 'categorical'}, 
            'creationTime': {'title': 'Creation time (days)', 
                             'format': ticker.StrMethodFormatter('{x:.0f}')},
            'length':       {'title': 'Length (m)', 
                             'format': ticker.StrMethodFormatter('{x:.3f}')},
            'SUF':          {'title': 'SUF (-)', 
                             'format': ticker.StrMethodFormatter('{x:.2e}')},
                             
             'rootId':      {'title': 'Root ID (-)',
                             'format': ticker.StrMethodFormatter('{x:.0f}'),
                             'type': 'categorical'},
        }
        self.camera_is_set = False
        self.current_file = None
        self.setWindowTitle("RSML Viewer by Daniel Leitner")
        self.resize(1600, 900)
        self.data = ViewerDataModel()
        self._create_actions()
        self._create_menus()
        self._create_central_widget()
        self._create_control_dock()

    def _create_actions(self):
        self.openAct = QtWidgets.QAction("Open RSML...", self, triggered=self.open_file)
        self.saveAct = QtWidgets.QAction("Save RSML...", self, triggered=self.save_file)
        self.saveVtpAct = QtWidgets.QAction("Save VTP...", self, triggered=self.save_vtp)
        self.exitAct = QtWidgets.QAction("Exit", self, triggered=self.close)

    def _create_menus(self):
        menubar = self.menuBar()
        fileMenu = menubar.addMenu("File")
        fileMenu.addAction(self.openAct)
        fileMenu.addAction(self.saveAct)
        fileMenu.addAction(self.saveVtpAct)
        fileMenu.addSeparator()
        fileMenu.addAction(self.exitAct)
    
    def _create_control_dock(self):
        # Create a fixed, left-side control panel without dock handles
        dock = QtWidgets.QDockWidget("Controls", self)
        dock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea)
        dock.setFeatures(QtWidgets.QDockWidget.NoDockWidgetFeatures)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)

        # Build the control widgets container
        container = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(container)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(12)

        # ——— A) Plot Data Modes ———
        plot_group = QtWidgets.QGroupBox("Plot Data")
        pg_layout = QtWidgets.QVBoxLayout(plot_group)
        for key, props in self.data_properties.items():
            btn = QtWidgets.QPushButton(props['title'])
            btn.clicked.connect(lambda checked, k=key: self.render_3d(k))
            pg_layout.addWidget(btn)
        layout.addWidget(plot_group)
        
        # ——— C) Camera Views ———
        cam_group = QtWidgets.QGroupBox("Camera Views")
        cg_layout = QtWidgets.QVBoxLayout(cam_group)
        for label, slot, tip in [
            ("X-View (X)", self._on_x_view, "X-axis aligned"),
            ("Y-View (Y)", self._on_y_view, "Y-axis aligned"),
            ("Z-View (Z)", self._on_z_view, "Z-axis aligned"),
            ("Tilted (V)", self._on_tilted_view, "Tilted view (3D)"),
            ("Center roots (R)", self._on_reset_view, "Fit to roots"),
        ]:
            btn = QtWidgets.QPushButton(label)
            btn.setToolTip(f"{tip} — shortcut: {label[-2]}")
            btn.clicked.connect(slot)
            cg_layout.addWidget(btn)
        layout.addWidget(cam_group)

        # ——— B) Screenshots ———
        save_group = QtWidgets.QGroupBox("Screenshots")
        sg_layout = QtWidgets.QVBoxLayout(save_group)
        btn_hires = QtWidgets.QPushButton("High-res (G)")
        btn_hires.clicked.connect(lambda: self._on_screenshot(mag=5))
        btn_quick = QtWidgets.QPushButton("Window-size (S)")
        btn_quick.clicked.connect(lambda: self._on_screenshot(mag=1))
        sg_layout.addWidget(btn_hires)
        sg_layout.addWidget(btn_quick)
        layout.addWidget(save_group)



        # Fill remaining space
        layout.addStretch(1)

        # Set the widget and hide the title bar for a clean look
        dock.setWidget(container)
        dock.setTitleBarWidget(QtWidgets.QWidget())
        self.control_dock = dock
 
    def _on_screenshot(self, mag):
        """Take a VTK+legend JPEG at given magnification."""
        fn = self.vtk_widget.GetRenderWindow().GetWindowName()
        vp.write_jpg(
            self.vtk_widget.GetRenderWindow(),
            fn,
            magnification=mag
        )

    def _on_x_view(self):
        cam = self.renderer.GetActiveCamera()
        cam.SetPosition([100, 0, 0.5 * (self.bounds[4] + self.bounds[5])])
        cam.SetViewUp(0, 0, 1)
        cam.OrthogonalizeViewUp()
        self.renderer.ResetCameraClippingRange()
        self.vtk_widget.GetRenderWindow().Render()

    def _on_y_view(self):
        cam = self.renderer.GetActiveCamera()
        cam.SetPosition([0, 100, 0.5 * (self.bounds[4] + self.bounds[5])])
        cam.SetViewUp(0, 0, 1)
        cam.OrthogonalizeViewUp()
        self.renderer.ResetCameraClippingRange()
        self.vtk_widget.GetRenderWindow().Render()

    def _on_z_view(self):
        cam = self.renderer.GetActiveCamera()
        cam.SetPosition([0, 0, 100])
        cam.SetViewUp(0, 1, 0)
        cam.OrthogonalizeViewUp()
        self.renderer.ResetCameraClippingRange()
        self.vtk_widget.GetRenderWindow().Render()

    def _on_tilted_view(self):
        cam = self.renderer.GetActiveCamera()
        cam.SetPosition([100, 0, 0.5 * (self.bounds[4] + self.bounds[5])])
        cam.SetViewUp(0, 0, 1)
        cam.Azimuth(30)
        cam.Elevation(30)
        cam.OrthogonalizeViewUp()
        self.renderer.ResetCameraClippingRange()
        self.vtk_widget.GetRenderWindow().Render()

    def _on_reset_view(self):
        """Reset the camera to fit the root system bounds."""
        self.renderer.ResetCamera()
        self.renderer.ResetCameraClippingRange()
        self.vtk_widget.GetRenderWindow().Render()
    
    
    def _create_colorbar_widget(self):
        qt_bg_color = self.palette().color(QtGui.QPalette.Window)
        self.colorbar_fig = Figure(figsize=(1, 8), facecolor='#f7f7f7')
        self.colorbar_canvas = FigureCanvas(self.colorbar_fig)
        self.colorbar_canvas.hide() # hide colorbar until RSML is loaded
        
        self.colorbar_canvas.setFixedWidth(220)
        self.cax = self.colorbar_fig.add_axes([0.25, 0.1, 0.15, 0.8])
        
        return self.colorbar_canvas

    def _vtk_to_mpl_colormap(self, vtk_lookup_table):
        """Converts a vtkLookupTable to a matplotlib colormap."""
        num_colors = vtk_lookup_table.GetNumberOfTableValues()
        colors = []
        for i in range(num_colors):
            r, g, b, a = vtk_lookup_table.GetTableValue(i)
            colors.append((r, g, b, a))
        return mpl.colors.ListedColormap(colors)

    def _update_colorbar(self, lookup_table, data_range, data_key):
        """Clears and redraws the Matplotlib colorbar."""
        self.cax.clear()

        props = self.data_properties.get(data_key, {'title': data_key})
        cmap = self._vtk_to_mpl_colormap(lookup_table)

        if props.get('type') == 'categorical':
            ticks = np.arange(int(data_range[0]), int(data_range[1]) + 1)
            boundaries = np.arange(int(data_range[0]) - 0.5, int(data_range[1]) + 1.5, 1)
            norm = mpl.colors.BoundaryNorm(boundaries, cmap.N)
        else: 
            ticks = np.linspace(data_range[0], data_range[1], 5)
            norm = mpl.colors.Normalize(vmin=data_range[0], vmax=data_range[1])

        cb = mpl.colorbar.ColorbarBase(self.cax,
                                       cmap=cmap,
                                       norm=norm,
                                       orientation='vertical',
                                       format=props.get('format'))
        
        cb.set_ticks(ticks)
        cb.set_label(props['title'], fontsize=14, labelpad=15)
        cb.ax.tick_params(labelsize=12)

        self.colorbar_canvas.draw()
        
    def _create_central_widget(self):
            # Main split: left panel (VTK + Colorbar) and right panel (Tabs)
            central = QtWidgets.QWidget()
            layout = QtWidgets.QHBoxLayout(central)

            # --- Create the Left Panel ---
            left_panel = QtWidgets.QWidget()
            left_layout = QtWidgets.QHBoxLayout(left_panel)
            left_layout.setContentsMargins(0, 0, 0, 0)
            left_layout.setSpacing(0)
            
            # VTK widget
            self.vtk_widget = QVTKRenderWindowInteractor()
            self.vtk_widget.setMinimumWidth(600)
            self.renderer = vtk.vtkRenderer()
            self.renderer.SetBackground(0.97, 0.97, 0.97)
            self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)
            self.vtk_widget.Initialize()
            self.vtk_widget.Start()

            self.vtk_widget.AddObserver(
                "KeyPressEvent",
                lambda obj, ev: keypress_callback_(obj, ev, self.bounds)
            )
            
            # Add the VTK widget and the new Colorbar widget to the left panel
            left_layout.addWidget(self.vtk_widget)
            left_layout.addWidget(self._create_colorbar_widget())

            # --- Create the Right Panel (Tabs) ---
            self.tabs = QtWidgets.QTabWidget()
            self._init_info_tab()
            self._init_depth_tab()
            self._init_dev_tab()
            self._init_suf_tab()
            self._init_krs_tab()

            # Add the left and right panels to the main layout
            layout.addWidget(left_panel, stretch=1)
            layout.addWidget(self.tabs, stretch=1)
            
            self.setCentralWidget(central)
            self.setStyleSheet("#buttonContainer { background-color: #f7f7f7;}")

    def _init_info_tab(self):
        info = QtWidgets.QWidget()
        self.tabs.addTab(info, "Information")

        # Outer layout for the tab
        vbox = QtWidgets.QVBoxLayout(info)

        # --- General group ---
        grp_general = QtWidgets.QGroupBox("General")
        grid_gen   = QtWidgets.QGridLayout(grp_general)
        vbox.addWidget(grp_general)

        # Create the two QLabel columns
        self.label_general_l = QtWidgets.QLabel()
        self.label_general_r = QtWidgets.QLabel()
        # Make left column fixed-width and align top
        self.label_general_l.setFixedWidth(150)
        self.label_general_l.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft)
        self.label_general_r.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft)
        self.label_general_l.setText(
            "Software\n"
            "Filename\n"
            "Number of plants (base nodes)\n"
            "Number of base roots (segments)\n"
            "Number of roots\n"
            "Number of nodes\n"
            "Bounding box\n"
            "Unit (length scale)\n"
            "Resolution\n"
        )
        grid_gen.addWidget(self.label_general_l, 0, 0)
        grid_gen.addWidget(self.label_general_r, 0, 1)

        # --- Properties group ---
        grp_prop = QtWidgets.QGroupBox("Properties (values per root)")
        grid_prop = QtWidgets.QGridLayout(grp_prop)
        vbox.addWidget(grp_prop)

        self.label_prop_l = QtWidgets.QLabel()
        self.label_prop_r = QtWidgets.QLabel()
        self.label_prop_l.setFixedWidth(150)
        self.label_prop_l.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft)
        self.label_prop_r.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft)
        self.label_prop_l.setText("Property\n")
        grid_prop.addWidget(self.label_prop_l, 0, 0)
        grid_prop.addWidget(self.label_prop_r, 0, 1)

        # --- Functions group ---
        grp_fun = QtWidgets.QGroupBox("Functions (values per node)")
        grid_fun = QtWidgets.QGridLayout(grp_fun)
        vbox.addWidget(grp_fun)

        self.label_fun_l = QtWidgets.QLabel()
        self.label_fun_r = QtWidgets.QLabel()
        self.label_fun_l.setFixedWidth(150)
        self.label_fun_l.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft)
        self.label_fun_r.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft)
        self.label_fun_l.setText("Function\n")
        grid_fun.addWidget(self.label_fun_l, 0, 0)
        grid_fun.addWidget(self.label_fun_r, 0, 1)

        # --- Using group ---
        grp_use = QtWidgets.QGroupBox("Using")
        grid_use = QtWidgets.QGridLayout(grp_use)
        vbox.addWidget(grp_use)

        self.label_use_l = QtWidgets.QLabel()
        self.label_use_r = QtWidgets.QLabel()
        self.label_use_l.setFixedWidth(150)
        self.label_use_l.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft)
        self.label_use_r.setAlignment(QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft)
        self.label_use_l.setText(
            "Radius\n"
            "Creation time\n"
            "Types\n"
        )
        grid_use.addWidget(self.label_use_l, 0, 0)
        grid_use.addWidget(self.label_use_r, 0, 1)

        # stretch so the info tab expands vertically
        vbox.addStretch(1)


    def _init_depth_tab(self):
        depth = QtWidgets.QWidget()
        self.tabs.addTab(depth, "Depth Profile")
        layout = QtWidgets.QVBoxLayout(depth)
        self.combo_depth = QtWidgets.QComboBox();
        self.combo_depth.addItems(["Length","Surface","Volume"])
        self.combo_depth.currentIndexChanged.connect(self.update_depth)
        layout.addWidget(self.combo_depth)
        fig, self.ax_depth = plt.subplots()
        self.canvas_depth = FigureCanvas(fig)
        layout.addWidget(self.canvas_depth)

    def _init_dev_tab(self):
        dev = QtWidgets.QWidget()
        self.tabs.addTab(dev, "Development")
        layout = QtWidgets.QVBoxLayout(dev)
        self.combo_dev = QtWidgets.QComboBox();
        self.combo_dev.addItems(["Length","Surface","Volume"])
        self.combo_dev.currentIndexChanged.connect(self.update_dev)
        layout.addWidget(self.combo_dev)
        fig, self.ax_dev = plt.subplots()
        self.canvas_dev = FigureCanvas(fig)
        layout.addWidget(self.canvas_dev)

    def _init_suf_tab(self):
        suf = QtWidgets.QWidget()
        self.tabs.addTab(suf, "Hydraulic Properties")
        layout = QtWidgets.QVBoxLayout(suf)
        self.combo_suf = QtWidgets.QComboBox();
        self.combo_suf.addItems(["Const1","Const2","Dyn1","Dyn2"])
        self.combo_suf.currentIndexChanged.connect(self.update_suf)
        layout.addWidget(self.combo_suf)
        fig, self.ax_suf = plt.subplots()
        self.canvas_suf = FigureCanvas(fig)
        layout.addWidget(self.canvas_suf)

    def _init_krs_tab(self):
        krs = QtWidgets.QWidget()
        self.tabs.addTab(krs, "Hydraulic Development")
        layout = QtWidgets.QVBoxLayout(krs)
        self.combo_krs = QtWidgets.QComboBox();
        self.combo_krs.addItems(["Const1","Const2","Dyn1","Dyn2"])
        self.combo_krs.currentIndexChanged.connect(self.update_krs)
        layout.addWidget(self.combo_krs)
        fig, self.ax_krs = plt.subplots()
        self.canvas_krs = FigureCanvas(fig)
        layout.addWidget(self.canvas_krs)

    def open_file(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open RSML", "", "RSML Files (*.rsml *.xml)")
        if path:
            self.camera_is_set = False
            self.data.open_rsml(path)
            self.current_file = path
            self.update_all()
            self.render_3d('subType')

    def save_file(self):
        if self.data.exists() and self.current_file:
            fn, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save RSML", "", "RSML Files (*.rsml)")
            if fn:
                pd = vp.segs_to_polydata(self.data.analyser)
                vt.write_rsml(fn, pd, 0, None, self.data.base_nodes)

    def save_vtp(self):
        if self.data.exists():
            fn, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save VTP", "", "VTP Files (*.vtp)")
            if fn:
                self.data.analyser.write(fn)

    def update_all(self):
        if self.data.exists():
            self.update_info()
            self.update_depth()
            self.update_dev()
            self.update_suf()
            self.update_krs()

    def update_info(self):
        """ Update the Information tab labels in the PyQt viewer """
        if getattr(self, 'current_file', None):
            fname = os.path.basename(self.current_file)
        else:
            try:
                fname = os.path.basename(self.data.filename)
            except Exception:
                fname = ''

        # --- GENERAL INFO ---
        total_nodes = sum(len(pl) for pl in self.data.polylines)
        meta        = self.data.metadata
        nop         = len(self.data.base_nodes)
        nobr        = len(self.data.base_segs)
        min_bb      = self.data.analyser.getMinBounds()
        max_bb      = self.data.analyser.getMaxBounds()
        res_str     = str(meta.resolution)

        # Left and right text blocks
        left_gen = (
            "Software\n"
            "Filename\n"
            "Number of plants (base nodes)\n"
            "Number of base roots (segments)\n"
            "Number of roots\n"
            "Number of nodes\n"
            "Bounding box\n"
            "Unit (length scale)\n"
            "Resolution\n"
        )
        right_gen = (
            f"{meta.software}\n"
            f"{fname}\n"
            f"{nop:g}\n"
            f"{nobr:g}\n"
            f"{len(self.data.polylines):g}\n"
            f"{total_nodes:g}\n"
            f"{min_bb} - {max_bb} (cm)\n"
            f"{meta.unit}\n"
            f"{res_str} (dots per {meta.unit})\n"
        )

        self.label_general_l.setText(left_gen)
        self.label_general_r.setText(right_gen)

        # --- PROPERTIES (per root) ---
        left_prop = ""
        right_prop = ""
        for key, vals in self.data.properties.items():
            arr = np.array(vals)
            left_prop += key + "\n"
            if arr.size:
                mn, mx = arr.min(), arr.max()
                unit   = getattr(meta.properties.get(key, None), 'unit', '')
                right_prop += f"[{mn:g}, {mx:g}] {unit}\n"
            else:
                right_prop += "[]\n"

        self.label_prop_l.setText(left_prop)
        self.label_prop_r.setText(right_prop)

        # --- FUNCTIONS (per node) ---
        left_fun = ""
        right_fun = ""
        for key, seqs in self.data.functions.items():
            flat = np.concatenate(seqs) if seqs else np.array([])
            left_fun += key + "\n"
            if flat.size:
                mn, mx = flat.min(), flat.max()
                unit   = getattr(meta.properties.get(key, None), 'unit', '')
                right_fun += f"[{mn:g}, {mx:g}] {unit}\n"
            else:
                right_fun += "[]\n"

        self.label_fun_l.setText(left_fun)
        self.label_fun_r.setText(right_fun)

        # --- USING tags ---
        tags     = self.data.tagnames
        left_use = (
            "Radius\n"
            "Creation time\n"
            "Types\n"
        )
        right_use = ""
        # Radii tag
        if tags[0] and hasattr(self.data, 'radii'):
            r = np.array(self.data.radii)
            right_use += (
                f"from tag '{tags[0]}' within [{r.min():g}, {r.max():g}] cm\n"
                if r.size else "no radii data\n"
            )
        else:
            right_use += "not found (default)\n"
        # Creation-time tag
        if tags[1] and hasattr(self.data, 'cts'):
            c = np.array(self.data.cts)
            right_use += (
                f"from tag '{tags[1]}' within [{c.min():g}, {c.max():g}] days\n"
                if c.size else "no creation-time data\n"
            )
        else:
            right_use += "not found\n"
        # Types tag
        if tags[2] and hasattr(self.data, 'types'):
            t = np.array(self.data.types)
            right_use += (
                f"from tag '{tags[2]}' within [{t.min():g}, {t.max():g}]\n"
                if t.size else "no types data\n"
            )
        else:
            right_use += "not found\n"

        self.label_use_l.setText(left_use)
        self.label_use_r.setText(right_use)


    def update_depth(self):
        if self.data.exists():
            viewer_plots.plot_depth_profile(self.data.analyser, self.ax_depth, self.combo_depth.currentIndex())
            self.canvas_depth.draw()

    def update_dev(self):
        if self.data.exists():
            viewer_plots.plot_rootsystem_development(self.data.analyser, self.ax_dev, self.combo_dev.currentIndex())
            self.canvas_dev.draw()

    def update_suf(self):
        if self.data.exists():
            viewer_conductivities.init_constant_scenario1(self.data.xylem_flux)
            viewer_plots.plot_suf(self.data, self.ax_suf, self.combo_suf.currentIndex())
            self.canvas_suf.draw()

    def update_krs(self):
        if self.data.exists():
            viewer_plots.plot_krs(self.data, self.ax_krs, self.combo_krs.currentIndex())
            self.canvas_krs.draw()
            
    def _build_rootid_actor(self):
        """
        Builds a VTK actor where each RSML polyline is colored by its Root ID (1..N).

        Notes:
        - Root IDs here are 1-based and correspond to the order of self.data.polylines.
        - This is meant as a visual assessment tool (quick QA), not a hydraulic/metric plot.
        """
        if not getattr(self.data, "polylines", None):
            return None, (0, 0)

        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()

        root_id_arr = vtk.vtkIntArray()
        root_id_arr.SetName("rootId")

        pid = 0
        n_roots = 0

        for rid, pl in enumerate(self.data.polylines, start=1):
            if pl is None or len(pl) < 2:
                continue
            n_roots = max(n_roots, rid)

            polyline = vtk.vtkPolyLine()
            polyline.GetPointIds().SetNumberOfIds(len(pl))

            for j, p in enumerate(pl):
                # p is expected to be (x, y, z) in RSML coordinates
                points.InsertNextPoint(float(p[0]), float(p[1]), float(p[2]))
                root_id_arr.InsertNextValue(int(rid))
                polyline.GetPointIds().SetId(j, pid)
                pid += 1

            lines.InsertNextCell(polyline)

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetLines(lines)
        polydata.GetPointData().SetScalars(root_id_arr)

        # Lookup table with N distinct-ish hues
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfTableValues(max(1, n_roots))
        lut.SetRange(1, max(1, n_roots))
        lut.Build()

        phi = 0.618033988749895  # golden ratio conjugate
        h = 0.0
        for i in range(max(1, n_roots)):
            h = (h + phi) % 1.0
            r, g, b = colorsys.hsv_to_rgb(h, 0.95, 0.95)
            lut.SetTableValue(i, r, g, b, 1.0)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polydata)
        mapper.SetLookupTable(lut)
        mapper.SetScalarRange(1, max(1, n_roots))
        mapper.SelectColorArray("rootId")
        mapper.SetScalarModeToUsePointData()

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        return actor, (1, max(1, n_roots))

    
    def render_3d(self, name):
        if not self.data.exists():
            QtWidgets.QMessageBox.warning(self, "Warning", "Open RSML file first")
            return

        self.colorbar_canvas.show() # show colorbar when RSML is loaded
        
        # If a root actor already exists, remove it
        if hasattr(self, 'root_actor') and self.root_actor:
            self.renderer.RemoveActor(self.root_actor)

        # Create the new 3D actor (plot_roots should now only return the actor)
        # Create the new 3D actor
        if name == "rootId":
            self.root_actor, data_range = self._build_rootid_actor()
        else:
            self.root_actor, data_range = vp.plot_roots(
                self.data.analyser, name,
                render=False, interactiveImage=False
            )

        if self.root_actor is None:
            QtWidgets.QMessageBox.warning(self, "Warning", f"Could not render mode: {name}")
            return

        self.renderer.AddActor(self.root_actor)

        # Optional but recommended: keep bounds updated for camera/key handlers that use self.bounds
        try:
            self.bounds = self.root_actor.GetBounds()
        except Exception:
            pass


        # Update the Matplotlib colorbar with info from the new actor
        mapper = self.root_actor.GetMapper()
        self._update_colorbar(mapper.GetLookupTable(), data_range, name)

        # Camera logic remains the same
        if not self.camera_is_set:
            self.renderer.ResetCamera()
            cam = self.renderer.GetActiveCamera()
            cam.ParallelProjectionOn()
            cam.Azimuth(30)
            cam.Elevation(60)
            cam.SetViewUp(0, 0, 1)
            cam.OrthogonalizeViewUp()
            self.renderer.ResetCameraClippingRange()
            self.camera_is_set = True

        self.vtk_widget.GetRenderWindow().Render()

    def multiple_suf(self):
        self.render_3d('SUF')

def create_rsml_window(rsml_path: str = None):
    
    win = MainWindow()
    if rsml_path:
        win.data.open_rsml(rsml_path)
        win.update_all()
    return win

if __name__ == "__main__":
    import sys
    from qtpy.QtWidgets import QApplication, QFileDialog

    app = QApplication(sys.argv)

    # Make the RSML path argument optional:
    if len(sys.argv) < 2:
        # no file provided => just open the viewer without loading
        win = create_rsml_window(None)
    else:
        # file provided => load it
        path = sys.argv[1]
        win = create_rsml_window(path)

    win.show()
    if len(sys.argv) >= 2:
        def _init_view():
            win.render_3d('subType')   # Type
            win._on_x_view()           # X-View
        QtCore.QTimer.singleShot(0, _init_view)
    sys.exit(app.exec_())