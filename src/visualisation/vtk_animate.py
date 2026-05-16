"""
vtk_animate.py  –  VTK-based animation of CPlantBox root systems

Provides :class:`AnimateRoots`, which renders a growing root system (and
optionally a soil grid / container) frame-by-frame in an interactive VTK
window or writes the frames to JPEG files for offline video assembly.

Typical usage::

    rs = pb.RootSystem()
    rs.openFile("params.xml")

    anim = AnimateRoots(rs)
    anim.root_name = "age"
    anim.start()           # open VTK window and configure camera

    for t in simulation_times:
        rs.simulate(dt)
        anim.update()      # refresh the scene after each time step

    anim.run()             # enter the VTK event loop (blocking)

Author: Daniel Leitner
"""

import vtk

import plantbox as pb
from plantbox.visualisation.vtk_plot import (
    plot_container,
    plot_mesh,
    plot_plant,
    plot_roots,
    render_window,
    segs_to_polydata,
    uniform_grid,
    write_jpg,
)


class AnimateRoots:
    """Interactive VTK animation of a growing root system.

    The class manages a VTK render window that can be refreshed on demand
    (e.g. driven by a simulation timer) to show the current state of a
    :class:`plantbox.RootSystem` or :class:`plantbox.Plant`.  Optionally a
    soil wireframe grid and / or a container surface can be overlaid.

    Attributes
    ----------
    rootsystem : pb.Organism or pb.MappedSegments, optional
        The plant / root system to visualise.  If *None* no root actors are
        added.
    root_name : str
        Name of the scalar field used to colour the root segments.
        Defaults to ``"subType"``.
    plant : bool
        If *True* use :func:`plot_plant` (renders leaf surfaces as polygons
        in addition to the tube plot); if *False* use :func:`plot_roots`
        (tube plot only).  Default *False*.
    container_sdf : plantbox SDF, optional
        Signed-distance function of the container shape.  When set the
        container surface is computed once in :meth:`start` and reused each
        frame.
    soil : vtkUniformGrid or dumux-rosi solver, optional
        Soil data source used by :meth:`create_soil_actors`.  Must be set
        when :attr:`soil_data` is *True*.
    soil_data : bool
        Whether to visualise the soil grid.  Defaults to *False*.
    cuts : bool
        If *True* render the soil as pressure-head YZ cut planes (requires
        ``plot_mesh_yz`` to be implemented); if *False* render as a wireframe
        bounding box.  Default *False*.
    min, max, res : array-like
        Bounding-box minimum / maximum and cell resolution for the soil grid.
        Only used when :attr:`soil_data` is *True*.
    avi_name : str, optional
        Base file name for JPEG frame export.  When set, each call to
        :meth:`update` writes ``<avi_name><frame_index>.jpg`` to disk.
    """

    def __init__(self, rootsystem=None, container_sdf=None):
        # Root system data source
        self.rootsystem = rootsystem
        self.root_name = "subType"
        self.plant = False  # True → include leaf polygons via plot_plant

        # Container (computed once in start, then reused)
        self.container_sdf = container_sdf
        self._container_actor = None

        # Soil overlay
        self.soil = None
        self.soil_data = False  # set to True together with self.soil to enable
        self.cuts = False  # True → YZ cut planes; False → wireframe box
        self.min = None
        self.max = None
        self.res = None

        # VTK state (populated by start())
        self.actors = []  # current list of vtkActors in the scene
        self.iren = None  # vtkRenderWindowInteractor
        self.color_bar = None  # active vtkScalarBarActor
        self.bounds = None  # [xmin, xmax, ymin, ymax, zmin, zmax]

        # Optional JPEG frame export
        self.avi_name = None
        self._frame_count = 0

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def start(self, axis="x"):
        """Open the VTK render window and position the camera.

        Builds the initial scene (root actors, soil actors, and the container
        if :attr:`container_sdf` is set) and creates the interactor.
        Call :meth:`run` afterwards to enter the interactive event loop, or
        call :meth:`update` repeatedly from a simulation loop first.

        Parameters
        ----------
        axis : str
            Initial camera orientation.  One of:

            * ``"x"`` – look along the X-axis (Y-Z plane visible)
            * ``"y"`` – look along the Y-axis (X-Z plane visible)
            * ``"z"`` – look along the Z-axis (X-Y plane visible)
            * ``"v"`` – oblique view (30° azimuth / 30° elevation)
        """
        if self.container_sdf:
            container_actors, _, _ = plot_container(self.container_sdf, render=False)
            self._container_actor = container_actors[0]

        self.create_root_actors()
        self.create_soil_actors()
        if self._container_actor:
            self.actors.append(self._container_actor)

        self.iren = render_window(self.actors, "AnimateRoots", self.color_bar, self.bounds)
        self._set_camera(axis)

    def run(self):
        """Enter the VTK event loop (blocking).

        Starts the :class:`vtkRenderWindowInteractor`.  The window remains
        open until the user presses *e* (exit).

        Raises
        ------
        RuntimeError
            If :meth:`start` has not been called yet.
        """
        if self.iren is None:
            raise RuntimeError("Call start() before run().")
        self.iren.Start()

    def update(self):
        """Refresh the scene with the current simulation state.

        Removes all existing root and soil actors, rebuilds them from the
        current :attr:`rootsystem` / :attr:`soil` data, and triggers a VTK
        re-render.  The container actor is preserved across frames (it does
        not change during a simulation).

        If :attr:`avi_name` is set the rendered frame is also written to a
        JPEG file named ``<avi_name><frame_index>.jpg``.

        Raises
        ------
        RuntimeError
            If :meth:`start` has not been called yet.
        """
        if self.iren is None:
            raise RuntimeError("Call start() before update().")

        renWin = self.iren.GetRenderWindow()
        ren = renWin.GetRenderers().GetFirstRenderer()

        # Remove root / soil actors from the previous frame
        for a in self.actors:
            ren.RemoveActor(a)
        # Remove 2-D overlays (scalar bar) so they can be re-added below
        for a in ren.GetActors2D():
            ren.RemoveActor2D(a)
        if self.color_bar:
            ren.AddActor2D(self.color_bar)

        self.actors = []
        self.create_root_actors()
        self.create_soil_actors()
        if self._container_actor:
            self.actors.append(self._container_actor)
        for a in self.actors:
            ren.AddActor(a)

        self.iren.Render()

        if self.avi_name:
            self._save_frame(renWin)

    # ------------------------------------------------------------------
    # Actor builders
    # ------------------------------------------------------------------

    def create_root_actors(self):
        """Build VTK actors for the current root system state.

        Appends one or more actors to :attr:`actors` and updates
        :attr:`color_bar` and :attr:`bounds`.  Does nothing when
        :attr:`rootsystem` is *None*.
        """
        if self.rootsystem is None:
            return

        pd = segs_to_polydata(self.rootsystem, 1.0, [self.root_name, "radius"])

        if self.plant:
            new_actors, root_cbar = plot_plant(self.rootsystem, self.root_name, render=False)
        else:
            new_actors, root_cbar = plot_roots(pd, self.root_name, "", render=False)

        if isinstance(new_actors, list):
            self.actors.extend(new_actors)
        else:
            self.actors.append(new_actors)

        self.color_bar = root_cbar
        self.bounds = pd.GetBounds()

    def create_soil_actors(self):
        """Build VTK actors for the soil grid overlay.

        Only executed when :attr:`soil_data` is *True* **and** :attr:`soil`
        is not *None*.  Appends actors to :attr:`actors` and updates
        :attr:`bounds`.

        Two display modes are available, controlled by :attr:`cuts`:

        * *False* (default) – render the grid as a wireframe bounding box
          using :func:`plot_mesh`.
        * *True* – render YZ cut planes of the pressure head.  This requires
          ``plot_mesh_yz`` (not yet implemented in vtk_plot.py).

        Raises
        ------
        NotImplementedError
            When :attr:`cuts` is *True* (``plot_mesh_yz`` is not yet available).
        """
        if not self.soil_data or self.soil is None:
            return

        if self.cuts:
            # TODO: implement plot_mesh_yz in vtk_plot.py and import it here.
            raise NotImplementedError("Cut-plane soil visualisation requires plot_mesh_yz, " "which is not yet available in vtk_plot.")

        grid = uniform_grid(self.min, self.max, self.res)
        mesh_actors, _ = plot_mesh(grid, "", "", render=False)
        if isinstance(mesh_actors, list):
            self.actors.extend(mesh_actors)
        else:
            self.actors.append(mesh_actors)
        self.bounds = grid.GetBounds()

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _set_camera(self, axis):
        """Position the active camera according to *axis*.

        Parameters
        ----------
        axis : str
            One of ``"x"``, ``"y"``, ``"z"``, ``"v"`` (see :meth:`start`).
        """
        if self.bounds is None or self.iren is None:
            return
        renWin = self.iren.GetRenderWindow()
        ren = renWin.GetRenderers().GetItemAsObject(0)
        camera = ren.GetActiveCamera()
        z_mid = 0.5 * (self.bounds[4] + self.bounds[5])
        if axis == "x":
            camera.SetPosition([100, 0, z_mid])
            camera.SetViewUp(0, 0, 1)
        elif axis == "y":
            camera.SetPosition([0, 100, z_mid])
            camera.SetViewUp(0, 0, 1)
        elif axis == "z":
            camera.SetPosition([0, 0, 100])
            camera.SetViewUp(0, 1, 0)
        elif axis == "v":
            camera.SetPosition([100, 0, z_mid])
            camera.SetViewUp(0, 0, 1)
            camera.Azimuth(30)
            camera.Elevation(30)

    def _save_frame(self, renWin):
        """Write the current render window to a JPEG file.

        The file is named ``<avi_name><frame_count>.jpg`` and the internal
        frame counter is incremented.

        Parameters
        ----------
        renWin : vtkRenderWindow
            The render window to capture.
        """
        filename = f"{self.avi_name}{self._frame_count}"
        write_jpg(renWin, filename)
        print(f"saved {filename}.jpg")
        self._frame_count += 1
