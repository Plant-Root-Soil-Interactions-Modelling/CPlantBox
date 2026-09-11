"""
Unit tests for the VTK visualisation helpers (vtk_plot, and the 2-D overlay
handling shared with vtk_animate).

The whole module is skipped when vtk is missing, as it is an optional
dependency (the "vis" extra).  Tests that need an OpenGL context are skipped
in addition when offscreen rendering is unavailable, e.g. on a headless
machine without Mesa, so the remaining tests still run there.
"""

import sys

sys.path.append("..")
sys.path.append("../src/")

import inspect
import os
import subprocess
import tempfile
import unittest

import plantbox as pb

try:
    import vtk

    from plantbox.visualisation import vtk_animate
    from plantbox.visualisation import vtk_plot as vp
except ImportError:
    vtk = None

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

_OFFSCREEN_PROBE = """
import vtk
renWin = vtk.vtkRenderWindow()
renWin.SetOffScreenRendering(1)
renWin.SetSize(10, 10)
renWin.AddRenderer(vtk.vtkRenderer())
renWin.Render()
"""


def offscreen_rendering_works():
    """
    True if VTK can render offscreen.

    Probed in a subprocess because a failing OpenGL context aborts the
    interpreter instead of raising, which would take the test run with it.
    """
    if vtk is None:
        return False
    return subprocess.run([sys.executable, "-c", _OFFSCREEN_PROBE], capture_output=True).returncode == 0


HAS_OFFSCREEN = offscreen_rendering_works()


def make_linear_sa(n_segs=4, radius=0.1, dz=-1.0):
    """
    SegmentAnalyser with n_segs segments along the z-axis.

    Uses the raw-vector constructor, so "radius" and "creationTime" are
    available as scalars ("age" and "subType" would need an owning organism).
    """
    nodes = [pb.Vector3d(0.0, 0.0, k * dz) for k in range(n_segs + 1)]
    segments = [pb.Vector2i(k, k + 1) for k in range(n_segs)]
    cts = [float(i) for i in range(n_segs)]
    radii = [radius] * n_segs
    return pb.SegmentAnalyser(nodes, segments, cts, radii)


def make_renderer():
    """
    A renderer holding what render_window() puts on stage: the root actor, the
    axes actor, and the scalar bar as a 2-D overlay.
    """
    actor, scalar_bar = vp.plot_roots(make_linear_sa(), "radius", render=False)
    ren = vtk.vtkRenderer()
    ren.AddActor(actor)
    ren.AddActor(vtk.vtkAxesActor())
    ren.AddViewProp(scalar_bar)
    return ren, scalar_bar


@unittest.skipIf(vtk is None, "vtk is not installed")
class TestVtkPlot(unittest.TestCase):
    """plot_roots and the polydata it builds"""

    def test_segs_to_polydata(self):
        """a SegmentAnalyser is converted to polydata with one cell per segment"""
        sa = make_linear_sa(n_segs=4)
        pd = vp.segs_to_polydata(sa, 1.0, ["radius", "creationTime"])
        self.assertEqual(pd.GetNumberOfPoints(), 5)
        self.assertEqual(pd.GetNumberOfCells(), 4)
        self.assertIsNotNone(pd.GetCellData().GetArray("radius"))

    def test_plot_roots_without_rendering(self):
        """render=False returns the actor and its scalar bar, and opens no window"""
        actor, scalar_bar = vp.plot_roots(make_linear_sa(), "radius", render=False)
        self.assertIsInstance(actor, vtk.vtkActor)
        self.assertIsInstance(scalar_bar, vtk.vtkScalarBarActor)
        self.assertEqual(scalar_bar.GetTitle(), "radius")
        self.assertIsNotNone(actor.GetMapper().GetLookupTable())

    def test_plot_roots_returns_lookup_table(self):
        """returnLut=True adds the lookup table, ranged over the plotted scalar"""
        _, _, lut = vp.plot_roots(make_linear_sa(radius=0.1), "radius", render=False, returnLut=True)
        self.assertIsInstance(lut, vtk.vtkLookupTable)
        mn, mx = lut.GetTableRange()
        self.assertLessEqual(mn, 0.1)
        self.assertGreaterEqual(mx, 0.1)


@unittest.skipIf(vtk is None, "vtk is not installed")
class TestScalarBarOverlay(unittest.TestCase):
    """
    The scalar bar is a 2-D prop that has to be added and removed without
    disturbing the 3-D actors.  VTK 9.7 dropped AddActor2D / RemoveActor2D and
    deprecated GetActors2D, so this is done through the ViewProp API.
    """

    def test_scalar_bar_is_added_as_2d_prop(self):
        """render_window's AddViewProp(scalar_bar) puts it among the 2-D props"""
        ren, scalar_bar = make_renderer()
        props_2d = [p for p in ren.GetViewProps() if isinstance(p, vtk.vtkActor2D)]
        self.assertEqual([id(p) for p in props_2d], [id(scalar_bar)])

    def test_overlay_sweep_keeps_3d_actors(self):
        """vtk_animate.update() removes the 2-D overlays only, not the plant"""
        ren, _ = make_renderer()
        n_actors = ren.GetActors().GetNumberOfItems()
        self.assertGreater(n_actors, 0)
        for a in [p for p in ren.GetViewProps() if isinstance(p, vtk.vtkActor2D)]:
            ren.RemoveViewProp(a)
        self.assertEqual([p for p in ren.GetViewProps() if isinstance(p, vtk.vtkActor2D)], [])
        self.assertEqual(ren.GetActors().GetNumberOfItems(), n_actors)

    def test_no_removed_or_deprecated_vtk_api(self):
        """
        Guard against the VTK 9.7 removals creeping back in: AddActor2D and
        RemoveActor2D no longer exist, GetActors2D is deprecated.  Their
        replacements are Add/RemoveViewProp and GetViewProps (the latter also
        returns the 3-D props, so it needs a vtkActor2D filter).
        """
        for module in (vp, vtk_animate):
            source = inspect.getsource(module)
            for removed in ("AddActor2D", "RemoveActor2D", "GetActors2D"):
                with self.subTest(module=module.__name__, api=removed):
                    self.assertNotIn(removed + "(", source)


@unittest.skipIf(vtk is None, "vtk is not installed")
@unittest.skipUnless(HAS_OFFSCREEN, "no offscreen OpenGL context available")
class TestWriteJpg(unittest.TestCase):
    """write_jpg renders offscreen, so it needs an OpenGL context"""

    def test_write_jpg(self):
        """a jpg is written, and the scalar bar title it pads is restored"""
        ren, scalar_bar = make_renderer()
        renWin = vtk.vtkRenderWindow()
        renWin.SetOffScreenRendering(1)
        renWin.SetSize(300, 300)
        renWin.AddRenderer(ren)
        with tempfile.TemporaryDirectory() as tmp:
            file_name = os.path.join(tmp, "test_write_jpg")
            vp.write_jpg(renWin, file_name, magnification=1)
            self.assertTrue(os.path.exists(file_name + ".jpg"))
            self.assertGreater(os.path.getsize(file_name + ".jpg"), 0)
        self.assertEqual(scalar_bar.GetTitle(), "radius")


if __name__ == "__main__":
    unittest.main()
