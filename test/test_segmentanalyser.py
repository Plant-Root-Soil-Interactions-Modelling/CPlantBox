"""
Unit tests for SegmentAnalyser.

Most tests use the raw-vector constructor (nodes, segments, CTs, radii) so
they are self-contained and fast.  The test_from_organism test requires a
modelparameter XML file and is skipped when it is not found.
"""

import sys

sys.path.append("..")
sys.path.append("../src/")

import math
import os
import tempfile
import unittest

import plantbox as pb

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

MODEL_FILE = "../modelparameter/structural/rootsystem/Anagallis_femina_Leitner_2010.xml"


def make_linear_sa(n_segs=4, radius=0.1, dz=-1.0):
    """
    Build a SegmentAnalyser with n_segs segments along the z-axis.

    Nodes: (0,0,0), (0,0,dz), (0,0,2*dz), ...
    Segment i runs from node i to node i+1.
    creationTime[i] = float(i); radius is constant.
    """
    nodes = [pb.Vector3d(0.0, 0.0, k * dz) for k in range(n_segs + 1)]
    segments = [pb.Vector2i(k, k + 1) for k in range(n_segs)]
    cts = [float(i) for i in range(n_segs)]
    radii = [radius] * n_segs
    return pb.SegmentAnalyser(nodes, segments, cts, radii)


def make_single_segment_sa(x1=0.0, y1=0.0, z1=0.0, x2=1.0, y2=0.0, z2=0.0, ct=0.0, radius=0.1):
    """Build a SegmentAnalyser containing exactly one segment."""
    nodes = [pb.Vector3d(x1, y1, z1), pb.Vector3d(x2, y2, z2)]
    segments = [pb.Vector2i(0, 1)]
    return pb.SegmentAnalyser(nodes, segments, [ct], [radius])


# ---------------------------------------------------------------------------
# Tests: constructors
# ---------------------------------------------------------------------------


class TestConstructors(unittest.TestCase):

    def test_empty_constructor(self):
        sa = pb.SegmentAnalyser()
        self.assertEqual(len(sa.nodes), 0)
        self.assertEqual(len(sa.segments), 0)

    def test_raw_constructor_node_count(self):
        sa = make_linear_sa(4)
        self.assertEqual(len(sa.nodes), 5)

    def test_raw_constructor_segment_count(self):
        sa = make_linear_sa(4)
        self.assertEqual(len(sa.segments), 4)

    def test_raw_constructor_creation_time(self):
        sa = make_linear_sa(3)
        cts = sa.getParameter("creationTime")
        self.assertEqual(len(cts), 3)
        for i, ct in enumerate(cts):
            self.assertAlmostEqual(ct, float(i))

    def test_raw_constructor_radius(self):
        sa = make_linear_sa(3, radius=0.05)
        for r in sa.getParameter("radius"):
            self.assertAlmostEqual(r, 0.05)

    def test_copy_constructor(self):
        sa = make_linear_sa(3)
        sa2 = pb.SegmentAnalyser(sa)
        self.assertEqual(len(sa2.segments), len(sa.segments))
        self.assertEqual(len(sa2.nodes), len(sa.nodes))

    @unittest.skipUnless(os.path.isfile(MODEL_FILE), "modelparameter file not found")
    def test_from_organism(self):
        rs = pb.RootSystem()
        rs.readParameters(MODEL_FILE, verbose=False)
        rs.initialize(False)
        rs.simulate(10.0)
        sa = pb.SegmentAnalyser(rs)
        self.assertGreater(len(sa.segments), 0)
        self.assertEqual(len(sa.segments), len(sa.getParameter("radius")))
        self.assertEqual(len(sa.segments), len(sa.getParameter("organType")))


# ---------------------------------------------------------------------------
# Tests: addSegments / addSegment
# ---------------------------------------------------------------------------


class TestAddSegments(unittest.TestCase):

    def test_addSegments_from_analyser_count(self):
        sa1 = make_linear_sa(3)
        sa2 = make_linear_sa(2)
        sa1.addSegments(sa2)
        self.assertEqual(len(sa1.segments), 5)

    def test_addSegments_node_index_offset(self):
        """Segment indices appended from sa2 must be offset by the original node count."""
        sa1 = make_linear_sa(1)  # nodes 0,1 → seg (0,1)
        n_nodes_before = len(sa1.nodes)
        sa2 = make_linear_sa(1)  # nodes 0,1 → seg (0,1)
        sa1.addSegments(sa2)
        last = sa1.segments[-1]
        self.assertGreaterEqual(last.x, n_nodes_before)
        self.assertGreaterEqual(last.y, n_nodes_before)

    def test_addSegment_append_increases_count(self):
        sa = make_linear_sa(2)
        n_before = len(sa.segments)
        sa.addSegment(pb.Vector2i(0, 1), 99.0, 0.2, False)
        self.assertEqual(len(sa.segments), n_before + 1)

    def test_addSegment_insert_at_front(self):
        sa = make_linear_sa(3)
        sa.addSegment(pb.Vector2i(0, 1), 99.0, 0.2, True)
        ct = sa.getParameter("creationTime")
        self.assertAlmostEqual(ct[0], 99.0)

    def test_addSegment_data_arrays_consistent(self):
        """All data arrays must stay the same length as segments after addSegment."""
        sa = make_linear_sa(3)
        sa.addData("custom", [1.0, 2.0, 3.0])
        sa.addSegment(pb.Vector2i(0, 1), 5.0, 0.1)
        for key, vals in sa.data.items():
            self.assertEqual(len(vals), len(sa.segments), f"data['{key}'] length mismatch after addSegment")


# ---------------------------------------------------------------------------
# Tests: getParameter
# ---------------------------------------------------------------------------


class TestGetParameter(unittest.TestCase):

    def test_creationTime(self):
        sa = make_linear_sa(4)
        cts = sa.getParameter("creationTime")
        for i, ct in enumerate(cts):
            self.assertAlmostEqual(ct, float(i))

    def test_radius(self):
        sa = make_linear_sa(4, radius=0.07)
        for r in sa.getParameter("radius"):
            self.assertAlmostEqual(r, 0.07)

    def test_length_unit_segments(self):
        sa = make_linear_sa(3, dz=-1.0)
        for l in sa.getParameter("length"):
            self.assertAlmostEqual(l, 1.0, places=10)

    def test_length_diagonal(self):
        sa = make_single_segment_sa(0.0, 0.0, 0.0, 1.0, 1.0, 0.0)
        self.assertAlmostEqual(sa.getParameter("length")[0], math.sqrt(2.0), places=10)

    def test_surface(self):
        r = 0.15
        sa = make_linear_sa(1, radius=r, dz=-1.0)
        expected = 2.0 * math.pi * r * 1.0
        self.assertAlmostEqual(sa.getParameter("surface")[0], expected, places=10)

    def test_volume(self):
        r = 0.2
        sa = make_linear_sa(1, radius=r, dz=-1.0)
        expected = math.pi * r * r * 1.0
        self.assertAlmostEqual(sa.getParameter("volume")[0], expected, places=10)

    def test_custom_data(self):
        sa = make_linear_sa(3)
        custom = [10.0, 20.0, 30.0]
        sa.addData("myval", custom)
        for exp, got in zip(custom, sa.getParameter("myval")):
            self.assertAlmostEqual(exp, got)


# ---------------------------------------------------------------------------
# Tests: getSegmentLength
# ---------------------------------------------------------------------------


class TestGetSegmentLength(unittest.TestCase):

    def test_unit_segment(self):
        sa = make_linear_sa(2, dz=-1.0)
        self.assertAlmostEqual(sa.getSegmentLength(0), 1.0, places=12)

    def test_3_4_5_triangle(self):
        sa = make_single_segment_sa(0.0, 0.0, 0.0, 3.0, 4.0, 0.0)
        self.assertAlmostEqual(sa.getSegmentLength(0), 5.0, places=12)


# ---------------------------------------------------------------------------
# Tests: getSummed
# ---------------------------------------------------------------------------


class TestGetSummed(unittest.TestCase):

    def test_total_length(self):
        n = 5
        sa = make_linear_sa(n, dz=-1.0)
        self.assertAlmostEqual(sa.getSummed("length"), float(n), places=10)

    def test_summed_after_filter(self):
        sa = make_linear_sa(4, dz=-1.0)  # CTs 0,1,2,3 → keep CT 0 and 1 → 2 segments of length 1
        sa.filter("creationTime", 0.0, 1.0)
        self.assertAlmostEqual(sa.getSummed("length"), 2.0, places=10)

    def test_summed_with_geometry(self):
        """Midpoint-based sum within a box of depth 2 should equal 2."""
        sa = make_linear_sa(4, dz=-1.0)
        geom = pb.SDF_PlantBox(100.0, 100.0, 2.0)
        # Segments 0 and 1 have midpoints at z=-0.5 and z=-1.5, both inside
        result = sa.getSummed("length", geom)
        self.assertAlmostEqual(result, 2.0, places=10)


# ---------------------------------------------------------------------------
# Tests: filter
# ---------------------------------------------------------------------------


class TestFilter(unittest.TestCase):

    def test_filter_range_count(self):
        sa = make_linear_sa(6)  # CTs 0..5
        sa.filter("creationTime", 1.0, 3.0)
        self.assertEqual(len(sa.segments), 3)

    def test_filter_range_values(self):
        sa = make_linear_sa(5)
        sa.filter("creationTime", 0.0, 2.0)
        for ct in sa.getParameter("creationTime"):
            self.assertGreaterEqual(ct, 0.0)
            self.assertLessEqual(ct, 2.0)

    def test_filter_exact_value_count(self):
        sa = make_linear_sa(5)
        sa.filter("creationTime", 2.0)
        self.assertEqual(len(sa.segments), 1)

    def test_filter_exact_value_content(self):
        sa = make_linear_sa(5)
        sa.filter("creationTime", 2.0)
        self.assertAlmostEqual(sa.getParameter("creationTime")[0], 2.0)

    def test_filter_removes_all(self):
        sa = make_linear_sa(3)
        sa.filter("creationTime", 100.0, 200.0)
        self.assertEqual(len(sa.segments), 0)

    def test_filter_data_arrays_consistent(self):
        sa = make_linear_sa(6)
        sa.filter("creationTime", 2.0, 4.0)
        n = len(sa.segments)
        for key, vals in sa.data.items():
            self.assertEqual(len(vals), n, f"data['{key}'] length mismatch after filter")


# ---------------------------------------------------------------------------
# Tests: crop / cropDomain
# ---------------------------------------------------------------------------


class TestCrop(unittest.TestCase):

    def test_cropDomain_total_length(self):
        """After cropping to depth 2, total length must be exactly 2."""
        sa = make_linear_sa(5, dz=-1.0)  # 5 segments, total length 5 cm
        sa.cropDomain(100.0, 100.0, 2.0)
        self.assertAlmostEqual(sa.getSummed("length"), 2.0, places=5)

    def test_cropDomain_reduces_segments(self):
        sa = make_linear_sa(5, dz=-1.0)
        sa.cropDomain(100.0, 100.0, 2.0)
        self.assertLess(len(sa.segments), 5)
        self.assertGreater(len(sa.segments), 0)

    def test_crop_with_sdf(self):
        sa = make_linear_sa(4, dz=-1.0)
        sa.crop(pb.SDF_PlantBox(100.0, 100.0, 1.5))
        self.assertAlmostEqual(sa.getSummed("length"), 1.5, places=5)

    def test_crop_data_arrays_consistent(self):
        sa = make_linear_sa(4, dz=-1.0)
        sa.cropDomain(100.0, 100.0, 2.0)
        n = len(sa.segments)
        for key, vals in sa.data.items():
            self.assertEqual(len(vals), n, f"data['{key}'] length mismatch after crop")


# ---------------------------------------------------------------------------
# Tests: pack
# ---------------------------------------------------------------------------


class TestPack(unittest.TestCase):

    def test_pack_reduces_nodes_after_filter(self):
        sa = make_linear_sa(8, dz=-1.0)
        n_nodes_before = len(sa.nodes)
        sa.filter("creationTime", 0.0, 2.0)  # 3 of 8 segments remain
        sa.pack()
        self.assertLess(len(sa.nodes), n_nodes_before)

    def test_pack_segment_indices_valid(self):
        sa = make_linear_sa(5)
        sa.filter("creationTime", 0.0, 1.0)
        sa.pack()
        n = len(sa.nodes)
        for s in sa.segments:
            self.assertGreaterEqual(s.x, 0)
            self.assertGreater(n, s.x)
            self.assertGreaterEqual(s.y, 0)
            self.assertGreater(n, s.y)

    def test_pack_preserves_length(self):
        sa = make_linear_sa(5, dz=-1.0)
        sa.filter("creationTime", 0.0, 2.0)
        total_before = sa.getSummed("length")
        sa.pack()
        self.assertAlmostEqual(sa.getSummed("length"), total_before, places=10)


# ---------------------------------------------------------------------------
# Tests: bounds
# ---------------------------------------------------------------------------


class TestBounds(unittest.TestCase):

    def test_min_bounds_z(self):
        sa = make_linear_sa(4, dz=-1.0)  # deepest node at z=-4
        mn = sa.getMinBounds()
        self.assertAlmostEqual(mn.z, -4.0)

    def test_max_bounds_z(self):
        sa = make_linear_sa(4, dz=-1.0)  # topmost node at z=0
        mx = sa.getMaxBounds()
        self.assertAlmostEqual(mx.z, 0.0)

    def test_min_bounds_x_y(self):
        sa = make_linear_sa(3, dz=-1.0)  # all nodes at x=y=0
        mn = sa.getMinBounds()
        self.assertAlmostEqual(mn.x, 0.0)
        self.assertAlmostEqual(mn.y, 0.0)


# ---------------------------------------------------------------------------
# Tests: addAge
# ---------------------------------------------------------------------------


class TestAddAge(unittest.TestCase):

    def test_age_equals_simtime_minus_ct(self):
        sa = make_linear_sa(4)  # CTs: 0,1,2,3
        sim_time = 5.0
        sa.addAge(sim_time)
        for age, ct in zip(sa.getParameter("age"), sa.getParameter("creationTime")):
            self.assertAlmostEqual(age, sim_time - ct, places=12)

    def test_age_clamped_to_zero(self):
        """Segments created after simTime should get age 0, not a negative value."""
        sa = make_linear_sa(3)  # CTs: 0,1,2
        sa.addAge(0.5)  # only segment with CT=0 is non-zero aged
        for age in sa.getParameter("age"):
            self.assertGreaterEqual(age, 0.0)


# ---------------------------------------------------------------------------
# Tests: addData
# ---------------------------------------------------------------------------


class TestAddData(unittest.TestCase):

    def test_segment_sized_data(self):
        sa = make_linear_sa(3)
        sa.addData("field", [1.0, 2.0, 3.0])
        self.assertEqual(list(sa.getParameter("field")), [1.0, 2.0, 3.0])

    def test_node_sized_data_converted(self):
        """Node-sized input is converted to segment data using the tip-node index."""
        sa = make_linear_sa(3)  # 4 nodes, 3 segments
        node_vals = [0.0, 1.0, 2.0, 3.0]
        sa.addData("node_field", node_vals)
        vals = sa.getParameter("node_field")
        self.assertEqual(len(vals), 3)
        for i, v in enumerate(vals):
            tip_idx = sa.segments[i].y
            self.assertAlmostEqual(v, node_vals[tip_idx])

    def test_wrong_size_raises(self):
        sa = make_linear_sa(3)
        with self.assertRaises(Exception):
            sa.addData("bad", [1.0, 2.0])  # neither n_segs (3) nor n_nodes (4)

    def test_overwrite_existing_field(self):
        sa = make_linear_sa(2)
        sa.addData("radius", [0.5, 0.5])
        for r in sa.getParameter("radius"):
            self.assertAlmostEqual(r, 0.5)


# ---------------------------------------------------------------------------
# Tests: distribution
# ---------------------------------------------------------------------------


class TestDistribution(unittest.TestCase):

    def test_scalar_layer_count(self):
        sa = make_linear_sa(4, dz=-1.0)
        d = sa.distribution("length", 0.0, -4.0, 4, False)
        self.assertEqual(len(d), 4)

    def test_scalar_sum_equals_total(self):
        sa = make_linear_sa(4, dz=-1.0)
        d = sa.distribution("length", 0.0, -4.0, 4, False)
        self.assertAlmostEqual(sum(d), 4.0, places=5)

    def test_scalar_midpoint_per_layer(self):
        """Midpoint-based (exact=False): each 1-cm layer should contain exactly 1 cm."""
        sa = make_linear_sa(4, dz=-1.0)
        d = sa.distribution("length", 0.0, -4.0, 4, False)
        for val in d:
            self.assertAlmostEqual(val, 1.0, places=5)

    def test_objects_layer_count(self):
        sa = make_linear_sa(4, dz=-1.0)
        layers = sa.distribution(0.0, -4.0, 4)
        self.assertEqual(len(layers), 4)

    def test_objects_type(self):
        sa = make_linear_sa(4, dz=-1.0)
        layers = sa.distribution(0.0, -4.0, 2)
        for layer in layers:
            self.assertIsInstance(layer, pb.SegmentAnalyser)

    def test_distribution2_shape(self):
        sa = make_linear_sa(4, dz=-1.0)
        d = sa.distribution2("length", 0.0, -4.0, -10.0, 10.0, 4, 2, False)
        self.assertEqual(len(d), 4)
        for row in d:
            self.assertEqual(len(row), 2)


# ---------------------------------------------------------------------------
# Tests: mapPeriodic
# ---------------------------------------------------------------------------


class TestMapPeriodic(unittest.TestCase):

    def _sa_crossing_x(self):
        """One segment from x=0.3 to x=0.8 that crosses the x=0.5 boundary."""
        return make_single_segment_sa(0.3, 0.0, 0.0, 0.8, 0.0, 0.0)

    def test_crossing_segment_is_split(self):
        sa = self._sa_crossing_x()
        sa.mapPeriodic(1.0, 1.0)
        self.assertEqual(len(sa.segments), 2)

    def test_total_length_conserved(self):
        sa = self._sa_crossing_x()
        length_before = sa.getSummed("length")
        sa.mapPeriodic(1.0, 1.0)
        self.assertAlmostEqual(sa.getSummed("length"), length_before, places=5)

    def test_nodes_in_periodic_range(self):
        """After mapping all x-coordinates must lie in [-xx/2, xx/2)."""
        sa = self._sa_crossing_x()
        sa.mapPeriodic(1.0, 1.0)
        for n in sa.nodes:
            self.assertGreaterEqual(n.x, -0.5 - 1e-9)
            self.assertLess(n.x, 0.5 + 1e-9)

    def test_non_crossing_segment_unchanged(self):
        """A segment wholly inside the domain should not be split."""
        sa = make_single_segment_sa(0.1, 0.0, 0.0, 0.4, 0.0, 0.0)
        sa.mapPeriodic(1.0, 1.0)
        self.assertEqual(len(sa.segments), 1)


# ---------------------------------------------------------------------------
# Tests: write (integration – file creation and basic format checks)
# ---------------------------------------------------------------------------


class TestWrite(unittest.TestCase):

    def test_write_vtp_creates_file(self):
        sa = make_linear_sa(3)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.vtp")
            sa.write(path)
            self.assertTrue(os.path.isfile(path))

    def test_write_vtp_is_valid_xml(self):
        sa = make_linear_sa(3)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.vtp")
            sa.write(path)
            with open(path) as f:
                content = f.read()
            self.assertIn("VTKFile", content)
            self.assertIn("PolyData", content)

    def test_write_dgf_creates_file(self):
        sa = make_linear_sa(3)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.dgf")
            sa.write(path)
            self.assertTrue(os.path.isfile(path))

    def test_write_dgf_starts_with_dgf(self):
        sa = make_linear_sa(3)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.dgf")
            sa.write(path)
            with open(path) as f:
                self.assertEqual(f.readline().strip(), "DGF")

    def test_write_vtp_custom_types(self):
        """Writing with a custom field should include it in the output."""
        sa = make_linear_sa(3)
        sa.addAge(5.0)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "out.vtp")
            sa.write(path, ["radius", "age"])
            with open(path) as f:
                content = f.read()
            self.assertIn("age", content)
            self.assertIn("radius", content)


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main(verbosity=2)
