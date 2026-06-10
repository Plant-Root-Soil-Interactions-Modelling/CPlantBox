"""
Tests for the GrassLeaf organ class using a full Poaceae plant simulation.

Parameters are chosen to keep runtime short:
  - stem:     lb=1 cm, r=1 cm/day  →  first leaf appears at t ≈ 1 day
  - GrassLeaf: sheathLength=1, bladeLength=2, leafGrowthDuration=5 (LinearGrowth)
               →  fully grown (3 cm) by t ≈ 6 days

setUpClass simulates 8 days so the first leaf is fully grown.

Run from this directory:
    python -m unittest test_grassleaf -v
"""

import sys

sys.path.append("..")
sys.path.append("../src/")
import math
import unittest

import numpy as np

import plantbox as pb

# ---------------------------------------------------------------------------
# Helper: Poaceae plant
# ---------------------------------------------------------------------------


class Poaceae(pb.Plant):
    """Plant subclass that substitutes GrassLeaf organs for the default Leaf."""

    def createLeaf(self, subType, delay, parent, pni):
        return pb.GrassLeaf(self, subType, delay, parent, pni)


def _build_plant(seed=42, sheath=1.0, blade=2.0, lgd=5.0, bending=0.0, blade_angle=0.3):
    """Return an initialised-but-not-yet-simulated Poaceae plant.

    Parameters
    ----------
    sheath, blade : target sheath / blade lengths [cm]
    lgd           : leafGrowthDuration [days]
    bending       : bladeBending [rad/cm]
    blade_angle   : bladeAngle [rad]
    """
    plant = Poaceae()
    plant.setSeed(seed)

    seed_rp = pb.SeedRandomParameter(plant)
    seed_rp.subType = 0
    plant.setOrganRandomParameter(seed_rp)

    root_rp = pb.RootRandomParameter(plant)
    root_rp.subType = 1
    root_rp.lmax = 0.1
    root_rp.r = 0.1
    plant.setOrganRandomParameter(root_rp)

    stem_rp = pb.StemRandomParameter(plant)
    stem_rp.subType = 1
    stem_rp.a = 0.1
    stem_rp.lmax = 20.0
    stem_rp.r = 1.0  # 1 cm/day → lb reached after 1 day
    stem_rp.la = 0.0
    stem_rp.lb = 1.0  # first lateral at 1 cm
    stem_rp.ln = 3.0  # inter-lateral spacing
    stem_rp.theta = 0.0  # upright
    stem_rp.dx = 0.5
    stem_rp.dxMin = 0.25
    stem_rp.successor = [[1]]
    stem_rp.successorP = [[1.0]]
    stem_rp.successorOT = [[pb.leaf]]
    stem_rp.successorNo = [1]
    plant.setOrganRandomParameter(stem_rp)

    gl_rp = pb.GrassLeafRandomParameter(plant)
    gl_rp.subType = 1
    gl_rp.a = 0.02
    gl_rp.bladeAngle = blade_angle
    gl_rp.bladeAngles = 0.0
    gl_rp.bladeLength = blade
    gl_rp.bladeLengths = 0.0
    gl_rp.bladeWidth = 0.8
    gl_rp.bladeWidths = 0.0
    gl_rp.sheathLength = sheath
    gl_rp.sheathLengths = 0.0
    gl_rp.bladeBending = bending
    gl_rp.bladeBendings = 0.0
    gl_rp.leafGrowthDuration = lgd
    gl_rp.leafGrowthDurations = 0.0
    gl_rp.f_gf = pb.LinearGrowth()
    plant.setOrganRandomParameter(gl_rp)

    return plant, gl_rp


# ---------------------------------------------------------------------------
# Test class
# ---------------------------------------------------------------------------


class TestGrassLeaf(unittest.TestCase):
    """Integration tests for GrassLeaf using a full plant simulation."""

    SHEATH = 1.0  # target sheathLength [cm]
    BLADE = 2.0  # target bladeLength  [cm]
    LGD = 5.0  # leafGrowthDuration  [days]

    @classmethod
    def setUpClass(cls):
        """Run the simulation once for the whole test class."""
        cls.plant, cls.gl_rp = _build_plant(
            seed=42,
            sheath=cls.SHEATH,
            blade=cls.BLADE,
            lgd=cls.LGD,
            bending=0.0,  # zero bending for clean geometry tests
            blade_angle=0.3,
        )
        cls.plant.initialize(verbose=False)
        cls.plant.simulate(8.0, verbose=False)  # leaf fully grown by day ≈6
        cls.leaves = cls.plant.getOrgans(pb.OrganTypes.leaf, True)

    # ------------------------------------------------------------------
    # Existence / type
    # ------------------------------------------------------------------

    def test_leaves_were_created(self):
        """At least one leaf organ must exist after 8 days."""
        self.assertGreater(len(self.leaves), 0, "No leaf organs found after simulation")

    def test_leaf_is_grassleaf_instance(self):
        """Every leaf organ returned for ot_leaf must be a pb.GrassLeaf instance
        (confirms that Poaceae.createLeaf override was effective).
        """
        for leaf in self.leaves:
            self.assertIsInstance(leaf, pb.GrassLeaf, "Expected pb.GrassLeaf; got " + type(leaf).__name__)

    # ------------------------------------------------------------------
    # param() / getGrassLeafRandomParameter()
    # ------------------------------------------------------------------

    def test_param_returns_grassleaf_specific_parameter(self):
        """param() must return a GrassLeafSpecificParameter."""
        gl = self.leaves[0]
        self.assertIsInstance(gl.param(), pb.GrassLeafSpecificParameter)

    def test_param_values_match_random_parameter(self):
        """With zero standard deviations realised param values equal the means."""
        gl = self.leaves[0]
        p = gl.param()
        rp = self.gl_rp
        self.assertAlmostEqual(p.sheathLength, rp.sheathLength)
        self.assertAlmostEqual(p.bladeLength, rp.bladeLength)
        self.assertAlmostEqual(p.bladeAngle, rp.bladeAngle)
        self.assertAlmostEqual(p.bladeBending, rp.bladeBending)
        self.assertAlmostEqual(p.leafGrowthDuration, rp.leafGrowthDuration)

    def test_getGrassLeafRandomParameter(self):
        """getGrassLeafRandomParameter() returns a GrassLeafRandomParameter."""
        gl = self.leaves[0]
        rp = gl.getGrassLeafRandomParameter()
        self.assertIsInstance(rp, pb.GrassLeafRandomParameter)

    # ------------------------------------------------------------------
    # Growth lengths
    # ------------------------------------------------------------------

    def test_total_length_positive(self):
        """getSheathLength() + getBladeLength() must be > 0 after growth."""
        gl = self.leaves[0]
        self.assertGreater(gl.getSheathLength() + gl.getBladeLength(), 0.0)

    def test_fully_grown_total_length(self):
        """Oldest leaf should be near full length (sheathLength + bladeLength)
        because leafGrowthDuration=5 days elapsed before t=8 days.
        """
        gl = self.leaves[0]
        p = gl.param()
        expected_total = p.sheathLength + p.bladeLength
        actual_total = gl.getSheathLength() + gl.getBladeLength()
        self.assertAlmostEqual(actual_total, expected_total, delta=0.2, msg="Fully grown leaf total length unexpected")

    def test_sheath_length_at_full_growth(self):
        """At full growth sheathLength must equal the parameter sheathLength."""
        gl = self.leaves[0]
        p = gl.param()
        self.assertAlmostEqual(gl.getSheathLength(), p.sheathLength, delta=0.1)

    def test_blade_length_at_full_growth(self):
        """At full growth bladeLength must equal the parameter bladeLength."""
        gl = self.leaves[0]
        p = gl.param()
        self.assertAlmostEqual(gl.getBladeLength(), p.bladeLength, delta=0.1)

    # ------------------------------------------------------------------
    # calcLength / calcAge (should be mutual inverses)
    # ------------------------------------------------------------------

    def test_calcLength_zero_age_is_zero(self):
        """calcLength(0) must equal 0 (no growth at zero age)."""
        gl = self.leaves[0]
        self.assertAlmostEqual(gl.calcLength(0.0), 0.0, places=9)

    def test_calcLength_at_full_duration(self):
        """calcLength(leafGrowthDuration) must equal sheathLength + bladeLength."""
        gl = self.leaves[0]
        p = gl.param()
        k = p.sheathLength + p.bladeLength
        self.assertAlmostEqual(gl.calcLength(p.leafGrowthDuration), k, delta=1e-6)

    def test_calcLength_calcAge_are_inverse(self):
        """calcLength(calcAge(l)) must recover l for several lengths in (0, k)."""
        gl = self.leaves[0]
        p = gl.param()
        k = p.sheathLength + p.bladeLength
        for frac in (0.05, 0.2, 0.5, 0.75, 0.95):
            l = frac * k
            age_ = gl.calcAge(l)
            l2 = gl.calcLength(age_)
            self.assertAlmostEqual(l2, l, delta=1e-6, msg=f"Round-trip failed at fraction {frac}: " f"calcLength(calcAge({l:.3f})) = {l2:.6f}")

    # ------------------------------------------------------------------
    # Node geometry
    # ------------------------------------------------------------------

    def test_node_count_greater_than_one(self):
        """A fully grown leaf must have more than one node
        (anchor + at least one turtle node).
        """
        gl = self.leaves[0]
        self.assertGreater(gl.getNumberOfNodes(), 1)

    def test_getNumberOfNodes_matches_node_access(self):
        """getNode(i) must succeed for every index 0..getNumberOfNodes()-1."""
        gl = self.leaves[0]
        n = gl.getNumberOfNodes()
        for i in range(n):
            pos = gl.getNode(i)  # must not raise
            self.assertIsInstance(pos, pb.Vector3d, f"getNode({i}) did not return a Vector3d")

    def test_getNode0_equals_parent_attachment(self):
        """getNode(0) must be exactly the parent node at parentNI
        (it is the anchor / leaf attachment point on the stem).
        """
        gl = self.leaves[0]
        parent = gl.getParent()
        pni = gl.getParentNI()
        parent_node = parent.getNode(pni)
        leaf_node0 = gl.getNode(0)
        self.assertAlmostEqual(leaf_node0.x, parent_node.x, places=9)
        self.assertAlmostEqual(leaf_node0.y, parent_node.y, places=9)
        self.assertAlmostEqual(leaf_node0.z, parent_node.z, places=9)

    def test_node_positions_differ(self):
        """Successive nodes must not all be at the same position
        (the leaf must actually have extended in space).
        """
        gl = self.leaves[0]
        n = gl.getNumberOfNodes()
        nodes = [gl.getNode(i) for i in range(n)]
        # Check that at least the first and last node differ
        p0 = nodes[0]
        pN = nodes[-1]
        dist = math.sqrt((pN.x - p0.x) ** 2 + (pN.y - p0.y) ** 2 + (pN.z - p0.z) ** 2)
        self.assertGreater(dist, 0.0, "First and last node are at the same position")

    # ------------------------------------------------------------------
    # getParameter()
    # ------------------------------------------------------------------

    def test_getParameter_sheathLength(self):
        """getParameter('sheathLength') returns the realised sheathLength."""
        gl = self.leaves[0]
        p = gl.param()
        self.assertAlmostEqual(gl.getParameter("sheathLength"), p.sheathLength)

    def test_getParameter_bladeLength(self):
        """getParameter('bladeLength') returns the realised bladeLength."""
        gl = self.leaves[0]
        p = gl.param()
        self.assertAlmostEqual(gl.getParameter("bladeLength"), p.bladeLength)

    def test_getParameter_bladeAngle(self):
        """getParameter('bladeAngle') returns the realised bladeAngle."""
        gl = self.leaves[0]
        p = gl.param()
        self.assertAlmostEqual(gl.getParameter("bladeAngle"), p.bladeAngle)

    def test_getParameter_bladeWidth(self):
        """getParameter('bladeWidth') returns the realised bladeWidth."""
        gl = self.leaves[0]
        p = gl.param()
        self.assertAlmostEqual(gl.getParameter("bladeWidth"), p.bladeWidth)

    def test_getParameter_leafGrowthDuration(self):
        """getParameter('leafGrowthDuration') returns the realised value."""
        gl = self.leaves[0]
        p = gl.param()
        self.assertAlmostEqual(gl.getParameter("leafGrowthDuration"), p.leafGrowthDuration)

    # ------------------------------------------------------------------
    # getLeafVis()
    # ------------------------------------------------------------------

    def test_getLeafVis_returns_two_points(self):
        """getLeafVis(i) must return exactly 2 Vector3d edge points for any
        valid node index.
        """
        gl = self.leaves[0]
        n = gl.getNumberOfNodes()
        for i in (0, n // 2, n - 1):
            vis = gl.getLeafVis(i)
            self.assertEqual(len(vis), 2, f"getLeafVis({i}) returned {len(vis)} points, expected 2")
            for pt in vis:
                self.assertIsInstance(pt, pb.Vector3d)

    # ------------------------------------------------------------------
    # toString()
    # ------------------------------------------------------------------

    def test_toString_contains_grassleaf(self):
        """toString() must produce a string containing 'GrassLeaf'."""
        gl = self.leaves[0]
        s = str(gl)
        self.assertIn("GrassLeaf", s)

    # ------------------------------------------------------------------
    # Multiple leaves / growth state
    # ------------------------------------------------------------------

    def test_inactive_leaf_fully_grown(self):
        """A leaf whose age >= leafGrowthDuration should be inactive and at
        maximum length.
        """
        gl = self.leaves[0]
        p = gl.param()
        if gl.getAge() >= p.leafGrowthDuration:
            self.assertFalse(gl.isActive(), "Fully grown leaf should be inactive")
            total = gl.getSheathLength() + gl.getBladeLength()
            expected = p.sheathLength + p.bladeLength
            self.assertAlmostEqual(total, expected, delta=0.05)

    def test_younger_leaves_still_growing(self):
        """If multiple leaves exist, later-emerging ones should be shorter
        than the first (oldest) leaf.
        """
        if len(self.leaves) < 2:
            self.skipTest("Only one leaf produced; skipping multi-leaf test")
        gl0 = self.leaves[0]
        gl1 = self.leaves[1]
        total0 = gl0.getSheathLength() + gl0.getBladeLength()
        total1 = gl1.getSheathLength() + gl1.getBladeLength()
        self.assertGreaterEqual(total0, total1, "First leaf should be longer than or equal to second")


if __name__ == "__main__":
    unittest.main(verbosity=2)
