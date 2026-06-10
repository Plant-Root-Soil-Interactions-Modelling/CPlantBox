"""
Tests for GrassLeafSpecificParameter and GrassLeafRandomParameter
exposed via the plantbox Python binding.

Run from this directory:
    python -m unittest test_grassleafparameter -v
"""

import sys

sys.path.append("..")
sys.path.append("../src/")
import unittest

import plantbox as pb


class TestGrassLeafSpecificParameter(unittest.TestCase):
    """Unit tests for GrassLeafSpecificParameter (per-instance realised parameters)."""

    def test_default_constructor(self):
        """Default-constructed instance has expected zero / C++ in-class defaults."""
        p = pb.GrassLeafSpecificParameter()
        self.assertEqual(p.subType, -1)
        self.assertAlmostEqual(p.a, 0.0)
        self.assertAlmostEqual(p.bladeAngle, 0.0)
        self.assertAlmostEqual(p.bladeWidth, 0.0)
        self.assertAlmostEqual(p.bladeLength, 0.0)
        self.assertAlmostEqual(p.sheathLength, 0.0)
        self.assertAlmostEqual(p.leafGrowthDuration, 0.0)
        self.assertAlmostEqual(p.bladeBending, 0.05)  # C++ in-class default

    def test_value_constructor(self):
        """Full constructor (subType, a, bladeAngle, bladeWidth, bladeLength,
        sheathLength, leafGrowthDuration, bladeBending) sets all fields.
        """
        p = pb.GrassLeafSpecificParameter(2, 0.03, 0.4, 0.9, 11.0, 5.5, 18.0, 0.08)
        self.assertEqual(p.subType, 2)
        self.assertAlmostEqual(p.a, 0.03)
        self.assertAlmostEqual(p.bladeAngle, 0.4)
        self.assertAlmostEqual(p.bladeWidth, 0.9)
        self.assertAlmostEqual(p.bladeLength, 11.0)
        self.assertAlmostEqual(p.sheathLength, 5.5)
        self.assertAlmostEqual(p.leafGrowthDuration, 18.0)
        self.assertAlmostEqual(p.bladeBending, 0.08)

    def test_fields_are_readwrite(self):
        """All exposed fields can be written and read back."""
        p = pb.GrassLeafSpecificParameter()
        p.bladeAngle = 1.1
        p.bladeWidth = 2.2
        p.bladeLength = 3.3
        p.sheathLength = 4.4
        p.leafGrowthDuration = 5.5
        p.bladeBending = 0.07
        self.assertAlmostEqual(p.bladeAngle, 1.1)
        self.assertAlmostEqual(p.bladeWidth, 2.2)
        self.assertAlmostEqual(p.bladeLength, 3.3)
        self.assertAlmostEqual(p.sheathLength, 4.4)
        self.assertAlmostEqual(p.leafGrowthDuration, 5.5)
        self.assertAlmostEqual(p.bladeBending, 0.07)

    def test_toString_contains_key_fields(self):
        """toString() output mentions all important parameter names."""
        p = pb.GrassLeafSpecificParameter(1, 0.02, 0.5, 0.8, 12.0, 6.0, 20.0, 0.1)
        s = str(p)
        for field in ("GrassLeaf", "bladeAngle", "bladeWidth", "bladeLength", "sheathLength", "leafGrowthDuration", "bladeBending"):
            self.assertIn(field, s, f"toString: field '{field}' not found in:\n{s}")


class TestGrassLeafRandomParameter(unittest.TestCase):
    """Unit tests for GrassLeafRandomParameter (stochastic type parameters)."""

    def _make(self):
        """Return (rp, org) with a fresh GrassLeafRandomParameter, subType=1.

        The Organism *org* must be kept alive so that the internal weak_ptr
        held by GrassLeafRandomParameter does not dangle.
        """
        org = pb.Organism()
        rp = pb.GrassLeafRandomParameter(org)
        rp.subType = 1
        return rp, org

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def test_default_values(self):
        """Constructor produces the default parameter values stated in the C++ header."""
        rp, _ = self._make()
        self.assertAlmostEqual(rp.bladeAngle, 0.3)
        self.assertAlmostEqual(rp.bladeAngles, 0.0)
        self.assertAlmostEqual(rp.bladeWidth, 0.5)
        self.assertAlmostEqual(rp.bladeWidths, 0.0)
        self.assertAlmostEqual(rp.bladeLength, 10.0)
        self.assertAlmostEqual(rp.bladeLengths, 0.0)
        self.assertAlmostEqual(rp.sheathLength, 5.0)
        self.assertAlmostEqual(rp.sheathLengths, 0.0)
        self.assertAlmostEqual(rp.bladeBending, 0.05)
        self.assertAlmostEqual(rp.bladeBendings, 0.0)
        self.assertAlmostEqual(rp.leafGrowthDuration, 20.0)
        self.assertAlmostEqual(rp.leafGrowthDurations, 0.0)

    def test_organ_type_is_leaf(self):
        """organType must equal ot_leaf (4)."""
        rp, _ = self._make()
        self.assertEqual(rp.organType, int(pb.OrganTypes.leaf))

    # ------------------------------------------------------------------
    # copy()
    # ------------------------------------------------------------------

    def test_copy_is_independent(self):
        """copy() returns a distinct object; mutations to the original do not
        propagate to the copy."""
        org = pb.Organism()
        rp = pb.GrassLeafRandomParameter(org)
        rp.bladeLength = 7.0
        rp2 = rp.copy(org)
        self.assertIsNot(rp, rp2)
        self.assertAlmostEqual(rp2.bladeLength, 7.0, msg="copy: value not carried over")
        rp.bladeLength = 99.0
        self.assertAlmostEqual(rp2.bladeLength, 7.0, msg="copy: mutation leaked into copy")

    def test_copy_preserves_all_grass_fields(self):
        """copy() transfers all GrassLeaf-specific parameter values."""
        org = pb.Organism()
        rp = pb.GrassLeafRandomParameter(org)
        rp.bladeAngle = 0.6
        rp.bladeAngles = 0.01
        rp.bladeWidth = 1.2
        rp.bladeWidths = 0.02
        rp.bladeLength = 8.0
        rp.bladeLengths = 0.3
        rp.sheathLength = 4.0
        rp.sheathLengths = 0.2
        rp.bladeBending = 0.03
        rp.bladeBendings = 0.001
        rp.leafGrowthDuration = 15.0
        rp.leafGrowthDurations = 0.5
        rp2 = rp.copy(org)
        self.assertAlmostEqual(rp2.bladeAngle, rp.bladeAngle)
        self.assertAlmostEqual(rp2.bladeAngles, rp.bladeAngles)
        self.assertAlmostEqual(rp2.bladeWidth, rp.bladeWidth)
        self.assertAlmostEqual(rp2.bladeWidths, rp.bladeWidths)
        self.assertAlmostEqual(rp2.bladeLength, rp.bladeLength)
        self.assertAlmostEqual(rp2.bladeLengths, rp.bladeLengths)
        self.assertAlmostEqual(rp2.sheathLength, rp.sheathLength)
        self.assertAlmostEqual(rp2.sheathLengths, rp.sheathLengths)
        self.assertAlmostEqual(rp2.bladeBending, rp.bladeBending)
        self.assertAlmostEqual(rp2.bladeBendings, rp.bladeBendings)
        self.assertAlmostEqual(rp2.leafGrowthDuration, rp.leafGrowthDuration)
        self.assertAlmostEqual(rp2.leafGrowthDurations, rp.leafGrowthDurations)

    # ------------------------------------------------------------------
    # realize()
    # ------------------------------------------------------------------

    def test_realize_returns_grassleaf_specific_parameter(self):
        """realize() produces a GrassLeafSpecificParameter instance."""
        rp, _ = self._make()
        p = rp.realize()
        self.assertIsInstance(p, pb.GrassLeafSpecificParameter)

    def test_realize_subtype_propagated(self):
        """Realised parameter carries the correct subType."""
        rp, _ = self._make()
        rp.subType = 3
        p = rp.realize()
        self.assertEqual(p.subType, 3)

    def test_realize_exact_means_with_zero_std(self):
        """With all standard deviations = 0, realize() returns exactly the
        mean values (no random perturbation can occur).
        """
        org = pb.Organism()
        rp = pb.GrassLeafRandomParameter(org)
        rp.subType = 1
        rp.bladeAngle = 0.4
        rp.bladeAngles = 0.0
        rp.bladeWidth = 0.9
        rp.bladeWidths = 0.0
        rp.bladeLength = 11.0
        rp.bladeLengths = 0.0
        rp.sheathLength = 5.5
        rp.sheathLengths = 0.0
        rp.bladeBending = 0.07
        rp.bladeBendings = 0.0
        rp.leafGrowthDuration = 18.0
        rp.leafGrowthDurations = 0.0
        p = rp.realize()
        self.assertAlmostEqual(p.bladeAngle, 0.4)
        self.assertAlmostEqual(p.bladeWidth, 0.9)
        self.assertAlmostEqual(p.bladeLength, 11.0)
        self.assertAlmostEqual(p.sheathLength, 5.5)
        self.assertAlmostEqual(p.bladeBending, 0.07)
        self.assertAlmostEqual(p.leafGrowthDuration, 18.0)

    def test_realize_values_nonneg_with_large_std(self):
        """Realised values must be >= 0 even when the standard deviation greatly
        exceeds the mean (the C++ implementation clamps with std::max(..., 0.0)).
        """
        org = pb.Organism()
        rp = pb.GrassLeafRandomParameter(org)
        rp.subType = 1
        rp.bladeLength = 1.0
        rp.bladeLengths = 20.0
        rp.sheathLength = 1.0
        rp.sheathLengths = 20.0
        rp.leafGrowthDuration = 1.0
        rp.leafGrowthDurations = 20.0
        for _ in range(30):
            p = rp.realize()
            self.assertGreaterEqual(p.bladeLength, 0.0, "realized bladeLength must be >= 0")
            self.assertGreaterEqual(p.sheathLength, 0.0, "realized sheathLength must be >= 0")
            self.assertGreaterEqual(p.leafGrowthDuration, 0.0, "realized leafGrowthDuration must be >= 0")

    # ------------------------------------------------------------------
    # toString / getParameter
    # ------------------------------------------------------------------

    def test_toString_works(self):
        """toString() produces a non-empty string without raising."""
        rp, _ = self._make()
        s = str(rp)
        self.assertIsInstance(s, str)
        self.assertGreater(len(s), 0)

    def test_getParameter_grass_specific_names(self):
        """getParameter() exposes all GrassLeaf-specific scalar parameters."""
        org = pb.Organism()
        rp = pb.GrassLeafRandomParameter(org)
        rp.subType = 1
        rp.bladeAngle = 0.35
        rp.bladeWidth = 0.6
        rp.bladeLength = 9.0
        rp.sheathLength = 4.5
        rp.bladeBending = 0.06
        rp.leafGrowthDuration = 17.0
        self.assertAlmostEqual(rp.getParameter("bladeAngle"), 0.35)
        self.assertAlmostEqual(rp.getParameter("bladeWidth"), 0.6)
        self.assertAlmostEqual(rp.getParameter("bladeLength"), 9.0)
        self.assertAlmostEqual(rp.getParameter("sheathLength"), 4.5)
        self.assertAlmostEqual(rp.getParameter("bladeBending"), 0.06)
        self.assertAlmostEqual(rp.getParameter("leafGrowthDuration"), 17.0)

    def test_getParameter_dev_suffix(self):
        """Standard deviations are accessible with the '_dev' suffix via getParameter()."""
        org = pb.Organism()
        rp = pb.GrassLeafRandomParameter(org)
        rp.bladeAngles = 0.01
        rp.bladeLengths = 0.5
        rp.sheathLengths = 0.2
        rp.leafGrowthDurations = 0.3
        self.assertAlmostEqual(rp.getParameter("bladeAngle_dev"), 0.01)
        self.assertAlmostEqual(rp.getParameter("bladeLength_dev"), 0.5)
        self.assertAlmostEqual(rp.getParameter("sheathLength_dev"), 0.2)
        self.assertAlmostEqual(rp.getParameter("leafGrowthDuration_dev"), 0.3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
