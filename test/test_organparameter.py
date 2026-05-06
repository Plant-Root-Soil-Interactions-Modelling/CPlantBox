import sys

sys.path.append("..")
sys.path.append("../src/")
import math
import os
import unittest

import plantbox as pb
from plantbox.rsml.rsml_reader import *


class TestOrganParameter(unittest.TestCase):
    """Unit tests for OrganRandomParameter and OrganSpecificParameter (base classes)."""

    def organ_example(self):
        """Returns an OrganRandomParameter with non-default scalar values used across tests."""
        plant = pb.Organism()
        orp = pb.OrganRandomParameter(plant)
        orp.name = "test organ"
        orp.subType = 3
        orp.a = 0.05
        orp.a_s = 0.01  # standard deviation of radius
        orp.dx = 0.5
        orp.dxMin = 1e-5
        orp.ldelay = 2.0
        orp.ldelays = 0.5
        return plant, orp

    def test_constructors(self):
        """Tests default values, copy creates a distinct object, and copy preserves all scalar fields."""
        plant = pb.Organism()
        orp = pb.OrganRandomParameter(plant)
        # check defaults defined in OrganRandomParameter
        self.assertEqual(orp.name, "organ", "default name unexpected")
        self.assertEqual(orp.organType, 0, "default organType unexpected")
        self.assertEqual(orp.subType, 0, "default subType unexpected")
        self.assertAlmostEqual(orp.a, 0.1, msg="default radius unexpected")
        self.assertAlmostEqual(orp.a_s, 0.0, msg="default radius std-dev unexpected")
        self.assertAlmostEqual(orp.dx, 0.25, msg="default dx unexpected")
        self.assertAlmostEqual(orp.dxMin, 1e-6, msg="default dxMin unexpected")
        self.assertAlmostEqual(orp.ldelay, -1.0, msg="default ldelay unexpected")
        self.assertAlmostEqual(orp.ldelays, 0.0, msg="default ldelays unexpected")
        # copy produces a separate object with identical values
        otp2 = orp.copy(plant)
        self.assertIsNot(orp, otp2, "copy: returned same object instead of a copy")
        self.assertEqual(otp2.name, orp.name, "copy: name not preserved")
        self.assertEqual(otp2.organType, orp.organType, "copy: organType not preserved")
        self.assertEqual(otp2.subType, orp.subType, "copy: subType not preserved")
        self.assertAlmostEqual(otp2.a, orp.a, msg="copy: radius not preserved")

    def test_parameter(self):
        """Tests getParameter() for plain names, _dev suffix, _mean suffix, and unknown names."""
        plant, orp = self.organ_example()
        # plain parameter look-up
        self.assertEqual(orp.getParameter("organType"), 0.0, "getParameter organType unexpected")
        self.assertEqual(orp.getParameter("subType"), 3.0, "getParameter subType unexpected")
        self.assertAlmostEqual(orp.getParameter("a"), 0.05, msg="getParameter a unexpected")
        self.assertAlmostEqual(orp.getParameter("dx"), 0.5, msg="getParameter dx unexpected")
        self.assertAlmostEqual(orp.getParameter("ldelay"), 2.0, msg="getParameter ldelay unexpected")
        # _dev suffix returns the standard deviation
        self.assertAlmostEqual(orp.getParameter("a_dev"), 0.01, msg="getParameter a_dev unexpected")
        self.assertAlmostEqual(orp.getParameter("ldelay_dev"), 0.5, msg="getParameter ldelay_dev unexpected")
        # _mean suffix is an alias for the mean value
        self.assertAlmostEqual(orp.getParameter("a_mean"), 0.05, msg="getParameter a_mean unexpected")
        # unknown parameter name raises RuntimeError
        with self.assertRaises(RuntimeError, msg="getParameter unknown name should raise RuntimeError"):
            orp.getParameter("nonexistent")

    def test_toString(self):
        """Tests __str__ in both non-verbose (one-liner) and verbose (full table) modes."""
        orp = pb.OrganRandomParameter(pb.Organism())
        orp.name = "the great organ"
        orp.organType = 1
        orp.subType = 2
        # non-verbose: compact one-liner
        self.assertEqual(orp.__str__(False), "name: the great organ, organType: 1, subType: 2.", "toString non-verbose: unexpected output")
        # verbose: should at minimum return a non-empty string without raising
        verbose_str = orp.__str__(True)
        self.assertIsInstance(verbose_str, str, "toString verbose: expected a string")
        self.assertGreater(len(verbose_str), 0, "toString verbose: returned empty string")

    def test_xml(self):
        """Writes parameters to XML and reads them back, checking all scalar fields round-trip correctly."""
        fname = "organ.xml"
        try:
            plant, orp = self.organ_example()
            orp.writeXML(fname)
            otp2 = pb.OrganRandomParameter(plant)
            otp2.readXML(fname)
            self.assertEqual(otp2.name, orp.name, "xml: name not preserved")
            self.assertEqual(otp2.subType, orp.subType, "xml: subType not preserved")
            self.assertAlmostEqual(otp2.a, orp.a, places=6, msg="xml: radius not preserved")
            self.assertAlmostEqual(otp2.a_s, orp.a_s, places=6, msg="xml: radius std-dev not preserved")
            self.assertAlmostEqual(otp2.dx, orp.dx, places=6, msg="xml: dx not preserved")
            self.assertAlmostEqual(otp2.ldelay, orp.ldelay, places=6, msg="xml: ldelay not preserved")
            self.assertAlmostEqual(otp2.ldelays, orp.ldelays, places=6, msg="xml: ldelays not preserved")
        finally:
            if os.path.exists(fname):
                os.remove(fname)

    def test_realize(self):
        """Tests realize(): correct class, correct subType, radius equals mean when std-dev is zero,
        and radius is always non-negative when std-dev is non-zero."""
        plant = pb.Organism()
        orp = pb.OrganRandomParameter(plant)
        orp.subType = 2
        orp.a = 0.08
        orp.a_s = 0.0  # deterministic: realised a must equal mean exactly
        p = orp.realize()
        self.assertEqual(p.__class__.__name__, "OrganSpecificParameter", "realize: unexpected class type")
        self.assertEqual(p.subType, 2, "realize: subType not preserved")
        self.assertAlmostEqual(p.a, 0.08, places=9, msg="realize: radius should equal mean when std-dev is zero")
        # with non-zero std-dev the radius must still be clamped to >= 0
        orp.a_s = 1e6  # extreme deviation: result is always max(a + noise, 0)
        for _ in range(20):
            p = orp.realize()
            self.assertGreaterEqual(p.a, 0.0, "realize: radius must be non-negative")


if __name__ == "__main__":
    unittest.main()
