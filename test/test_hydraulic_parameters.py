import sys

sys.path.append("..")
sys.path.append("../src")

import io
import json
import os
import tempfile
import unittest
from contextlib import redirect_stdout

import numpy as np

import plantbox as pb
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters


class TestPlantHydraulicParameters(unittest.TestCase):

    def setUp(self):
        self.params = PlantHydraulicParameters()
        self.ot_root = int(pb.OrganTypes.root)
        self.ot_stem = int(pb.OrganTypes.stem)
        self.ot_leaf = int(pb.OrganTypes.leaf)

    # --- Default state ---

    def test_default_kr_mode(self):
        self.assertEqual(self.params.krMode, "const")

    def test_default_kx_mode(self):
        self.assertEqual(self.params.kxMode, "const")

    def test_default_kr_root_nonzero(self):
        self.assertGreater(self.params.kr_f(0, 5.0, 0, self.ot_root), 0.0)

    def test_default_kr_stem_zero(self):
        self.assertAlmostEqual(self.params.kr_f(0, 5.0, 0, self.ot_stem), 0.0)

    def test_default_kr_leaf_zero(self):
        self.assertAlmostEqual(self.params.kr_f(0, 5.0, 0, self.ot_leaf), 0.0)

    # --- Constant conductivities ---

    def test_kr_const_all_subtypes(self):
        self.params.set_kr_const(0.01)
        for st in range(0, 10):
            self.assertAlmostEqual(self.params.kr_f(0, 5.0, st, self.ot_root), 0.01, msg=f"subType={st}")

    def test_kx_const_all_subtypes(self):
        self.params.set_kx_const(0.5)
        for st in range(0, 10):
            self.assertAlmostEqual(self.params.kx_f(0, 5.0, st, self.ot_root), 0.5, msg=f"subType={st}")

    def test_kr_const_single_subtype(self):
        self.params.set_kr_const(1e-4)
        self.params.set_kr_const(0.02, subType=2)
        self.assertAlmostEqual(self.params.kr_f(0, 5.0, 2, self.ot_root), 0.02)
        self.assertAlmostEqual(self.params.kr_f(0, 5.0, 1, self.ot_root), 1e-4)

    def test_kx_const_list_subtypes(self):
        self.params.set_kx_const(1.0)
        self.params.set_kx_const(2.0, subType=[1, 3])
        self.assertAlmostEqual(self.params.kx_f(0, 5.0, 1, self.ot_root), 2.0)
        self.assertAlmostEqual(self.params.kx_f(0, 5.0, 3, self.ot_root), 2.0)
        self.assertAlmostEqual(self.params.kx_f(0, 5.0, 2, self.ot_root), 1.0)

    def test_kr_const_mode_string(self):
        self.params.set_kr_const(0.01)
        self.assertEqual(self.params.krMode, "const")

    def test_kx_const_mode_string(self):
        self.params.set_kx_const(0.5)
        self.assertEqual(self.params.kxMode, "const")

    def test_kr_const_age_independent(self):
        self.params.set_kr_const(0.03)
        val_young = self.params.kr_f(0, 1.0, 0, self.ot_root)
        val_old = self.params.kr_f(0, 100.0, 0, self.ot_root)
        self.assertAlmostEqual(val_young, val_old)

    # --- Organ-type scoping ---

    def test_kr_const_default_organtype_is_root(self):
        initial_stem = self.params.kr_f(0, 5.0, 0, self.ot_stem)
        self.params.set_kr_const(0.99)
        self.assertAlmostEqual(self.params.kr_f(0, 5.0, 0, self.ot_root), 0.99)
        self.assertAlmostEqual(self.params.kr_f(0, 5.0, 0, self.ot_stem), initial_stem)

    def test_kr_const_all_organ_types(self):
        self.params.set_kr_const(0.1, organType=-1)
        for ot in [self.ot_root, self.ot_stem, self.ot_leaf]:
            self.assertAlmostEqual(self.params.kr_f(0, 5.0, 0, ot), 0.1, msg=f"organType={ot}")

    def test_kr_const_stem_only(self):
        self.params.set_kr_const(0.007, organType=self.ot_stem)
        self.assertAlmostEqual(self.params.kr_f(0, 5.0, 0, self.ot_stem), 0.007)

    # --- Age-dependent conductivities ---

    def test_kr_age_exact_point(self):
        self.params.set_kr_age_dependent([0.0, 10.0, 20.0], [0.001, 0.002, 0.001])
        self.assertAlmostEqual(self.params.kr_f(0, 10.0, 0, self.ot_root), 0.002)

    def test_kr_age_interpolation(self):
        self.params.set_kr_age_dependent([0.0, 10.0], [0.0, 1.0])
        self.assertAlmostEqual(self.params.kr_f(0, 5.0, 0, self.ot_root), 0.5)

    def test_kx_age_interpolation(self):
        self.params.set_kx_age_dependent([0.0, 20.0], [0.0, 2.0])
        self.assertAlmostEqual(self.params.kx_f(0, 10.0, 0, self.ot_root), 1.0)

    def test_kr_age_extrapolation_high(self):
        self.params.set_kr_age_dependent([0.0, 10.0], [0.001, 0.002])
        self.assertAlmostEqual(self.params.kr_f(0, 50.0, 0, self.ot_root), 0.002)

    def test_kr_age_extrapolation_low(self):
        self.params.set_kr_age_dependent([0.0, 10.0], [0.001, 0.002])
        self.assertAlmostEqual(self.params.kr_f(0, -1.0, 0, self.ot_root), 0.001)

    def test_kr_age_mode_string(self):
        self.params.set_kr_age_dependent([0.0, 10.0], [0.001, 0.002])
        self.assertEqual(self.params.krMode, "age dependent")

    def test_kx_age_mode_string(self):
        self.params.set_kx_age_dependent([0.0, 10.0], [0.01, 0.1])
        self.assertEqual(self.params.kxMode, "age dependent")

    def test_kr_age_per_subtype(self):
        self.params.set_kr_age_dependent([0.0, 10.0], [0.0, 1.0])
        self.params.set_kr_age_dependent([0.0, 10.0], [0.0, 2.0], subType=1)
        self.assertAlmostEqual(self.params.kr_f(0, 10.0, 1, self.ot_root), 2.0)
        self.assertAlmostEqual(self.params.kr_f(0, 10.0, 0, self.ot_root), 1.0)

    # --- Per-segment conductivities ---

    def test_kr_per_segment_mode_string(self):
        self.params.setKrValues([0.01, 0.02, 0.03])
        self.assertEqual(self.params.krMode, "per segment")

    def test_kx_per_segment_mode_string(self):
        self.params.setKxValues([0.1, 0.2])
        self.assertEqual(self.params.kxMode, "per segment")

    def test_kr_per_segment_returns_correct_value(self):
        self.params.setKrValues([0.01, 0.05, 0.03])
        self.assertAlmostEqual(self.params.kr_f(0, 0.0, 0, self.ot_root), 0.01)
        self.assertAlmostEqual(self.params.kr_f(1, 0.0, 0, self.ot_root), 0.05)
        self.assertAlmostEqual(self.params.kr_f(2, 0.0, 0, self.ot_root), 0.03)

    def test_kx_per_segment_returns_correct_value(self):
        self.params.setKxValues([0.1, 0.3])
        self.assertAlmostEqual(self.params.kx_f(0, 0.0, 0, self.ot_root), 0.1)
        self.assertAlmostEqual(self.params.kx_f(1, 0.0, 0, self.ot_root), 0.3)

    # --- setMode ---

    def test_set_mode_restores_const(self):
        self.params.setKrAgeDependent([0.0, 10.0], [0.001, 0.002], 0, self.ot_root)
        self.params.setMode("const", "const")
        self.assertEqual(self.params.krMode, "const")
        self.assertEqual(self.params.kxMode, "const")

    def test_set_mode_age_dependent(self):
        self.params.setMode("age dependent", "age dependent")
        self.assertEqual(self.params.krMode, "age dependent")
        self.assertEqual(self.params.kxMode, "age dependent")

    def test_set_mode_per_segment(self):
        self.params.setMode("per segment", "per segment")
        self.assertEqual(self.params.krMode, "per segment")
        self.assertEqual(self.params.kxMode, "per segment")

    def test_set_mode_distance_dependent(self):
        self.params.setMode("distance dependent", "distance dependent")
        self.assertEqual(self.params.krMode, "distance dependent")
        self.assertEqual(self.params.kxMode, "distance dependent")

    def test_set_mode_mixed(self):
        self.params.setMode("age dependent", "const")
        self.assertEqual(self.params.krMode, "age dependent")
        self.assertEqual(self.params.kxMode, "const")

    # --- JSON write / read round-trip ---

    def test_json_file_created(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "out")
            self.params.write_parameters(path)
            self.assertTrue(os.path.isfile(path + ".json"))

    def test_json_contains_expected_keys(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "out")
            self.params.write_parameters(path)
            with open(path + ".json") as f:
                data = json.load(f)
        for key in ["mode", "kr_ages", "kr_values", "kx_ages", "kx_values", "krPerSegment", "kxPerSegment"]:
            self.assertIn(key, data)

    def test_write_read_roundtrip_const(self):
        self.params.set_kr_const(0.012)
        self.params.set_kx_const(0.34)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "params")
            self.params.write_parameters(path)
            params2 = PlantHydraulicParameters()
            params2.read_parameters(path)
        self.assertAlmostEqual(params2.kr_f(0, 5.0, 0, self.ot_root), 0.012)
        self.assertAlmostEqual(params2.kx_f(0, 5.0, 0, self.ot_root), 0.34)

    def test_write_read_roundtrip_age_dependent(self):
        ages = [0.0, 10.0, 20.0]
        values = [0.001, 0.002, 0.0015]
        self.params.set_kr_age_dependent(ages, values)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "params")
            self.params.write_parameters(path)
            params2 = PlantHydraulicParameters()
            params2.read_parameters(path)
        self.assertAlmostEqual(params2.kr_f(0, 0.0, 0, self.ot_root), 0.001)
        self.assertAlmostEqual(params2.kr_f(0, 10.0, 0, self.ot_root), 0.002)
        self.assertAlmostEqual(params2.kr_f(0, 15.0, 0, self.ot_root), 0.00175)

    def test_write_read_restores_mode(self):
        self.params.set_kr_age_dependent([0.0, 10.0], [0.001, 0.002])
        self.params.set_kx_const(0.5)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "params")
            self.params.write_parameters(path)
            params2 = PlantHydraulicParameters()
            params2.read_parameters(path)
        self.assertEqual(params2.krMode, "age dependent")
        self.assertEqual(params2.kxMode, "const")

    def test_write_per_segment_prints_warning(self):
        self.params.setKrValues([0.01, 0.02])
        buf = io.StringIO()
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "params")
            with redirect_stdout(buf):
                self.params.write_parameters(path)
        self.assertIn("Warning", buf.getvalue())

    def test_write_read_multiple_subtypes(self):
        self.params.set_kr_age_dependent([0.0, 10.0], [0.001, 0.002], subType=0)
        self.params.set_kr_age_dependent([0.0, 10.0], [0.003, 0.006], subType=1)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "params")
            self.params.write_parameters(path)
            params2 = PlantHydraulicParameters()
            params2.read_parameters(path)
        self.assertAlmostEqual(params2.kr_f(0, 10.0, 0, self.ot_root), 0.002)
        self.assertAlmostEqual(params2.kr_f(0, 10.0, 1, self.ot_root), 0.006)


if __name__ == "__main__":
    unittest.main()
    print("done.")
