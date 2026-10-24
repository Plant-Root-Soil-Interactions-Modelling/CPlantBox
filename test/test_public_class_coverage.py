"""Behavioral tests for core classes that lacked focused unit coverage."""
# ruff: noqa: E402, I001

import math
import os
import sys
import unittest
import importlib.util

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "../src")))

import plantbox as pb

if not hasattr(pb, "__path__"):
    pb.__path__ = [os.path.abspath(os.path.join(_HERE, "../src"))]

from functional.RootSoilSoluteInteraction import RootSoilSoluteInteraction
from functional.van_genuchten import Parameters as VanGenuchtenParameters
from rsml.rsml_data import RsmlData
from rsml.rsml_writer import LinkedPolylines, Metadata, Property
from plantbox.functional.Perirhizal import PerirhizalPython
from plantbox.functional.PerirhizalHeterogeneous import PerirhizalHetereogeneous
from plantbox.functional.Photosynthesis import PhotosynthesisPython
from plantbox.functional.phloem_flux import PhloemFluxPython
from plantbox.structural.MappedOrganism import MappedPlantPython
from plantbox.structural.StaticPlant import StaticPlant
from plantbox.visualisation.vtk_animate import AnimatePlant, AnimateRoots


class TestQuaternion(unittest.TestCase):

    def test_rotation_and_inverse(self):
        angle = math.pi / 2
        q = pb.Quaternion(math.cos(angle / 2), 0.0, 0.0, math.sin(angle / 2))
        rotated = q.Rotate(pb.Vector3d(1.0, 0.0, 0.0))
        self.assertAlmostEqual(rotated.x, 0.0, places=12)
        self.assertAlmostEqual(rotated.y, 1.0, places=12)
        self.assertAlmostEqual(rotated.z, 0.0, places=12)

        restored = q.inverse().Rotate(rotated)
        self.assertAlmostEqual(restored.x, 1.0, places=12)
        self.assertAlmostEqual(restored.y, 0.0, places=12)
        self.assertAlmostEqual(restored.z, 0.0, places=12)


class TestSignedDistanceFunctions(unittest.TestCase):

    def setUp(self):
        self.origin = pb.Vector3d(0.0, 0.0, 0.0)
        self.box = pb.SDF_Cuboid(pb.Vector3d(-1.0, -1.0, -1.0), pb.Vector3d(1.0, 1.0, 1.0))

    def test_base_distance_is_unconstrained(self):
        sdf = pb.SignedDistanceFunction()
        self.assertLess(sdf.getDist(self.origin), -1e90)
        self.assertEqual(sdf.dist([0.0, 0.0, 0.0]), sdf.getDist(self.origin))

    def test_cuboid_and_container_distances(self):
        self.assertLess(self.box.getDist(self.origin), 0.0)
        self.assertGreater(self.box.getDist(pb.Vector3d(2.0, 0.0, 0.0)), 0.0)

        container = pb.SDF_PlantContainer(2.0, 2.0, 4.0, 0.0)
        self.assertLess(container.getDist(pb.Vector3d(0.0, 0.0, -1.0)), 0.0)
        self.assertGreater(container.getDist(pb.Vector3d(10.0, 0.0, 0.0)), 0.0)

    def test_transform_and_boolean_operations(self):
        translated = pb.SDF_RotateTranslate(self.box, pb.Vector3d(3.0, 0.0, 0.0))
        self.assertLess(translated.getDist(pb.Vector3d(3.0, 0.0, 0.0)), 0.0)
        self.assertGreater(translated.getDist(self.origin), 0.0)

        other = pb.SDF_RotateTranslate(self.box, pb.Vector3d(3.0, 0.0, 0.0))
        intersection = pb.SDF_Intersection(self.box, other)
        union = pb.SDF_Union(self.box, other)
        difference = pb.SDF_Difference(self.box, other)
        complement = pb.SDF_Complement(self.box)

        self.assertGreater(intersection.getDist(self.origin), 0.0)
        self.assertLess(union.getDist(self.origin), 0.0)
        self.assertLess(difference.getDist(self.origin), 0.0)
        self.assertGreater(difference.getDist(pb.Vector3d(3.0, 0.0, 0.0)), 0.0)
        self.assertGreater(complement.getDist(self.origin), 0.0)
        self.assertLess(complement.getDist(pb.Vector3d(2.0, 0.0, 0.0)), 0.0)

    def test_half_plane_distance_and_gradient(self):
        plane = pb.SDF_HalfPlane(self.origin, pb.Vector3d(0.0, 0.0, 1.0))
        point_above = pb.Vector3d(0.0, 0.0, 2.0)
        self.assertAlmostEqual(plane.getDist(point_above), 2.0)
        gradient = plane.getGradient(self.origin)
        self.assertAlmostEqual(gradient.x, 0.0, places=8)
        self.assertAlmostEqual(gradient.y, 0.0, places=8)
        self.assertAlmostEqual(gradient.z, 1.0, places=8)


class TestSoilLookupClasses(unittest.TestCase):

    def test_base_and_sdf_lookup(self):
        base = pb.SoilLookUp()
        self.assertEqual(base.getValue(pb.Vector3d(5.0, 0.0, 0.0)), 1.0)

        box = pb.SDF_PlantBox(2.0, 2.0, 2.0)
        lookup = pb.SoilLookUpSDF(box, 4.0, 0.0, 0.5)
        self.assertEqual(lookup.getValue(pb.Vector3d(0.0, 0.0, -1.0)), 4.0)
        self.assertEqual(lookup.getValue(pb.Vector3d(20.0, 0.0, -1.0)), 0.0)

    def test_lookup_composition(self):
        first = pb.ProportionalElongation()
        first.setScale(2.0)
        second = pb.ProportionalElongation()
        second.setScale(3.0)
        product = pb.MultiplySoilLookUps(first, second)
        self.assertEqual(product.getValue(pb.Vector3d()), 6.0)

        scaled = pb.ProportionalElongation()
        scaled.setBaseLookUp(first)
        scaled.setScale(0.5)
        self.assertEqual(scaled.getValue(pb.Vector3d()), 1.0)

    def test_one_dimensional_grid(self):
        grid = pb.Grid1D(3, [0.0, 1.0, 2.0], [10.0, 20.0, 30.0])
        self.assertEqual(grid.map(1.5), 1)
        self.assertEqual(grid.getValue(pb.Vector3d(0.0, 0.0, 1.5)), 20.0)

        equidistant = pb.EquidistantGrid1D(0.0, 2.0, [10.0, 20.0, 30.0])
        self.assertEqual(equidistant.map(1.5), 1)
        self.assertEqual(equidistant.getValue(pb.Vector3d(0.0, 0.0, 1.5)), 20.0)

    def test_three_dimensional_grids(self):
        x = pb.Grid1D(3, [0.0, 1.0, 2.0], [0.0, 0.0, 0.0])
        y = pb.Grid1D(3, [0.0, 1.0, 2.0], [0.0, 0.0, 0.0])
        z = pb.Grid1D(3, [0.0, 1.0, 2.0], [0.0, 0.0, 0.0])
        grid = pb.RectilinearGrid3D(x, y, z)
        grid.setData(1, 0, 1, 12.5)
        self.assertEqual(grid.map(1.5, 0.5, 1.5), 10)
        self.assertEqual(grid.getValue(pb.Vector3d(1.5, 0.5, 1.5)), 12.5)
        self.assertEqual(grid.getGridPoint(1, 0, 1).z, 1.0)

        equidistant = pb.EquidistantGrid3D(2.0, 2.0, 2.0, 2, 2, 2)
        equidistant.setData(0, 0, 0, 7.0)
        self.assertEqual(equidistant.getValue(pb.Vector3d(0.0, 0.0, 0.0)), 7.0)


class TestGrowthClasses(unittest.TestCase):

    def test_growth_function_base_reports_unimplemented_operation(self):
        growth = pb.GrowthFunction()
        with self.assertRaises(RuntimeError):
            growth.getLength(1.0, 1.0, 10.0, None)

    def test_exponential_growth_length_and_inverse(self):
        growth = pb.ExponentialGrowth()
        age, rate, maximum = 2.0, 3.0, 10.0
        length = maximum * (1.0 - math.exp(-(rate / maximum) * age))
        self.assertAlmostEqual(growth.getLength(age, rate, maximum, None), length)
        self.assertAlmostEqual(growth.getAge(length, rate, maximum, None), age)

    def test_cw_limited_growth_falls_back_to_exponential(self):
        growth = pb.CWLimitedGrowth()
        exponential = pb.ExponentialGrowth()
        self.assertAlmostEqual(
            growth.getLength(2.0, 3.0, 10.0, None),
            exponential.getLength(2.0, 3.0, 10.0, None),
        )
        self.assertAlmostEqual(
            growth.getAge(2.0, 3.0, 10.0, None),
            exponential.getAge(2.0, 3.0, 10.0, None),
        )


class TestTropismClasses(unittest.TestCase):

    def test_concrete_tropisms_copy_to_the_requested_plant(self):
        plant = pb.Plant()
        target = pb.Plant()
        soil = pb.ProportionalElongation()
        soil.setScale(1.0)
        tropisms = [
            pb.Gravitropism(plant, 1.0, 0.1),
            pb.Plagiotropism(plant, 1.0, 0.1),
            pb.Exotropism(plant, 1.0, 0.1),
            pb.Hydrotropism(plant, 1.0, 0.1, soil),
        ]
        expected_types = [pb.Gravitropism, pb.Plagiotropism, pb.Exotropism, pb.Hydrotropism]

        for tropism, expected_type in zip(tropisms, expected_types):
            with self.subTest(tropism=expected_type.__name__):
                self.assertIsInstance(tropism.copy(target), expected_type)
                self.assertIs(tropism.getPlant(), plant)


class TestPlantVisualiser(unittest.TestCase):

    def test_empty_visualiser_geometry_lifecycle(self):
        visualiser = pb.PlantVisualiser()
        self.assertFalse(visualiser.HasGeometry())
        visualiser.SetGeometryResolution(8)
        visualiser.SetLeafResolution(4)
        visualiser.SetComputeMidlineInLeaf(True)
        visualiser.ResetGeometry()
        self.assertFalse(visualiser.HasGeometry())
        self.assertEqual(len(visualiser.GetGeometry()), 0)
        self.assertEqual(len(visualiser.GetGeometryIndices()), 0)


class TestAdditionalCoreClasses(unittest.TestCase):

    def test_static_root_construction(self):
        plant = pb.Plant()
        random_parameter = pb.RootRandomParameter(plant)
        root = pb.StaticRoot(0, random_parameter.realize(), 0.0, 0)
        self.assertIsInstance(root, pb.StaticRoot)

    def test_perirhizal_default_state(self):
        model = pb.Perirhizal()
        self.assertIsNone(model.ms)

    def test_photosynthesis_construction(self):
        plant = pb.MappedPlant()
        parameters = pb.PlantHydraulicParameters()
        model = pb.Photosynthesis(plant, parameters)
        self.assertIs(model.plant, plant)

    def test_phloem_flux_is_bound_when_optional_backend_is_enabled(self):
        if not hasattr(pb, "PhloemFlux"):
            self.skipTest("PhloemFlux requires the optional PIAFMUNCH backend")
        plant = pb.MappedPlant()
        parameters = pb.PlantHydraulicParameters()
        model = pb.PhloemFlux(plant, parameters, -500.0, 350e-6)
        self.assertIs(model.plant, plant)

    def test_unconstructible_exudation_model_is_registered(self):
        self.assertTrue(callable(pb.ExudationModel))


class TestPythonHelperClasses(unittest.TestCase):

    def test_van_genuchten_parameters(self):
        parameters = VanGenuchtenParameters([0.1, 0.4, 0.02, 1.6, 10.0])
        self.assertAlmostEqual(parameters.theta_R, 0.1)
        self.assertAlmostEqual(parameters.theta_S, 0.4)
        self.assertAlmostEqual(parameters.m, 1.0 - 1.0 / 1.6)
        self.assertEqual(list(parameters), [0.1, 0.4, 0.02, 1.6, 10.0])

    def test_linked_polylines_defaults_and_metadata(self):
        polyline = LinkedPolylines()
        self.assertEqual(polyline.branchnumber, -1)
        self.assertEqual(polyline.polyline, [])
        self.assertEqual(polyline.laterals, [])

        metadata = Metadata()
        metadata.add_property(Property("radius", "float", "cm", [0.1]))
        LinkedPolylines.set_metadata(metadata)
        self.assertIs(LinkedPolylines.metadata, metadata)

    def test_rsml_data_setters(self):
        data = RsmlData()
        self.assertFalse(data.exists())
        data.set_rsml([[[0.0, 0.0, 0.0]]], {}, {}, Metadata())
        self.assertTrue(data.exists())
        data.set_selected([0.1], [float("nan")], [2], ["radius"])
        self.assertEqual(data.cts.tolist(), [0.0])
        self.assertEqual(data.max_ct, 0.0)

    def test_root_soil_solute_flux(self):
        class Segments:
            segments = [0]
            radii = [0.1]

            @staticmethod
            def segLength():
                return [1.0]

        model = RootSoilSoluteInteraction(Segments(), 2.0, 1.0, 0.5, 0.0)
        flux = model.solute_fluxes([1.5])
        expected = -2.0 * math.pi * 0.1 * 1.0 * 1e-4 * ((1.5 - 0.5) * 2.0 / (1.0 + (1.5 - 0.5))) * 1e3
        self.assertAlmostEqual(flux[0], expected)
        self.assertEqual(model.solute_fluxes([0.1])[0], 0.0)


class TestPythonModelWrappers(unittest.TestCase):

    def test_mapped_plant_adapter_delegates_and_converts_arrays(self):
        adapter = MappedPlantPython()
        self.assertIsInstance(adapter._base, pb.MappedPlant)
        self.assertEqual(adapter.get_nodes().shape, (0,))
        adapter.bladeLength = [3.0]
        self.assertEqual(adapter._base.bladeLength, [3.0])

    def test_static_plant_initial_state(self):
        plant = StaticPlant()
        self.assertIsInstance(plant, pb.Plant)
        self.assertEqual(plant.static_organs, {})
        self.assertIsNone(plant.data)
        self.assertAlmostEqual(plant.polyline_length_(0, 2, [[0, 0, 0], [0, 3, 0], [0, 3, 4]]), 7.0)

    def test_perirhizal_python_defaults_and_heterogeneous_layers(self):
        model = PerirhizalPython()
        self.assertIsNone(model.lookup_table)
        self.assertEqual(model.alpha_0, 0.3)

        heterogeneous = PerirhizalHetereogeneous()
        heterogeneous.addLayers([-20.0, 0.0, -10.0], [2, 0, 1])
        self.assertEqual(heterogeneous.corners.tolist(), [0.0, -10.0, -20.0])
        self.assertEqual(heterogeneous.layer_indices.tolist(), [0, 1, 2])
        heterogeneous.addVanGenuchtenParamters(VanGenuchtenParameters([0.1, 0.4, 0.02, 1.6, 10.0]))
        self.assertEqual(len(heterogeneous.sp), 1)
        self.assertEqual(heterogeneous.lookup_tables, [None])

    def test_photosynthesis_python_wrapper_initial_state(self):
        plant = pb.MappedPlant()
        parameters = pb.PlantHydraulicParameters()
        model = PhotosynthesisPython(plant, parameters)
        self.assertIs(model.plant, plant)
        self.assertAlmostEqual(model.pCO2, 350e-6)

    def test_phloem_flux_python_wrapper_initial_state(self):
        if not hasattr(pb, "PhloemFlux"):
            self.skipTest("PhloemFlux requires the optional PIAFMUNCH backend")
        plant = pb.MappedPlant()
        parameters = pb.PlantHydraulicParameters()
        model = PhloemFluxPython(plant, parameters)
        self.assertIs(model.plant, plant)

    def test_animation_classes_can_be_constructed_without_rendering(self):
        root_animation = AnimateRoots(pb.RootSystem())
        plant_animation = AnimatePlant(pb.Plant())
        self.assertIsInstance(root_animation, AnimateRoots)
        self.assertIsInstance(plant_animation, AnimateRoots)


class TestOptionalFiPyClasses(unittest.TestCase):

    @unittest.skipUnless(importlib.util.find_spec("fipy"), "FiPy is an optional dependency")
    def test_fipy_variable_adapters_and_mesh(self):
        from plantbox.functional.CellVariablemod import CellVariablemod, FaceVariablemod
        from plantbox.functional.Mesh1Dmod import Mesh1Dmod

        mesh = Mesh1Dmod(
            radiiVertices=[0.1, 0.1, 0.1],
            radiiCells=[0.1, 0.1],
            length=[0.5, 0.5],
            vertexCoords=[[0.0, 0.5, 1.0]],
            faceVertexIDs=[[0, 1, 2]],
            cellFaceIDs=[[0, 1], [1, 2]],
        )
        cell_variable = CellVariablemod(mesh=mesh, value=1.0)
        face_variable = FaceVariablemod(mesh=mesh, value=2.0)
        self.assertIs(cell_variable.faceValue._variableClass, FaceVariablemod)
        self.assertIs(cell_variable.faceGrad._variableClass, FaceVariablemod)
        self.assertIsInstance(face_variable, FaceVariablemod)
        self.assertEqual(face_variable.shape[-1], 3)


if __name__ == "__main__":
    unittest.main()
