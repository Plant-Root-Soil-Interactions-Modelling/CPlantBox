import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
import numpy as np


class TestHyphaeParameter(unittest.TestCase):

    def hyphae_example(self):
        """ hyphae parameter set with non-default values for all parameters introduced by HyphaeRandomParameter """
        self.plant = pb.Organism()
        self.hrp = pb.HyphaeRandomParameter(self.plant)
        self.hrp.name = "runnerhyphae"
        self.hrp.subType = 1
        self.hrp.dx = 0.1
        self.hrp.v = 0.12
        self.hrp.vs = 0.02
        self.hrp.b = 1.
        self.hrp.bs = 0.3
        self.hrp.b_prob = 0.4
        self.hrp.b_prob_s = 0.05
        self.hrp.hlt = 5.
        self.hrp.hlts = 0.5
        self.hrp.theta = 15. / 180. * np.pi
        self.hrp.thetas = 2. / 180. * np.pi
        self.hrp.distTH = 0.2
        self.hrp.ana = 0.5
        self.hrp.la = 0.5
        self.hrp.las = 0.05
        self.hrp.lb = 0.5
        self.hrp.lbs = 0.05
        self.hrp.ln = 0.5
        self.hrp.lns = 0.05
        self.hrp.lmax = 5.
        self.hrp.lmaxs = 0.5
        self.hrp.tropismT = 1
        self.hrp.tropismN = 2.
        self.hrp.tropismS = 0.1
        self.hrp.successor = [[1]]
        self.hrp.successorP = [[1.]]

    # all parameters introduced by HyphaeRandomParameter (set to non-default values in hyphae_example)
    hyphae_parameters = ["v", "vs", "b", "bs", "b_prob", "b_prob_s", "hlt", "hlts", "theta", "thetas", "distTH", "ana",
                         "la", "las", "lb", "lbs", "ln", "lns", "lmax", "lmaxs", "tropismT", "tropismN", "tropismS"]

    def set_deviations_zero(self):
        """ sets all standard deviations of hyphae_example to zero """
        for name in ["vs", "bs", "b_prob_s", "hlts", "thetas", "las", "lbs", "lns", "lmaxs"]:
            setattr(self.hrp, name, 0.)

    def test_constructors(self):
        """ tests constructor and copy """
        self.hyphae_example()
        hrp = self.hrp  # rename
        plant2 = pb.Organism()
        hrp2 = hrp.copy(plant2)
        self.assertIsNot(hrp, hrp2, "copy: organ type parameter set is not copied")
        self.assertEqual(hrp2.__class__.__name__, "HyphaeRandomParameter", "copy: unexpected class type")
        self.assertEqual(hrp2.name, hrp.name, "copy: value unexpected")
        self.assertEqual(hrp2.organType, hrp.organType, "copy: value unexpected")
        self.assertEqual(hrp2.subType, hrp.subType, "copy: value unexpected")
        self.assertEqual(hrp2.dx, hrp.dx, "copy: value unexpected")
        self.assertEqual(hrp2.successor, hrp.successor, "copy: value unexpected")
        for name in self.hyphae_parameters:
            self.assertEqual(getattr(hrp2, name), getattr(hrp, name), "copy: value unexpected for " + name)
        hrp.v = 0.9  # copy must be independent of the original
        self.assertEqual(hrp2.v, 0.12, "copy: copied parameter set changes with the original")
        self.assertEqual(hrp2.getParameter("v"), 0.12, "copy: class introspection points to the original")

    def test_parameter(self):
        """ tests default values via getParameter() """
        hrp = pb.HyphaeRandomParameter(pb.Organism())
        self.assertEqual(hrp.getParameter("organType"), pb.hyphae, "getParameter: value unexpected")
        self.assertEqual(hrp.getParameter("subType"), -1, "getParameter: value unexpected")
        defaults = {"v": 0.13, "v_dev": 0.01, "b": 2., "b_dev": 0.5, "b_prob": 0., "b_prob_dev": 0., "hlt": 10.,
                    "hlt_dev": 0.1, "theta": 30. / 180. * np.pi, "theta_dev": 0., "distTH": 0., "ana": 1.,
                    "la": 0.003, "la_dev": 0., "lb": 0.001, "lb_dev": 0., "ln": 0.005, "ln_dev": 0., "lmax": 10.,
                    "lmax_dev": 0., "tropismT": 2, "tropismN": 1., "tropismS": 0.3}
        for name, value in defaults.items():
            self.assertAlmostEqual(hrp.getParameter(name), value, 12, "getParameter: unexpected default for " + name)

    def test_parameter_changed(self):
        """ tests that changing the parameters manually is seen by getParameter() """
        self.hyphae_example()
        expected = {"subType": 1, "v": 0.12, "v_mean": 0.12, "v_dev": 0.02, "b": 1., "b_dev": 0.3, "b_prob": 0.4,
                    "b_prob_dev": 0.05, "hlt": 5., "hlt_dev": 0.5, "theta": 15. / 180. * np.pi,
                    "theta_dev": 2. / 180. * np.pi, "distTH": 0.2, "ana": 0.5, "la": 0.5, "la_dev": 0.05, "lb": 0.5,
                    "lb_dev": 0.05, "ln": 0.5, "ln_dev": 0.05, "lmax": 5., "lmax_dev": 0.5, "tropismT": 1,
                    "tropismN": 2., "tropismS": 0.1}
        for name, value in expected.items():
            self.assertAlmostEqual(self.hrp.getParameter(name), value, 12, "getParameter: value unexpected for " + name)

    def test_toString(self):
        """ tests __str__ output """
        hrp = pb.HyphaeRandomParameter(pb.Organism())
        hrp.name = "the hyphae"
        self.assertEqual(hrp.__str__(False), "name: the hyphae, organType: " + str(int(pb.hyphae)) + ", subType: -1.",
                         "toString: value unexpected")
        verbose = hrp.__str__(True)
        for name in ["v", "b", "b_prob", "hlt", "theta", "distTH", "ana", "la", "lb", "ln", "lmax", "tropismT",
                     "tropismN", "tropismS"]:
            self.assertIn(name, verbose, "toString: parameter " + name + " missing in verbose output")

    def test_xml(self):
        """ tests reading and writing xml parameter file """
        self.hyphae_example()
        hrp = self.hrp  # rename
        hrp.writeXML("hyphae.xml")
        hrp2 = pb.HyphaeRandomParameter(self.plant)
        hrp2.readXML("hyphae.xml")
        self.assertEqual(hrp2.name, hrp.name, "xml read/write: value unexpected")
        self.assertEqual(hrp2.organType, hrp.organType, "xml read/write: value unexpected")
        self.assertEqual(hrp2.subType, hrp.subType, "xml read/write: value unexpected")
        self.assertAlmostEqual(hrp2.dx, hrp.dx, 7, "xml read/write: value unexpected")
        for name in self.hyphae_parameters:
            self.assertAlmostEqual(getattr(hrp2, name), getattr(hrp, name), 7, "xml read/write: value unexpected for " + name)
        self.assertEqual(hrp2.successor[0][0], hrp.successor[0][0], "xml read/write: value unexpected")
        self.assertAlmostEqual(hrp2.successorP[0][0], hrp.successorP[0][0], 7, "xml read/write: value unexpected")

    def test_realize(self):
        """ tests realize() with laterals and zero standard deviations """
        self.hyphae_example()
        self.set_deviations_zero()
        p = self.hrp.realize()
        self.assertEqual(p.__class__.__name__, "HyphaeSpecificParameter", "realize: unexpected class type")
        self.assertEqual(p.subType, 1, "realize: unexpected sub type")
        self.assertAlmostEqual(p.a, self.hrp.a, 12, "realize: unexpected radius")
        self.assertAlmostEqual(p.v, 0.12, 12, "realize: unexpected tip elongation rate")
        self.assertAlmostEqual(p.b, 1., 12, "realize: unexpected branching rate")
        self.assertAlmostEqual(p.hlt, 5., 12, "realize: unexpected hyphal lifetime")
        self.assertAlmostEqual(p.theta, 15. / 180. * np.pi, 12, "realize: unexpected branching angle")
        self.assertTrue(p.laterals, "realize: hyphae with successors should have laterals")
        self.assertAlmostEqual(p.la, 0.5, 12, "realize: unexpected apical zone")
        self.assertAlmostEqual(p.lb, 0.5, 12, "realize: unexpected basal zone")
        nob = (self.hrp.lmax - self.hrp.la - self.hrp.lb) / self.hrp.ln + 1  # mean number of branching nodes
        self.assertEqual(len(p.ln) + 1, round(nob), "realize: internodal distances +1 should be number of branching nodes")
        for ln in p.ln:
            self.assertAlmostEqual(ln, 0.5, 12, "realize: unexpected inter-lateral distance")
        self.assertAlmostEqual(p.getK(), self.hrp.lmax, 10, "realize: maximal length getK() should equal lmax")

    def test_realize_no_laterals(self):
        """ tests realize() without successors: hyphae consists of the apical zone only """
        self.hyphae_example()
        self.set_deviations_zero()
        self.hrp.successor = []
        self.hrp.successorP = []
        p = self.hrp.realize()
        self.assertFalse(p.laterals, "realize: hyphae without successors should not have laterals")
        self.assertEqual(p.lb, 0., "realize: basal zone should be zero without laterals")
        self.assertAlmostEqual(p.la, self.hrp.lmax, 12, "realize: apical zone should be lmax without laterals")
        self.assertEqual(len(p.ln), 0, "realize: no inter-lateral distances expected without laterals")
        self.assertAlmostEqual(p.getK(), self.hrp.lmax, 12, "realize: maximal length getK() should equal lmax")

    def test_realize_non_negative(self):
        """ tests that negative realized values are cut off at 0 """
        self.hyphae_example()
        self.set_deviations_zero()
        self.hrp.v, self.hrp.b, self.hrp.hlt, self.hrp.theta = -1., -1., -1., -1.
        p = self.hrp.realize()
        self.assertEqual(p.v, 0., "realize: v should not be negative")
        self.assertEqual(p.b, 0., "realize: b should not be negative")
        self.assertEqual(p.hlt, 0., "realize: hlt should not be negative")
        self.assertEqual(p.theta, 0., "realize: theta should not be negative")

    def test_realize_deviation(self):
        """ tests that v, b, and hlt are drawn with the given mean and standard deviation """
        self.hyphae_example()
        self.plant.setSeed(1)
        self.hrp.v, self.hrp.vs = 0.4, 0.04
        self.hrp.b, self.hrp.bs = 2., 0.2
        self.hrp.hlt, self.hrp.hlts = 10., 1.
        v, b, hlt = [], [], []
        for i in range(0, 2000):
            p = self.hrp.realize()
            v.append(p.v)
            b.append(p.b)
            hlt.append(p.hlt)
        self.assertAlmostEqual(np.mean(v) / 0.4, 1., 1, "realize: unexpected mean of v")
        self.assertAlmostEqual(np.std(v) / 0.04, 1., 1, "realize: unexpected standard deviation of v")
        self.assertAlmostEqual(np.mean(b) / 2., 1., 1, "realize: unexpected mean of b")
        self.assertAlmostEqual(np.std(b) / 0.2, 1., 1, "realize: unexpected standard deviation of b")
        self.assertAlmostEqual(np.mean(hlt) / 10., 1., 1, "realize: unexpected mean of hlt")
        self.assertAlmostEqual(np.std(hlt) / 1., 1., 1, "realize: unexpected standard deviation of hlt")

    def test_realize_lmax(self):
        """ tests that the realized lmax equals the maximal length getK(), with and without laterals """
        self.hyphae_example()
        self.set_deviations_zero()
        p = self.hrp.realize()
        self.assertAlmostEqual(p.lmax, self.hrp.lmax, 10, "realize: unexpected lmax with laterals")
        self.assertAlmostEqual(p.lmax, p.getK(), 12, "realize: lmax should equal getK()")
        self.hrp.successor = []
        self.hrp.successorP = []
        p = self.hrp.realize()
        self.assertAlmostEqual(p.lmax, self.hrp.lmax, 12, "realize: unexpected lmax without laterals")

    def test_set_in_organism(self):
        """ tests if the hyphae parameter can be set and retrieved from an organism """
        self.hyphae_example()
        self.plant.setOrganRandomParameter(self.hrp)
        hrp_plant = self.plant.getOrganRandomParameter(pb.hyphae, self.hrp.subType)
        self.assertEqual(hrp_plant.__class__.__name__, "HyphaeRandomParameter", "organism: unexpected class type")
        for name in self.hyphae_parameters:
            self.assertEqual(getattr(hrp_plant, name), getattr(self.hrp, name), "organism: value unexpected for " + name)
        self.assertEqual(len(self.plant.getOrganRandomParameter(pb.hyphae)), 1, "organism: expected exactly one hyphae type")


if __name__ == '__main__':
    unittest.main()
