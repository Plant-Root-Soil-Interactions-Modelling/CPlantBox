import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
import numpy as np
# from rsml.rsml_reader import * QUESTION? do we need this it causes a missing import flag

class TestMycParameter(unittest.TestCase):

    def mycroot_example(self):
        self.plant = pb.Organism()
        self.mrrp = pb.MycorrhizalRootRandomParameter(self.plant)
        self.mrrp.la = 1.5
        self.mrrp.lb = 5.5
        self.mrrp.ln = 1.25
        self.mrrp.lns = 0.12
        self.mrrp.lmax = 7 * self.mrrp.ln + self.mrrp.la + self.mrrp.lb
        self.mrrp.subType = 1
        self.mrrp.successor = [[4, 5, 6]]
        self.mrrp.successorP = [[0.4, 0.1, 0.5]]
        self.mrrp.ldelay = 5
        self.mrrp.ldelays = 2
        self.mrrp.vi = 0.5
        self.mrrp.vis = 0.05
        self.mrrp.minAge = 0.1
        self.mrrp.maxAge = 10
        self.mrrp.lmbd = 0.5
        self.mrrp.lmbds = 0.02
        self.mrrp.maxColonization = 0.5
        self.mrrp.hyphalDelay = 0.5
        self.mrrp.highresolution = 0
        self.mrrp.dx_inf = 0.2
        self.mrrp.hyphalEmergenceDensity = 5

    # all parameters introduced by MycorrhizalRootRandomParameter (set to non-default values in mycroot_example)
    myc_parameters = ["lmbd", "lmbds", "minAge", "maxAge", "vi", "vis", "maxColonization", "hyphalDelay",
                      "highresolution", "dx_inf", "hyphalEmergenceDensity"]

    def test_constructors(self):
        """ tests constructor and copy """
        self.mycroot_example()
        mrrp = self.mrrp  # rename
        mrrp.theta = 123
        mrrp.thetas = 456
        mrrp.gf = 789
        plant2 = pb.Organism()
        otp2 = mrrp.copy(plant2)
        self.assertIsNot(mrrp, otp2, "copy: organ type parameter set is not copied")
        self.assertEqual(otp2.__class__.__name__, "MycorrhizalRootRandomParameter", "copy: unexpected class type")
        self.assertEqual(otp2.name, mrrp.name, "copy: value unexpected")
        self.assertEqual(otp2.organType, mrrp.organType, "copy: value unexpected")
        self.assertEqual(otp2.subType, mrrp.subType, "copy: value unexpected")
        self.assertEqual(otp2.a, mrrp.a, "copy: value unexpected")
        self.assertEqual(otp2.theta, mrrp.theta, "copy: value unexpected")
        self.assertEqual(otp2.thetas, mrrp.thetas, "copy: value unexpected")
        self.assertEqual(otp2.gf, mrrp.gf, "copy: value unexpected")
        for name in self.myc_parameters:
            self.assertEqual(getattr(otp2, name), getattr(mrrp, name), "copy: value unexpected for " + name)
        mrrp.lmbd = 0.9  # copy must be independent of the original
        self.assertEqual(otp2.lmbd, 0.5, "copy: copied parameter set changes with the original")
        self.assertEqual(otp2.getParameter("lmbd"), 0.5, "copy: class introspection points to the original")

    def test_parameter(self):
        """ tests getParameter() """
        self.plant = pb.Organism()
        mrrp = pb.MycorrhizalRootRandomParameter(self.plant)
        mrrp.lns = 0.123
        mrrp.la = 12
        mrrp.vi = 0.5
        ot = mrrp.getParameter("organType")  # test defaults
        st = mrrp.getParameter("subType")
        gf = mrrp.getParameter("gf")
        ln = mrrp.getParameter("ln")
        lns = mrrp.getParameter("ln_dev")
        la = mrrp.getParameter("la_mean")  # we can always add "_mean" to avoid naming conflicts
        lmbd = mrrp.getParameter("lmbd_mean")
        lmbds = mrrp.getParameter("lmbd_dev")
        vi = mrrp.getParameter("vi")
        vis = mrrp.getParameter("vi_dev")
        self.assertEqual(ot, 2., "getParameter: value unexpected")
        self.assertEqual(st, -1., "getParameter: value unexpected")
        self.assertEqual(gf, 1., "getParameter: value unexpected")
        self.assertEqual(ln, 1., "getParameter: value unexpected")
        self.assertEqual(lns, 0.123, "getParameter: value unexpected")
        self.assertEqual(la, 12, "getParameter: value unexpected")
        self.assertEqual(lmbd, 0.15, "getParameter: value unexpected")
        self.assertEqual(lmbds, 0., "getParameter: value unexpected")
        self.assertEqual(vi, 0.5, "getParameter: value unexpected")
        self.assertEqual(vis, 0., "getParameter: value unexpected")
        mrrp.theta = 123  # change values
        mrrp.thetas = 456
        theta = mrrp.getParameter("theta")
        thetas = mrrp.getParameter("theta_dev")
        self.assertEqual(theta, 123, "getParameter: value unexpected")
        self.assertEqual(thetas, 456, "getParameter: value unexpected")

    def test_parameter_defaults(self):
        """ tests default values of the mycorrhizal parameters via getParameter() """
        mrrp = pb.MycorrhizalRootRandomParameter(pb.Organism())
        defaults = {"lmbd": 0.15, "lmbd_dev": 0., "minAge": 0., "maxAge": 32., "vi": 0.13, "vi_dev": 0.,
                    "maxColonization": 1., "hyphalDelay": 0., "highresolution": 1., "dx_inf": 0.1,
                    "hyphalEmergenceDensity": 1.}
        for name, value in defaults.items():
            self.assertAlmostEqual(mrrp.getParameter(name), value, 12, "getParameter: unexpected default for " + name)

    def test_parameter_changed(self):
        """ tests that changing the mycorrhizal parameters manually is seen by getParameter() """
        self.mycroot_example()
        mrrp = self.mrrp  # rename
        expected = {"lmbd": 0.5, "lmbd_mean": 0.5, "lmbd_dev": 0.02, "minAge": 0.1, "maxAge": 10., "vi": 0.5,
                    "vi_mean": 0.5, "vi_dev": 0.05, "maxColonization": 0.5, "hyphalDelay": 0.5,
                    "highresolution": 0., "dx_inf": 0.2, "hyphalEmergenceDensity": 5.}
        for name, value in expected.items():
            self.assertAlmostEqual(mrrp.getParameter(name), value, 12, "getParameter: value unexpected for " + name)

    def test_toString(self):
        """ tests __str__ output """
        self.mrrp = pb.MycorrhizalRootRandomParameter(pb.Organism())
        mrrp = self.mrrp  # rename
        mrrp.name = "the root"
        self.assertEqual(mrrp.__str__(False), "name: the root, organType: 2, subType: -1.", "toString: value unexpected")
        verbose = mrrp.__str__(True)
        for name in ["lmbd", "minAge", "maxAge", "vi", "maxColonization", "hyphalDelay", "highresolution", "dx_inf",
                     "hyphalEmergenceDensity"]:
            self.assertIn(name, verbose, "toString: parameter " + name + " missing in verbose output")

    def test_xml(self):
        """ write the organ as xml, and rereads it """
        self.mycroot_example()
        mrrp = self.mrrp  # rename
        mrrp.name = "lateral"
        mrrp.subType = 2
        mrrp.writeXML("mycroot.xml")
        otp2 = pb.MycorrhizalRootRandomParameter(self.plant)
        otp2.readXML("mycroot.xml")
        self.assertEqual(otp2.ldelay, mrrp.ldelay, "xml: value unexpected")
        self.assertEqual(otp2.ldelays, mrrp.ldelays, "xml: value unexpected")
        self.assertEqual(otp2.name, mrrp.name, "xml: value unexpected")
        self.assertEqual(otp2.organType, mrrp.organType, "xml: value unexpected")
        self.assertEqual(otp2.subType, mrrp.subType, "xml: value unexpected")
        self.assertEqual(otp2.nob(), mrrp.nob(), "xml: value unexpected")  # value
        self.assertEqual(otp2.lns, mrrp.lns, "xml: value unexpected")  # dev
        for name in self.myc_parameters:
            self.assertAlmostEqual(getattr(otp2, name), getattr(mrrp, name), 7, "xml: value unexpected for " + name)
        for i in range(0, 3):
            self.assertEqual(otp2.successor[0][i], mrrp.successor[0][i], "xml: value unexpected")
        for i in range(0, 3):
            self.assertAlmostEqual(otp2.successorP[0][i], mrrp.successorP[0][i], 7, "xml: value unexpected")

    def test_realize(self):
        """ calls realize """
        self.mycroot_example()
        p = self.mrrp.realize()
        self.assertEqual(p.__class__.__name__, "MycorrhizalRootSpecificParameter", "realize: unexpected class type")
        self.assertEqual(p.subType, 1, "realize: unexpected sub type")
        self.assertEqual(p.a, 0.1, "realize: unexpected value")
        self.assertEqual(len(p.ln) + 1, self.mrrp.nob(), "realize: internodal distances +1 should be  number of laterals")

    def test_realize_colonization(self):
        """ tests the realized colonization parameters lmbd and vi """
        self.mycroot_example()
        self.mrrp.lmbds = 0.
        self.mrrp.vis = 0.
        p = self.mrrp.realize()
        self.assertEqual(p.lmbd, self.mrrp.lmbd, "realize: lmbd unexpected for zero deviation")
        self.assertEqual(p.vi, self.mrrp.vi, "realize: vi unexpected for zero deviation")
        self.mrrp.lmbd, self.mrrp.lmbds = -1., 0.1  # negative values are cut off at 0
        self.mrrp.vi, self.mrrp.vis = -1., 0.1
        p = self.mrrp.realize()
        self.assertEqual(p.lmbd, 0., "realize: lmbd should not be negative")
        self.assertEqual(p.vi, 0., "realize: vi should not be negative")

    def test_realize_deviation(self):
        """ tests that lmbd and vi are drawn with the given mean and standard deviation """
        self.mycroot_example()
        self.plant.setSeed(1)
        self.mrrp.lmbd, self.mrrp.lmbds = 0.5, 0.05
        self.mrrp.vi, self.mrrp.vis = 0.4, 0.04
        lmbd, vi = [], []
        for i in range(0, 2000):
            p = self.mrrp.realize()
            lmbd.append(p.lmbd)
            vi.append(p.vi)
        self.assertAlmostEqual(np.mean(lmbd), 0.5, 2, "realize: unexpected mean of lmbd")
        self.assertAlmostEqual(np.std(lmbd), 0.05, 2, "realize: unexpected standard deviation of lmbd")
        self.assertAlmostEqual(np.mean(vi), 0.4, 2, "realize: unexpected mean of vi")
        self.assertAlmostEqual(np.std(vi), 0.04, 2, "realize: unexpected standard deviation of vi")


if __name__ == '__main__':
    unittest.main()
