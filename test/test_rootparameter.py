import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
from rsml.rsml_reader import *


class TestStemParameter(unittest.TestCase):

    def root_example(self):
        """ example root parameters used below """
        self.plant = pb.Organism()
        self.rrp = pb.RootRandomParameter(self.plant)
        self.rrp.la = 1.5
        self.rrp.lb = 5.5
        self.rrp.ln = 1.25
        self.rrp.lns = 0.12
        self.rrp.lmax = 7 * self.rrp.ln + self.rrp.la + self.rrp.lb
        self.rrp.subType = 1
        self.rrp.successor = [[4, 5, 6]]
        self.rrp.successorP = [[0.4, 0.1, 0.5]]
        self.rrp.ldelay = 5
        self.rrp.ldelays = 2

    def test_constructors(self):
        """ tests constructor and copy """
        plant = pb.Organism()
        rrp = pb.RootRandomParameter(plant)
        rrp.theta = 123
        rrp.thetas = 456
        rrp.gf = 789
        otp2 = rrp.copy(plant)
        self.assertIsNot(rrp, otp2, "copy: organ type parameter set is not copied")
        self.assertEqual(otp2.name, rrp.name, "copy: value unexpected")
        self.assertEqual(otp2.organType, rrp.organType, "copy: value unexpected")
        self.assertEqual(otp2.subType, rrp.subType, "copy: value unexpected")
        self.assertEqual(otp2.a, rrp.a, "copy: value unexpected")
        self.assertEqual(otp2.theta, rrp.theta, "copy: value unexpected")
        self.assertEqual(otp2.thetas, rrp.thetas, "copy: value unexpected")
        self.assertEqual(otp2.gf, rrp.gf, "copy: value unexpected")

    def test_parameter(self):
        """ tests getParameter() """
        self.plant = pb.Organism()
        rrp = pb.RootRandomParameter(self.plant)
        rrp.lns = 0.123
        rrp.la = 12
        ot = rrp.getParameter("organType")  # test defaults
        st = rrp.getParameter("subType")
        gf = rrp.getParameter("gf")
        ln = rrp.getParameter("ln")
        lns = rrp.getParameter("ln_dev")
        la = rrp.getParameter("la_mean")  # we can always add "_mean" to avoid naming conflicts
        self.assertEqual(ot, 2., "getParameter: value unexpected")
        self.assertEqual(st, -1., "getParameter: value unexpected")
        self.assertEqual(gf, 1., "getParameter: value unexpected")
        self.assertEqual(ln, 1., "getParameter: value unexpected")
        self.assertEqual(lns, 0.123, "getParameter: value unexpected")
        self.assertEqual(la, 12, "getParameter: value unexpected")
        rrp.theta = 123  # change values
        rrp.thetas = 456
        theta = rrp.getParameter("theta")
        thetas = rrp.getParameter("theta_dev")
        self.assertEqual(theta, 123, "getParameter: value unexpected")
        self.assertEqual(thetas, 456, "getParameter: value unexpected")

    def test_toString(self):
        """ tests __str__ output """
        self.rrp = pb.RootRandomParameter(pb.Organism())
        rrp = self.rrp  # rename
        rrp.name = "the root"
        self.assertEqual(rrp.__str__(False), "name: the root, organType: 2, subType: -1.", "toString: value unexpected")
        # print(rrp)

    def test_xml(self):
        """ write the organ as xml, and rereads it """
        self.root_example()
        rrp = self.rrp  # rename
        rrp.name = "lateral"
        rrp.subType = 2
        rrp.writeXML("root.xml")
        otp2 = pb.RootRandomParameter(self.plant)
        otp2.readXML("root.xml")
        self.assertEqual(otp2.ldelay, rrp.ldelay, "xml: value unexpected")
        self.assertEqual(otp2.ldelays, rrp.ldelays, "xml: value unexpected")
        self.assertEqual(otp2.name, rrp.name, "xml: value unexpected")
        self.assertEqual(otp2.organType, rrp.organType, "xml: value unexpected")
        self.assertEqual(otp2.subType, rrp.subType, "xml: value unexpected")
        self.assertEqual(otp2.nob(), rrp.nob(), "xml: value unexpected")  # value
        self.assertEqual(otp2.lns, rrp.lns, "xml: value unexpected")  # dev
        for i in range(0, 3):
            self.assertEqual(otp2.successor[0][i], rrp.successor[0][i], "xml: value unexpected")
        for i in range(0, 3):
            self.assertAlmostEqual(otp2.successorP[0][i], rrp.successorP[0][i], 7, "xml: value unexpected")

    def test_realize(self):
        """ calls realize """
        self.root_example()
        p = self.rrp.realize()
        self.assertEqual(p.__class__.__name__, "RootSpecificParameter", "realize: unexpected class type")
        self.assertEqual(p.subType, 1, "realize: unexpected sub type")
        self.assertEqual(p.a, 0.1, "realize: unexpected value")
        self.assertEqual(len(p.ln) + 1, self.rrp.nob(), "realize: internodal distances +1 should be  number of laterals")


if __name__ == '__main__':
    unittest.main()
