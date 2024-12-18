import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
from rsml.rsml_reader import *

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

    def test_constructors(self):
        """ tests constructor and copy """
        plant = pb.Organism()
        mrrp = pb.MycorrhizalRootRandomParameter(plant)
        mrrp.theta = 123
        mrrp.thetas = 456
        mrrp.gf = 789
        otp2 = mrrp.copy(plant)
        self.assertIsNot(mrrp, otp2, "copy: organ type parameter set is not copied")
        self.assertEqual(otp2.name, mrrp.name, "copy: value unexpected")
        self.assertEqual(otp2.organType, mrrp.organType, "copy: value unexpected")
        self.assertEqual(otp2.subType, mrrp.subType, "copy: value unexpected")
        self.assertEqual(otp2.a, mrrp.a, "copy: value unexpected")
        self.assertEqual(otp2.theta, mrrp.theta, "copy: value unexpected")
        self.assertEqual(otp2.thetas, mrrp.thetas, "copy: value unexpected")
        self.assertEqual(otp2.gf, mrrp.gf, "copy: value unexpected")
        self.assertEqual(otp2.p, mrrp.p,"copy: value unexpected")
        self.assertEqual(otp2.minAge, mrrp.minAge,"copy: value unexpected")
        self.assertEqual(otp2.maxAge, mrrp.maxAge,"copy: value unexpected")
        self.assertEqual(otp2.vi, mrrp.vi,"copy: value unexpected")
        self.assertEqual(otp2.infected, mrrp.infected,"copy: value unexpected")

    def test_parameter(self):
        """ tests getParameter() """
        self.plant = pb.Organism()
        mrrp = pb.MycorrhizalRootRandomParameter(self.plant)
        mrrp.lns = 0.123
        mrrp.la = 12
        ot = mrrp.getParameter("organType")  # test defaults
        st = mrrp.getParameter("subType")
        gf = mrrp.getParameter("gf")
        ln = mrrp.getParameter("ln")
        lns = mrrp.getParameter("ln_dev")
        la = mrrp.getParameter("la_mean")  # we can always add "_mean" to avoid naming conflicts
        self.assertEqual(ot, 2., "getParameter: value unexpected")
        self.assertEqual(st, -1., "getParameter: value unexpected")
        self.assertEqual(gf, 1., "getParameter: value unexpected")
        self.assertEqual(ln, 1., "getParameter: value unexpected")
        self.assertEqual(lns, 0.123, "getParameter: value unexpected")
        self.assertEqual(la, 12, "getParameter: value unexpected")
        mrrp.theta = 123  # change values
        mrrp.thetas = 456
        theta = mrrp.getParameter("theta")
        thetas = mrrp.getParameter("theta_dev")
        self.assertEqual(theta, 123, "getParameter: value unexpected")
        self.assertEqual(thetas, 456, "getParameter: value unexpected")

    def test_toString(self):
        """ tests __str__ output """
        self.mrrp = pb.MycorrhizalRootRandomParameter(pb.Organism())
        mrrp = self.mrrp  # rename
        mrrp.name = "the root"
        self.assertEqual(mrrp.__str__(False), "name: the root, organType: 2, subType: -1.", "toString: value unexpected")
        # print(mrrp)

    def test_xml(self):
        """ write the organ as xml, and rereads it """
        self.mycroot_example()
        mrrp = self.mrrp  # rename
        mrrp.name = "lateral"
        mrrp.subType = 2
        mrrp.writeXML("root.xml")
        otp2 = pb.MycorrhizalRootRandomParameter(self.plant)
        otp2.readXML("root.xml")
        self.assertEqual(otp2.ldelay, mrrp.ldelay, "xml: value unexpected")
        self.assertEqual(otp2.ldelays, mrrp.ldelays, "xml: value unexpected")
        self.assertEqual(otp2.name, mrrp.name, "xml: value unexpected")
        self.assertEqual(otp2.organType, mrrp.organType, "xml: value unexpected")
        self.assertEqual(otp2.subType, mrrp.subType, "xml: value unexpected")
        self.assertEqual(otp2.nob(), mrrp.nob(), "xml: value unexpected")  # value
        self.assertEqual(otp2.lns, mrrp.lns, "xml: value unexpected")  # dev
        self.assertEqual(otp2.p, mrrp.p,"xml: value unexpected")
        self.assertEqual(otp2.vi, mrrp.vi,"xml: value unexpected")
        self.assertEqual(otp2.maxAge, mrrp.maxAge,"xml: value unexpected")
        for i in range(0, 3):
            self.assertEqual(otp2.successor[0][i], mrrp.successor[0][i], "xml: value unexpected")
        for i in range(0, 3):
            self.assertAlmostEqual(otp2.successorP[0][i], mrrp.successorP[0][i], 7, "xml: value unexpected")

    def test_realize(self):
        """ calls realize """
        self.mycroot_example()
        p = self.mrrp.realize()
        self.assertEqual(p.__class__.__name__, "RootSpecificParameter", "realize: unexpected class type")
        self.assertEqual(p.subType, 1, "realize: unexpected sub type")
        self.assertEqual(p.a, 0.1, "realize: unexpected value")
        self.assertEqual(len(p.ln) + 1, self.mrrp.nob(), "realize: internodal distances +1 should be  number of laterals")


if __name__ == '__main__':
    unittest.main()
