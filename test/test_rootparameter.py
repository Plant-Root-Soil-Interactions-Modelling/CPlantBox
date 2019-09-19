import unittest
import sys
sys.path.append("..")
import rootbox as rb
from rsml import *


class TestStemParameter(unittest.TestCase):

    def root_example(self):
        """ example root parameters used below """
        self.plant = rb.Organism()
        self.rtp = rb.StemRandomParameter(self.plant)
        self.rtp.la = 1.5
        self.rtp.lb = 5.5
        self.rtp.ln = 1.25
        self.rtp.lns = 0.12
        self.rtp.nob = 8
        self.rtp.subType = 1

    def add_successors(self):
        """ add successor sub types to the example"""
        l = rb.std_vector_int_()
        l.append(4)
        l.append(5)
        l.append(6)
        self.rtp.successor = l
        l = rb.std_vector_double_()
        l.append(0.4)
        l.append(0.1)
        l.append(0.5)
        self.rtp.successorP = l

    def test_constructors(self):
        """ tests constructor and copy """
        plant = rb.Organism()
        otp = rb.StemRandomParameter(plant)
        otp.theta = 123
        otp.thetas = 456
        otp.gf = 789
        otp2 = otp.copy(plant)
        self.assertIsNot(otp, otp2, "copy: organ type parameter set is not copied")
        self.assertEqual(otp2.name, otp.name, "copy: value unexpected")
        self.assertEqual(otp2.organType, otp.organType, "copy: value unexpected")
        self.assertEqual(otp2.subType, otp.subType, "copy: value unexpected")
        self.assertEqual(otp2.a, otp.a, "copy: value unexpected")
        self.assertEqual(otp2.theta, otp.theta, "copy: value unexpected")
        self.assertEqual(otp2.thetas, otp.thetas, "copy: value unexpected")
        self.assertEqual(otp2.gf, otp.gf, "copy: value unexpected")

    def test_parameter(self):
        """ tests getParameter() """
        rtp = rb.StemRandomParameter(rb.Organism())
        rtp.lns = 0.123
        rtp.la = 12
        ot = rtp.getParameter("organType")  # test defaults
        st = rtp.getParameter("subType")
        gf = rtp.getParameter("gf")
        ln = rtp.getParameter("ln")
        lns = rtp.getParameter("ln_dev")
        la = rtp.getParameter("la_mean")  # we can always add "_mean" to avoid naming conflicts
        self.assertEqual(ot, 3., "getParameter: value unexpected")
        self.assertEqual(st, -1., "getParameter: value unexpected")
        self.assertEqual(gf, 1., "getParameter: value unexpected")
        self.assertEqual(ln, 1., "getParameter: value unexpected")
        self.assertEqual(lns, 0.123, "getParameter: value unexpected")
        self.assertEqual(la, 12, "getParameter: value unexpected")
        rtp.theta = 123  # change values
        rtp.thetas = 456
        theta = rtp.getParameter("theta")
        thetas = rtp.getParameter("theta_dev")
        self.assertEqual(theta, 123, "getParameter: value unexpected")
        self.assertEqual(thetas, 456, "getParameter: value unexpected")

    def test_toString(self):
        """ tests __str__ output """
        self.rtp = rb.StemRandomParameter(rb.Organism())
        self.add_successors()
        rtp = self.rtp  # rename
        rtp.name = "the root"
        self.assertEqual(rtp.__str__(False), "Name: the root, organType: 3, subType, -1", "toString: value unexpected")
        # print(rtp)

    def test_xml(self):
        """ write the organ as xml, and rereads it """
        self.root_example()
        self.add_successors()
        otp = self.rtp  # rename
        otp.name = "lateral"
        otp.subType = 2
        otp.writeXML("root.xml")
        otp2 = rb.StemRandomParameter(self.plant)
        otp2.readXML("root.xml")
        self.assertEqual(otp2.name, otp.name, "xml: value unexpected")
        self.assertEqual(otp2.organType, otp.organType, "xml: value unexpected")
        self.assertEqual(otp2.subType, otp.subType, "xml: value unexpected")
        self.assertEqual(otp2.nob, otp.nob, "xml: value unexpected")  # value
        self.assertEqual(otp2.lns, otp.lns, "xml: value unexpected")  # dev
        for i in range(0, 3):
            self.assertEqual(otp2.successor[i], otp.successor[i], "xml: value unexpected")
        for i in range(0, 3):
            self.assertAlmostEqual(otp2.successorP[i], otp.successorP[i], 7, "xml: value unexpected")

    def test_realize(self):
        """ calls realize """
        self.root_example()
        p = self.rtp.realize()
        self.assertEqual(p.__class__.__name__, "StemSpecificParameter", "realize: unexpected class type")
        self.assertEqual(p.subType, 1, "realize: unexpected sub type")
        self.assertEqual(p.a, 0.1, "realize: unexpected value")
        self.assertEqual(len(p.ln) + 1, self.rtp.nob, "realize: internodal distances +1 should be  number of laterals")


if __name__ == '__main__':
    unittest.main()
