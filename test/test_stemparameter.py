import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
from rsml.rsml_reader import *


class TestStemParameter(unittest.TestCase):

    def stem_example(self):
        self.plant = pb.Organism()
        self.srp = pb.StemRandomParameter(self.plant)
        self.srp.la = 1.5
        self.srp.lb = 5.5
        self.srp.ln = 1.25
        self.srp.lns = 0.12
        self.srp.lmax = 7 * self.srp.ln + self.srp.la + self.srp.lb
        # print(self.srp.las)

        self.srp.successor = [[4, 5, 6]]
        self.srp.successorP = [[0.4, 0.1, 0.5]]

        # print(self.srp.successorP[0])

    def test_constructors(self):
        """ tests constructor and copy """
        plant = pb.Organism()
        otp = pb.StemRandomParameter(plant)
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
        # print(self.assertEqual(otp2.name, otp.name, "copy: value unexpected"))

    def test_parameter(self):
        """ tests getParameter() """
        self.plant = pb.Organism()
        srp = pb.StemRandomParameter(self.plant)
        srp.lns = 0.123
        srp.la = 12
        ot = srp.getParameter("organType")  # test defaults
        st = srp.getParameter("subType")
        gf = srp.getParameter("gf")
        ln = srp.getParameter("ln")
        lns = srp.getParameter("ln_dev")
        la = srp.getParameter("la_mean")  # we can always add "_mean" to avoid naming conflicts
        self.assertEqual(ot, 3., "getParameter: value unexpected")
        self.assertEqual(st, -1., "getParameter: value unexpected")
        self.assertEqual(gf, 1., "getParameter: value unexpected")
        self.assertEqual(ln, 1., "getParameter: value unexpected")
        self.assertEqual(lns, 0.123, "getParameter: value unexpected")
        self.assertEqual(la, 12, "getParameter: value unexpected")
        srp.theta = 123  # change values
        srp.thetas = 456
        theta = srp.getParameter("theta")
        thetas = srp.getParameter("theta_dev")
        self.assertEqual(theta, 123, "getParameter: value unexpected")
        self.assertEqual(thetas, 456, "getParameter: value unexpected")
        # print([theta,thetas])

    def test_toString(self):
        self.srp = pb.StemRandomParameter(pb.Organism())

        srp = self.srp  # rename
        srp.name = "the stem"

        self.assertEqual(srp.__str__(False), "name: the stem, organType: 3, subType: -1.", "toString: value unexpected")

    def test_xml(self):
        """ write the organ as xml, and rereads it """
        self.stem_example()

        otp = self.srp  # rename
        otp.name = "lateral"
        otp.subType = 1

        otp.writeXML("stem.xml")
        otp2 = pb.StemRandomParameter(self.plant)
        otp2.readXML("stem.xml")
        self.assertEqual(otp2.name, otp.name, "xml: value unexpected")
        self.assertEqual(otp2.organType, otp.organType, "xml: value unexpected")

        self.assertEqual(otp2.subType, otp.subType, "xml: value unexpected")
        self.assertEqual(otp2.lmax, otp.lmax, "xml: value unexpected")  # value
        self.assertEqual(otp2.nob(), otp.nob(), "xml: value unexpected")  # value
        self.assertEqual(otp2.lns, otp.lns, "xml: value unexpected")  # dev
        for i in range(0, 3):
            self.assertEqual(otp2.successor[0][i], otp.successor[0][i], "xml: value unexpected")
        for i in range(0, 3):
            self.assertAlmostEqual(otp2.successorP[0][i], otp.successorP[0][i], 7, "xml: value unexpected")

    def test_realize(self):
        """ calls realize """
        self.stem_example()
        p = self.srp.realize()
        self.assertEqual(p.__class__.__name__, "StemSpecificParameter", "realize: unexpected class type")
        self.assertEqual(p.subType, -1, "realize: unexpected sub type")
        self.assertEqual(p.a, 0.1, "realize: unexpected value")
        self.assertEqual(len(p.ln) + 1 , self.srp.nob(), "realize: internodal distances +1 should be  number of laterals")
        # print(p)


if __name__ == '__main__':
    unittest.main()

