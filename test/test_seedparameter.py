import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
from rsml.rsml_reader import *

import matplotlib.pyplot as plt


class TestShootParameter(unittest.TestCase):

    def shoot_example(self):
        """ example seed parameters used below """
        self.plant = pb.Organism()
        self.srp = pb.SeedRandomParameter(self.plant)
        self.srp.firstB = 10
        self.srp.delayB = 7
        self.srp.nC = 5
        self.srp.firstSB = 21
        self.srp.delaySB = 4
        self.srp.delayRC = 15
        self.srp.nz = 1.
        self.srp.subType = 0

    def test_constructors(self):
        """ tests constructor and copy """
        plant = pb.Organism()
        srp = pb.SeedRandomParameter(plant)
        srp.firstB = 123
        srp.firstBs = 456
        srp.nz = 789
        srp2 = srp.copy(plant)
        self.assertIsNot(srp, srp2, "copy: organ type parameter set is not copied")
        self.assertEqual(srp2.name, srp.name, "copy: value unexpected")
        self.assertEqual(srp2.organType, srp.organType, "copy: value unexpected")
        self.assertEqual(srp2.subType, srp.subType, "copy: value unexpected")
        self.assertEqual(srp2.firstB, srp.firstB, "copy: value unexpected")
        self.assertEqual(srp2.firstBs, srp.firstBs, "copy: value unexpected")
        self.assertEqual(srp2.nz, srp.nz, "copy: value unexpected")

    def test_parameter(self):
        """ tests getParameter() """
        plant = pb.Organism()
        srp = pb.SeedRandomParameter(plant)
        srp.firstB = 0.123
        srp.delaySBs = 12
        ot = srp.getParameter("organType")  # test defaults
        st = srp.getParameter("subType")
        nz = srp.getParameter("nz")
        firstB = srp.getParameter("firstB")
        delaySBs = srp.getParameter("delaySB_dev")
        firstB2 = srp.getParameter("firstB_mean")  # we can always add "_mean" to avoid naming conflicts
        self.assertEqual(ot, 1., "getParameter: value unexpected")
        self.assertEqual(st, 0, "getParameter: value unexpected")
        self.assertEqual(nz, 1., "getParameter: value unexpected")
        self.assertEqual(firstB, 0.123, "getParameter: value unexpected")
        self.assertEqual(firstB, firstB2, "getParameter: value unexpected")
        self.assertEqual(delaySBs, 12, "getParameter: value unexpected")
        srp.firstB = 123  # change values
        srp.delaySBs = 456
        firstB = srp.getParameter("firstB")
        delaySBs = srp.getParameter("delaySB_dev")
        self.assertEqual(firstB, 123, "getParameter: value unexpected")
        self.assertEqual(delaySBs, 456, "getParameter: value unexpected")

    def test_toString(self):
        """ tests __str__ output """
        self.srp = pb.SeedRandomParameter(pb.Organism())
        srp = self.srp  # rename
        srp.name = "the seed"
        # print(srp.__str__(False))
        # print(srp)
        self.assertEqual(srp.__str__(False), "name: the seed, organType: 1, subType: 0.", "toString: value unexpected")

    def test_xml(self):
        """ write the organ as xml, and rereads it """
        self.shoot_example()
        otp = self.srp  # rename
        otp.name = "best_seed"
        otp.subType = 0
        otp.writeXML("seed.xml")
        otp2 = pb.SeedRandomParameter(self.plant)
        otp2.readXML("seed.xml")
        self.assertEqual(otp2.name, otp.name, "xml: value unexpected")
        self.assertEqual(otp2.organType, otp.organType, "xml: value unexpected")
        self.assertEqual(otp2.subType, otp.subType, "xml: value unexpected")
        self.assertEqual(otp2.firstB, otp.firstB, "xml: value unexpected")  # value
        self.assertEqual(otp2.delayRC, otp.delayRC, "xml: value unexpected")  # value
        self.assertEqual(otp2.firstBs, otp.firstBs, "xml: value unexpected")  # dev
        self.assertEqual(otp2.delayRCs, otp.delayRCs, "xml: value unexpected")  # dev

    def test_realize(self):
        """ calls realize """
        self.shoot_example()
        p = self.srp.realize()
        self.assertEqual(p.__class__.__name__, "SeedSpecificParameter", "realize: unexpected class type")
        self.assertEqual(p.subType, 0, "realize: unexpected sub type")
        self.assertEqual(p.firstB, 10, "realize: unexpected value")


if __name__ == '__main__':
    unittest.main()
