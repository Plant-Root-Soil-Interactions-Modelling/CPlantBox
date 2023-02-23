import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
from rsml.rsml_reader import *


class TestOrganParameter(unittest.TestCase):

    def test_constructors(self):
        """ tests constructor and copy """
        plant = pb.Organism()
        orp = pb.OrganRandomParameter(plant)
        self.assertEqual(orp.name, "organ", "copy: default value unexpected")
        self.assertEqual(orp.organType, 0, "copy: default value unexpected")
        self.assertEqual(orp.subType, 0, "copy: default value unexpected")
        otp2 = orp.copy(plant)
        self.assertIsNot(orp, otp2, "copy: organ type parameter set is not copied")
        self.assertEqual(otp2.name, orp.name, "copy: value unexpected")
        self.assertEqual(otp2.organType, orp.organType, "copy: value unexpected")
        self.assertEqual(otp2.subType, orp.subType, "copy: value unexpected")

    def test_parameter(self):
        """ tests getParameter() """
        orp = pb.OrganRandomParameter(pb.Organism())
        orp.organType = 1
        orp.subType = 2
        ot = orp.getParameter("organType")
        st = orp.getParameter("subType")
        self.assertEqual(ot, 1., "getParameter: value unexpected")
        self.assertEqual(st, 2., "getParameter: value unexpected")

    def test_toString(self):
        """ tests __str__ output """
        orp = pb.OrganRandomParameter(pb.Organism())
        orp.name = "the great organ"
        orp.organType = 1
        orp.subType = 2
        self.assertEqual(orp.__str__(False), "name: the great organ, organType: 1, subType: 2.", "toString: value unexpected")
        # print(orp)

    def test_xml(self):
        """ write the organ type parameter as xml, and rereads it """
        plant = pb.Organism()
        orp = pb.OrganRandomParameter(plant)
        orp.name = "the great organ"
        orp.subType = 2
        orp.writeXML("organ.xml")
        otp2 = pb.OrganRandomParameter(plant)
        otp2.readXML("organ.xml")
        self.assertEqual(otp2.name, orp.name, "xml: value unexpected")
        self.assertEqual(otp2.subType, orp.subType, "xml: value unexpected")

    def test_realize(self):
        """ calls realize """
        o = pb.Organism()
        orp = pb.OrganRandomParameter(o)
        orp.subType = 2
        p = orp.realize()
        self.assertEqual(p.__class__.__name__, "OrganSpecificParameter", "realize: unexpected class type")
        self.assertEqual(p.subType, 2, "realize: unexpected sub type")


if __name__ == '__main__':
    unittest.main()
