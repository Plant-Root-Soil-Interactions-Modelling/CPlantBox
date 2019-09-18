import unittest
import py_rootbox as rb
from rsml import *


class TestOrganParameter(unittest.TestCase):

    def test_constructors(self):
        """ tests constructor and copy """
        plant = rb.Organism()
        otp = rb.OrganRandomParameter(plant)
        self.assertEqual(otp.name, "organ", "copy: default value unexpected")
        self.assertEqual(otp.organType, 0, "copy: default value unexpected")
        self.assertEqual(otp.subType, 0, "copy: default value unexpected")
        otp2 = otp.copy(plant)
        self.assertIsNot(otp, otp2, "copy: organ type parameter set is not copied")
        self.assertEqual(otp2.name, otp.name, "copy: value unexpected")
        self.assertEqual(otp2.organType, otp.organType, "copy: value unexpected")
        self.assertEqual(otp2.subType, otp.subType, "copy: value unexpected")

    def test_parameter(self):
        """ tests getParameter() """
        otp = rb.OrganRandomParameter(rb.Organism())
        otp.organType = 1
        otp.subType = 2
        ot = otp.getParameter("organType")
        st = otp.getParameter("subType")
        self.assertEqual(ot, 1., "getParameter: value unexpected")
        self.assertEqual(st, 2., "getParameter: value unexpected")

    def test_toString(self):
        """ tests __str__ output """
        otp = rb.OrganRandomParameter(rb.Organism())
        otp.name = "the great organ"
        otp.organType = 1
        otp.subType = 2
        self.assertEqual(otp.__str__(False), "Name: the great organ, organType: 1, subType, 2", "toString: value unexpected")
        # print(otp)

    def test_xml(self):
        """ write the organ type parameter as xml, and rereads it """
        plant = rb.Organism()
        otp = rb.OrganRandomParameter(plant)
        otp.name = "the great organ"
        otp.organType = 1  # seed
        otp.subType = 2
        otp.writeXML("organ.xml")
        otp2 = rb.OrganRandomParameter(plant)
        otp2.readXML("organ.xml")
        self.assertEqual(otp2.name, otp.name, "xml: value unexpected")
        self.assertEqual(otp2.organType, otp.organType, "xml: value unexpected")
        self.assertEqual(otp2.subType, otp.subType, "xml: value unexpected")

    def test_realize(self):
        """ calls realize """
        otp = rb.OrganRandomParameter(rb.Organism())
        otp.subType = 2
        p = otp.realize();
        self.assertEqual(p.__class__.__name__, "OrganSpecificParameter", "realize: unexpected class type")
        self.assertEqual(p.subType, 2, "realize: unexpected sub type")


if __name__ == '__main__':
    unittest.main()
