import unittest
import ../rootbox as rb
from rsml import *


class TestStemParameter(unittest.TestCase):

    def stem_example(self):
        self.plant = rb.Organism()
        self.srp = rb.StemRandomParameter(self.plant)
        self.srp.la = 1.5
        self.srp.lb = 5.5
        self.srp.ln = 1.25
        self.srp.lns = 0.12
        self.srp.nob = 8
        #print(self.srp.las)

    def add_successors(self):
        """ add successor sub types to the example"""
        l = rb.std_vector_int_()
        l.append(4)
        l.append(5)
        l.append(6)
        self.srp.successor = l
        l = rb.std_vector_double_()
        l.append(0.4)
        l.append(0.1)
        l.append(0.5)
        self.srp.successorP = l
        #print(self.srp.successorP[0])

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
        #print(self.assertEqual(otp2.name, otp.name, "copy: value unexpected"))

    def test_parameter(self):
        """ tests getParameter() """
        srp = rb.StemRandomParameter(rb.Organism())
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
        #print([theta,thetas])

    def test_toString(self):
        self.srp = rb.StemRandomParameter(rb.Organism())
        self.add_successors()
        srp = self.srp  # rename
        srp.name = "the stem"
        self.assertEqual(srp.__str__(False), "Name: the stem, organType: 3, subType, -1", "toString: value unexpected")
        print(srp.name)

    def test_xml(self):
        """ write the organ as xml, and rereads it """
        self.stem_example()
        self.add_successors()
        otp = self.srp  # rename
        otp.name = "lateral"
        otp.subType = 1
        otp.nob = (otp.k-otp.la-otp.lb)/otp.ln + 1;
        otp.writeXML("stem.xml")
        otp2 = rb.StemRandomParameter(self.plant)
        otp2.readXML("stem.xml")
        otp2.nob = (otp2.k-otp2.la-otp2.lb)/otp2.ln + 1;
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
        self.stem_example()
        p = self.srp.realize()
        self.assertEqual(p.__class__.__name__, "StemSpecificParameter", "realize: unexpected class type")
        self.assertEqual(p.subType, -1, "realize: unexpected sub type")
        self.assertEqual(p.a, 0.1, "realize: unexpected value")
        self.assertEqual(len(p.ln) + 1, self.srp.nob, "realize: internodal distances +1 should be  number of laterals")
        print(p)


if __name__ == '__main__':
    unittest.main()



