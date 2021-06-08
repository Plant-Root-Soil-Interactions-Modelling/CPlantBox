import unittest
import sys; sys.path.append("..")import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb


class TestLeafParameter(unittest.TestCase):

    def leaf_example(self):
        self.plant = pb.Organism()
        self.lrp = pb.LeafRandomParameter(self.plant)
        self.lrp.la = 1.5
        self.lrp.lb = 6
        self.lrp.ln = 1.25
        self.lrp.lns = 0.12
        self.lrp.lmax = 3 * self.lrp.ln + self.lrp.la + self.lrp.lb
        self.lrp.subType = 1

    def add_successors(self):
        """ add successor sub types to the example"""
        self.lrp.successor = [4, 5, 6]
        self.lrp.successorP = [0.4, 0.1, 0.5]
        # print(self.srp.successorP[0])

    def test_constructors(self):
        """ tests constructor and copy """
        plant = pb.Organism()
        otp = pb.LeafRandomParameter(plant)
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
        lrp = pb.LeafRandomParameter(pb.Organism())
        lrp.lns = 0.123
        lrp.la = 12
        ot = lrp.getParameter("organType")  # test defaults
        st = lrp.getParameter("subType")
        gf = lrp.getParameter("gf")
        ln = lrp.getParameter("ln")
        lns = lrp.getParameter("ln_dev")
        la = lrp.getParameter("la_mean")  # we can always add "_mean" to avoid naming conflicts
        self.assertEqual(ot, 4., "getParameter: value unexpected")
        self.assertEqual(st, -1., "getParameter: value unexpected")
        self.assertEqual(gf, 1., "getParameter: value unexpected")
        self.assertEqual(ln, 1., "getParameter: value unexpected")
        self.assertEqual(lns, 0.123, "getParameter: value unexpected")
        self.assertEqual(la, 12, "getParameter: value unexpected")
        lrp.theta = 123  # change values
        lrp.thetas = 456
        theta = lrp.getParameter("theta")
        thetas = lrp.getParameter("theta_dev")
        self.assertEqual(theta, 123, "getParameter: value unexpected")
        self.assertEqual(thetas, 456, "getParameter: value unexpected")
        # print([theta,thetas])

    def test_toString(self):
        self.lrp = pb.LeafRandomParameter(pb.Organism())
        self.add_successors()
        lrp = self.lrp  # rename
        lrp.name = "the leaf"
        self.assertEqual(lrp.__str__(False), "name: the leaf, organType: 4, subType: -1.", "toString: value unexpected")
        # print(lrp.name)

    def test_xml(self):
        """ write the organ as xml, and rereads it """
        self.leaf_example()
        self.add_successors()
        otp = self.lrp  # rename
        otp.name = "lateral"
        otp.subType = 1
        otp.lmax = 17;
        otp.writeXML("leaf.xml")
        otp2 = pb.LeafRandomParameter(self.plant)
        otp2.readXML("leaf.xml")
        self.assertEqual(otp2.name, otp.name, "xml: value unexpected")
        self.assertEqual(otp2.organType, otp.organType, "xml: value unexpected")
        self.assertEqual(otp2.subType, otp.subType, "xml: value unexpected")
        self.assertEqual(otp2.lmax, otp.lmax, "xml: value unexpected")  # value
        self.assertEqual(otp2.lns, otp.lns, "xml: value unexpected")  # dev
        for i in range(0, 3):
            self.assertEqual(otp2.successor[i], otp.successor[i], "xml: value unexpected")
        for i in range(0, 3):
            self.assertAlmostEqual(otp2.successorP[i], otp.successorP[i], 7, "xml: value unexpected")

    def test_realize(self):
        """ calls realize """
        self.leaf_example()
        p = self.lrp.realize()
        self.assertEqual(p.__class__.__name__, "LeafSpecificParameter", "realize: unexpected class type")
        self.assertEqual(p.subType, 1, "realize: unexpected sub type")
        self.assertEqual(p.a, 0.1, "realize: unexpected value")
        self.assertEqual(len(p.ln) + 1, self.lrp.nob(), "realize: internodal distances +1 should be  number of laterals")
        # print(p)

    def test_radial_leaf_geometry(self):
        self.plant = pb.Organism()
        lrp = pb.LeafRandomParameter(self.plant)
        lrp.la = 3.5
        lrp.lb = 1.
        lrp.ln = 3.
        lrp.lmax = lrp.la + lrp.lb + lrp.ln  
        """ radial geometry """
        phi = np.array([-90, -45, 0., 45, 90]) / 180. * np.pi
        l = np.array([3, 2.2, 1.7, 2, 3.5])
        N = 15  # N is rather high for testing
        lrp.createLeafRadialGeometry(phi, l, N)  
        self.assertEqual(lrp.leafMid(), 3., "unexpected leaf mid")
        self.assertEqual(lrp.leafLength(), 6.5, "unexpected leaf length")        
        y_ = np.linspace(0, lrp.leafLength(), N)
        geom = lrp.getLeafGeometry()
        x_ = np.array([ x[-1] for x in geom])
        self.assertEqual(x_.shape[0], y_.shape[0], "leaf geometry has wrong size");
#         plt.plot(x_, y_, "g-*")
#         plt.plot(-x_, y_, "g-*")
#         plt.show()
    
    def test_leaf_geometry(self):
        self.plant = pb.Organism()
        lrp = pb.LeafRandomParameter(self.plant)
        lrp.la = 3.5
        lrp.lb = 1.
        lrp.ln = 3.
        lrp.lmax = lrp.la + lrp.lb + lrp.ln  
        """ radial geometry """
        phi = np.array([-3, -3 * 0.7, 0., 3.5 * 0.7, 3.5]) / 180. * np.pi
        l = np.array([0., 2.2 * 0.7, 1.7, 2.*0.7, 0.])
        N = 15  # N is rather high for testing
        lrp.createLeafGeometry(phi, l, N)  
        self.assertEqual(lrp.leafMid(), 3., "unexpected leaf mid")
        self.assertEqual(lrp.leafLength(), 6.5, "unexpected leaf length")        
        y_ = np.linspace(0, lrp.leafLength(), N)
        geom = lrp.getLeafGeometry()
        x_ = np.array([ x[-1] for x in geom])
        self.assertEqual(x_.shape[0], y_.shape[0], "leaf geometry has wrong size");
        plt.plot(x_, y_, "g-*")
        plt.plot(-x_, y_, "g-*")
        plt.show()        


if __name__ == '__main__':
    unittest.main()

