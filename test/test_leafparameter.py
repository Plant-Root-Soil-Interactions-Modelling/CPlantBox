import sys; sys.path.append("..")
import unittest

import plantbox as pb

import numpy as np
import matplotlib.pyplot as plt


class TestLeafParameter(unittest.TestCase):

    def leaf_example(self):
        self.plant = pb.Organism()
        self.lrp = pb.LeafRandomParameter(self.plant)
        self.lrp.la = 1.5
        self.lrp.lb = 6
        self.lrp.ln = 1.25
        self.lrp.lns = 0.12
        self.lrp.lmax = 3 * self.lrp.ln + self.lrp.la + self.lrp.lb
        self.lrp.areaMax = 10
        self.lrp.subType = 1

    def add_successors(self):
        """ add successor sub types to the example"""
        self.lrp.successor = [[4, 5, 6]]
        self.lrp.successorP = [[0.4, 0.1, 0.5]]
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
        """ tests toString() """
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
            self.assertEqual(otp2.successor[0][i], otp.successor[0][i], "xml: value unexpected")
        for i in range(0, 3):
            self.assertAlmostEqual(otp2.successorP[0][i], otp.successorP[0][i], 7, "xml: value unexpected")

    def test_realize(self):
        """ calls realize """
        self.leaf_example()
        self.add_successors()
        p = self.lrp.realize()
        self.assertEqual(p.__class__.__name__, "LeafSpecificParameter", "realize: unexpected class type")
        self.assertEqual(p.subType, 1, "realize: unexpected sub type")
        self.assertEqual(p.a, 0.1, "realize: unexpected value for radius")
        self.assertEqual(len(p.ln) + 1, self.lrp.nob(), "realize: internodal distances +1 should be  number of laterals")
        self.assertEqual(p.leafArea, 10., "realize: unexpected value for leaf area")
        # print(p)

    def test_radial_leaf_geometry(self):
        """ tests if a radial leaf geometry can be set"""

        self.plant = pb.Organism()
        lrp = pb.LeafRandomParameter(self.plant)
        lrp.la = 3.5
        lrp.lb = 1.
        lrp.ln = 3.
        lrp.lmax = lrp.la + lrp.lb + lrp.ln
        lrp.areaMax = 50
#         """ radial geometry """
#         phi = np.array([-90, -45, 0., 45, 90]) / 180. * np.pi
#         l = np.array([3, 2.2, 1.7, 2, 3.5])
#         N = 105  # N is rather high for testing
#         lrp.createLeafRadialGeometry(phi, l, N)

        lrp.la, lrp.lb, lrp.lmax, lrp.ln, lrp.r, lrp.dx = 5, 1, 11, 5, 1, 0.5
        phi = np.array([-90., -67.5, -45, -22.5, 0, 22.5, 45, 67.5, 90]) / 180. * np.pi
        l = np.array([2 for x_i in range(len(phi))])  # l = np.array([5., 1, 5, 1, 5, 1, 5, 1, 5])
        assert(l.shape == phi.shape)
        N = 500  # N is rather high for testing
        lrp.createLeafRadialGeometry(phi, l, N)

        # self.assertEqual(lrp.leafMid(), 3., "unexpected leaf mid")
        # self.assertEqual(lrp.leafLength(), 6.5, "unexpected leaf length")
        yy = np.linspace(0, lrp.leafLength(), N)
        geom = lrp.leafGeometry
        x_, y_ = [], []
        for i, x in enumerate(geom):
#             if len(x) > 1:
#                 print(len(x), x, [yy[i]] * len(x))
#             if len(x) == 1:
#                 print(len(x), x, [yy[i]] * len(x))
            x_.extend(x)
            y_.extend([yy[i]] * len(x))
        x_ = np.array(x_)
        y_ = np.array(y_)

        # self.assertEqual(x_.shape[0], y_.shape[0], "leaf geometry has wrong size");
        a = lrp.areaMax / lrp.leafLength()
        # self.assertAlmostEqual(2 * np.sum(x_) / N * a * lrp.leafLength(), lrp.areaMax , 2, "xml: value unexpected")
        # plt.plot(x_ * a, y_, "g*")
        # plt.plot(-x_ * a, y_, "g*")
        # plt.ylim([0, 10])
        # plt.xlim([-7.5, 7.5])
        # plt.show()

    def test_leaf_geometry(self):
        """ tests if a leaf geometry can be set"""
        self.plant = pb.Organism()
        lrp = pb.LeafRandomParameter(self.plant)
        lrp.la = 3.5
        lrp.lb = 1.
        lrp.ln = 3.
        lrp.lmax = lrp.la + lrp.lb + lrp.ln
        lrp.areaMax = 10
        """ radial geometry """
        y = np.array([-3, -3 * 0.7, 0., 3.5 * 0.7, 3.5])
        l = np.array([0., 2.2 * 0.7, 1.7, 1.8 * 0.7, 0.])
        N = 105  # N is rather high for testing
        lrp.createLeafGeometry(y, l, N)
        self.assertEqual(lrp.leafMid(), 3., "unexpected leaf mid")
        self.assertEqual(lrp.leafLength(), 6.5, "unexpected leaf length")
        y_ = np.linspace(0, lrp.leafLength(), N)
        geom = lrp.leafGeometry
        x_ = np.array([ x[-1] for x in geom])
        self.assertEqual(x_.shape[0], y_.shape[0], "leaf geometry has wrong size");
        a = lrp.areaMax / lrp.leafLength()
        self.assertAlmostEqual(2 * np.sum(x_) / N * a * lrp.leafLength(), lrp.areaMax , 2, "xml: value unexpected")
        # plt.plot(x_ * a, y_, "g-*")
        # plt.plot(-x_ * a, y_, "g-*")
        # plt.ylim([0, 7])
        # plt.xlim([-2, 2])
        # plt.show()


if __name__ == '__main__':
    unittest.main()

