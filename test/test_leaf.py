"""
Copyright 2019, Forschungszentrum Jülich GmbH, licensed under GNU GPLv3
@authors Daniel Leitner, Andrea Schnepf, Xiaoran Zhou
"""

"""
TODO length, area (with no laterals)
TODO lenght, area (with laterals)
"""
import sys; sys.path.append("..")
import unittest

import plantbox as pb

import numpy as np
import matplotlib.pyplot as plt


def leafAge(l, r, k):  # leaf age at a certain length
    return -np.log(1 - l / k) * k / r


def leafLength(t, r, k):  # leaf length at a certain age
    return k * (1 - np.exp(-r * t / k))


def leafLateralLength(t, et, r, k):  # length of first order laterals (without second order laterals)
    i, l = 0, 0
    while et[i] < t:
        age = t - et[i]
        l += leafLength(age, r, k)
        i += 1
    return l


class TestLeaf(unittest.TestCase):

    def test_leaf(self):
        """ leaf without lateral leafs """
        ons = pb.Vector3d(0., 0., 1.)
        plant = pb.Plant()  # store organism (not owned by Organ, or OrganRandomParameter)
        p0 = pb.LeafRandomParameter(plant)
        p0.name, p0.subType, p0.la, p0.lb, p0.lmax, p0.ln, p0.r, p0.dx = "leaf", 1, 3.5, 1., 7.5, 3, 1, 0.1

        phi = np.array([-90, -45, 0., 45, 90]) / 180. * np.pi
        l = np.array([3, 2.2, 1.7, 2, 3.5])
        N = 100  # N is rather high for testing
        p0.createLeafRadialGeometry(phi, l, N)
#         y = np.array([-3, -3 * 0.7, 0., 3.5 * 0.7, 3.5])
#         l = np.array([0., 2.2 * 0.7, 1.7, 1.8 * 0.7, 0.])
#         N = 105  # N is rather high for testing
#         p0.createLeafGeometry(y, l, N)

        plant.setOrganRandomParameter(p0)  # the organism manages the type parameters and takes ownership
        plant.setOrganRandomParameter(pb.SeedRandomParameter(plant))
        # because we cannot pass a nullptr to pb.Leaf(...) L48
        param0 = p0.realize()  # set up leaf by hand (without a leaf syleaf)
        param0.la, param0.lb = 0, 0  # its important parent has zero length, otherwise creation times are messed up
        parentleaf = pb.Leaf(1, param0, True, True, 0., 0., ons, 0, False, 0)  # takes ownership of param0
        parentleaf.setOrganism(plant)
        parentleaf.addNode(pb.Vector3d(0, 0, -3), 0)  # there is no nullptr in Python

        leaf = pb.Leaf(plant, p0.subType, 0, parentleaf , 0)
        leaf.setOrganism(plant)
        leaf.simulate(7)
        # vp.plot_leaf(leaf)

#         nodes = []
#         nodes_leaf = []
#         for i in range(0, leaf.getNumberOfNodes()):
#             n = leaf.getNode(i)
#             nodes.append(np.array([n.x, n.y, n.z]))
#             n2 = leaf.getLeafVis(i)
#             for n_ in n2:
#                 print(n_)
#                 nodes_leaf.append(np.array([n_.x, n_.y, n_.z]))
#         nodes_leaf = np.array(nodes_leaf)
#         nodes = np.array(nodes)
#         print(nodes)
#         fig = plt.figure()
#         ax = fig.add_subplot(projection='3d')
#         ax.scatter(nodes_leaf[:, 0], nodes_leaf[:, 1], nodes_leaf[:, 2], marker='v')
#         ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2], marker='o')
#         ax.set_xlabel('X Label')
#         ax.set_ylabel('Y Label')
#         ax.set_zlabel('Z Label')
#         plt.show()

#     def leaf_example_rtp(self):
#         """ an example used in the tests below, a main leaf with laterals """
#         self.plant = pb.Organism()  # store organism (not owned by Organ, or OrganRandomParameter)
#         p0 = pb.LeafRandomParameter(self.plant)
#         p0.name, p0.subType, p0.la, p0.lb, p0.lmax, p0.ln, p0.r, p0.dx = "main", 1, 1, 10, 100, (89. / 19.), 1, 0.5
#         p0.successor = [3]
#         p0.successorP = [1.]
#         p1 = pb.LeafRandomParameter(self.plant)
#         p1.name, p1.subType, p1.la, p1.ln, p1.r, p1.dx = "lateral", 3, 25, 0, 2, 0.1
#         self.p0, self.p1 = p0, p1  # needed at later point
#         self.plant.setOrganRandomParameter(p0)  # the organism manages the type parameters and takes ownership
#         self.plant.setOrganRandomParameter(p1)
#         # TODO (first node is not set, if seed is used)
#         self.plant.setOrganRandomParameter(pb.SeedRandomParameter(self.plant))
#         self.seed = pb.Seed(self.plant)  # store parent (not owned by child Organ)
#         #
#         param0 = p0.realize()  # set up leaf by hand (without a leaf syleaf)
#         param0.la, param0.lb = 0, 0  # its important parent has zero length, otherwise creation times are messed up
#         parentleaf = pb.Leaf(1, param0, True, True, 0., 0., pb.Vector3d(0, 0, -1), 0, 0, False, 0)  # takes ownership of param0
#         parentleaf.setOrganism(self.plant)
#         parentleaf.addNode(pb.Vector3d(0, 0, -3), 0)  # there is no nullptr in Python
#         self.parentleaf = parentleaf  # store parent (not owned by child Organ)
#         #
#         self.leaf = pb.Leaf(self.plant, p0.subType, pb.Vector3d(0, 0, -1), 0, self.parentleaf , 0, 0)
#         self.leaf.setOrganism(self.plant)
#
#     def leaf_length_test(self, dt, l, subDt):
#         """ simulates a single leaf and checks length against analytic length """
#         nl, nl2, non, meanDX = [], [], [], []
#         for t in dt:
#             for i in range(0, subDt):
#
#                 self.leaf.simulate(t / subDt)
#             nl.append(self.leaf.getParameter("length"))
#             non.append(self.leaf.getNumberOfNodes())
#             meanDX.append(nl[-1] / non[-1])
#
#             # length from geometry
#             poly = np.zeros((non[-1], 3))  #
#             for i in range(0, non[-1]):
#                 v = self.leaf.getNode(i)
#                 poly[i, 0] = v.x
#                 poly[i, 1] = v.y
#                 poly[i, 2] = v.z
#             d = np.diff(poly, axis = 0)
#             sd = np.sqrt((d ** 2).sum(axis = 1))
#             nl2.append(sum(sd))
#         for i in range(0, len(dt)):
#             self.assertAlmostEqual(l[i], nl[i], 10, "numeric and analytic lengths do not agree in time step " + str(i + 1))
#             self.assertAlmostEqual(l[i], nl2[i], 10, "numeric and analytic lengths do not agree in time step " + str(i + 1))
#             self.assertLessEqual(meanDX[i], 0.5, "axial resolution dx is too large")
#             self.assertLessEqual(0.25, meanDX[i], "axial resolution dx is unexpected small")
#
#     def test_constructors(self):
#         """ tests two kinds of constructors and copy"""
#         self.leaf_example_rtp()
#         # 1. constructor from scratch
#         param = self.p0.realize()
#         leaf = pb.Leaf(1, param, True, True, 0., 0., pb.Vector3d(0, 0, -1), 0, 0, False, 0)
#         leaf.setOrganism(self.plant)
#         leaf.addNode(pb.Vector3d(0, 0, -3), 0)  # parent must have at least one nodes
#         # 2. used in simulation (must have parent, since there is no nullptr in Pyhton)
#         leaf2 = pb.Leaf(self.plant, self.p1.subType, pb.Vector3d(0, 0, -1), 0, leaf, 0, 0)
#         leaf.addChild(leaf2)
#         # 3. deep copy (with a factory function)
#         plant2 = pb.Organism()
#         leaf3 = leaf.copy(plant2)
#         self.assertEqual(str(leaf), str(leaf3), "deep copy: the organs shold be equal")
#         self.assertIsNot(leaf.getParam(), leaf3.getParam(), "deep copy: organs have same parameter set")
#         # TODO check if OTP were copied
#
#     def test_leaf_length(self):
#         """ tests if numerical leaf length agrees with analytic solutions at 4 points in time with two scales of dt"""
#         self.leaf_example_rtp()
#         times = np.array([0., 7., 15., 30., 60.])
#         dt = np.diff(times)
#         k = self.leaf.param().getK()  # maximal leaf length
#         self.assertAlmostEqual(k, 100, 12, "example leaf has wrong maximal length")
#         l = leafLength(times[1:], self.p0.r, k)  # analytical leaf length
#         leaf = self.leaf.copy(self.plant)
#         self.leaf_length_test(dt, l, 1)  # large dt
#         self.leaf = leaf
#         self.leaf_length_test(dt, l, 1000)  # very fine dt
#
#     def test_leaf_length_including_laterals(self):
#         """ tests if numerical leaf length agrees with analytic solution including laterals """
#         self.leaf_example_rtp()
#         times = np.array([0., 7., 15., 30., 60.])
#         dt = np.diff(times)
#         p = self.leaf.param()  # rename
#         k = p.getK()
#         et = np.zeros((p.nob()))
#         l = 0
#         et[0] = leafAge(p.la + p.lb + l, p.r, k)
#         for i in range(0, p.nob() - 1):  # calculate lateral emergence times
#             l += p.ln[i]
#             et[i + 1] = leafAge(p.la + p.lb + l, p.r, k + 1e-12)
#         l = leafLength(times[1:], p.r, k)  # zero order lengths
#         l1 = []
#         r2 = self.p1.r
#         k2 = self.p1.la  # consists of lateral zone only
#         for t in times[1:]:
#             l1.append(leafLateralLength(t, et, r2, k2))
#         analytic_total = l + l1
#
#         for subDX in [1, 1000]:
#             numeric_total = []
#             for t in times[1:]:
#                 leaf = self.leaf.copy(self.plant)
#                 self.leaf_length_test([t], [leafLength(t, p.r, k)], subDX)
#                 organs = self.leaf.getOrgans()
#                 nl = 0
#                 for o in organs:
#                     nl += o.getParameter("length")
#                 numeric_total.append(nl);
#                 self.leaf = leaf
#             for i in range(0, len(times[1:])):
#                 self.assertAlmostEqual(numeric_total[i], analytic_total[i], 10, "numeric and analytic total lengths do not agree in time step " + str(i + 1))
#
#     def test_geometry(self):
#         """ tests if nodes can be retrieved from the organ """
#         # TODO make plot for plausibility
#
#     def test_parameter(self):
#         """ tests some parameters on sequential organ list """
#         self.leaf_example_rtp()
#         self.leaf.simulate(30)
#         organs = self.leaf.getOrgans()
#         type, age, radius, order, ct = [], [], [], [], []
#         for o in organs:
#             type.append(o.getParameter("subType"))
#             age.append(o.getParameter("age"))
#             ct.append(o.getParameter("creationTime"))
#             radius.append(o.getParameter("radius"))
#             order.append(o.getParameter("order"))
#         self.assertEqual(type, [1.0, 3.0, 3.0, 3.0, 3.0], "getParameter: unexpected leaf sub types")
#         self.assertEqual(order, [1.0, 2.0, 2.0, 2.0, 2.0], "getParameter: unexpected leaf sub types")  # +1, because of artificial parent leaf
#         for i in range(0, 5):
#             self.assertEqual(age[i], 30 - ct[i], "getParameter: unexpected leaf sub types")  # +1, because of artificial parent leaf
#
#     def test_dynamics(self):
#         """ tests if nodes created in last time step are correct """  #
#         self.leaf_example_rtp()
#         r = self.leaf
#         r.simulate(.5, True)
#         self.assertEqual(r.hasMoved(), False, "dynamics: node movement during first step")
#         r.simulate(1e-1, True)
#         self.assertEqual(r.hasMoved(), False, "dynamics: movement, but previous node at axial resolution")
#         r.simulate(1e-1, True)
#         self.assertEqual(r.hasMoved(), True, "dynamics: node was expected to move, but did not")
#         r.simulate(2.4, True)
#         self.assertEqual(r.getNumberOfNodes() - r.getOldNumberOfNodes(), 5, "dynamics: unexcpected number of new nodes")
#
#     def test_leafgrow(self):
#         """ tests if the leaf can create leaf """  #
#         self.leaf_example_rtp()
#         r = self.leaf


if __name__ == '__main__':
    unittest.main()
