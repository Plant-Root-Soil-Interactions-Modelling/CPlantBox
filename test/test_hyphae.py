import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb

import numpy as np
from scipy.linalg import norm


def arc_length(organ, i):
    """ length along the organ from its first node to node i """
    nodes = np.array([[n.x, n.y, n.z] for n in organ.getNodes()[:i + 1]])
    return np.sum(norm(np.diff(nodes, axis = 0), axis = 1)) if len(nodes) > 1 else 0.


class TestHyphae(unittest.TestCase):

    def hyphae_example(self, b_prob = 0., distTH = 0., ana = 1., hlt = 10., laterals = True, hed = 1., seed = 1):
        """ a mycorrhizal plant with a single colonized tap root, hyphae emerge from the colonized root nodes """
        self.plant = pb.MycorrhizalPlant(seed)
        srp = pb.SeedRandomParameter(self.plant)
        srp.subType = 0
        self.plant.setOrganRandomParameter(srp)
        self.rrp = pb.MycorrhizalRootRandomParameter(self.plant)
        self.rrp.name, self.rrp.subType, self.rrp.lmax, self.rrp.r, self.rrp.dx = "taproot", 1, 10., 2., 0.5
        self.rrp.lmbd, self.rrp.highresolution, self.rrp.hyphalEmergenceDensity = 2., 0., hed
        self.plant.setOrganRandomParameter(self.rrp)
        self.hrp = pb.HyphaeRandomParameter(self.plant)  # root hyphae are always of sub type 1
        self.hrp.name, self.hrp.subType, self.hrp.dx = "hyphae", 1, 0.05
        self.hrp.v, self.hrp.vs, self.hrp.hlt, self.hrp.hlts = 1., 0., hlt, 0.
        self.hrp.lmax, self.hrp.lmaxs, self.hrp.la, self.hrp.lb, self.hrp.ln = 1., 0., 0.2, 0.2, 0.2
        self.hrp.b, self.hrp.bs, self.hrp.b_prob = 2., 0., b_prob
        self.hrp.distTH, self.hrp.ana = distTH, ana
        if laterals:
            self.hrp.successor = [[1]]
            self.hrp.successorP = [[1.]]
        self.plant.setOrganRandomParameter(self.hrp)
        self.plant.initialize(False)

    def simulate(self, time, dt = 0.25):
        for i in range(0, int(round(time / dt))):
            self.plant.simulate(dt, False)
        return self.plant.getOrgans(pb.hyphae)

    def test_emergence(self):
        """ tests that hyphae emerge from colonized nodes of the mycorrhizal root """
        self.hyphae_example()
        hyphae = self.simulate(3.)
        self.assertGreater(len(hyphae), 0, "emergence: no hyphae emerged")
        root_hyphae = [h for h in hyphae if h.getParent().organType() == pb.root]
        self.assertGreater(len(root_hyphae), 0, "emergence: no hyphae emerged from the root")
        for h in root_hyphae:
            root = h.getParent()
            self.assertGreater(root.getNodeColonization(h.parentNI), 0, "emergence: hypha emerged from a node that is not colonized")
            self.assertEqual(h.getNodeId(0), root.getNodeId(h.parentNI), "emergence: hypha does not start at its root node")
            self.assertGreaterEqual(h.getNodeCT(0), root.getNodeColonizationTime(h.parentNI) - 1.e-9,
                                    "emergence: hypha emerged before its root node was colonized")
        root = self.plant.getOrgans(pb.root)[0]
        self.assertLessEqual(len(root_hyphae), self.rrp.hyphalEmergenceDensity * root.getLength() + 1,
                             "emergence: more hyphae than given by the hyphal emergence density")

    def test_constructors(self):
        """ tests deep copy of a hypha """
        self.hyphae_example()
        hyphae = self.simulate(2.)
        h = [h for h in hyphae if h.getNumberOfChildren() > 0][0]
        plant2 = pb.MycorrhizalPlant()
        h2 = h.copy(plant2)
        self.assertEqual(h2.__class__.__name__, "Hyphae", "copy: unexpected class type")
        self.assertIsNot(h2.param(), h.param(), "copy: hyphae share the same specific parameter set")
        self.assertEqual(str(h2.param()), str(h.param()), "copy: hyphae have different parameter values")
        self.assertEqual(h2.getNumberOfNodes(), h.getNumberOfNodes(), "copy: unexpected number of nodes")
        self.assertAlmostEqual(h2.getLength(), h.getLength(), 12, "copy: unexpected length")
        self.assertEqual(h2.getNumberOfChildren(), h.getNumberOfChildren(), "copy: unexpected number of children")
        self.assertEqual(h2.getParameter("hyphalTreeIndex"), h.getParameter("hyphalTreeIndex"), "copy: unexpected hyphal tree index")
        self.assertIsNot(h2.getChild(0), h.getChild(0), "copy: children are not deep copied")

    def test_parameter(self):
        """ tests getParameter() of the hyphae organ """
        self.hyphae_example()
        hyphae = self.simulate(1.)
        for h in hyphae:
            self.assertEqual(h.getParameter("organType"), pb.hyphae, "getParameter: value unexpected")
            self.assertEqual(h.getParameter("subType"), 1, "getParameter: value unexpected")
            self.assertAlmostEqual(h.getParameter("v"), self.hrp.v, 12, "getParameter: value unexpected")
            self.assertAlmostEqual(h.getParameter("hlt"), self.hrp.hlt, 12, "getParameter: value unexpected")
            self.assertAlmostEqual(h.getParameter("theta"), self.hrp.theta, 12, "getParameter: value unexpected")
            self.assertAlmostEqual(h.getParameter("lmax"), self.hrp.lmax, 10, "getParameter: value unexpected")
            self.assertGreaterEqual(h.getParameter("hyphalTreeIndex"), 0, "getParameter: hyphal tree index is not set")
            self.assertEqual(h.getParameter("mergePointID"), -1, "getParameter: no anastomosis expected (distTH = 0)")

    def test_length(self):
        """ tests that hyphae grow linearly with v up to lmax """
        self.hyphae_example()
        hyphae = self.simulate(3.)
        for h in hyphae:
            l = min(self.hrp.v * h.getAge(), self.hrp.lmax)
            self.assertAlmostEqual(h.getLength(), l, 6, "length: numeric and analytic length do not agree")
            self.assertAlmostEqual(arc_length(h, h.getNumberOfNodes() - 1), h.getLength(), 6, "length: length and geometry do not agree")
            self.assertLessEqual(np.max(np.diff([arc_length(h, i) for i in range(h.getNumberOfNodes())])), self.hrp.dx + 1.e-9,
                                 "length: segment longer than dx")
        self.assertTrue(any(abs(h.getLength() - self.hrp.lmax) < 1.e-6 for h in hyphae), "length: no hypha reached lmax")
        for h in hyphae:
            if abs(h.getLength() - self.hrp.lmax) < 1.e-6:
                self.assertFalse(h.isActive(), "length: hypha at lmax should be inactive")

    def test_lifetime(self):
        """ tests that hyphae die after the hyphal lifetime hlt """
        self.hyphae_example(hlt = 1.)
        hyphae = self.simulate(3.)
        dead = 0
        for h in hyphae:
            self.assertLessEqual(h.getAge(), self.hrp.hlt + 1.e-9, "lifetime: hypha older than hlt")
            self.assertLessEqual(h.getLength(), self.hrp.v * self.hrp.hlt + 1.e-6, "lifetime: hypha grew after its death")
            if not h.isAlive():
                dead += 1
                self.assertAlmostEqual(h.getAge(), self.hrp.hlt, 9, "lifetime: dead hypha has unexpected age")
        self.assertGreater(dead, 0, "lifetime: no hypha died")

    def test_lateral_branching(self):
        """ tests position of lateral hyphae along the parent hypha (b_prob = 0, no tip splitting) """
        self.hyphae_example(b_prob = 0.)
        hyphae = self.simulate(3.)
        p = 0
        for h in hyphae:
            pnis = sorted([h.getChild(i).parentNI for i in range(h.getNumberOfChildren())])
            self.assertEqual(len(pnis), len(set(pnis)), "lateral branching: two laterals at the same node")
            for k, pni in enumerate(pnis):  # laterals at lb, lb + ln, lb + 2*ln, ...
                p += 1
                self.assertAlmostEqual(arc_length(h, pni), self.hrp.lb + k * self.hrp.ln, 6, "lateral branching: unexpected lateral position")
            nob = round((self.hrp.lmax - self.hrp.la - self.hrp.lb) / self.hrp.ln + 1)
            self.assertLessEqual(len(pnis), nob, "lateral branching: too many laterals")
        self.assertGreater(p, 0, "lateral branching: no laterals")

    def test_no_laterals(self):
        """ tests that hyphae without successors do not branch """
        self.hyphae_example(laterals = False)
        hyphae = self.simulate(3.)
        self.assertGreater(len(hyphae), 0, "no laterals: no hyphae emerged")
        for h in hyphae:
            self.assertEqual(h.getNumberOfChildren(), 0, "no laterals: hypha without successor has laterals")

    def test_tip_splitting(self):
        """ tests tip splitting (b_prob = 1): the tip splits into two hyphae after age*b > 0.5 """
        self.hyphae_example(b_prob = 1.)
        hyphae = self.simulate(3.)
        splits = 0
        for h in hyphae:
            n = h.getNumberOfChildren()
            if n > 0:
                splits += 1
                self.assertEqual(n, 2, "tip splitting: tip should split into exactly two hyphae")
                self.assertEqual(h.getChild(0).parentNI, h.getNumberOfNodes() - 1, "tip splitting: split not at the tip")
                self.assertEqual(h.getChild(1).parentNI, h.getNumberOfNodes() - 1, "tip splitting: split not at the tip")
                self.assertFalse(h.isActive(), "tip splitting: split hypha should be inactive")
                self.assertGreater(h.getAge() * self.hrp.b, 0.5, "tip splitting: split too early")
        self.assertGreater(splits, 0, "tip splitting: no tip splits")

    def test_hyphal_tree_index(self):
        """ tests that each root hypha starts a new hyphal tree and laterals inherit the index """
        self.hyphae_example()
        hyphae = self.simulate(3.)
        root_indices = [h.getParameter("hyphalTreeIndex") for h in hyphae if h.getParent().organType() == pb.root]
        self.assertEqual(len(root_indices), len(set(root_indices)), "hyphal tree index: root hyphae share an index")
        for h in hyphae:
            self.assertGreaterEqual(h.getParameter("hyphalTreeIndex"), 0, "hyphal tree index: index is negative")
            for i in range(h.getNumberOfChildren()):
                self.assertEqual(h.getChild(i).getParameter("hyphalTreeIndex"), h.getParameter("hyphalTreeIndex"),
                                 "hyphal tree index: lateral has a different index than its parent")

    def test_anastomosis(self):
        """ tests that merged hyphae are inactive and end at a node of a hypha of another hyphal tree """
        self.hyphae_example(distTH = 0.3, hed = 3.)  # dense hyphae, so that hyphae of different trees meet
        hyphae = self.simulate(3.)
        merged = [h for h in hyphae if h.getParameter("mergePointID") > -1]
        self.assertGreater(len(merged), 0, "anastomosis: no anastomosis occurred")
        for h in merged:
            self.assertFalse(h.isActive(), "anastomosis: merged hypha should be inactive")
            other = h.getMergedHyphae()
            self.assertIsNotNone(other, "anastomosis: merged hypha is not set")
            self.assertNotEqual(other.getParameter("hyphalTreeIndex"), h.getParameter("hyphalTreeIndex"),
                                "anastomosis: hypha merged with its own hyphal tree")
            mpid = int(h.getParameter("mergePointID"))
            self.assertIn(mpid, list(other.getNodeIds()), "anastomosis: merge point is not a node of the merged hypha")
            merge_node = other.getNode(list(other.getNodeIds()).index(mpid))
            last = h.getNode(h.getNumberOfNodes() - 1)
            self.assertAlmostEqual(norm(np.array([last.x, last.y, last.z]) - np.array([merge_node.x, merge_node.y, merge_node.z])), 0., 12,
                                   "anastomosis: hypha does not end at the merge point")
        self.assertEqual(sum(self.plant.getAnastomosisPoints(pb.hyphae)), len(merged), "anastomosis: unexpected number of anastomosis points")

    def test_no_anastomosis(self):
        """ tests that no anastomosis occurs for distTH = 0 or ana = 0 """
        for distTH, ana in [(0., 1.), (0.1, 0.)]:
            self.hyphae_example(distTH = distTH, ana = ana)
            hyphae = self.simulate(3.)
            merged = [h for h in hyphae if h.getParameter("mergePointID") > -1]
            self.assertEqual(len(merged), 0, "no anastomosis: anastomosis for distTH = {:g}, ana = {:g}".format(distTH, ana))


if __name__ == '__main__':
    unittest.main()
