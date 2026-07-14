import math
import unittest

import plantbox as pb


class TestOrgan(unittest.TestCase):
    """Unit tests for the Organ base class.

    The 'hand with two fingers' example is a minimal organ tree used throughout:
      hand  (created from explicit state, length=15, two children added manually)
      ├── thumb        (delayed 4 days, pni=0)
      └── little_finger (delayed 3 days, pni=0)
    """

    def hand_example(self):
        """Creates the hand/thumb/little_finger organ tree and stores it on self."""
        self.ons = pb.Vector3d(0.0, 0.0, 1.0)
        self.human1 = pb.Organism()
        otp = pb.OrganRandomParameter(self.human1)
        self.human1.setOrganRandomParameter(otp)
        op = otp.realize()
        self.hand = pb.Organ(self.human1.getOrganIndex(), op, True, True, 0, 15.0, self.ons, 0, False, 0)
        self.hand.setOrganism(self.human1)
        self.thumb = pb.Organ(self.human1, self.hand, 0, 0, 4, 0)  # delayed 4 days
        self.little_finger = pb.Organ(self.human1, self.hand, 0, 0, 3, 0)  # delayed 3 days
        self.hand.addChild(self.thumb)
        self.hand.addChild(self.little_finger)

    def add_nodes(self):
        """Adds geometry to the hand example.

        Two getNodeIndex() calls make the global node counter disagree with the
        organ node counter so tests can verify id tracking.
        Resulting node layout:
          hand:          (0,0,0), (0,0,1.5), (0,-1,1.6), (0,1,1.6)
          thumb:         (0,-1,1.6) [shared id], (0,-2,2.5)
          little_finger: (0,1,1.6)  [shared id], (0,1.7,2.5)
        """
        self.human1.getNodeIndex()
        self.human1.getNodeIndex()  # to make global and organ index disagree
        self.hand.addNode(pb.Vector3d(0, 0, 0), 0)
        self.hand.addNode(pb.Vector3d(0, 0, 1.5), 0)
        self.hand.addNode(pb.Vector3d(0, -1, 1.6), 0)  # thumb base
        self.hand.addNode(pb.Vector3d(0, 1, 1.6), 0)  # little finger base
        thumb = self.hand.getNodeId(2)
        lf = self.hand.getNodeId(3)
        self.thumb.addNode(pb.Vector3d(0, -1, 1.6), thumb, 4)
        self.thumb.addNode(pb.Vector3d(0, -2, 2.5), 4)
        self.little_finger.addNode(pb.Vector3d(0, 1, 1.6), lf, 3)
        self.little_finger.addNode(pb.Vector3d(0, 1.7, 2.5), 3)

    def test_constructors(self):
        """Tests three different construction modes and deep copy."""
        ons = pb.Vector3d(0.0, 0.0, 1.0)
        human1 = pb.Organism()
        otp = pb.OrganRandomParameter(human1)
        human1.setOrganRandomParameter(otp)
        op = otp.realize()
        # 1. explicit-state constructor
        hand = pb.Organ(human1.getOrganIndex(), op, True, True, 0.0, 15.0, ons, 0, False, 0)
        hand.setOrganism(human1)
        # 2. simulation constructor (requires a parent; no nullptr in Python)
        thumb = pb.Organ(human1, hand, 0, 0, 4, 0)
        little_finger = pb.Organ(human1, hand, 0, 0, 3, 0)
        hand.addChild(thumb)
        hand.addChild(little_finger)
        # 3. deep copy via factory function
        human2 = pb.Organism()
        human2.setOrganRandomParameter(otp.copy(human2))
        hand2 = hand.copy(human2)
        self.assertEqual(str(hand), str(hand2), "deep copy: organs should compare equal")
        self.assertIsNot(hand.getParam(), hand2.getParam(), "deep copy: parameter sets must be distinct objects")
        self.assertEqual(str(hand.getParam()), str(hand2.getParam()), "deep copy: parameter values should match")
        self.assertEqual(str(hand.getOrganRandomParameter()), str(hand2.getOrganRandomParameter()), "deep copy: random parameter values should match")
        self.assertIsNot(hand.getOrganism(), hand2.getOrganism(), "deep copy: organisms must be distinct")

    def test_simulation(self):
        """Tests that ages are correct after simulate(10)."""
        self.hand_example()
        self.hand.simulate(10)
        self.assertEqual(self.hand.getAge(), 10, "wrong organ age for hand")
        self.assertEqual(self.thumb.getAge(), 6, "wrong organ age for thumb (delay=4)")
        self.assertEqual(self.little_finger.getAge(), 7, "wrong organ age for little_finger (delay=3)")

    def test_sequential(self):
        """Tests that getOrgans() returns the correct flat list of organs with >1 node."""
        self.hand_example()
        self.add_nodes()
        ring = pb.Organ(self.human1, self.thumb, 0, 0, 4, 0)  # ring on the thumb
        self.thumb.addChild(ring)
        ring.addNode(pb.Vector3d(0, -1, 1.6), self.thumb.getNodeId(1), 4)
        ring.addNode(pb.Vector3d(0, -1, 1.6), 4)
        organs = self.hand.getOrgans()
        self.assertEqual(len(organs), 4, "wrong number of organs in sequential list")

    def test_geometry(self):
        """Tests node counts, segment counts, and node coordinates."""
        self.hand_example()
        self.add_nodes()
        organs = [self.hand, self.thumb, self.little_finger]
        expected_nodes = [4, 2, 2]  # number of nodes per organ
        for oi, o in enumerate(organs):
            self.assertEqual(o.getNumberOfNodes(), expected_nodes[oi], "wrong number of nodes")
            self.assertEqual(len(o.getSegments()), expected_nodes[oi] - 1, "wrong number of segments")
            for i in range(o.getNumberOfNodes()):
                self.assertEqual(o.getNode(i).x, 0, "x coordinate should be zero")

    def test_parameter(self):
        """Tests getParameter() for 'age' and topological 'order'."""
        self.hand_example()
        self.add_nodes()
        self.assertEqual(self.little_finger.getParameter("age"), -3, "wrong age for little_finger")
        self.assertEqual(self.thumb.getParameter("age"), -4, "wrong age for thumb")
        self.assertEqual(self.hand.getParameter("order"), 0, "wrong order for hand")
        self.assertEqual(self.little_finger.getParameter("order"), 1, "wrong order for little_finger")
        self.assertEqual(self.thumb.getParameter("order"), 1, "wrong order for thumb")

    def test_dynamics(self):
        """Tests new-node detection via getOldNumberOfNodes() across simulate() calls."""
        self.hand_example()
        self.hand.simulate(1)
        self.add_nodes()
        # nodes added manually after simulate count as new (oldNumberOfNodes was 0)
        self.assertEqual(self.hand.getNumberOfNodes() - self.hand.getOldNumberOfNodes(), 4, "wrong number of new nodes for hand")
        self.assertEqual(self.thumb.getNumberOfNodes() - self.thumb.getOldNumberOfNodes(), 2, "wrong number of new nodes for thumb")
        self.assertEqual(self.little_finger.getNumberOfNodes() - self.little_finger.getOldNumberOfNodes(), 2, "wrong number of new nodes for little_finger")
        # second simulate: no nodes added, so delta = 0 for all
        self.hand.simulate(1)
        self.assertEqual(self.hand.getNumberOfNodes() - self.hand.getOldNumberOfNodes(), 0, "wrong number of new nodes after second simulate")
        self.assertEqual(self.little_finger.getNumberOfNodes() - self.little_finger.getOldNumberOfNodes(), 0, "wrong number of new nodes for little_finger after second simulate")
        self.assertEqual(self.thumb.getNumberOfNodes() - self.thumb.getOldNumberOfNodes(), 0, "wrong number of new nodes for thumb after second simulate")
        # adding a node manually between simulates is detected
        self.hand.simulate(1)
        self.little_finger.addNode(pb.Vector3d(0, 1, 1.6), 6)
        self.assertEqual(self.little_finger.getNumberOfNodes() - self.little_finger.getOldNumberOfNodes(), 1, "wrong number of new nodes after manual addNode")

    # -------------------------------------------------------------------------
    # New tests
    # -------------------------------------------------------------------------

    def test_tree_accessors(self):
        """Tests getId, getNumberOfChildren, getChild, getParent, parentNI, isAlive, isActive."""
        self.hand_example()
        # id is a non-negative integer
        self.assertIsInstance(self.hand.getId(), int, "getId: expected int")
        self.assertGreaterEqual(self.hand.getId(), 0, "getId: expected non-negative")
        # each child gets a different id
        self.assertNotEqual(self.thumb.getId(), self.little_finger.getId(), "ids must be unique")
        # child count and retrieval
        self.assertEqual(self.hand.getNumberOfChildren(), 2, "getNumberOfChildren: expected 2")
        self.assertIs(self.hand.getChild(0), self.thumb, "getChild(0): expected thumb")
        self.assertIs(self.hand.getChild(1), self.little_finger, "getChild(1): expected little_finger")
        # parent link
        self.assertIs(self.thumb.getParent(), self.hand, "getParent: thumb's parent should be hand")
        self.assertEqual(self.thumb.parentNI, 0, "parentNI: expected 0")
        # alive / active flags set at construction
        self.assertTrue(self.hand.isAlive(), "isAlive: hand should be alive")
        self.assertTrue(self.hand.isActive(), "isActive: hand should be active")

    def test_getNumberOfLaterals(self):
        """Tests getNumberOfLaterals() before and after children emerge."""
        self.hand_example()
        # before simulate: both children have negative age (not yet emerged)
        self.assertEqual(self.hand.getNumberOfLaterals(), 0, "getNumberOfLaterals: expected 0 before simulate")
        # after simulate(10): thumb (age=6) and little_finger (age=7) have both emerged
        self.hand.simulate(10)
        self.assertEqual(self.hand.getNumberOfLaterals(), 2, "getNumberOfLaterals: expected 2 after simulate(10)")

    def test_node_accessors(self):
        """Tests getNodeId, getNodeIds, getNodeCT, getNodes, getLength, dx, dxMin."""
        self.hand_example()
        self.add_nodes()
        # getNodeIds returns a list of the same length as the node count
        node_ids = self.hand.getNodeIds()
        self.assertEqual(len(node_ids), 4, "getNodeIds: expected 4 ids")
        # thumb's first node is inserted at hand's node index 2 and shares the global id
        self.assertEqual(self.thumb.getNodeId(0), self.hand.getNodeId(2), "getNodeId: thumb base should share global id with hand node 2")
        # getNodes() returns a list of Vector3d
        nodes = self.hand.getNodes()
        self.assertEqual(len(nodes), 4, "getNodes: expected 4 nodes")
        # node creation times — all hand nodes added with t=0
        self.assertAlmostEqual(self.hand.getNodeCT(0), 0.0, msg="getNodeCT hand[0]: expected 0")
        # thumb's first node was added with t=4
        self.assertAlmostEqual(self.thumb.getNodeCT(0), 4.0, msg="getNodeCT thumb[0]: expected 4")
        # little_finger's first node was added with t=3
        self.assertAlmostEqual(self.little_finger.getNodeCT(0), 3.0, msg="getNodeCT little_finger[0]: expected 3")
        # getLength(bool): hand was constructed with length=15.0 and epsilonDx=0
        self.assertAlmostEqual(self.hand.getLength(False), 15.0, msg="getLength(False): expected 15.0")
        self.assertAlmostEqual(self.hand.getLength(True), 15.0, msg="getLength(True): expected 15.0")
        # dx / dxMin from default OrganRandomParameter
        self.assertAlmostEqual(self.hand.dx(), 0.25, msg="dx: expected default 0.25")
        self.assertAlmostEqual(self.hand.dxMin(), 1e-6, msg="dxMin: expected default 1e-6")

    def test_hasMoved_epsilon(self):
        """Tests hasMoved() and getEpsilon() at construction and after simulate()."""
        self.hand_example()
        # at construction both flags are at their initial values
        self.assertFalse(self.hand.hasMoved(), "hasMoved: should be False at construction")
        self.assertAlmostEqual(self.hand.getEpsilon(), 0.0, msg="getEpsilon: should be 0 at construction")
        # base Organ::simulate() resets moved=False; no createSegments is called
        self.hand.simulate(1)
        self.assertFalse(self.hand.hasMoved(), "hasMoved: should be False after base simulate")
        self.assertAlmostEqual(self.hand.getEpsilon(), 0.0, msg="getEpsilon: should remain 0 after simulate")

    def test_getOrgans_filter(self):
        """Tests getOrgans() with organ-type filtering."""
        self.hand_example()
        self.add_nodes()
        # all organs in this example are type 0 (ot_organ)
        all_organs = self.hand.getOrgans(-1)
        type0_organs = self.hand.getOrgans(0)
        self.assertEqual(len(all_organs), len(type0_organs), "getOrgans: ot=0 and ot=-1 should return the same count")
        self.assertGreater(len(all_organs), 0, "getOrgans: expected at least one organ")
        # filtering for a type not present in this tree returns empty
        roots = self.hand.getOrgans(2)  # ot_root = 2
        self.assertEqual(len(roots), 0, "getOrgans(2): no roots in hand example")

    def test_parameter_extended(self):
        """Tests getParameter() for structural, state, and geometry parameter names."""
        self.hand_example()
        self.add_nodes()
        # identity / topology
        self.assertEqual(self.hand.getParameter("id"), self.hand.getId(), "getParameter id: should match getId()")
        self.assertEqual(self.hand.getParameter("organType"), 0.0, "getParameter organType: expected 0")
        self.assertEqual(self.hand.getParameter("numberOfChildren"), 2.0, "getParameter numberOfChildren: expected 2")
        # state flags (returned as float)
        self.assertEqual(self.hand.getParameter("alive"), 1.0, "getParameter alive: expected 1")
        self.assertEqual(self.hand.getParameter("active"), 1.0, "getParameter active: expected 1")
        # geometry
        self.assertEqual(self.hand.getParameter("numberOfNodes"), 4.0, "getParameter numberOfNodes: expected 4")
        self.assertEqual(self.hand.getParameter("numberOfSegments"), 3.0, "getParameter numberOfSegments: expected 3")
        # radius-related parameters from the realised OrganSpecificParameter
        a = self.hand.param().a
        self.assertAlmostEqual(self.hand.getParameter("a"), a, msg="getParameter a: should match param().a")
        self.assertAlmostEqual(self.hand.getParameter("radius"), a, msg="getParameter radius: should equal a")
        self.assertAlmostEqual(self.hand.getParameter("diameter"), 2.0 * a, msg="getParameter diameter: should be 2*a")

    def test_toString(self):
        """Tests that __str__ returns a non-empty string without raising."""
        self.hand_example()
        self.add_nodes()
        s = str(self.hand)
        self.assertIsInstance(s, str, "toString: expected a string")
        self.assertGreater(len(s), 0, "toString: returned empty string")

    def test_orgVolume(self):
        """Tests orgVolume() for the default cylinder approximation (V = π·a²·L)."""
        self.hand_example()
        a = self.hand.param().a
        length = self.hand.getLength(False)  # theoretical length = 15.0
        expected = math.pi * a * a * length
        self.assertAlmostEqual(self.hand.orgVolume(), expected, places=9, msg="orgVolume: expected π·a²·length")


if __name__ == "__main__":
    unittest.main()
    print("done.")
