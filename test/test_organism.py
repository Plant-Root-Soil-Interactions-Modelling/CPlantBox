import sys; sys.path.append(".."); sys.path.append("../src/python_modules")

import unittest
import sys; sys.path.append("..")
import plantbox as pb

import matplotlib.pyplot as plt

import plantbox as pb
from rsml_reader import *


class TestOrganism(unittest.TestCase):

    def hand_example(self):
        """ an example used in the tests below (same as test_organ), a hand with two fingers """
        self.ons = pb.Matrix3d(pb.Vector3d(0., 0., 1.), pb.Vector3d(0., 1., 0.), pb.Vector3d(1., 0., 0.))
        self.human1 = pb.Organism()  # same example as in test_constructor ...
        otp = pb.OrganRandomParameter(self.human1)
        self.human1.setOrganRandomParameter(otp)
        op = otp.realize()
        self.hand = pb.Organ(self.human1.getOrganIndex(), op, True, True, 0, 15., self.ons, 0, False, 0)
        self.hand.setOrganism(self.human1)
        self.thumb = pb.Organ(self.human1, self.hand, 0, 0, 4, self.ons, 2)  # delayedfor 4 days
        self.little_finger = pb.Organ(self.human1, self.hand, 0, 0, 3, self.ons, 3)  # delayed for 3 days
        self.hand.addChild(self.thumb)
        self.hand.addChild(self.little_finger)
        self.human1.addOrgan(self.hand)

    def add_nodes(self):
        """ used in the tests below, adds nodes to the hand example """
        self.hand.addNode(pb.Vector3d(0, 0, 0), 0)

        self.hand.addNode(pb.Vector3d(0, 0, 1.5), 0)
        self.hand.addNode(pb.Vector3d(0, -1, 1.6), 0)  # thumb
        self.hand.addNode(pb.Vector3d(0, 1, 1.6), 0)  # little finger
        nodes = np.array((list(map(np.array, self.human1.getNodes()))))
        print(nodes)
        thumb = self.hand.getNodeId(2)
        lf = self.hand.getNodeId(3)
        self.thumb.addNode(pb.Vector3d(0, -1, 1.6), thumb, 4)
        nodes = np.array((list(map(np.array, self.human1.getNodes()))))
        print("add 0, -1, 1.6\n", nodes)
        self.thumb.addNode(pb.Vector3d(0, -2, 2.5), 4)
        nodes = np.array((list(map(np.array, self.human1.getNodes()))))
        print("add 0, -2, 2.5\n",nodes)
        self.little_finger.addNode(pb.Vector3d(0, 1, 1.6), lf, 3)
        self.little_finger.addNode(pb.Vector3d(0, 1.7, 2.5), 3)
##old version
        #self.hand.addNode(pb.Vector3d(0, 0, -1.5), 0)
        #self.hand.addNode(pb.Vector3d(0, -1, -1.6), 0)  # thumb
        #self.hand.addNode(pb.Vector3d(0, 1, -1.6), 0)  # little finger
        #thumb = self.hand.getNodeId(2)
        #lf = self.hand.getNodeId(3)
        #self.thumb.addNode(pb.Vector3d(0, -1, -1.6), thumb, 4)
        #self.thumb.addNode(pb.Vector3d(0, -2, -2.5), 4)
        #self.little_finger.addNode(pb.Vector3d(0, 1, -1.6), lf, 3)
        #self.little_finger.addNode(pb.Vector3d(0, 1.7, -2.5), 3)
##oldversion end

    def test_geometry(self):
        """ tests ability to retrieve geometry """
        self.hand_example()
        self.add_nodes()
        nodes = np.array((list(map(np.array, self.human1.getNodes()))))
        print(nodes)
        self.assertEqual(nodes.shape, (6, 3), "geometry: number of nodes unexpected")
        segs = np.array((list(map(np.array, self.human1.getSegments()))))
        self.assertEqual(np.sum(np.sum(segs.flat != np.array([[0, 1], [1, 2], [2, 3], [2, 4], [3, 5]]).flat)), 0, "geometry: segments ids are unexcpected")

##new tests
    def test_parameter(self):
        """ test if getParameter works """
        self.hand_example()
        self.add_nodes()
        age = self.human1.getParameter("age")
        simtime = 10;
        self.human1.simulate(simtime)
        organs = self.human1.getOrgans()
        age = self.human1.getParameter("age")
        self.assertEqual(len(organs), len(age) , "parameter: organ size unequal to parameter size")
        self.assertEqual(age[2], -3 + simtime, "parameter: wrong ages")

    def test_rsml(self):
        """ checks rmsl functionality with Python rsml reader """
        self.hand_example()
        self.add_nodes()
        self.human1.writeRSML("organism.rsml")
        pl, props, funcs, _ = read_rsml("organism.rsml")
        pl2 = [[[0.0, 0.0, 0.0], [0.0, 0.0, -1.5], [0.0, -1.0, -1.6], [0.0, 1.0, -1.6]], [[0.0, -2.0, -2.5]], [[0.0, 1.7, -2.5]]]
        self.assertEqual(pl, pl2, "rsml: polylines are not equal")
        self.assertEqual(props["age"], [0, -4, -3] , "rsml: polylines are not equal")



if __name__ == '__main__':
    # todo test XML ?
    unittest.main()
