import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
import rsml.rsml_reader as rsml_reader
import rsml.rsml_writer as rsml_writer
import visualisation.vtk_tools as vt
import visualisation.vtk_plot as vp

import matplotlib.pyplot as plt
import numpy as np


class TestOrganism(unittest.TestCase):

    def hand_example(self):
        """ an example used in the tests below (same as test_organ), a hand with two fingers """
        self.ons = pb.Matrix3d(pb.Vector3d(0., 0., 1.), pb.Vector3d(0., 1., 0.), pb.Vector3d(1., 0., 0.))
        self.human1 = pb.Organism()  # same example as in test_constructor ...
        otp = pb.OrganRandomParameter(self.human1)
        self.human1.setOrganRandomParameter(otp)
        op = otp.realize()
        self.hand = pb.Organ(self.human1.getOrganIndex(), op, True, True, 0, 15., pb.Vector3d(0, 0, -1), 0, False, 0)
        self.hand.setOrganism(self.human1)
        self.thumb = pb.Organ(self.human1, self.hand, 0, 0, 4,  0)  # delayedfor 4 days
        self.little_finger = pb.Organ(self.human1, self.hand, 0, 0, 3, 0)  # delayed for 3 days
        self.hand.addChild(self.thumb)
        self.hand.addChild(self.little_finger)
        self.human1.addOrgan(self.hand)

    def add_nodes(self):
        """ used in the tests below, adds nodes to the hand example """
        self.hand.addNode(pb.Vector3d(0, 0, 0), 0)
        self.hand.addNode(pb.Vector3d(0, 0, 1.5), 0)
        self.hand.addNode(pb.Vector3d(0, -1, 1.6), 0)  # thumb
        self.hand.addNode(pb.Vector3d(0, 1, 1.6), 0)  # little finger
        thumb = self.hand.getNodeId(2)
        lf = self.hand.getNodeId(3)
        self.thumb.addNode(pb.Vector3d(0, -1, 1.6), thumb, 4)
        self.thumb.addNode(pb.Vector3d(0, -2, 2.5), 4)
        self.little_finger.addNode(pb.Vector3d(0, 1, 1.6), lf, 3)
        self.little_finger.addNode(pb.Vector3d(0, 1.7, 2.5), 3)

    def test_rsml_cplantbox(self):
        """ checks cplantbox rmsl writer functionality with Python rsml reader """
        self.hand_example()
        self.add_nodes()
        self.human1.writeRSML("organism.rsml")
        pl, props, funcs, _ = rsml_reader.read_rsml("organism.rsml")
        pl2 = [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.5], [0.0, -1.0, 1.6], [0.0, 1.0, 1.6]], [[0.0, -2.0, 2.5]], [[0.0, 1.7, 2.5]]]
        self.assertEqual(pl, pl2, "rsml: polylines are not equal")
        self.assertEqual(props["age"], [0, -4, -3] , "rsml: polylines are not equal")

    def test_rsml_vt(self):
        """ checks vtk_tools rmsl writer functionality, first reading with the Python rsml reader 
        TODO this is not a test, yet
        """
        pd = vt.read_vtp("Dumux-VTP.vtp")
        meta = rsml_writer.Metadata()
        meta.unit = "m"
        meta.add_property(rsml_writer.Property("radius [m]", "float", "m", None))
        order_id = 5
        vt.write_rsml("organism2.rsml", pd, order_id, meta)  # meta is optional now


if __name__ == '__main__':
    # todo test XML ?
    unittest.main()
