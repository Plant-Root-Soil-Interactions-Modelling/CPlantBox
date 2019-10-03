import unittest
import sys
import numpy as np
sys.path.append("..")
import plantbox as pb


class TestRoot(unittest.TestCase):

    def seed_example_srp(self):
        """ an example used in the tests below, a main root with laterals """
        pass

    def test_constructors(self):
        """ tests two kinds of constructors and copy"""
        self.seed_example_srp()
        # 1. constructor from scratch
        param = self.p0.realize()
        seed = pb.Organ(self.plant.getOrganIndex(), param, True, True, 0, 15., False, 0)
        seed.setOrganism(self.plant)
        seed.addNode(pb.Vector3d(0, 0, -3), 0)
        # 2. used in simulation (must have parent, since there is no nullptr in Pyhton)
        seed2 = pb.Seed(self.plant)
        # 3. deep copy (with a factory function)
        plant2 = pb.Organism()
        seed3 = root.copy(plant2)
        print("root")
        self.assertEqual(str(seed), str(seed3), "deep copy: the seed string representations shold be equal")
        self.assertIsNot(seed.getParam(), seed3.getParam(), "deep copy: organs have same parameter set")
        self.assertEqual(str(seed.param()), str(seed3.param()), "deep copy: organs have different parameter values")  # type RootSpecificParameter

    def test_parameter(self):
        """ tests some parameters on the seed """
        self.root_example_rtp()
        self.root.simulate(30)
        organs = self.root.getOrgans()
        type, age, radius, order, ct = [], [], [], [], []
        for o in organs:
            type.append(o.getParameter("subType"))
            age.append(o.getParameter("age"))
            ct.append(o.getParameter("creationTime"))
            radius.append(o.getParameter("radius"))
            order.append(o.getParameter("order"))
        self.assertEqual(type, [1.0, 2.0, 2.0, 2.0, 2.0], "getParameter: unexpected root sub types")
        self.assertEqual(order, [1.0, 2.0, 2.0, 2.0, 2.0], "getParameter: unexpected root sub types")  # +1, because of artificial parent root
        for i in range(0, 5):
            self.assertEqual(age[i], 30 - ct[i], "getParameter: unexpected root sub types")  # +1, because of artificial parent root

    def test_initialize(self):
        """test initialization """
        pass


if __name__ == '__main__':
    unittest.main()
