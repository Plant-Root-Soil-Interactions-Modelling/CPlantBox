import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
from rsml.rsml_reader import *

path = "../modelparameter/structural/plant/"


def rootLength(t, r, k):  # root length at a certain age
    return k * (1 - np.exp(-r * t / k))


class TestPlant(unittest.TestCase):

    def root_example_rrp(self):
        """ an example used in the tests below, a main root with laterals """
        self.plant = pb.Organism()  # store organism (not owned by Organ, or OrganRandomParameter)
        p0 = pb.RootRandomParameter(self.plant)
        p0.name, p0.subType, p0.la, p0.lb, p0.lmax, p0.ln, p0.r, p0.dx = "taproot", 1, 10., 1., 100., 1., 1.5, 0.5
        p0.successor = [[2]]
        p0.successorP = [[1.]]
        p1 = pb.RootRandomParameter(self.plant)
        p1.name, p1.subType, p1.lmax, p1.r, p1.dx = "lateral", 2, 25., 2., 0.1
        self.p0, self.p1 = p0, p1  # needed at later point
        self.plant.setOrganRandomParameter(p0)  # the organism manages the type parameters and takes ownership
        self.plant.setOrganRandomParameter(p1)
        srp = pb.SeedRandomParameter(self.plant)
        self.plant.setOrganRandomParameter(srp)

        param0 = p0.realize()  # set up root by hand (without a root system)
        param0.la, param0.lb = 0, 0  # its important parent has zero length, otherwise creation times are messed up
        parentroot = pb.Root(1, param0, True, True, 0., 0., pb.Vector3d(0, 0, -1), 0, False, 0)  # takes ownership of param0
        parentroot.setOrganism(self.plant)
        parentroot.addNode(pb.Vector3d(0, 0, -3), 0)  # there is no nullptr in Python

        self.parentroot = parentroot  # store parent (not owned by child Organ)
        self.root = pb.Root(self.plant, p0.subType,  0, self.parentroot , 0)
        self.root.setOrganism(self.plant)

    def test_constructors(self):
        """ tests two kinds of constructors and copy"""
        self.root_example_rrp()
        # 1. constructor from scratch
        param = self.p0.realize()
        root = pb.Root(1, param, True, True, 0., 0., pb.Vector3d(0, 0, -1), 0, False, 0)
        root.setOrganism(self.plant)
        root.addNode(pb.Vector3d(0, 0, -3), 0)  # parent must have at least one nodes
        # 2. used in simulation (must have parent, since there is no nullptr in Pyhton)
        root2 = pb.Root(self.plant, self.p1.subType, 0, root, 0)
        root.addChild(root2)
        # 3. deep copy (with a factory function)
        plant2 = pb.Organism()
        root3 = root.copy(plant2)
        self.assertEqual(str(root), str(root3), "deep copy: the root string representations shold be equal")
        self.assertIsNot(root.getParam(), root3.getParam(), "deep copy: roots have same specific parameter set")  # type OrganSpecificParameter
        self.assertEqual(str(root.param()), str(root3.param()), "deep copy: roots have different parameter values")  # type RootSpecificParameter


    def test_copy(self):
        """ checks if the root system can be copied, and if randomness works """
        seed = 110  # random seed
        name = "Brassica_oleracea_Vansteenkiste_2014"
        rs = pb.Plant()  # the original
        rs.readParameters("../modelparameter/structural/rootsystem/" + name + ".xml", verbose = False)
        rs.setSeed(seed)
        rs.initialize(False)
        rs2 = rs.copy()  # copy root system
        n1 = rs.rand()
        self.assertIsNot(rs2, rs, "copy: not a copy")
        self.assertEqual(str(rs), str(rs2), "copy: the organisms should be equal")
        self.assertEqual(rs2.rand(), n1, "copy: random generator seed was not copied")
        rs.simulate(10)
        rs2.simulate(10)
        n2 = rs.rand()
        self.assertEqual(rs2.rand(), n2, "copy: simulation is not deterministic")
        rs3 = pb.Plant()  # rebuild same
        rs3.readParameters("../modelparameter/structural/rootsystem/" + name + ".xml", verbose = False)
        rs3.setSeed(seed)
        rs3.initialize(False)
        self.assertEqual(rs3.rand(), n1, "copy: random generator seed was not copied")
        rs3.simulate(10)
        self.assertEqual(rs3.rand(), n2, "copy: simulation is not deterministic")

if __name__ == '__main__':
    # MANY tests missing !!!

    unittest.main()
