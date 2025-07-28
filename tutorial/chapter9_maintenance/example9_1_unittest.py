import sys; sys.path.append("../.."); sys.path.append("../../src/")
import unittest

import plantbox as pb
from rsml.rsml_reader import *

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

    def test_parameter(self):
        """ tests some parameters on sequential organ list """
        self.root_example_rrp()
        simtime = 30.
        self.root.simulate(simtime, False)
        organs = self.root.getOrgans()
        age, ct = [], []
        for o in organs:
            age.append(o.getParameter("age"))
            ct.append(o.getParameter("creationTime"))

        nol = round(self.root.getParameter("numberOfLaterals"))
        for i in range(0, nol):
            self.assertAlmostEqual(age[i], simtime - ct[i], 10, "getParameter: age and creation time does not agree")
    
    def test_copy(self):
        """ checks if the root system can be copied, and if randomness works """
        seed = 110  # random seed
        path = "../../modelparameter/structural/plant/"
        name = "Brassica_oleracea_Vansteenkiste_2014"
        rs = pb.Plant()  # the original
        rs.readParameters(path + name + ".xml", verbose = False)
        rs.setSeed(seed)
        rs.initialize(False)
        rs2 = rs.copy()  # copy root system
        self.assertIsNot(rs2, rs, "copy: not a copy")
        self.assertEqual(str(rs), str(rs2), "copy: the organisms should be equal")

if __name__ == '__main__':
    unittest.main()