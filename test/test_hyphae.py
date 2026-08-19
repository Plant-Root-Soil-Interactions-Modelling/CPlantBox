import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb

import numpy as np
from scipy.linalg import norm

class TestHyphae(unittest.TestCase):
    def hyphae_example(self): # TODO check which parameters are needed and how structure should be
        self.plant = pb.Organism()
        self.hrp = pb.HyphaeRandomParameter(self.plant)
        self.hrp.dx = 0.1
        self.hrp.distTH = 0.2
        self.hrp.subType = 1

    def test_constructors(self): # TODO check again all parameters needed
        """ tests constructor and copy """
        plant = pb.Organism()
        hrp = pb.HyphaeRandomParameter(plant)
        hrp.dx = 0.15
        hrp.distTH = 0.25
        hrp.subType = 2
        hrp2 = hrp.copy(plant)
        self.assertIsNot(hrp, hrp2, "copy: organ type parameter set is not copied")
        self.assertEqual(hrp2.name, hrp.name, "copy: value unexpected")
        self.assertEqual(hrp2.organType, hrp.organType, "copy: value unexpected")
        self.assertEqual(hrp2.subType, hrp.subType, "copy: value unexpected")
        self.assertEqual(hrp2.dx, hrp.dx, "copy: value unexpected")
        self.assertEqual(hrp2.distTH, hrp.distTH, "copy: value unexpected")
    
    def test_parameter(self):
        """ tests getParameter() """
        self.plant = pb.Organism()
        hrp = pb.HyphaeRandomParameter(self.plant)
        hrp.dx = 0.2
        hrp.distTH = 0.3
        ot = hrp.getParameter("organType")  # test defaults but this works only on RandomPArameter but need the Organ
        st = hrp.getParameter("subType")
        hTI = hrp.getParameter("hyphalTreeIndex")
        mPID = hrp.getParameter("mergePointID")
        self.assertEqual(ot, pb.hyphae, "getParameter: value unexpected")
        self.assertEqual(st, 0, "getParameter: value unexpected")
        self.assertEqual(hTI, -1, "getParameter: value unexpected")
        self.assertEqual(mPID, -1, "getParameter: value unexpected")

    def test_dynamics(self):
        """ tests hyphae growth simulation """
        self.hyphae_example()
        self.plant.setOrganRandomParameter(self.hrp)
        root_rp = pb.RootRandomParameter(self.plant)
        root_rp.hyphalEmergenceDensity = 10.0  # high density to have hyphae emerging
        self.plant.setOrganRandomParameter(root_rp)
        self.plant.initialize(True)
        initial_hyphae_count = len(self.plant.getOrgansOfType(pb.hyphae))
        sim_time = 1.0
        dt = 0.1
        steps = int(sim_time / dt)
        for _ in range(steps):
            self.plant.simulate(dt, False)
        final_hyphae_count = len(self.plant.getOrgansOfType(pb.hyphae))
        self.assertGreater(final_hyphae_count, initial_hyphae_count, "Hyphae did not grow as expected")

    def test_branching(self):
        """ tests hyphae branching behavior """
        self.hyphae_example()
        self.hrp.subType = 1  # assuming subType 1 enables branching
        self.plant.setOrganRandomParameter(self.hrp)
        root_rp = pb.RootRandomParameter(self.plant)
        root_rp.hyphalEmergenceDensity = 5.0
        self.plant.setOrganRandomParameter(root_rp)
        self.plant.initialize(True)
        sim_time = 2.0
        dt = 0.1
        steps = int(sim_time / dt)
        for _ in range(steps):
            self.plant.simulate(dt, False)
        hyphae = self.plant.getOrgansOfType(pb.hyphae)
        lateral_branch_points = 0
        tip_split_branch_points = 0
        for h in hyphae:
            if len(h.getChildren()) > 2: # if a hypha has more than 2 children it can only be a hypha with lateral branches, not a tip split
                lateral_branch_points += len(h.getChildren())  # count lateral branches beyond the main two
            elif len(h.getChildren()) == 2: # if a hypha has exactly 2 children, it can be either a tip split or a lateral branch
                if h.getChildren()[0].getOrigin() == h.getChildren()[1].getOrigin(): # if both children originate from the same point, it's a tip split
                    tip_split_branch_points += 1  # count tips that split into two
                else:
                    lateral_branch_points += 1  # count lateral branches
            elif len(h.getChildren()) == 1: # if a hypha has exactly 1 child, it can only be a lateral branch
                lateral_branch_points += 1  # count single lateral branches
        self.assertGreater(tip_split_branch_points, 0, "No tip splitting occurred in hyphae")
        self.assertGreater(lateral_branch_points, 0, "No lateral branching occurred in hyphae")

    def test_tip_splitting_frequency(self):
        """ tests hyphae tip splitting frequency """
        self.hyphae_example()
        self.hrp.subType = 2  # assuming subType 2 enables tip splitting
        self.hrp.b = 2  # set branching parameter to promote tip splitting
        self.hrp.b_prob = 1.0  # set branching probability to 1 to ensure tip splitting occurs
        self.plant.setOrganRandomParameter(self.hrp)
        self.plant.initialize(True)
        sim_time = 2.0
        dt = 0.1
        steps = int(sim_time / dt)
        for _ in range(steps):
            self.plant.simulate(dt, False)
        hyphae = self.plant.getOrgansOfType(pb.hyphae)
        tip_split_count = sum(1 for h in hyphae if len(h.getChildren()) == 2 and h.getChildren()[0].getOrigin() == h.getChildren()[1].getOrigin())
        self.assertEqual(tip_split_count, sim_time*self.hrp.b, "No tip splitting occurred in hyphae") # tip splitting only dependent on time and branching parameter b, so total number of tip splits should be equal to sim_time*b

    def test_lateral_branching_distance(self):
        self.hyphae_example()
        self.hrp.subType = 1  # assuming subType 1 enables lateral branching
        self.hrp.b_prob = 0
        self.plant.setOrganRandomParameter(self.hrp)
        self.plant.initialize(True)
        sim_time = 2.0
        dt = 0.1
        steps = int(sim_time / dt)
        for _ in range(steps):
            self.plant.simulate(dt, False)
        hyphae = self.plant.getOrgansOfType(pb.hyphae)
        lateral_branch_distances = []
        for h in hyphae:
            if len(h.getChildren()) > 1:  # if a hypha has more than one child, it can have lateral branches
                for i in range(1, len(h.getChildren())):  # start from 1 to skip the main child
                    lateral_branch_distances.append(norm(h.getChildren()[i].getOrigin() - h.getChildren()[i-1].getOrigin()))
        for d in lateral_branch_distances:
            self.assertGreaterEqual(d, self.hrp.ln - self.hrp.lns, "Lateral branch distance is greater than minimum distance threshold")
            self.assertLessEqual(d, self.hrp.ln + self.hrp.lns, "Lateral branch distance is less than maximum distance threshold")
    

    def test_hyphal_tree_index(self):
        """ tests hyphal tree index """
        self.hyphae_example()
        self.plant.setOrganRandomParameter(self.hrp)
        self.plant.initialize(True)
        sim_time = 2.0
        dt = 0.1
        steps = int(sim_time / dt)
        for _ in range(steps):
            self.plant.simulate(dt, False)
        hyphae = self.plant.getOrgansOfType(pb.hyphae)
        for h in hyphae:
            if len(h.getChildren()) > 0:  # only check hyphae that have children
                tree_index = [h.getHyphalTreeIndex()]*len(h.getChildren())
                self.assertEqual(tree_index, h.getChildren().getHyphalTreeIndex(), "Hyphal tree index of parent does not match that of child")
            self.assertGreaterEqual(tree_index, 0, "Hyphal tree index is negative")

    def test_merged_hyphae_index(self):
        """ tests merged hyphae index """
        self.hyphae_example()
        self.plant.setOrganRandomParameter(self.hrp)
        self.plant.initialize(True)
        sim_time = 2.0
        dt = 0.1
        steps = int(sim_time / dt)
        for _ in range(steps):
            self.plant.simulate(dt, False)
        hyphae = self.plant.getOrgansOfType(pb.hyphae)
        for h in hyphae:
            ###### This is only a basic check if the merged hyphae index is not negative
            ###### Need to check if the index actually points to a node that exists
            ###### and then need to check if the second to last node is actually within the distance threshold
            ##### this implies a double for loop but this is very expensive and should be avoided if possible
            ###### then should check that the hyphal tree index of merged hypha and original hypha are different
            if h.getParameter("mergePointID") > 0:  # only check hyphae that have merged
                self.assertGreaterEqual(h.getParameter("mergePointID"), 0, "Merged hyphae index is negative")
                self.assertNotEqual(h.getParameter("mergePointID"), h.getHyphalTreeIndex(), "Merged hyphae index should not be equal to hyphal tree index")
  
    def test_anastomosis(self):
        """ tests hyphae anastomosis behavior """
        self.hyphae_example()
        self.hrp.distTH = 0.1  # set distance threshold for anastomosis
        self.plant.setOrganRandomParameter(self.hrp)
        self.hrp.hyphalEmergenceDensity = 20.0  # high density to promote anastomosis
        self.plant.setOrganRandomParameter(self.hrp)
        self.plant.initialize(True)
        sim_time = 3.0
        dt = 0.1
        steps = int(sim_time / dt)
        for _ in range(steps):
            self.plant.simulate(dt, False)
        anastomosis_points = self.plant.getAnastomosisPoints(5)
        self.assertGreater(len(anastomosis_points), 0, "No anastomosis points detected in hyphae")

if __name__ == '__main__':
    unittest.main()       