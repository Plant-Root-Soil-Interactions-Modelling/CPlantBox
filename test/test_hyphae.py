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
        self.hrp.hyphalEmergenceDensity = 2.0
        self.hrp.highresolution = 1

    def test_constructors(self): # TODO check again all parameters needed
        """ tests constructor and copy """
        plant = pb.Organism()
        hrp = pb.HyphaeRandomParameter(plant)
        hrp.dx = 0.15
        hrp.distTH = 0.25
        hrp.subType = 2
        hrp.hyphalEmergenceDensity = 3.0
        hrp.highresolution = 0
        hrp2 = hrp.copy(plant)
        self.assertIsNot(hrp, hrp2, "copy: organ type parameter set is not copied")
        self.assertEqual(hrp2.name, hrp.name, "copy: value unexpected")
        self.assertEqual(hrp2.organType, hrp.organType, "copy: value unexpected")
        self.assertEqual(hrp2.subType, hrp.subType, "copy: value unexpected")
        self.assertEqual(hrp2.dx, hrp.dx, "copy: value unexpected")
        self.assertEqual(hrp2.distTH, hrp.distTH, "copy: value unexpected")
        self.assertEqual(hrp2.hyphalEmergenceDensity, hrp.hyphalEmergenceDensity, "copy: value unexpected")
        self.assertEqual(hrp2.highresolution, hrp.highresolution, "copy: value unexpected")
    
    def test_parameter(self):
        """ tests getParameter() """
        self.plant = pb.Organism()
        hrp = pb.HyphaeRandomParameter(self.plant)
        hrp.dx = 0.2
        hrp.distTH = 0.3
        ot = hrp.getParameter("organType")  # test defaults
        st = hrp.getParameter("subType")
        dx = hrp.getParameter("dx")
        dth = hrp.getParameter("distTH")
        hed = hrp.getParameter("hyphalEmergenceDensity")
        hr = hrp.getParameter("highresolution")
        self.assertEqual(ot, pb.hyphae, "getParameter: value unexpected")
        self.assertEqual(st, 0, "getParameter: value unexpected")
        self.assertEqual(dx, 0.2, "getParameter: value unexpected")
        self.assertEqual(dth, 0.3, "getParameter: value unexpected")
        self.assertEqual(hed, 1.0, "getParameter: value unexpected")
        self.assertEqual(hr, 1, "getParameter: value unexpected")

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
        branch_points = 0
        for h in hyphae:
            if len(h.getChildren()) > 1:
                branch_points += 1
        self.assertGreater(branch_points, 0, "No branching occurred in hyphae")

    def test_anastomosis(self):
        """ tests hyphae anastomosis behavior """
        self.hyphae_example()
        self.hrp.distTH = 0.1  # set distance threshold for anastomosis
        self.plant.setOrganRandomParameter(self.hrp)
        root_rp = pb.RootRandomParameter(self.plant)
        root_rp.hyphalEmergenceDensity = 20.0  # high density to promote anastomosis
        self.plant.setOrganRandomParameter(root_rp)
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