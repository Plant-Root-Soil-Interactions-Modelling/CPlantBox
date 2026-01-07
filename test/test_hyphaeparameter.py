import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
from rsml.rsml_reader import *

class TestHyphaeParameter(unittest.TestCase):

    def hyphae_lowres_example(self): # TODO check which parameters are needed
        self.plant = pb.Organism()
        self.hrp = pb.HyphaeRandomParameter(self.plant)
        self.hrp.dx = 0.1
        self.hrp.distTH = 0.2
        self.hrp.subType = 1
    
    def hyphae_highres_example(self): # TODO check which parameters are needed
        self.plant = pb.Organism()
        self.hrp = pb.HyphaeRandomParameter(self.plant)
        self.hrp.dx = 0.1
        self.hrp.distTH = 0.15
        self.hrp.subType = 2

    def test_constructors(self): # TODO check again all parameters needed
        """ tests constructor and copy """
        plant = pb.Organism()
        hrp = pb.HyphaeRandomParameter(plant)
        hrp.dx = 0.15
        hrp.distTH = 0.25
        hrp.subType = 1
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
        self.assertEqual(ot, pb.hyphae, "getParameter: value unexpected")
        self.assertEqual(st, 0, "getParameter: value unexpected")
        self.assertEqual(dx, 0.2, "getParameter: value unexpected")
        self.assertEqual(dth, 0.3, "getParameter: value unexpected")
    
    def test_xml(self): # TODO check again all parameters needed
        """ tests reading and writing xml parameter file """
        self.hyphae_lowres_example()
        # write parameters to xml
        self.plant.writeParameters("hyphae_test_parameters.xml", 'plant', True)
        # read parameters from xml
        plant2 = pb.Organism()
        plant2.readParameters("hyphae_test_parameters.xml", fromFile = True, verbose = True)
        hrp2 = plant2.getOrganRandomParameter(pb.hyphae)[0]
        self.assertEqual(hrp2.dx, self.hrp.dx, "xml read/write: value unexpected")
        self.assertEqual(hrp2.distTH, self.hrp.distTH, "xml read/write: value unexpected")
        self.assertEqual(hrp2.subType, self.hrp.subType, "xml read/write: value unexpected")

    def test_realize(self):
        """ tests if hyphae parameter can be set and realized in a plant """
        self.hyphae_lowres_example()
        self.plant.initialize(True)
        # set the hyphae parameter to the plant
        self.plant.setOrganRandomParameter(self.hrp)
        # check if the parameter is realized in the plant
        hrp_plant = self.plant.getOrganRandomParameter(pb.hyphae)[0]
        self.assertEqual(hrp_plant.dx, self.hrp.dx, "realize: value unexpected")
        self.assertEqual(hrp_plant.distTH, self.hrp.distTH, "realize: value unexpected")
        self.assertEqual(hrp_plant.subType, self.hrp.subType, "realize: value unexpected")

    if __name__ == '__main__':
        unittest.main()