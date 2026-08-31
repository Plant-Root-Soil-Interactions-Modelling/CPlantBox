import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
import numpy as np

# from rsml.rsml_reader import *

class TestHyphaeParameter(unittest.TestCase):

    def hyphae_example(self): # TODO check which parameters are needed
        self.plant = pb.Organism()
        self.hrp = pb.HyphaeRandomParameter(self.plant)
        self.hrp.dx = 0.1
        self.hrp.subType = 1
        self.hrp.v = 0.12
        self.hrp.hlt = 5
        self.theta = 15
        self.hrp.distTH = 0.2
        self.hrp.ana = 0.5
        self.hrp.lb = 0.5
        self.hrp.la = 0.5
        self.hrp.ln = 0.5
        self.hrp.lmax = 5.0
        self.hrp.b = 1
        self.hrp.b_prob = 1.0

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
        self.assertEqual(hrp2.v, hrp.v, "copy: value unexpected")
        self.assertEqual(hrp2.hlt, hrp.hlt, "copy: value unexpected")
        self.assertEqual(hrp2.theta, hrp.theta, "copy: value unexpected")
        self.assertEqual(hrp2.ana, hrp.ana, "copy: value unexpected")
        self.assertEqual(hrp2.lb, hrp.lb, "copy: value unexpected")
        self.assertEqual(hrp2.la, hrp.la, "copy: value unexpected")
        self.assertEqual(hrp2.b, hrp.b, "copy: value unexpected")
        self.assertEqual(hrp2.b_prob, hrp.b_prob, "copy: value unexpected")

    def test_parameter(self):
        """ tests getParameter() """
        self.plant = pb.Organism()
        hrp = pb.HyphaeRandomParameter(self.plant)
        hrp.dx = 0.2
        hrp.distTH = 0.3
        hrp.subType = 1
        ot = hrp.getParameter("organType")  # test defaults
        st = hrp.getParameter("subType")
        dth = hrp.getParameter("distTH")
        v = hrp.getParameter("v")
        hlt = hrp.getParameter("hlt")
        theta = hrp.getParameter("theta")
        ana = hrp.getParameter("ana")
        lb = hrp.getParameter("lb")
        la = hrp.getParameter("la")
        ln = hrp.getParameter("ln")
        lmax = hrp.getParameter("lmax")
        b = hrp.getParameter("b")
        b_prob = hrp.getParameter("b_prob")
        self.assertEqual(ot, pb.hyphae, "getParameter: value unexpected")
        self.assertEqual(st, 1, "getParameter: value unexpected")
        self.assertEqual(dth, 0.3, "getParameter: value unexpected")
        self.assertEqual(v, 0.13, "getParameter: value unexpected")
        self.assertEqual(hlt, 10, "getParameter: value unexpected")
        self.assertEqual(theta, 30/180*np.pi, "getParameter: value unexpected")
        self.assertEqual(ana, 1.0, "getParameter: value unexpected")
        self.assertEqual(lb, 0.001, "getParameter: value unexpected")
        self.assertEqual(la, 0.003, "getParameter: value unexpected")
        self.assertEqual(ln, 0.005, "getParameter: value unexpected")
        self.assertEqual(lmax, 10.0, "getParameter: value unexpected")
        self.assertEqual(b, 2, "getParameter: value unexpected")
        self.assertEqual(b_prob, 0, "getParameter: value unexpected")
    
    def test_xml(self): # TODO check again all parameters needed
        """ tests reading and writing xml parameter file """
        self.hyphae_example()
        hrp = self.hrp
        # write parameters to xml
        hrp.writeXML("hyphae.xml")
        # read parameters from xml
        hrp2 = pb.HyphaeRandomParameter(self.plant)
        hrp2.readXML("hyphae.xml")
        # hrp2 = plant2.getOrganRandomParameter(pb.hyphae)
        self.assertEqual(hrp2.distTH, self.hrp.distTH, "xml read/write: value unexpected")
        self.assertEqual(hrp2.subType, self.hrp.subType, "xml read/write: value unexpected")
        self.assertEqual(hrp2.nob(), self.hrp.nob(), "xml read/write: value unexpected")
        self.assertEqual(hrp2.lns, self.hrp.lns, "xml read/write: value unexpected")
        self.assertEqual(hrp2.v, self.hrp.v, "xml read/write: value unexpected")
        self.assertEqual(hrp2.hlt, self.hrp.hlt, "xml read/write: value unexpected")
        self.assertEqual(hrp2.successor[0][0], self.hrp.successor[0][0], "xml read/write: value unexpected")


    def test_realize(self):
        """ tests if hyphae parameter can be set and realized in a plant """
        self.hyphae_example()
        self.plant.initialize(True)
        # set the hyphae parameter to the plant
        self.plant.setOrganRandomParameter(self.hrp)
        # check if the parameter is realized in the plant
        hrp_plant = self.plant.getOrganRandomParameter(pb.hyphae)
        self.assertEqual(hrp_plant.distTH, self.hrp.distTH, "realize: value unexpected")
        self.assertEqual(hrp_plant.subType, self.hrp.subType, "realize: value unexpected")

if __name__ == '__main__':
    unittest.main()