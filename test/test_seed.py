import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
from rsml.rsml_reader import *

import matplotlib.pyplot as plt


class TestRoot(unittest.TestCase):

    def get_srp(self):
        self.plant = pb.Organism()  # need to store this
        srp = pb.SeedRandomParameter(self.plant)
        srp.seedPos = pb.Vector3d(0, 0, -2)
        srp.firstB = 3  # basal
        srp.delayB = 5
        srp.maxB = 7
        srp.nC = 7
        srp.firstSB = 30  # shoot borne
        srp.delaySB = 10
        srp.delayRC = 70
        srp.nz = 0.7
        srp.maxTil = 3  # stem
        srp.simtime = 10
        self.plant.setOrganRandomParameter(srp)
        return srp

    def test_constructors(self):
        """ tests two kinds of constructors and copy """
        srp = self.get_srp()
        # print(srp)
        # 1. constructor from scratch
        param = srp.realize()
        seed = pb.Seed(self.plant.getOrganIndex(), param, True, True, 0, 0, False, 0)
        seed.setOrganism(self.plant)
        # 2. used in simulation (must have parent, since there is no nullptr in Pyhton)
        seed2 = pb.Seed(self.plant)
        # 3. deep copy (with a factory function)
        plant2 = pb.Organism()
        seed3 = seed.copy(plant2)
        self.assertEqual(str(seed), str(seed3), "deep copy: the seed string representations shold be equal")
        self.assertIsNot(seed.getParam(), seed3.getParam(), "deep copy: organs have same parameter set")
        self.assertEqual(str(seed.param()), str(seed3.param()), "deep copy: organs have different parameter values")  # type RootSpecificParameter

    def test_initialize(self):
        """ test initialization (! most important) """
        srp = self.get_srp()
        seed = pb.Seed(self.plant)
        # print(seed)
        tap_rp = pb.RootRandomParameter(self.plant)
        tap_rp.name = "taproot"
        tap_rp.subType = seed.tapType  # 1
        basal_rp = pb.RootRandomParameter(self.plant)
        basal_rp.name = "basalroot"
        basal_rp.subType = seed.basalType  # 4
        shootborne_rp = pb.RootRandomParameter(self.plant)
        shootborne_rp.name = "shootborneroot"
        shootborne_rp.subType = seed.shootborneType  # 5
        self.plant.setOrganRandomParameter(tap_rp)
        self.plant.setOrganRandomParameter(basal_rp)
        self.plant.setOrganRandomParameter(shootborne_rp)
        seed.initialize(False)  # A CRootBox initializations (now StemRandomParameter set)

        # print("CRootBox initialization (including shootborne)")
        organs = seed.baseOrgans()
        st_ = []
        ot_ = []
        for o in organs:
            st_ .append(o.getParameter("subType"))
            ot_ .append(o.getParameter("organType"))
        norc = int(((90 - srp.firstSB) / srp.delayRC) + 0.5)
        self.assertEqual(norc, seed.getNumberOfRootCrowns(), "wrong number of root crowns")
        maxB = min(srp.maxB, int(((365 - srp.firstB) / srp.delayB) + 0.5))
        self.assertEqual(len(st_), norc * srp.nC + 1 + maxB, "wrong number of roots created")

        seed = pb.Seed(self.plant)
        stem_rp = pb.StemRandomParameter(self.plant)
        stem_rp.name = "stem"
        stem_rp.subType = seed.mainStemType  # 1
        tiller_rp = pb.StemRandomParameter(self.plant)
        tiller_rp.name = "tiller"
        tiller_rp.subType = seed.tillerType  # 4

        self.plant.setOrganRandomParameter(stem_rp)
        self.plant.setOrganRandomParameter(tiller_rp)
        seed.initialize(False)  # A CPlantBox initializations (now StemRandomParameter set)
        organs = seed.baseOrgans()

        # print("CPlantBox initialization (no shootborne, since they emerge from the stem)")
        st_ = []
        ot_ = []
        for o in organs:
            # print(o)
            st_ .append(o.getParameter("subType"))
            ot_ .append(o.getParameter("organType"))
        st = [1.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 1.0, 4.0, 4.0, 4.0]
        ot = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0]

        for i in range(0, len(organs)):
            self.assertEqual(st_[i], st[i], "no tap root produced")
            self.assertEqual(ot_[i], ot[i], "tap root produced wrong organ type")

    def test_advanced_initialize(self):
        """ tests if factory functions can be overwritten in Pyhton """
        # todo
        pass


if __name__ == '__main__':
    unittest.main()
