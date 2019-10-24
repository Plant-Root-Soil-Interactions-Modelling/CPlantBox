import unittest
import sys
sys.path.append("..")
import plantbox as pb
from rsml import *

path = "../modelparameter/plant/"


class TestPlant(unittest.TestCase):

    def test_CPlantBox(self):
        """tests the CPlantBox function defined in CPlantBox_PiafMunch.py"""
        p = pb.Plant()
        p.openXML(path + "Heliantus_Pagès_2013_new.xml")

        seeds = p.getOrganRandomParameter(pb.OrganTypes.seed)
        roots = p.getOrganRandomParameter(pb.OrganTypes.root)
        stems = p.getOrganRandomParameter(pb.OrganTypes.stem)
        leafs = p.getOrganRandomParameter(pb.OrganTypes.leaf)
#         for p_ in seeds:
#             print(p_)
#         for p_ in roots[1:]:
#             print(p_)
#         for p_ in stems[1:]:
#             print(p_)
#         for p_ in leafs[1:]:
#             print(p_)

        self.assertEqual([len(seeds), len(roots[1:]), len(stems[1:]), len(leafs[1:])], [1, 3, 3, 1],
                         "test_CPlantBox: read wrong number of random parameter from xml")

        p.initialize(True)
        p.simulate(76)
        p.write("morningglory.vtp")

#     def test_CPlantBox_analysis(self):
#         """tests the CPlantBox_analysis function defined in CPlantBox_PiafMunch.py"""
#         p = pb.Plant()
#         p.openXML(path + "Heliantus_Pagès_2013_new.xml")
#         p.initialize()
#         p.simulate(76)
#         ana = pb.SegmentAnalyser(p)
#         ana.write(write("morningglory_ama.vtp"))

    def test_convert(self):
        """tests the functions needed by the convert function of CPlantBox_PiafMunch.py"""
        pass


if __name__ == '__main__':
    # MANY tests missing !!!

    unittest.main()
