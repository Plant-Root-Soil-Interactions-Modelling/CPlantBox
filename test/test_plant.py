import unittest
import sys
sys.path.append("..")
import plantbox as pb
from rsml import *

path = "../modelparameter/plant/"


class TestPlant(unittest.TestCase):

    def test_CPlantBox(self):
        """tests the functions needed by CPlantBox defined in CPlantBox_PiafMunch.py"""
        p = pb.Plant()
        p.openXML(path + "Heliantus_Pagès_2013.xml")

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

    def test_CPlantBox_analysis(self):
        """tests the functions needed by CPlantBox_analysis defined in CPlantBox_PiafMunch.py"""
        p = pb.Plant()
        p.openXML(path + "Heliantus_Pagès_2013.xml")
        p.initialize()
        p.simulate(76)
        ana = pb.SegmentAnalyser(p)
        ana.write("morningglory_ama.vtp")

    def test_convert(self):
        """tests the functions needed by the convert function of CPlantBox_PiafMunch.py"""
        p = pb.Plant()
        p.openXML(path + "Heliantus_Pagès_2013.xml")
        p.initialize()
        p.simulate(76)
        nodes = np.array([np.array(a)/100 for a in p.getNodes()]) # convert to numpy array, and from cm to m 
        print(nodes.shape)         
        rseg = np.array([np.array(s) for s in p.getSegments(pb.OrganTypes.root)]) # root system segments
        print(rseg.shape)
        sseg = np.array([np.array(s) for s in p.getSegments(pb.OrganTypes.stem)]) # stem system segments
        print(sseg.shape)
#         lseg = v2ai(plant.getNodesOrganType())
        l = np.array([ o.getParameter("organType") for o in p.getSegmentOrigins()])        
        print(l.shape)
#         plant_ana = pb.SegmentAnalyser(p) 
#         node_connection_o = seg2a(p.getSegments(15)) # plant segments        
        pass


if __name__ == '__main__':
    # MANY tests missing !!!

    unittest.main()
