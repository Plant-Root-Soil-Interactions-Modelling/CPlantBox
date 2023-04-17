import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
from rsml.rsml_reader import *

path = "../modelparameter/structural/plant/"


def rootLength(t, r, k):  # root length at a certain age
    return k * (1 - np.exp(-r * t / k))


class TestPlant(unittest.TestCase):

    def test_CPlantBox(self):
        """tests the functions needed by CPlantBox defined in CPlantBox_PiafMunch.py"""
        p = pb.Plant()
        p.readParameters(path + "Heliantus_Pagès_2013.xml",fromFile = True, verbose = False)

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

        p.initialize(False)
        p.simulate(76, False)
        p.write("morningglory.vtp")

    def test_CPlantBox_analysis(self):
        """tests the functions needed by CPlantBox_analysis defined in CPlantBox_PiafMunch.py"""
        p = pb.Plant()
        p.readParameters(path + "Heliantus_Pagès_2013.xml", fromFile = True,verbose = False)
        p.initialize(False)
        p.simulate(76)
        ana = pb.SegmentAnalyser(p)
        ana.write("morningglory_ama.vtp")

    def test_convert(self):
        """tests the functions needed by the convert function of CPlantBox_PiafMunch.py"""
        p = pb.Plant()
        p.readParameters(path + "Heliantus_Pagès_2013.xml",fromFile = True, verbose = False)
        p.initialize(False)
        p.simulate(76)
        nodes = np.array([np.array(a) / 100 for a in p.getNodes()])  # convert to numpy array, and from cm to m
        #print(nodes.shape)
        rseg = np.array([np.array(s) for s in p.getSegments(pb.OrganTypes.root)])  # root system segments
        #print(rseg.shape)
        sseg = np.array([np.array(s) for s in p.getSegments(pb.OrganTypes.stem)])  # stem system segments
        #print(sseg.shape)
#         lseg = v2ai(plant.getNodesOrganType())
        l = np.array([ o.getParameter("organType") for o in p.getSegmentOrigins()])
        #print(l.shape)
#         plant_ana = pb.SegmentAnalyser(p)
#         node_connection_o = seg2a(p.getSegments(15)) # plant segments
        pass
        
    def test_CPlantBox_step(self):
        """tests the functions needed by CPlantBox defined in CPlantBox_PiafMunch.py"""
        p1 = pb.MappedPlant(2)
        p1.readParameters(path + "Heliantus_Pagès_2013.xml",fromFile = True, verbose = False)
        p1.initialize(verbose = False, stochastic = False)
        p1.simulate(76, False)
        p1.write("test_CPlantBox_1step.vtp")
        
        p2 = pb.MappedPlant(2)
        p2.readParameters(path + "Heliantus_Pagès_2013.xml",fromFile = True, verbose = False)
        p2.initialize(verbose = False, stochastic = False)
        for i in range(100) :
            p2.simulate(76/100, False)
        p2.write("test_CPlantBox_100step.vtp")
        root1 = p1.getOrgans(2)
        root2 = p2.getOrgans(2)
        self.assertAlmostEqual(len(root1), len(root2), 10, "number of root organs do not agree")
        root1 = p1.getOrgans(3)
        root2 = p2.getOrgans(3)
        self.assertAlmostEqual(len(root1), len(root2), 10, "number of stem organs do not agree")
        root1 = p1.getOrgans(4)
        root2 = p2.getOrgans(4)
        self.assertAlmostEqual(len(root1), len(root2), 10, "number of leaf organs do not agree")                   

    def test_DB_delay(self):
        p = pb.MappedPlant(2)
        p.readParameters(path + "Heliantus_Pagès_2013.xml",fromFile = True, verbose = False)
        
        p.initializeDB()
        time = 76
        
        rrp = p.getOrganRandomParameter(pb.root)[1]
        rrp.ldelay = 3
        rrp.ldelays = 0.
        #print(rrp)
        p.setOrganRandomParameter(rrp)
        
        rrp = p.getOrganRandomParameter(pb.root)[2]
        rrp.ldelay = 5
        rrp.ldelays = 0.
        #print(rrp)
        p.setOrganRandomParameter(rrp)
        
        rrp = p.getOrganRandomParameter(pb.root)[3]
        rrp.ldelay = 5
        rrp.ldelays = 0.
        #print(rrp)
        p.setOrganRandomParameter(rrp)
        tl, rl = [], []
        p.simulate(time, False)
        for i, r in enumerate(p.getOrgans(pb.root)):
            rl.append(r.getLength())
            et, dl = 0, 0  # no delay for basal roots
            rsp = r.getParam()
            if rsp.subType > 1:
                dl = r.getParent().getOrganRandomParameter().ldelay  # only works because deviation == 0
                et = r.getParent().getNodeCT(r.parentNI) + dl
            self.assertAlmostEqual(r.getAge(), (time - et), 10, "numeric and analytic age of root n#" + str(i + 1) + " do not agree")
    def test_missing_laterals(self):
        """can CPB handle it when no lateral appear at a branching point
            because of probabilistic branching?
            currently, just tests that cpb does not crash
        """
        p = pb.MappedPlant(2)
        p.readParameters(path + "test_missing_laterals.xml", fromFile = True,verbose = False)        
        p.initialize(False)
        time = 76
        p.simulate(time, False)
        p.write("test_missing_laterals.vtp")


if __name__ == '__main__':
    # MANY tests missing !!!

    unittest.main()
