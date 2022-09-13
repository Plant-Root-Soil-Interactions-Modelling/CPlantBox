import unittest
import sys; sys.path.append(".."); sys.path.append("../src/python_modules")
sys.path.append("src/python_modules")
import plantbox as pb
from importlib import reload  # Python 3.4+
pb = reload(pb)
from rsml_reader import *
import vtk_plot as vp
import os
currentpath = os.getcwd()
if currentpath[-4:] == "tBox":
    path = "/modelparameter/plant/"
else:
    path = "/../modelparameter/plant/"



def rootLength(t, r, k):  # root length at a certain age
    return k * (1 - np.exp(-r * t / k))   
class TestPlant(unittest.TestCase):

    def test_CPlantBox(self):
        """tests the functions needed by CPlantBox defined in CPlantBox_PiafMunch.py"""
        p = pb.MappedPlant()
        print(os.getcwd() +path + "Heliantus_Pagès_2013.xml")
        p.openXML(os.getcwd() +path + "Heliantus_Pagès_2013.xml")

        p.initialize(True)
        #p.simulate(100, True)

        #for p_ in p.getOrganRandomParameter(pb.stem):
        #    if (p_.subType > 0):
        #        print(p_.subType, "radius", p_.a, "lmax", p_.lmax, p_.ln, p_.lb,  p_.nob())

        # nodes = p.getNodes()
        # segseg = p.getSegments()
        # orgs = p.getOrgans()
        # p.getNumberOfOrgans()
        #print(p.segLength())
        # for org in orgs:
        #     print(org.getId(), org.organType(), org.getOrigin(), org.getNumberOfNodes(),
        #           org.getNumberOfChildren(), org.getNodeId(0))
        #     print("     ",org.getNode(org.getNumberOfNodes() - 1))
        #
        #     nns = [ np.array(org.getNode(n_))
        #             for n_ in range(org.getNumberOfNodes())]
        #     print(nns)
        #     for nkid in range(org.getNumberOfChildren()):
        #         orgkid = org.getChild(nkid)
        #         print("kid :",orgkid.getId(), orgkid.organType(), orgkid.getOrigin(), orgkid.getNumberOfNodes(),
        #               orgkid.getNumberOfChildren(),orgkid.getNodeId(0))
        #         print("     ",orgkid.getNode(orgkid.getNumberOfNodes() - 1))
        #     #print(nns)
        # p.write("morningglory_" + str(0) + ".vtp")
        # raise Exception
        for i in range(76):
            p.simulate(1, True)
            p.write("morningglory_"+ str(i) +".vtp")
        vp.plot_plant(p,p_name = "organType")

    def _test_CPlantBox_analysis(self):
        """tests the functions needed by CPlantBox_analysis defined in CPlantBox_PiafMunch.py"""
        p = pb.Plant()
        p.openXML(path + "Heliantus_Pagès_2013.xml")
        p.initialize()
        p.simulate(76)
        ana = pb.SegmentAnalyser(p)
        ana.write("morningglory_ama.vtp")

    def _test_convert(self):
        """tests the functions needed by the convert function of CPlantBox_PiafMunch.py"""
        p = pb.Plant()
        p.openXML(path + "Heliantus_Pagès_2013.xml")
        p.initialize()
        p.simulate(76)
        nodes = np.array([np.array(a) / 100 for a in p.getNodes()])  # convert to numpy array, and from cm to m 
        print(nodes.shape)         
        rseg = np.array([np.array(s) for s in p.getSegments(pb.OrganTypes.root)])  # root system segments
        print(rseg.shape)
        sseg = np.array([np.array(s) for s in p.getSegments(pb.OrganTypes.stem)])  # stem system segments
        print(sseg.shape)
#         lseg = v2ai(plant.getNodesOrganType())
        l = np.array([ o.getParameter("organType") for o in p.getSegmentOrigins()])        
        print(l.shape)
#         plant_ana = pb.SegmentAnalyser(p) 
#         node_connection_o = seg2a(p.getSegments(15)) # plant segments        
        pass

    def _test_DB_delay(self):
        p = pb.MappedPlant()
        p.readParameters(path + "Heliantus_Pagès_2013.xml")
        rrp = p.getOrganRandomParameter(pb.root)[1]
        rrp.ldelay = 3
        rrp = p.getOrganRandomParameter(pb.root)[2]
        rrp.ldelay = 5
        p.setOrganRandomParameter(rrp)
        p.initializeDB()
        time = 76
        p.simulate(time)
        tl, rl = [], []
        for i, r in enumerate(p.getOrgans(pb.root)):
            rl.append(r.getLength())
            et, dl = 0, 0 #no delay for basal roots
            rsp = r.getParam()
            if rsp.subType > 1:
                dl = r.getParent().getOrganRandomParameter().ldelay #only works because deviation == 0
                et = r.getParent().getNodeCT(r.parentNI) + dl
            self.assertAlmostEqual(r.getAge(), (time -et), 10, "numeric and analytic age of root n#" + str(i + 1) +" do not agree")

if __name__ == '__main__':
    # MANY tests missing !!!

    unittest.main()
