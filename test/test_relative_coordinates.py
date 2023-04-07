"""
Copyright 2019, Forschungszentrum Jülich GmbH, licensed under GNU GPLv3
"""

"""
Objective: 
Make a plant grow where the root system (absolute coordinates) is supposed to be symetric
to the aboveground part (relative coordinates during Plant::simulate()) 
according to the ground (surface where (x, y, 0))
it is assumed that if node at the tip of leaves and 3rd root lateral have symetric coordinates, 
then the rest of the root system is symetric to the aboveground organs.

"""
import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt


class TestRelCoord(unittest.TestCase):

    def test_coord_100steps(self):
        pl = pb.MappedPlant()
        path = "../modelparameter/structural/plant/"
        name = "test_relcoord"
        pl.readParameters(path + name + ".xml")
        # stochastic = False => thus rand() always give 0.5 for Tropism. => result only change according to sigma
        pl.initialize(stochastic = False)
        
        dt = 1
        steps = 100
        
        for step in range(steps):
            pl.simulate(1, False)
                               
        pl.write("test_relcoord_100steps.vtp")
        params = pl.organParam
        seedPosx = params[1][0].seedPos.x
        seedPosy = params[1][0].seedPos.y
        seedPosz = params[1][0].seedPos.z
        # print(params[1][0].seedPos)
        roots = pl.getOrgans(2)
        leaves = pl.getOrgans(4)
        allStems=pl.getOrgans(3)
        mainStem = allStems[0]
        stems = allStems[1:]          

        rootSubtypes = [ o.param().subType for o in roots]
                                                            

        roots1 = np.array(roots)[np.where([st == 1 for st in rootSubtypes])[0]]
        roots2 = np.array(roots)[np.where([st == 2 for st in rootSubtypes])[0]]
        roots3 = np.array(roots)[np.where([st == 3 for st in rootSubtypes])[0]]

                                                                                 
                                                                                 

        rootTipsX = [np.array(r.getNode(r.getNumberOfNodes() - 1))[0] for r in roots1]
        rootTipsY = [np.array(r.getNode(r.getNumberOfNodes() - 1))[1] for r in roots1]
        rootTipsZ = [np.array(r.getNode(r.getNumberOfNodes() - 1))[2] for r in roots1]
        stemTipsX = [np.array(r.getNode(r.getNumberOfNodes() - 1))[0] for r in [mainStem]]
        stemTipsY = [np.array(r.getNode(r.getNumberOfNodes() - 1))[1] for r in [mainStem]]
        stemTipsZ = [np.array(r.getNode(r.getNumberOfNodes() - 1))[2] for r in [mainStem]]
        for i in range(0, len(rootTipsX)):
            self.assertAlmostEqual(rootTipsX[i] - seedPosx, -(stemTipsX[i] - seedPosx), 10, "coord X for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsY)):
            # print(rootTipsY[i]*10**16, leafTipsY[i]*10**16)
            self.assertAlmostEqual(rootTipsY[i] * 10 ** 16, stemTipsY[i] * 10 ** 16, 10, "coord Y for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsZ)):
            # print(rootTipsZ[i]-seedPosz, leafTipsZ[i]-seedPosz, seedPosz)
            self.assertAlmostEqual(rootTipsZ[i] - seedPosz, -(stemTipsZ[i] - seedPosz), 10, "coord Z for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        
        rootTipsX = [np.array(r.getNode(0))[0] for r in roots2]
        rootTipsY = [np.array(r.getNode(0))[1] for r in roots2]
        rootTipsZ = [np.array(r.getNode(0))[2] for r in roots2]
        stemTipsX = [np.array(r.getNode(0))[0] for r in stems]
        stemTipsY = [np.array(r.getNode(0))[1] for r in stems]
        stemTipsZ = [np.array(r.getNode(0))[2] for r in stems]
        for i in range(0, len(rootTipsX)):
            self.assertAlmostEqual(rootTipsX[i] - seedPosx, -(stemTipsX[i] - seedPosx),5, "coord X for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsY)):
            # print(rootTipsY[i]*10**16, leafTipsY[i]*10**16)
            self.assertAlmostEqual(rootTipsY[i] , stemTipsY[i], 5, "coord Y for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsZ)):
            # print(rootTipsZ[i]-seedPosz, leafTipsZ[i]-seedPosz, seedPosz)
            self.assertAlmostEqual(rootTipsZ[i] - seedPosz, -(stemTipsZ[i] - seedPosz), 10, "coord Z for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        
        rootTipsX = [np.array(r.getNode(r.getNumberOfNodes() - 1))[0] for r in roots2]
        rootTipsY = [np.array(r.getNode(r.getNumberOfNodes() - 1))[1] for r in roots2]
        rootTipsZ = [np.array(r.getNode(r.getNumberOfNodes() - 1))[2] for r in roots2]
        stemTipsX = [np.array(r.getNode(r.getNumberOfNodes() - 1))[0] for r in stems]
        stemTipsY = [np.array(r.getNode(r.getNumberOfNodes() - 1))[1] for r in stems]
        stemTipsZ = [np.array(r.getNode(r.getNumberOfNodes() - 1))[2] for r in stems]
        for i in range(0, len(rootTipsX)):
            # print(rootTipsX[i], leafTipsX[i])
            self.assertAlmostEqual(rootTipsX[i] - seedPosx, -(stemTipsX[i] - seedPosx),5, "coord X for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsY)):
            # print(rootTipsY[i]*10**16, leafTipsY[i]*10**16)
            self.assertAlmostEqual(rootTipsY[i] , stemTipsY[i], 5, "coord Y for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsZ)):
            # print(rootTipsZ[i]-seedPosz, leafTipsZ[i]-seedPosz, seedPosz)
            self.assertAlmostEqual(rootTipsZ[i] - seedPosz, -(stemTipsZ[i] - seedPosz), 5, "coord Z for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        
        
        rootTipsX = [np.array(r.getNode(r.getNumberOfNodes() - 1))[0] for r in roots3]
        rootTipsY = [np.array(r.getNode(r.getNumberOfNodes() - 1))[1] for r in roots3]
        rootTipsZ = [np.array(r.getNode(r.getNumberOfNodes() - 1))[2] for r in roots3]

        leafTipsX = [np.array(r.getNode(r.getNumberOfNodes() - 1))[0] for r in leaves]
        leafTipsY = [np.array(r.getNode(r.getNumberOfNodes() - 1))[1] for r in leaves]
        leafTipsZ = [np.array(r.getNode(r.getNumberOfNodes() - 1))[2] for r in leaves]

        for i in range(0, len(rootTipsX)-1):
                                               
            self.assertAlmostEqual(rootTipsX[i] - seedPosx, -(leafTipsX[i] - seedPosx), 10, "coord X for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsY)):
            # print(rootTipsY[i]*10**16, leafTipsY[i]*10**16)
            self.assertAlmostEqual(rootTipsY[i] * 10 ** 16, leafTipsY[i] * 10 ** 16, 10, "coord Y for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsZ)):
            # print(rootTipsZ[i]-seedPosz, leafTipsZ[i]-seedPosz, seedPosz)
            self.assertAlmostEqual(rootTipsZ[i] - seedPosz, -(leafTipsZ[i] - seedPosz), 10, "coord Z for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")

def test_coord_1step(self):
        pl = pb.MappedPlant()
        path = "../modelparameter/structural/plant/"
        name = "test_relcoord"
        pl.readParameters(path + name + ".xml")
        # stochastic = False => thus rand() always give 0.5 for Tropism. => result only change according to sigma
        pl.initialize(stochastic = False)
        pl.simulate(100, False)
        pl.write("test_relcoord_1step.vtp")
        params = pl.organParam
        seedPosx = params[1][0].seedPos.x
        seedPosy = params[1][0].seedPos.y
        seedPosz = params[1][0].seedPos.z
        # print(params[1][0].seedPos)
        roots = pl.getOrgans(2)
        leaves = pl.getOrgans(4)
        allStems=pl.getOrgans(3)
        mainStem = allStems[0]
        stems = allStems[1:]
        
        #print("stems non",[r.getNumberOfNodes() for r in allStems])
        #print("stems length",[r.getLength() for r in allStems])
        #print("stems segLength",[np.array(stems[0].getSegments())])
        segNodeId = np.array(list(map(lambda x: np.array(x), stems[0].getSegments())), dtype = np.int64)
        nodesId = stems[0].getNodeIds()
        nodes =  np.array(list(map(lambda x: np.array(x), [stems[0].getNode(i) for i in range(stems[0].getNumberOfNodes())])))
        #print("nodes",nodes)
        poly = np.zeros((stems[0].getNumberOfNodes(), 3))  #
        for i in range( stems[0].getNumberOfNodes()):
            v = stems[0].getNode(i)
            poly[i, 0] = v.x
            poly[i, 1] = v.y
            poly[i, 2] = v.z
        d = np.diff(poly, axis = 0)
        sd = np.sqrt((d ** 2).sum(axis = 1))
        #print("stems segLength",sd)
        
        #print("leaves pni",[r.parentNI for r in leaves])

        rootSubtypes = [ o.param().subType for o in roots]
                                                            

        roots1 = np.array(roots)[np.where([st == 1 for st in rootSubtypes])[0]]
        roots2 = np.array(roots)[np.where([st == 2 for st in rootSubtypes])[0]]
        roots3 = np.array(roots)[np.where([st == 3 for st in rootSubtypes])[0]]

                                                                                 
                                                                                 

        rootTipsX = [np.array(r.getNode(r.getNumberOfNodes() - 1))[0] for r in roots1]
        rootTipsY = [np.array(r.getNode(r.getNumberOfNodes() - 1))[1] for r in roots1]
        rootTipsZ = [np.array(r.getNode(r.getNumberOfNodes() - 1))[2] for r in roots1]
        stemTipsX = [np.array(r.getNode(r.getNumberOfNodes() - 1))[0] for r in [mainStem]]
        stemTipsY = [np.array(r.getNode(r.getNumberOfNodes() - 1))[1] for r in [mainStem]]
        stemTipsZ = [np.array(r.getNode(r.getNumberOfNodes() - 1))[2] for r in [mainStem]]
        for i in range(0, len(rootTipsX)):
            self.assertAlmostEqual(rootTipsX[i] - seedPosx, -(stemTipsX[i] - seedPosx), 10, "coord X for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsY)):
            # print(rootTipsY[i]*10**16, leafTipsY[i]*10**16)
            self.assertAlmostEqual(rootTipsY[i] * 10 ** 16, stemTipsY[i] * 10 ** 16, 10, "coord Y for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsZ)):
            #print(rootTipsZ[i], stemTipsZ[i],rootTipsZ[i]-seedPosz, stemTipsZ[i]-seedPosz, seedPosz)
            self.assertAlmostEqual(rootTipsZ[i] - seedPosz, -(stemTipsZ[i] - seedPosz), 10, "coord Z for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        
        rootTipsX = [np.array(r.getNode(0))[0] for r in roots2]
        rootTipsY = [np.array(r.getNode(0))[1] for r in roots2]
        rootTipsZ = [np.array(r.getNode(0))[2] for r in roots2]
        stemTipsX = [np.array(r.getNode(0))[0] for r in stems]
        stemTipsY = [np.array(r.getNode(0))[1] for r in stems]
        stemTipsZ = [np.array(r.getNode(0))[2] for r in stems]
        for i in range(0, len(rootTipsX)):
            self.assertAlmostEqual(rootTipsX[i] - seedPosx, -(stemTipsX[i] - seedPosx),5, "coord X for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsY)):
            # print(rootTipsY[i]*10**16, leafTipsY[i]*10**16)
            self.assertAlmostEqual(rootTipsY[i] , stemTipsY[i], 5, "coord Y for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsZ)):
            #print(rootTipsZ[i], stemTipsZ[i],rootTipsZ[i]-seedPosz, stemTipsZ[i]-seedPosz, seedPosz)
            self.assertAlmostEqual(rootTipsZ[i] - seedPosz, -(stemTipsZ[i] - seedPosz), 10, "coord Z for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        
        rootTipsX = [np.array(r.getNode(r.getNumberOfNodes() - 1))[0] for r in roots2]
        rootTipsY = [np.array(r.getNode(r.getNumberOfNodes() - 1))[1] for r in roots2]
        rootTipsZ = [np.array(r.getNode(r.getNumberOfNodes() - 1))[2] for r in roots2]
        stemTipsX = [np.array(r.getNode(r.getNumberOfNodes() - 1))[0] for r in stems]
        stemTipsY = [np.array(r.getNode(r.getNumberOfNodes() - 1))[1] for r in stems]
        stemTipsZ = [np.array(r.getNode(r.getNumberOfNodes() - 1))[2] for r in stems]
        for i in range(0, len(rootTipsX)):
            # print(rootTipsX[i], leafTipsX[i])
            self.assertAlmostEqual(rootTipsX[i] - seedPosx, -(stemTipsX[i] - seedPosx),5, "coord X for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsY)):
            # print(rootTipsY[i]*10**16, leafTipsY[i]*10**16)
            self.assertAlmostEqual(rootTipsY[i] , stemTipsY[i], 5, "coord Y for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsZ)):
            # print(rootTipsZ[i]-seedPosz, leafTipsZ[i]-seedPosz, seedPosz)
            self.assertAlmostEqual(rootTipsZ[i] - seedPosz, -(stemTipsZ[i] - seedPosz), 5, "coord Z for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        
        
        rootTipsX = [np.array(r.getNode(r.getNumberOfNodes() - 1))[0] for r in roots3]
        rootTipsY = [np.array(r.getNode(r.getNumberOfNodes() - 1))[1] for r in roots3]
        rootTipsZ = [np.array(r.getNode(r.getNumberOfNodes() - 1))[2] for r in roots3]

        leafTipsX = [np.array(r.getNode(r.getNumberOfNodes() - 1))[0] for r in leaves]
        leafTipsY = [np.array(r.getNode(r.getNumberOfNodes() - 1))[1] for r in leaves]
        leafTipsZ = [np.array(r.getNode(r.getNumberOfNodes() - 1))[2] for r in leaves]

        for i in range(0, len(rootTipsX)-1):
                                               
            self.assertAlmostEqual(rootTipsX[i] - seedPosx, -(leafTipsX[i] - seedPosx), 10, "coord X for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsY)):
            # print(rootTipsY[i]*10**16, leafTipsY[i]*10**16)
            self.assertAlmostEqual(rootTipsY[i] * 10 ** 16, leafTipsY[i] * 10 ** 16, 10, "coord Y for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsZ)):
            # print(rootTipsZ[i]-seedPosz, leafTipsZ[i]-seedPosz, seedPosz)
            self.assertAlmostEqual(rootTipsZ[i] - seedPosz, -(leafTipsZ[i] - seedPosz), 10, "coord Z for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
if __name__ == '__main__':
    unittest.main()

