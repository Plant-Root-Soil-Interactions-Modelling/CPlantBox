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
import plantbox.visualisation.vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt


def angle_between_segments(points, i):
    p = np.asarray(points, dtype=float)

    v1 = p[i]   - p[i-1]
    v2 = p[i+1] - p[i]

    # Normalize vectors
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)

    # Dot product and angle
    cos_theta = np.clip(np.dot(v1, v2), -1.0, 1.0)
    angle_rad = np.arccos(cos_theta)

    return np.degrees(angle_rad)

def segment_lengths(points):
    pts = np.asarray(points, dtype=float)
    diffs = pts[1:] - pts[:-1]    
    return np.linalg.norm(diffs, axis=1)
    
class TestRelCoord(unittest.TestCase):

    def plant_example(self, steps):
        pl = pb.MappedPlant()
        path = "../modelparameter/structural/plant/"
        name = "test_relcoord"
        pl.readParameters(path + name + ".xml")
        pl.initialize()
        # stochastic = False => thus rand() always give 0.5 for Tropism. => result only change according to sigma
        pl.setStochastic(False)
        
        dt = 100/steps
        
        for step in range(steps):
            pl.simulate(dt, False)
        
        #pl.write("test_relcoord_" + str(steps) +"steps.vtp")
        #vp.plot_plant(pl, "subType")
        params = pl.getOrganRandomParameter(1)[0]
        seedPosx = params.seedPos.x
        seedPosy = params.seedPos.y
        seedPosz = params.seedPos.z
        
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
            self.assertAlmostEqual(rootTipsX[i] - seedPosx, -(stemTipsX[i] - seedPosx), 5, "coord X for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsY)):
            # print(rootTipsY[i]*10**16, leafTipsY[i]*10**16)
            self.assertAlmostEqual(rootTipsY[i] * 10 ** 16, stemTipsY[i] * 10 ** 16, 5, "coord Y for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsZ)):
            # print(rootTipsZ[i]-seedPosz, leafTipsZ[i]-seedPosz, seedPosz)
            self.assertAlmostEqual(rootTipsZ[i] - seedPosz, -(stemTipsZ[i] - seedPosz), 5, "coord Z for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        
        
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
            self.assertAlmostEqual(rootTipsZ[i] - seedPosz, -(stemTipsZ[i] - seedPosz), 5, "coord Z for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        
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
                                               
            self.assertAlmostEqual(rootTipsX[i] - seedPosx, -(leafTipsX[i] - seedPosx), 5, "coord X for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsY)):
            # print(rootTipsY[i]*10**16, leafTipsY[i]*10**16)
            self.assertAlmostEqual(rootTipsY[i] * 10 ** 16, leafTipsY[i] * 10 ** 16, 5, "coord Y for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")
        for i in range(0, len(rootTipsZ)):
            # print(rootTipsZ[i]-seedPosz, leafTipsZ[i]-seedPosz, seedPosz)
            self.assertAlmostEqual(rootTipsZ[i] - seedPosz, -(leafTipsZ[i] - seedPosz), 5, "coord Z for tip of 3rd lat root and leaf n°" + str(i) + " not symetric")

    def test_coord_100steps(self):
        self.plant_example(100)
        
    def test_coord_1step(self):
        self.plant_example(1)
        
if __name__ == '__main__':
    unittest.main()

