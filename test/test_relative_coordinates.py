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

import unittest
import sys; sys.path.append(".."); sys.path.append("../src/python_modules")
import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import vtk_plot as vp


class TestRelCoord(unittest.TestCase):

    def test_coord(self):
        pl = pb.MappedPlant()  
        path = "../modelparameter/plant/"  
        name = "test_relcoord"
        pl.readParameters(path + name + ".xml") 
        #stochastic = False => thus rand() always give 0.5 for Tropism. => result only change according to sigma        
        pl.initialize(stochastic = False)
        dt = 1
        steps = 100
        for step in range(steps):
            pl.simulate(1, False)
        pl.write("test_relcoord.vtp")
        params = pl.organParam
        seedPosx = params[1][0].seedPos.x
        seedPosy = params[1][0].seedPos.y
        seedPosz = params[1][0].seedPos.z
        #print(params[1][0].seedPos)
        roots = pl.getOrgans(2)
        leaves = pl.getOrgans(4)
        #stems = pl.getOrgans(3)
        
        rootSubtypes = [ o.param().subType for o in roots]
        #stemSubtypes = [ o.param().subType for o in stems]
        
        #roots1 = np.array(roots)[np.where([st == 1 for st in rootSubtypes])[0]]
        #roots2 = np.array(roots)[np.where([st == 2 for st in rootSubtypes])[0]]
        roots3 = np.array(roots)[np.where([st == 3 for st in rootSubtypes])[0]]
        
        #stems1 = np.array(roots)[np.where([st == 1 for st in stemSubtypes])[0]]
        #stems3 = np.array(roots)[np.where([st == 3 for st in stemSubtypes])[0]]
        
        
        rootTipsX = [np.array(r.getNode(r.getNumberOfNodes()-1))[0] for r in roots3]
        rootTipsY = [np.array(r.getNode(r.getNumberOfNodes()-1))[1] for r in roots3]
        rootTipsZ = [np.array(r.getNode(r.getNumberOfNodes()-1))[2] for r in roots3]
        
        leafTipsX = [np.array(r.getNode(r.getNumberOfNodes()-1))[0] for r in leaves]
        leafTipsY = [np.array(r.getNode(r.getNumberOfNodes()-1))[1] for r in leaves]
        leafTipsZ = [np.array(r.getNode(r.getNumberOfNodes()-1))[2] for r in leaves]
        
        for i in range(0, len(rootTipsX)):
            #print(rootTipsX[i], leafTipsX[i])
            self.assertAlmostEqual(rootTipsX[i]-seedPosx, -(leafTipsX[i]-seedPosx), 10, "tip of 3rd lat root and leaf n°"+ str(i) +" not symetric")
        for i in range(0, len(rootTipsY)):
            #print(rootTipsY[i]*10**16, leafTipsY[i]*10**16)
            self.assertAlmostEqual(rootTipsY[i]*10**16, leafTipsY[i]*10**16, 10,  "tip of 3rd lat root and leaf n°"+ str(i) +" not symetric")
        for i in range(0, len(rootTipsZ)):
            #print(rootTipsZ[i]-seedPosz, leafTipsZ[i]-seedPosz, seedPosz)
            self.assertAlmostEqual(rootTipsZ[i]-seedPosz, -(leafTipsZ[i]-seedPosz), 10,  "tip of 3rd lat root and leaf n°"+ str(i) +" not symetric")
        

if __name__ == '__main__':
    unittest.main()


