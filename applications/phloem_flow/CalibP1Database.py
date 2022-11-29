""" water movement within the root (static soil) """


import sys; 
import math
import os
import numpy as np


def toTry():
    
    Qs =np.array([ 120.,  230.,  340.,  450.,  560., ])*1e-6
    #Qs =np.array([230.,240,250,  260,270,280,290,300,310,320,330])*1e-6
    #MulimSuc =np.array([0.3,0.6,0.9,1.2,1.5])
    MulimSuc =np.array([0.4,0.9,1.4])
    GrRatio = np.array([3, 6, 10]) 
    CarbonCost = np.array([1, 6, 10]) 
    nodeD = np.array([3,4,5,6,7,8])
    maxrun = len(Qs) * len(MulimSuc) * len(nodeD) * len(GrRatio)*len(CarbonCost)


    Qsv,MulimSucv,GrRatiov, CarbonCostv, nodeDv = np.meshgrid(Qs,MulimSuc,GrRatio, CarbonCost, nodeD)
    Qsv=Qsv.flatten()
    MulimSucv=MulimSucv.flatten()
    nodeDv=nodeDv.flatten()
    GrRatiov = GrRatiov.flatten()
    CarbonCostv = CarbonCostv.flatten()
    dictPara= {'Qsv' : Qsv,
                    'nodeDv':nodeDv, 
               'MulimSucv':MulimSucv,
              'GrRatiov':GrRatiov,
              'CarbonCostv':CarbonCostv}
    return dictPara

