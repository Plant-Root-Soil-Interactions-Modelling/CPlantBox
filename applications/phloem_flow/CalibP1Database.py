""" water movement within the root (static soil) """


import sys; 
import math
import os
import numpy as np


def toTry():
    
    #Qs =np.array([ 120.,  230.,  340.,  450.,  560., ])*1e-6
    Qs =np.array([230.,240,250,  260,270,280,290,300,310,320,330])*1e-6
    #MulimSuc =np.array([0.3,0.6,0.9,1.2,1.5])
    MulimSuc =np.array([0.4,0.5,0.6,0.9,1.2,1.3,1.4])
    GrRatio = np.array([3]) 
    nodeD = np.array([3,4,5,6,7,8])
    maxrun = len(Qs) * len(MulimSuc) * len(nodeD) * len(GrRatio)


    Qsv,MulimSucv, nodeDv,GrRatiov = np.meshgrid(Qs,MulimSuc, nodeD,GrRatio)
    Qsv=Qsv.flatten()
    MulimSucv=MulimSucv.flatten()
    nodeDv=nodeDv.flatten()
    GrRatiov = GrRatiov.flatten()
    dictPara= {'Qsv' : Qsv,
                    'nodeDv':nodeDv, 
               'MulimSucv':MulimSucv,
              'GrRatiov':GrRatiov}
    return dictPara

