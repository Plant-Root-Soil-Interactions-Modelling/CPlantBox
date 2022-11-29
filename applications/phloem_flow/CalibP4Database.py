""" water movement within the root (static soil) """


import sys; 
import math
import os
import numpy as np


def toTry():
    
    #Qs =np.array([ 120.,  230.,  340.,  450.,  560., ])*1e-6
    Qs =np.array([230.])*1e-6
    #MulimSuc =np.array([0.3,0.6,0.9,1.2,1.5])
    MulimSuc =np.array([0.6,0.9, 1.2])
    GrRatio = np.array([3]) 
    nodeD = np.array([3,4,5,6,7,8])
    kss = np.array([0.01, 0.100, 0.2, 0.3, 0.4])
    kaa = np.array([ 1. , 3.])
    maxrun = len(Qs) * len(MulimSuc) * len(nodeD) * len(GrRatio)*len(kaa)*len(kss)


    Qsv,MulimSucv,kssv,kaav, nodeDv,GrRatiov = np.meshgrid(Qs,MulimSuc,kss,kaa,nodeD,GrRatio)
    Qsv=Qsv.flatten()
    MulimSucv=MulimSucv.flatten()
    kaav=kaav.flatten()
    kssv = kssv.flatten()
    nodeDv=nodeDv.flatten()
    GrRatiov = GrRatiov.flatten()
    dictPara= {'Qsv' : Qsv,
                    'nodeDv':nodeDv, 
               'MulimSucv':MulimSucv,
                    'kaav':kaav, 
               'kssv':kssv,
              'GrRatiov':GrRatiov}
    return dictPara

