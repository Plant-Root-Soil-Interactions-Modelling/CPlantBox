""" water movement within the root (static soil) """


import sys; 
import math
import os
import numpy as np


def toTry():
    
    Qs =np.array([ 120.,  230.,  340.,  450.,  560., ])*1e-6
    MulimSuc =np.array([0.1,0.2,0.3])
    GrRatio = np.array([3]) 
    nodeD = np.array([3,4,5,6,7,8])#0
    kss_ = np.array([0.02,0.2,2.0])
    kaa_ = np.array([0.1,1.0,10.])
    maxrun = len(Qs) * len(MulimSuc) * len(nodeD) * len(GrRatio)* len(kss_) *len(kaa_)


    Qsv,MulimSucv, nodeDv,GrRatiov, kss_v, kaa_v = np.meshgrid(Qs,MulimSuc, nodeD,GrRatio, kss_, kaa_)
    Qsv=Qsv.flatten()
    MulimSucv=MulimSucv.flatten()
    nodeDv=nodeDv.flatten()
    GrRatiov = GrRatiov.flatten()
    kss_v = kss_v.flatten()
    kaa_v = kaa_v.flatten()
    dictPara= {'Qsv' : Qsv,
                    'nodeDv':nodeDv, 
               'MulimSucv':MulimSucv,
              'GrRatiov':GrRatiov,
              'kss_v': kss_v,
              'kaa_v': kaa_v}
    return dictPara

