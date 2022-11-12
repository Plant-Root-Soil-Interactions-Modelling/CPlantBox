""" water movement within the root (static soil) """


import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np

from datetime import datetime, timedelta

def toTry():
    Qs =np.array([ 120.,  230.,  340.,  450.,  560., ])*1e-6# 670.,  780.,  890.,1000.
    MulimSuc = np.array([0.1*i for i in range(4,20)])
    nodeD = np.array([0,3,4,5,6,7,8])
    maxrun = len(Qs) * len(MulimSuc) * len(nodeD)

    Qsv,MulimSucv, nodeDv = np.meshgrid(Qs,MulimSuc, nodeD)
    Qsv=Qsv.flatten()
    MulimSucv=MulimSucv.flatten()
    nodeDv=nodeDv.flatten()
    assert len(Qsv) == maxrun,"len(Qsv) != maxrun"
    dictPara= {'Qsv' : Qsv,
                    'nodeDv':nodeDv, 'MulimSucv':MulimSucv}
    return dictPara
