""" water movement within the root (static soil) """


import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np

from datetime import datetime, timedelta


def toTry():
    Qs =np.array([ 120., 340., 560.,   780., 1000.])*1e-6
    PRate =np.array([1+i*2 for i in range(3)])

    Mulimaux = np.array([10*(1+i*3) for i in range(3)])
    nodeD = np.array([0,3,4,5,6,7,8])

    maxrun = len(PRate) * len(Mulimaux) * len(nodeD)*len(Qs)


    PRatev,Qsv,Mulimauxv, nodeDv = np.meshgrid(PRate,Qs,Mulimaux, nodeD)
    Qsv=Qsv.flatten()
    PRatev=PRatev.flatten()
    Mulimauxv=Mulimauxv.flatten()
    nodeDv=nodeDv.flatten()

    assert len(PRatev) == (len(PRate) * len(Mulimaux) * len(nodeD)* len(Qs)),"len(Qsv) != maxrun"
    dictPara= {'Qsv' : Qsv,
               'PRatev':PRatev,
                    'Mulimauxv': Mulimauxv,
                    'nodeDv':nodeDv}
    return dictPara

