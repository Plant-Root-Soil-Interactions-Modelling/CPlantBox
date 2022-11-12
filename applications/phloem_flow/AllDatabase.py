""" water movement within the root (static soil) """


import sys; 
import math
import os
import numpy as np


def toTry():
    
    Qs =np.array([ 120.,  230.,  340.,  450.,  560., ])*1e-6
    MulimSuc =np.linspace(0.4,2.0,6)#np.array([0.1*i for i in range(4,20)])
    Ratiothrshold = np.linspace(0.1,1,6) #np.array([1/10*i for i in range(10)])
    nodeD = np.array([0,3,4,5,6,7,8])
    maxrun = len(Qs) * len(MulimSuc) * len(nodeD) * len(Ratiothrshold)


    Qsv,MulimSucv, nodeDv,Ratiothrsholdv = np.meshgrid(Qs,MulimSuc, nodeD,Ratiothrshold)
    Qsv=Qsv.flatten()
    MulimSucv=MulimSucv.flatten()
    nodeDv=nodeDv.flatten()
    Ratiothrsholdv = Ratiothrsholdv.flatten()
    dictPara= {'Qsv' : Qsv,
                    'nodeDv':nodeDv, 
               'MulimSucv':MulimSucv,
              'Ratiothrshold':Ratiothrsholdv}
    return dictPara

