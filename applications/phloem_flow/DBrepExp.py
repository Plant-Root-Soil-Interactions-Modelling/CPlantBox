""" water movement within the root (static soil) """


import sys; 
import math
import os
import numpy as np
import pandas as pd

def toTry():
    
    Qs =np.array([200.,300.,400.,500.])#,  500.,700 ])*1e-6
    MulimSuc =np.array([0.4,0.9,1.25,1.5,2.])#np.array([0.4,0.9,1.4])
    GrRatio = np.array([3]) 
    CarbonCost = np.array([1]) 
    nodeD = np.array([0,7]) #3,4,5,6,7,
    kss = np.array([0.02])#0.01, 0.4,0.6])
    kaa = np.array([ 1. ])#, 3.,  5.,  10.])
    Klight = np.array([0])#0.001,0.005,0.01,0.02])
    Berthlim = np.array([ 1,2,3,4,6,8 ])


    Qsv,MulimSucv,GrRatiov, CarbonCostv, kssv,kaav,Klightv,Berthlimv, nodeDv = np.meshgrid(Qs,MulimSuc,GrRatio, CarbonCost, kss,kaa,Klight,Berthlim,nodeD)
    maxrun = len(Qsv)
    Qsv=Qsv.flatten()
    MulimSucv=MulimSucv.flatten()
    nodeDv=nodeDv.flatten()
    GrRatiov = GrRatiov.flatten()
    CarbonCostv = CarbonCostv.flatten()
    kssv = kssv.flatten()
    kaav = kaav.flatten()
    Klightv = Klightv.flatten()
    Berthlimv = Berthlimv.flatten()
    dictPara= {'Qsv' : Qsv,
               'nodeDv':nodeDv, 
               'MulimSucv':MulimSucv,
               'GrRatiov':GrRatiov,
               'CarbonCostv':CarbonCostv,
               'kssv':kssv,
               'kaav':kaav,
              'Klightv':Klightv,
              'Berthlimv':Berthlimv}
    return dictPara


#-1 failur, 1 succes, 0 wait
def doCondition_(rinput, timeSinceDecap_, tt, simDuration, memory,nodeD_):
    orgs = rinput.plant.getOrgans(3, True)
    toKeep = np.array([org.getParameter("subType") <= 2 for org in orgs])
    orgs = np.array(orgs)[toKeep]
    budStage = np.array([org.budStage for org in orgs]) 
    budLength = np.array([org.getLength(False) for org in orgs])[1:] 
    maxbudL = np.where(budLength == max(budLength))[0][0] +1
    sumBr = sum(budStage[1:]>1)
    sumactiv = sum(budStage[1:]>0)
    msbs = budStage[0]
    #if (msbs == 2) and (sumactiv > 0): #thrshold too low
    #    print(tt,nodeD_,"fail (msbs == 2) and (sumactiv > 0)",  budStage, timeSinceDecap_)
    #    return -1
    if memory >=0:
        if simDuration > 2:
            print(tt,nodeD_,"too slow", budStage, timeSinceDecap_)
            return -1
        if (msbs == 2) and (sumBr > 0): #thrshold too low
            print(tt,nodeD_,"fail (msbs == 2) and (sumBr > 0)", budStage, timeSinceDecap_)
            return -1
        if(timeSinceDecap_ >=2) and (sumactiv ==0): #thrshold too high
            print(tt,nodeD_,"fail (timeSinceDecap_ >=2) and (sumactiv ==0)",budStage, timeSinceDecap_)
            return -1
        if(timeSinceDecap_ >= 6.9) and (sumBr ==0): #thrshold too high
            print(tt,nodeD_,"fail (timeSinceDecap_ >= 6.9) and (sumBr ==0)",budStage, timeSinceDecap_)
            return -1
        
    return memory
