""" water movement within the root (static soil) """


import sys; 
import math
import os
import numpy as np


def toTry():
    
    # Qs = np.array([500 ])*1e-6 #np.array([200.,  340.,  450.,  560., ])*1e-6
    # #Qs =np.array([230.,240,250,  260,270,280,290,300,310,320,330])*1e-6
    # #MulimSuc =np.array([0.3,0.6,0.9,1.2,1.5])
    # MulimSuc =np.array([0.4,0.9,1.4])
    # GrRatio = np.array([3])# 3, 6, 10 np.array([3, 6, 10]) 
    # CarbonCost = np.array([1, 6, 10]) 
    # nodeD = np.array([4,8])
    # Klight = np.array([0.01,0.02,0.1])
    # maxrun = len(Qs) * len(MulimSuc) * len(nodeD) * len(GrRatio)*len(CarbonCost)
    
    Qs = np.array([300,400,500])*1e-6 
    MulimSuc =np.array([0.8,1.2,1.6,2,2.4])
    GrRatio = np.array([3])# 3, 6, 10 np.array([3, 6, 10]) 
    CarbonCost = np.array([1]) 
    nodeD = np.array([4,7])
    BerthLim = np.array([0.5])
    
    kss = np.array([ 1. ])
    kaa = np.array([ 0.2 ])
    Klight = np.array([0,0.001,0.005,0.01,0.02])
    maxrun = len(Qs) * len(MulimSuc) * len(nodeD) * len(GrRatio)*len(CarbonCost) * len(Klight)


    Qsv,MulimSucv,GrRatiov, CarbonCostv,Klightv, nodeDv = np.meshgrid(Qs,MulimSuc,GrRatio, CarbonCost,Klight, nodeD)
    Qsv=Qsv.flatten()
    MulimSucv=MulimSucv.flatten()
    nodeDv=nodeDv.flatten()
    GrRatiov = GrRatiov.flatten()
    CarbonCostv = CarbonCostv.flatten()
    Klightv = Klightv.flatten()
    dictPara= {'Qsv' : Qsv,
                    'nodeDv':nodeDv, 
               'MulimSucv':MulimSucv,
              'GrRatiov':GrRatiov,
              'CarbonCostv':CarbonCostv,
              'Klightv':Klightv}
    return dictPara

def toTryold():
    
    # Qs = np.array([500 ])*1e-6 #np.array([200.,  340.,  450.,  560., ])*1e-6
    # #Qs =np.array([230.,240,250,  260,270,280,290,300,310,320,330])*1e-6
    # #MulimSuc =np.array([0.3,0.6,0.9,1.2,1.5])
    # MulimSuc =np.array([0.4,0.9,1.4])
    # GrRatio = np.array([3])# 3, 6, 10 np.array([3, 6, 10]) 
    # CarbonCost = np.array([1, 6, 10]) 
    # nodeD = np.array([4,8])
    # Klight = np.array([0.01,0.02,0.1])
    # maxrun = len(Qs) * len(MulimSuc) * len(nodeD) * len(GrRatio)*len(CarbonCost)
    
    Qs = np.array([300,500 ,700])*1e-6 
    MulimSuc =np.array([0.4,0.75,1.1,1.4])
    GrRatio = np.array([1,3,6,9])# 3, 6, 10 np.array([3, 6, 10]) 
    CarbonCost = np.array([1/3,1, 3]) 
    nodeD = np.array([4,7])
    Klight = np.array([0.001,0.005,0.01,0.02])
    maxrun = len(Qs) * len(MulimSuc) * len(nodeD) * len(GrRatio)*len(CarbonCost) * len(Klight)


    Qsv,MulimSucv,GrRatiov, CarbonCostv,Klightv, nodeDv = np.meshgrid(Qs,MulimSuc,GrRatio, CarbonCost,Klight, nodeD)
    Qsv=Qsv.flatten()
    MulimSucv=MulimSucv.flatten()
    nodeDv=nodeDv.flatten()
    GrRatiov = GrRatiov.flatten()
    CarbonCostv = CarbonCostv.flatten()
    Klightv = Klightv.flatten()
    dictPara= {'Qsv' : Qsv,
                    'nodeDv':nodeDv, 
               'MulimSucv':MulimSucv,
              'GrRatiov':GrRatiov,
              'CarbonCostv':CarbonCostv,
              'Klightv':Klightv}
    return dictPara
def doCondition_(rinput, timeSinceDecap_, tt, simDuration, memory,nodeD_):
    orgs = rinput.plant.getOrgans(3, True)
    toKeep = np.array([org.getParameter("subType") <= 2 for org in orgs])
    orgs = np.array(orgs)[toKeep]
    budStage = np.array([org.budStage for org in orgs]) 
    sumactiv = sum(budStage[1:]>0)
    msbs = budStage[0]
    if memory >=0:
        if simDuration > 2:
            print(tt,nodeD_,"too slow", budStage, timeSinceDecap_)
            return -1
        if (msbs == 2) and (sumactiv > 0): #thrshold too low
            print(tt,nodeD_,"fail (msbs == 2) and (sumactiv > 0)", budStage, timeSinceDecap_)
            return -1
        if(timeSinceDecap_ >=2) and (sumactiv ==0): #thrshold too high
            print(tt,nodeD_,"fail (timeSinceDecap_ >=2) and (sumactiv ==0)",budStage, timeSinceDecap_)
            return -1
        if (timeSinceDecap_ >= 2) and (nodeD_ >= 6) and (budStage[nodeD_-1] < 1):#last branch did not grow
            print(tt,nodeD_,"fail (timeSinceDecap_ >= 2) and (nodeD_ >= 6) and (budStage[nodeD_-1] < 1)",budStage, timeSinceDecap_)
         #   return -1
        if (timeSinceDecap_ >= 2) and (nodeD_ < 6) and (budStage[2] < 1):#last branch did not grow
            print(tt,nodeD_,"fail (timeSinceDecap_ >= 2) and (nodeD_ < 6) and (budStage[2] < 1)",budStage, timeSinceDecap_)
          #  return -1
        #if (nodeD_ == 4) and (budStage[3] > 1):#last branch of four grew #might be too difficult
         #   print(tt,nodeD_,"fail (nodeD_ == 4) and (budStage[3] > 0)",budStage, timeSinceDecap_)
          #  return -1
        else:
            return 0
    return memory