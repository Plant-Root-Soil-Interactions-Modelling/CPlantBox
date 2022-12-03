""" water movement within the root (static soil) """


import sys; 
import math
import os
import numpy as np
import pandas as pd

def toTry_old():
    
    dfP2 = pd.read_csv('resultsConsolidated/forP2.csv')


    # dictPara= {'Qsv' : dfP2['Qsv'],
    #                 'nodeDv':dfP2['nodeDv'], 
    #            'MulimSucv':dfP2['MulimSucv'],
    #                 'kaav':dfP2['kaav'], 
    #            'kssv':dfP2['kssv'],
    #           'GrRatiov':dfP2['GrRatiov'],
    #           'CarbonCostv':dfP2['CarbonCostv']}
    # df = pd.DataFrame(data=dictPara)
    if sum(dfP2.duplicated()) > 0 :
        print("CalibP2Database::toTry(): sum(df.duplicated()) > 0")
        raise Exception
    return dfP2 #dictPara



def toTry():
    
    Qs =np.array([200,400.,700])#,  500.,700 ])*1e-6
    MulimSuc =np.array([0.4,0.9,1.4])#np.array([0.4,0.9,1.4])
    GrRatio = np.array([3, 6]) 
    CarbonCost = np.array([1]) 
    nodeD = np.array([8]) #3,4,5,6,7,
    kss = np.array([0.01, 0.4,0.6])
    kaa = np.array([ 1. , 3.,  5.,  10.])
    Klight = np.array([0.001,0.005,0.01,0.02])
    maxrun = len(Qs) * len(MulimSuc) * len(nodeD) * len(GrRatio)*len(CarbonCost)*len(kss)*len(kaa)*len(Klight)


    Qsv,MulimSucv,GrRatiov, CarbonCostv, kssv,kaav,Klightv, nodeDv = np.meshgrid(Qs,MulimSuc,GrRatio, CarbonCost, kss,kaa,Klight,nodeD)
    Qsv=Qsv.flatten()
    MulimSucv=MulimSucv.flatten()
    nodeDv=nodeDv.flatten()
    GrRatiov = GrRatiov.flatten()
    CarbonCostv = CarbonCostv.flatten()
    kssv = kssv.flatten()
    kaav = kaav.flatten()
    Klightv = Klightv.flatten()
    dictPara= {'Qsv' : Qsv,
               'nodeDv':nodeDv, 
               'MulimSucv':MulimSucv,
               'GrRatiov':GrRatiov,
               'CarbonCostv':CarbonCostv,
               'kssv':kssv,
               'kaav':kaav,
              'Klightv':Klightv}
    return dictPara


def toTry_old2():
    
    Qs =np.array([400.,700 ])*1e-6
    MulimSuc =np.array([0.4,0.9,1.4])
    GrRatio = np.array([3]) 
    CarbonCost = np.array([1]) 
    nodeD = np.array([4,6,8])
    kss = np.array([0.01, 0.4,0.6])
    kaa = np.array([ 1. , 3.,  5.,10.])
    Klight = np.array([0.001,0.005,0.01,0.02])
    maxrun = len(Qs) * len(MulimSuc) * len(nodeD) * len(GrRatio)*len(CarbonCost)*len(kss)*len(kaa)*len(Klight)


    Qsv,MulimSucv,GrRatiov, CarbonCostv, kssv,kaav,Klightv, nodeDv = np.meshgrid(Qs,MulimSuc,GrRatio, CarbonCost, kss,kaa,Klight,nodeD)
    Qsv=Qsv.flatten()
    MulimSucv=MulimSucv.flatten()
    nodeDv=nodeDv.flatten()
    GrRatiov = GrRatiov.flatten()
    CarbonCostv = CarbonCostv.flatten()
    kssv = kssv.flatten()
    kaav = kaav.flatten()
    Klightv = Klightv.flatten()
    dictPara= {'Qsv' : Qsv,
               'nodeDv':nodeDv, 
               'MulimSucv':MulimSucv,
               'GrRatiov':GrRatiov,
               'CarbonCostv':CarbonCostv,
               'kssv':kssv,
               'kaav':kaav,
              'Klightv':Klightv}
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
        if((sumBr >3) and (nodeD_ != 5) and (nodeD_ != 6)): #thrshold too low
            print(tt,nodeD_,"(fail sumBr >3) and (nodeD_ != 5) and (nodeD_ != 6)",budStage, timeSinceDecap_)
            return -1
        if (timeSinceDecap_ >= 6.9) and (nodeD_ >= 6) and (budStage[nodeD_-1] < 2):#last branch did not grow
            print(tt,nodeD_,"fail  (timeSinceDecap_ >= 6.9) and (nodeD >= 6) and (budStage[nodeD-2] < 2)",budStage, timeSinceDecap_)
            return -1
        #if (timeSinceDecap_ >= 6.9) and (nodeD_ >= 6) and (maxbudL != nodeD_-1):#last branch did not dominate
         #   print(tt,nodeD_,"fail  (timeSinceDecap_ >= 6.9) and (nodeD >= 6) and (maxbudL != nodeD_-1)",maxbudL, budStage, timeSinceDecap_)
          #  return -1
        if (nodeD_ == 8) and (budStage[1:5] > 1):#branch 2 did not grow 
            print(tt,nodeD_,"(nodeD_ == 8) and (budStage[1:5] > 1)",budStage, timeSinceDecap_)
            return -1
        if (timeSinceDecap_ >= 6.9) and (nodeD_ < 6) and (budStage[2] < 2):#branch 2 did not grow 
            print(tt,nodeD_,"fail (timeSinceDecap_ >= 6.9) and (nodeD < 6) and (budStage[2] < 2)",budStage, timeSinceDecap_)
            return -1
        if (timeSinceDecap_ >= 6.9) and (nodeD_ < 6) and (maxbudL != 2):#branch 2 did not dominate 
            print(tt,nodeD_,"fail (timeSinceDecap_ >= 6.9) and (nodeD < 6) and (maxbudL != 2)",budStage, timeSinceDecap_)
            return -1
        if (timeSinceDecap_ >= 6.9) and (nodeD_ == 4) and (budStage[3] ==2):#last branch of four grew #might be too difficult
            print(tt,nodeD_,"fail (nodeD == 4) and (budStage[3] ==2)",budStage, timeSinceDecap_)
            return -1
    return memory
