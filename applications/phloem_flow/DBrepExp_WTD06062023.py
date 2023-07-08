""" water movement within the root (static soil) """


import sys; 
import math
import os
import numpy as np
import pandas as pd

def toTry():
    
    Qs =1e-6*np.array([500.*0.64])#,300.,400.,500.])#,  500.,700 ])*1e-6
    MulimSuc = np.array([0.07,0.1,0.2,0.3,0.4,0.6])#np.array([0.4,0.9,1.25,1.5,2.])#np.array([0.4,0.9,1.4])
    GrRatio = np.array([1]) 
    CarbonCost = np.array([1]) 
    nodeD = np.array([4,5,6,7,8,9])
    kss = np.array([0.02])#0.01, 0.4,0.6])
    kaa = np.array([ 1. ])#, 3.,  5.,  10.])
    Klight = np.array([0])#0.001,0.005,0.01,0.02])
    Berthlim = np.array([ 1,2,3,4,6,8 ])
    betaMeso = np.array([0.5, 1, 1.5, 2, 2.5]) 


    Qsv, MulimSucv, GrRatiov, CarbonCostv, kssv, kaav, Klightv, Berthlimv, nodeDv,betaMesov \
            = np.meshgrid(Qs, MulimSuc, GrRatio, CarbonCost, kss, kaa, Klight, Berthlim, nodeD,betaMeso)
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
    betaMesov = betaMesov.flatten()
    dictPara= {'Qsv' : Qsv,
               'nodeDv':nodeDv, 
               'MulimSucv':MulimSucv,
               'GrRatiov':GrRatiov,
               'CarbonCostv':CarbonCostv,
               'kssv':kssv,
               'kaav':kaav,
              'Klightv':Klightv,
              'Berthlimv':Berthlimv,
              'betaMesov':betaMesov}
    return dictPara

        

#-1 failur, 1 succes, 0 wait
def doCondition_(rinput, timeSinceDecap_, tt, simDuration, memory,nodeD_,allInputs):
    Mstem = rinput.plant.getOrgans(3, False)[0]
    orgs = np.array([Mstem.getChild(nkdb) for nkdb in range(Mstem.getNumberOfChildren()) if Mstem.getChild(nkdb).organType() == 3])
    budStage = np.array([org.budStage for org in orgs]) 
    budLength = np.array([org.getLength(False) for org in orgs])
    maxbudL = np.where(budLength == max(budLength))[0][0] +1
    sumBr = sum(budStage>1)
    sumactiv = sum(budStage>0)
    msbs = Mstem.budStage
    if memory >=0:
        if simDuration > 2/24:#too slow
            print(tt,nodeD_,"fail simDuration > 2/24", budStage, timeSinceDecap_)
            return -1
        if (msbs == 2) and (sumactiv > 0): #thrshold too low
            print(tt,nodeD_,"fail (msbs == 2) and (sumactiv > 0)", budStage, timeSinceDecap_,allInputs["thresholdSuc"] ,
                  allInputs["BerthLim"],allInputs["GrRatioLats"])
            return -1
        if(timeSinceDecap_ >=1) and (sumactiv ==0) and (msbs != 2): #thrshold too high
            print(tt,nodeD_,"fail (timeSinceDecap_ >=1) and (sumactiv ==0) and (msbs != 2)",budStage,np.round( timeSinceDecap_*24*10)/10,
                 allInputs["thresholdSuc"] ,
                  allInputs["BerthLim"],allInputs["GrRatioLats"])
            return -1
        if(timeSinceDecap_ >= 8.9) and (sumBr ==0) and (msbs != 2): #thrshold too high
            print(tt,nodeD_,"fail (timeSinceDecap_ >= 8.9) and (sumBr ==0) and (msbs != 2)",budStage, timeSinceDecap_,allInputs["thresholdSuc"] ,
                  allInputs["BerthLim"],allInputs["GrRatioLats"])
            return -1
        if(budStage[1] < 2) and (sumBr >0) : #thrshold wrong
            print(tt,nodeD_,"fail (budStage[1] < 2) and (sumBr >0) ",budStage, timeSinceDecap_,allInputs["thresholdSuc"] ,
                  allInputs["BerthLim"],allInputs["GrRatioLats"])
            return -1
        
    return memory

#evaluate quality of output
def successPoints_(rinput, timeSinceDecap_, tt, simDuration, memory,nodeD_,allInputs):
    Mstem = rinput.plant.getOrgans(3, False)[0]
    orgs = np.array([Mstem.getChild(nkdb) for nkdb in range(Mstem.getNumberOfChildren()) if Mstem.getChild(nkdb).organType() == 3])
    budStage = np.array([org.budStage for org in orgs]) 
    budLength = np.array([org.getLength(False) for org in orgs])
    maxbudL = np.where(budLength == max(budLength))[0][0] +1
    secMaxBudL = np.where(budLength == sorted(budLength)[-2] )[0][0] +1
    sumBr = sum(budStage>1)
    sumactiv = sum(budStage>0)
    msbs = Mstem.budStage
    score = 0
    if memory < 0:
        score = memory
    else:        
        
        if nodeD_ >= 7: #later cut:
            if maxbudL ==( nodeD_ -1):
                score += 1
            if secMaxBudL ==( nodeD_ -1):
                score += 1
        if maxbudL == 2: #all cuts
            score += 1
        if nodeD_ == 4: 
            if sumBr == 1 :
                score += 1
        if nodeD_ > 4:
            if sumactiv == (nodeD_ -1) :
                score += 1
    print(tt,nodeD_,"successPoints",score, budStage,budLength,allInputs["thresholdSuc"] , allInputs["BerthLim"],allInputs["GrRatioLats"])    
    return score