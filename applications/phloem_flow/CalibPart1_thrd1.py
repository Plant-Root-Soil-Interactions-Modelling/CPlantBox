import sys; 
import math
import os
import numpy as np
import time
import psutil
start_time_ = time.time()

from masterI import runSim
from CalibP1Database import toTry


#def AllAuxCmasterFunc(N):
    
N = int(sys.argv[1])
directoryN = "/"+os.path.basename(__file__)[:-3]+str(N)+"/"
print("N",N)
main_dir=os.environ['PWD']#dir of the file
results_dir = main_dir +"/results"+directoryN
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    import shutil
    shutil.rmtree(results_dir)
    os.makedirs(results_dir)

isCluster = (os.environ['HOME'] == '/home/m.giraud')
maxcore =  os.cpu_count()

totrun = (256 - 1) * (N -1)#how much done
params = toTry()
Qsv=params['Qsv']
MulimSucv=params['MulimSucv']
nodeDv=params['nodeDv']
GrRatiov= params['GrRatiov']
CarbonCostv= params['CarbonCostv']
Klightv= params['Klightv']
maxrun = len(Qsv) #tot to do
n_jobs = min((maxrun - totrun),#left to do
             maxcore - 1) #can do


print("isCluster:",isCluster,"maxcore",maxcore,"to do",maxrun,"already did:", totrun, "will do",n_jobs, "directoryN",directoryN)

#-1 failur, 1 succes, 0 wait
def doCondition_(rinput, timeSinceDecap_, tt, simDuration):
    orgs = rinput.plant.getOrgans(3, True)
    toKeep = np.array([org.getParameter("subType") <= 2 for org in orgs])
    orgs = np.array(orgs)[toKeep]
    budStage = np.array([org.budStage for org in orgs]) 
    sumactiv = sum(budStage[1:]>0)
    msbs = budStage[0]
    nodeD_ = len(budStage) - 1 #decapitated under node
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
        print(tt,nodeD_,"fail  (timeSinceDecap_ >= 6.9) and (nodeD >= 6) and (budStage[nodeD-2] < 2)",budStage, timeSinceDecap_)
        return -1
    if (timeSinceDecap_ >= 2) and (nodeD_ < 6) and (budStage[2] < 1):#last branch did not grow
        print(tt,nodeD_,"fail  (timeSinceDecap_ >= 6.9) and (nodeD >= 6) and (budStage[nodeD-2] < 2)",budStage, timeSinceDecap_)
        return -1
    if (nodeD_ == 4) and (budStage[3] > 0):#last branch of four grew #might be too difficult
        print(tt,nodeD_,"fail (nodeD_ == 4) and (budStage[3] > 0)",budStage, timeSinceDecap_)
        return -1
    else:
        return 0

runSim(directoryN_ = directoryN, doVTP = 1,verbosebase = False, 
        PRate_ = 6.8e-3,thresholdAux = 0, RatiothresholdAux =0.461,
        UseRatiothresholdAux = True,
        Qmax_ = 700e-6, thresholdSuc = 1.4,  
        GrRatio = 9, CarbonCost = 1, maxLBud = 1,
        maxLBudDormant = 0.1,maxLBudDormant_1 = 0.15,
        budGR = 0.1,L_dead_threshold=100,
        nodeD = 5, thread = 0,
        testTime=4, dtBefore = 1/24, dtAfter= 30/(24*60),
        start_time = start_time_,
        dt_write = 0, dtSIM_write = 1/(24*60), 
        doPrint = True, doDict = False,auxin_D=0.,doDiffLights=True)

        

