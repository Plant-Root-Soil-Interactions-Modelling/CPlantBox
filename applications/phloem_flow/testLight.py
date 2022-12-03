import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np
import time
import psutil
start_time_ = time.time()

from masterI import runSim
from CalibP2Database import toTry


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
maxcore = os.cpu_count()

totrun = (256 - 1) * (N -1)#how much done
params = toTry()
Qsv= params['Qsv']
MulimSucv= params['MulimSucv']
nodeDv=params['nodeDv']
GrRatiov= params['GrRatiov']
CarbonCostv= params['CarbonCostv']
kssv=params['kssv']
kaav= params['kaav']
maxrun = len(Qsv) #tot to do
n_jobs = min((maxrun - totrun),#left to do
             maxcore - 1) #can do

print("current_process = psutil.Process()")
current_process = psutil.Process()
subproc_before = set([p.pid for p in current_process.children(recursive=True)])
print("parallelizer = Parallel(n_jobs=n_jobs)")
parallelizer = Parallel(n_jobs=n_jobs)

print("isCluster:",isCluster,"maxcore",maxcore,"to do",maxrun,"already did:", totrun, "will do",n_jobs, "directoryN",directoryN)

#-1 failur, 1 succes, 0 wait
def doCondition_(rinput, timeSinceDecap_, tt):
    orgs = rinput.plant.getOrgans(3, True)
    toKeep = np.array([org.getParameter("subType") <= 2 for org in orgs])
    orgs = np.array(orgs)[toKeep]
    budStage = np.array([org.budStage for org in orgs]) 
    budLength = np.array([org.getLength(False) for org in orgs])[1:] 
    maxbudL = np.where(budLength == max(budLength))[0][0] +1
    nodeD_ = len(budStage) - 1 #decapitated under node
    sumBr = sum(budStage[1:]>1)
    sumactiv = sum(budStage[1:]>0)
    msbs = budStage[0]
    #if (msbs == 2) and (sumactiv > 0): #thrshold too low
    #    print(tt,nodeD_,"fail (msbs == 2) and (sumactiv > 0)",  budStage, timeSinceDecap_)
    #    return -1
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
    if (timeSinceDecap_ >= 6.9) and (nodeD_ >= 6) and (maxbudL != nodeD_-1):#last branch did not dominate
        print(tt,nodeD_,"fail  (timeSinceDecap_ >= 6.9) and (nodeD >= 6) and (maxbudL != nodeD_-1)",maxbudL, budStage, timeSinceDecap_)
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
    return 0

tasks_iterator = (delayed(runSim)
                        (directoryN_ = directoryN, doVTP = False, verbosebase = False,
                         PRate_ = 6.8e-3, thresholdAux = 0, 
                         RatiothresholdAux = 1,
                         Qmax_ = Qsv[i+totrun], thresholdSuc = MulimSucv[i+totrun], 
                         GrRatio = GrRatiov[i+totrun], CarbonCost = CarbonCostv[i+totrun],
                         maxLBud = 1.,maxLBudDormant = 0.1,maxLBudDormant_1 = 0.3,
                         budGR = 0.1,L_dead_threshold=100.,
                         kss=kssv[i+totrun],kaa=kaav[i+totrun],
                         UseRatiothresholdAux = True,
                         nodeD = nodeDv[i+totrun], thread = i,
                         testTime=7, dtBefore = 1/24, dtAfter= 30/(60*24),
                         start_time = start_time_,
                         doPrint = True, doDict = False,
                         dt_write = 0, dtSIM_write = 20/(60*24),auxin_D=0.,
                        doCondition = doCondition_,
                        doDiffLights = True)
                    for i in range(n_jobs))

        
results = parallelizer(tasks_iterator)

def write_file_array(name, data):
    name2 = 'results'+ directoryN+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(','.join([num for num in map(str, data)])  +'\n')
        
write_file_array("successThreads", results)
print("DONE", directoryN)

print("successThreads",np.array(results))
successThreads = np.array(results)[results!= -1] + totrun
print("successThreadsbis",successThreads)
write_file_array("successThreadsbis", successThreads)
