import sys; 
from joblib import Parallel, delayed
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
maxcore = 8
if isCluster:
    maxcore = 256

totrun = (maxcore - 1) * (N -1)#how much done
params = toTry()
Qsv=params['Qsv']
MulimSucv=params['MulimSucv']
nodeDv=params['nodeDv']
GrRatiov= params['GrRatiov']
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
    sumactiv = sum(budStage[1:]>0)
    msbs = budStage[0]
    if (msbs == 2) and (sumactiv > 0): #thrshold too low
        print(tt,"fail", budStage, timeSinceDecap_)
        return -1
    if(timeSinceDecap_ >= 6.5) and (sumactiv ==0): #thrshold too high
        print(tt,"fail",budStage, timeSinceDecap_)
        return -1
    if(timeSinceDecap_ >= 0) and (sumactiv > 0): #thrshold too high
        print(tt,"success",budStage, timeSinceDecap_)
        return 1
    else:
        return 0

tasks_iterator = (delayed(runSim)
                        (directoryN_ = directoryN, doVTP = False, verbosebase = False,
                         PRate_ = 6.8e-3, thresholdAux = 0, 
                         RatiothresholdAux = 1,
       Qmax_ = Qsv[i+totrun], thresholdSuc = MulimSucv[i+totrun], 
                         GrRatio = GrRatiov[i+totrun], 
                         maxLBud = 1., budGR = 0.1,L_dead_threshold=100.,
                         #kss=kss_v[i+totrun],kaa=kaa_v[i+totrun],
       UseRatiothresholdAux = True,
                         nodeD = nodeDv[i+totrun], thread = i,
                         testTime=7, dtBefore = 1/24, dtAfter= 30/(60*24),
                        start_time = start_time_,
                         doPrint = True, doDict = False,
                         dt_write = 0, dtSIM_write = 30/(60*24),auxin_D=0.,
                        doCondition = doCondition_)
                    for i in range(n_jobs))

        
results = parallelizer(tasks_iterator)

def write_file_array(name, data):
    name2 = 'results'+ directoryN+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(','.join([num for num in map(str, data)])  +'\n')
        
write_file_array("successThreads", results)
print("DONE", directoryN)
subproc_after = set([p.pid for p in current_process.children(recursive=True)]) - subproc_before

# for subproc in subproc_after:
#     kid_process = psutil.Process(subproc)
#     print('Killing n1 process with pid {}'.format(subproc))
#     kid_process.terminate()
print("successThreads",results)
    
#     import signal
# import platform
# # get the current PID for safe terminate server if needed:
# PID = os.getpid()
# if platform.system() is not 'Windows':
#     PGID = os.getpgid(PID)
# if platform.system() is not 'Windows':
#     os.killpg(PGID, signal.SIGKILL)
# else:
#     os.kill(PID, signal.SIGTERM)