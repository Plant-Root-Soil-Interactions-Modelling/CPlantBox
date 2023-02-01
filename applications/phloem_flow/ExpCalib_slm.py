import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np
import time
import psutil
start_time_ = time.time()

from masterI_0412 import runSim
from DBrepExp_DW import toTry
from DBrepExp_DW import doCondition_


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
#GrRatiov= params['GrRatiov']
Berthlim= params['Berthlimv']
Klightv= params['Klightv']
maxrun = len(Qsv) #tot to do
n_jobs = min((maxrun - totrun),#left to do
             maxcore - 1) #can do

#kssv= params['kssv']
#kaav= params['kaav']

print("current_process = psutil.Process()")
current_process = psutil.Process()
subproc_before = set([p.pid for p in current_process.children(recursive=True)])
print("parallelizer = Parallel(n_jobs=n_jobs)")
parallelizer = Parallel(n_jobs=n_jobs)

print("isCluster:",isCluster,"maxcore",maxcore,"to do",maxrun,"already did:", totrun, "will do",n_jobs, "directoryN",directoryN)

#-1 failur, 1 succes, 0 wait


tasks_iterator = (delayed(runSim)
                        (directoryN_ = directoryN, doVTP = 0, verbosebase = False,
             PRate_ = 6.8e-3, PRBA = 1,  PRBD=1, thresholdAux = 0, 
                         RatiothresholdAux = 1,useLength = 1,
                         Qmax_ = Qsv[i+totrun], Klight = Klightv[i+totrun],
                         thresholdSuc = 0., 
                         GrRatio = 3, CarbonCost =1,#BerthLim = 0.5,
                         maxLBud = np.array([1.]),  
                         maxLBudDormant = np.array([0.1,0.15,0.05]),
                         budGR = 0.1,L_dead_threshold=200.,
                         kss=0.2,kaa=1,
                        BerthLim = 100,#so high, as if no threshold 
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