import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np
import time
import psutil
start_time_ = time.time()

from AllDatabase import toTry


#def AllAuxCmasterFunc(N):
N = ""
if len(sys.argv)>1:    
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

totrun = 0#(maxcore - 1) * (N -1)
params = toTry()
Qsv=params['Qsv']
MulimSucv=params['MulimSucv']
nodeDv=params['nodeDv']
GrRatiov= params['GrRatiov']
kss_v= params['kss_v']
kaa_v= params['kaa_v']
maxrun = len(Qsv)
n_jobs =2# min((maxrun - totrun),#left to do
          #   maxcore - 1) #can do
print(os.getpid())
print("psutil.Process().pid",psutil.Process())
current_process = psutil.Process()
subproc_before = set([p.pid for p in current_process.children(recursive=True)])
parallelizer = Parallel(n_jobs=n_jobs)

print("isCluster:",isCluster,"maxcore",maxcore,"to do",maxrun,"already did:", totrun, "will do",n_jobs, "directoryN",directoryN)

def runSim(thread):
    print("high, I am n#",thread)

tasks_iterator = (delayed(runSim)
                        (thread = i)
                    for i in range(n_jobs))
print("test 1",set([p.pid for p in current_process.children(recursive=False)]))
print("test2", [i for i in range(n_jobs)])
parallelizer(tasks_iterator)
print("in masterfile")
#raise Exception #overwise hangs?
#os.abort()
subproc_after = set([p.pid for p in current_process.children(recursive=True)]) - subproc_before
print("the processes",subproc_before,subproc_after)
for subproc in subproc_after:
    kid_process = psutil.Process(subproc)
    print('Killing n1 process with pid {}'.format(subproc), "memory_info",kid_process.memory_info())
    #kid_process.terminate()
#for subproc in subproc_after:
 #   kid_process = psutil.Process(subproc)
  #  print('Killing n2 process with pid {}'.format(subproc), "memory_info",kid_process.memory_info())
#slurmstepd: error: *** JOB 6441 ON node03 CANCELLED AT 2022-11-19T10:00:54 ***
#/home/m.giraud/.local/lib/python3.8/site-packages/joblib/externals/loky/backend/resource_tracker.py:318: UserWarning: resource_tracker: There appear to be 8 leaked semlock objects to clean up at shutdown
 # warnings.warn('resource_tracker: There appear to be %d '