import sys; 
import math
import os
import numpy as np
import time
import psutil
start_time_ = time.time()
import gc
from masterI import runSim
from AllDatabase import toTry
import multiprocessing as mulp

from mpi4py import MPI
comm = MPI.COMM_WORLD
maxcore = comm.Get_size()
rank = comm.Get_rank()

N = ""
if len(sys.argv)>1:    
    N = int(sys.argv[1])
directoryN = "/"+os.path.basename(__file__)[:-3]+str(N)+"/"
print("N",N)
main_dir=os.environ['PWD']#dir of the file
results_dir = main_dir +"/results"+directoryN


isCluster = (os.environ['HOME'] == '/home/m.giraud')
params = toTry()
Qsv=params['Qsv']
MulimSucv=params['MulimSucv']
nodeDv=params['nodeDv']
GrRatiov= params['GrRatiov']
kss_v= params['kss_v']
kaa_v= params['kaa_v']
maxrun = len(Qsv)
n_jobs =2

#current_process = psutil.Process()

#"pid",current_process.pid, 
print("rank",rank,"isCluster:",isCluster,"maxcore",maxcore,"to do",maxrun, "directoryN",directoryN)



runSim(directoryN_ = directoryN, doVTP = False, 
        verbosebase = False,
        PRate_ = 6.8e-3, thresholdAux = 0, 
        RatiothresholdAux = 1,
        Qmax_ = Qsv[rank], thresholdSuc = MulimSucv[rank], 
        GrRatio = 10,#GrRatiov[rank], 
        maxLBud = 1., budGR = 0.1,L_dead_threshold=100.,
        kss=kss_v[rank],kaa=kaa_v[rank],
        UseRatiothresholdAux = True,
        nodeD = 3,#nodeDv[rank], 
        thread = rank,
        testTime=0.5, dtBefore = 1/24, dtAfter= 1/24,
        start_time = start_time_,
        doPrint = True, doDict = False,
        dt_write = 0, dtSIM_write = 1,auxin_D=0.)

gc.collect()

# subproc_after = set([p.pid for p in  psutil.Process().children(recursive=True)])# - subproc_before
# print(os.getpid())
# for subproc in subproc_after:
#     kid_process = psutil.Process(subproc)
#     print('Killing n1 process with pid {}'.format(subproc))#, "memory_info",kid_process.memory_info())
    #kid_process.terminate()
    
#slurmstepd: error: *** JOB 6441 ON node03 CANCELLED AT 2022-11-19T10:00:54 ***
#/home/m.giraud/.local/lib/python3.8/site-packages/joblib/externals/loky/backend/resource_tracker.py:318: UserWarning: resource_tracker: There appear to be 8 leaked semlock objects to clean up at shutdown
 # warnings.warn('resource_tracker: There appear to be %d '