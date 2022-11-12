""" water movement within the root (static soil) """


import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np
import time
start_time_ = time.time()

from masterI import runSim
from ConlyDatabase import toTry

main_dir=os.environ['PWD']#dir of the file
directoryN = "/"+os.path.basename(__file__)[:-3]+"/"#"/a_threshold/"
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
    maxcore = 110
parallelizer = Parallel(n_jobs=maxcore - 1)


Qs =np.array([ 120.,  230.,  340.,  450.,  560., ])*1e-6# 670.,  780.,  890.,1000.
MulimSuc = np.array([0.1*i for i in range(1,4)])
nodeD = np.array([0,3,4,5,6,7,8])
maxrun = len(Qs) * len(MulimSuc) * len(nodeD)

Qsv,MulimSucv, nodeDv = np.meshgrid(Qs,MulimSuc, nodeD)
Qsv=Qsv.flatten()
MulimSucv=MulimSucv.flatten()
nodeDv=nodeDv.flatten()
totrun = 0
maxrun = len(Qsv)

n_jobs = min((maxrun - totrun),
                 maxcore - 1)
print("leftToDo:", maxrun - totrun)
tasks_iterator = (delayed(runSim)
                        (directoryN_ = directoryN, doVTP = False, verbosebase = False,
                         PRate_ = 0, thresholdAux = 0, RatiothresholdAux =0,
       Qmax_ = Qsv[i+totrun], thresholdSuc = MulimSucv[i+totrun], 
       useCWGr = True, UseRatiothresholdAux = False,
                         nodeD = nodeDv[i+totrun], thread = i,
       activeAtThreshold_auxin = False, activeAtThreshold_suc = True,
                         testTime=7, dtBefore = 1/24, dtAfter= 1/(60*24),
                        start_time = start_time_,
                         doPrint = True, doDict = False,
                         dt_write = 0, dtSIM_write = 10/(60*24), auxin_D=0.)
                    for i in range(n_jobs))

parallelizer(tasks_iterator)
totrun += n_jobs