import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np
import time
start_time_ = time.time()

from masterI import runSim
from AllDatabase import toTry

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
    maxcore = 256
parallelizer = Parallel(n_jobs=maxcore - 1)

Qs =np.array([ 120.,  230.,  340.,  450.,  560.,  670.,  780.,  890.,
         1000.])*1e-6
nodeD = np.array([0,3,4,5,6,7,8])
maxrun = len(Qs)  * len(nodeD)

Qsv,nodeDv = np.meshgrid(Qs, nodeD)
Qsv=Qsv.flatten()
nodeDv=nodeDv.flatten()
assert len(Qsv) == maxrun,"len(Qsv) != maxrun"

totrun = 0

n_jobs = min((maxrun - totrun),
                 maxcore - 1)
print("leftToDo:", maxrun - totrun)
tasks_iterator = (delayed(runSim)
                        (directoryN_ = directoryN, doVTP = False, verbosebase = False,
                         PRate_ = 0., thresholdAux = 0, 
                         RatiothresholdAux =0.,
       Qmax_ = Qsv[i+totrun], thresholdSuc = 0., 
       useCWGr = True, UseRatiothresholdAux = False,
                         nodeD = nodeDv[i+totrun], thread = i,
       activeAtThreshold_auxin = False, activeAtThreshold_suc = False,
                         testTime=7, dtBefore = 1/24, dtAfter= 1/(60*24),
                        start_time = start_time_,
                         doPrint = True, doDict = False,
                         dt_write = 0, dtSIM_write = 0, auxin_D=0.)
                    for i in range(n_jobs))

parallelizer(tasks_iterator)
totrun += n_jobs