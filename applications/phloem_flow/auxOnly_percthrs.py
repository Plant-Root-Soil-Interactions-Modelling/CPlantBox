""" water movement within the root (static soil) """


import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np

import time
from masterI import runSim
start_time_ = time.time()

main_dir=os.environ['PWD']#dir of the file
directoryN = "/"+os.path.basename(__file__)[:-3]+"/"#"/a_threshold/"
results_dir = main_dir +"/results"+directoryN



if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    print(len(test))
    for item in test:
        if item != '.ipynb_checkpoints':
            try:
                os.remove(item)
            except:
                pass
                    
    test = os.listdir(results_dir)
    print(len(test))
    if len(test) >1:
        print("dir not empty")
        raise Exception

isCluster = (os.environ['HOME'] == '/home/m.giraud')
maxcore = 256
parallelizer = Parallel(n_jobs=maxcore - 1)

PRate =np.linspace(0.1,10,6)
RatioTh =np.linspace(1,90,6)
nodeD = np.array([0,3,4,5,6,7,8])

maxrun = len(PRate) * len(nodeD) * len(RatioTh)


PRatev,nodeDv, RatioThv = np.meshgrid(PRate, nodeD, RatioTh)
PRatev=PRatev.flatten()
nodeDv=nodeDv.flatten()
RatioThv=RatioThv.flatten()
assert len(PRatev) == (len(PRate) *  len(nodeD)* len(RatioTh)),"len(Qsv) != maxrun"

totrun = 0
while totrun < maxrun:
    n_jobs = min((maxrun - totrun),
                     maxcore - 1)
    print("leftToDo:", maxrun - totrun,maxrun , totrun, n_jobs)
    tasks_iterator = (delayed(runSim)
                            (directoryN_ = directoryN, doVTP = False, verbosebase = False,
                             PRate_ = PRatev[i+totrun], thresholdAux =0, 
                             RatiothresholdAux =RatioThv[i+totrun],
           Qmax_ =0, thresholdSuc = 0.8,
           useCWGr = False, UseRatiothresholdAux = True,
                             nodeD = nodeDv[i+totrun], thread = i,
           activeAtThreshold_auxin = True, activeAtThreshold_suc = False,
                             testTime=7, dtBefore = 1/24, dtAfter= 1/(24*60),
                            start_time = start_time_)
                        for i in range(n_jobs))
    outputs = parallelizer(tasks_iterator)
    totrun += n_jobs