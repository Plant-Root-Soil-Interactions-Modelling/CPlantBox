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
    import shutil
    shutil.rmtree(results_dir)
    os.makedirs(results_dir)
        

isCluster = (os.environ['HOME'] == '/home/m.giraud')
maxcore = 8
if isCluster:
    maxcore = 71
parallelizer = Parallel(n_jobs=maxcore - 1)



thrshold = np.array([1/10*i for i in range(10)])
nodeD = np.array([0,3,4,5,6,7,8])

maxrun = len(thrshold) * len(nodeD);maxrun



thrsholdv, nodeDv = np.meshgrid(thrshold,nodeD)
thrsholdv=thrsholdv.flatten()
nodeDv=nodeDv.flatten()

assert len(thrsholdv) == (len(thrshold) * len(nodeD)),"len(Qsv) != maxrun"

totrun = 0
#while totrun < maxrun:
n_jobs = min((maxrun - totrun),
                 maxcore - 1)
temp_time = time.time()
print("leftToDo:", maxrun - totrun,temp_time , start_time_)
print("This sampling run took %5.4f seconds." % (temp_time - start_time_))
tasks_iterator = (delayed(runSim)
                        (directoryN_ = directoryN, doVTP = False, verbosebase = False,
                         PRate_ = 1, thresholdAux = 0, 
                         RatiothresholdAux = thrsholdv[i + totrun] ,maxLBud = 1., budGR = 0.1,L_dead_threshold=100.,
                         kss=1,kaa=0.2,
       Qmax_ = 0, thresholdSuc = 0.8,GrRatio = 10,
       useCWGr = False, UseRatiothresholdAux = True,
                         nodeD = nodeDv[i+totrun], thread = i,
       activeAtThreshold_auxin = True, activeAtThreshold_suc = False,
                         testTime=7, dtBefore = 1/24, dtAfter= 1/(60*24), 
                         start_time = start_time_, 
                         doPrint = True, doDict = False, 
      dt_write = 0,dtSIM_write = 1, #in s
                         auxin_D=0.)
                    for i in range(n_jobs))
parallelizer(tasks_iterator)
totrun += n_jobs