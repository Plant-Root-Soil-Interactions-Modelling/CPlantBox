""" water movement within the root (static soil) """


import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np

import time
from masterI import runSim
start_time_ = time.time()
from AuxTrophCDatabase import toTry

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



params = toTry()
PRatev = params['PRatev']
Qsv=params['Qsv']
Mulimauxv=params['Mulimauxv']
nodeDv=params['nodeDv']
totrun = 0
maxrun = len(PRatev)

n_jobs = min((maxrun - totrun),
                 maxcore - 1)
print("leftToDo:", maxrun - totrun)
tasks_iterator = (delayed(runSim)
                        (directoryN_ = directoryN, 
                         PRate_ = PRatev[i+totrun],
                         threshold = Mulimauxv[i+totrun], 
                         Qmax_ = Qsv[i+totrun],
                         doVTP = False,
                         nodeD = nodeDv[i+totrun],
                        thread = i,
                        verbosebase = False,
                         testTime=7) 
                    for i in range(n_jobs))
parallelizer(tasks_iterator)
totrun += n_jobs