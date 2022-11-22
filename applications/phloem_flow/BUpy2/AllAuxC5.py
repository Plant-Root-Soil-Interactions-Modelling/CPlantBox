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



totrun = (maxcore - 1) * (int(directoryN[-1]) -1)
AllAuxCmasterFunc(totrun,maxcore, directoryN)