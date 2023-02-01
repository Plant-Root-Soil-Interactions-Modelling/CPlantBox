import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np
import time
import psutil
start_time_ = time.time()

from masterI_1701 import runAllSim



nameSimS = np.array(["DW","SLM","WT","WTD"])
lightLevelS = np.array(["low","medium","high"])

nameSimS, lightLevelS = np.meshgrid(nameSimS,lightLevelS)
nameSimS=nameSimS.flatten()
lightLevelS=lightLevelS.flatten()
n_jobs = len(nameSimS)

parallelizer = Parallel(n_jobs=n_jobs)

print("n_jobs",n_jobs)

tasks_iterator = (delayed(runAllSim)
                        (i = i_,arg1 = nameSimS[i_],arg2 = lightLevelS[i_])
                    for i_ in range(n_jobs))

results = parallelizer(tasks_iterator)