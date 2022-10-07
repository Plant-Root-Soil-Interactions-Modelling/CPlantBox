""" water movement within the root (static soil) """


import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np

from datetime import datetime, timedelta
from a_threshold import runSim

main_dir=os.environ['PWD']#dir of the file
directoryN = "/"+os.path.basename(__file__)[:-3]+"/"#"/a_threshold/"
results_dir = main_dir +"/results"+directoryN

def runSimTest(directoryN_,Qmax_ = 1000e-6, threshold = 0.8, doVTP = False, doDecapitation = False):
    directoryN = directoryN_
    strQ = str(int(np.round(Qmax_*1e6)))
    strTh = str(int(threshold*10))
    strDecap = str(doDecapitation)[0]
    print(directoryN_,Qmax_ , threshold , doVTP , doDecapitation )
    name = "test"
    name2 = 'results'+ directoryN+ name+ "_"+ strQ + "_"+strTh+"_"+strDecap+ '.txt'
    print(name2)

if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
   test = os.listdir(results_dir)
   for item in test:
       os.remove(results_dir+item)

isCluster = (os.environ['HOME'] == '/home/m.giraud')
maxcore = 8
if isCluster:
    maxcore = 256
parallelizer = Parallel(n_jobs=maxcore - 1)

Qs = np.logspace(1,3,num=10)*1e-6
MulimSuc = np.array([0.1*i for i in range(1,9)])
decapitate = np.array([True,False])
maxrun = len(Qs) * len(MulimSuc) * len(decapitate)

# QQ_ = np.concatenate((Qs,Qs))
# deca_ = np.concatenate((np.full(True, len(Qs)),np.full(False, len(Qs))))

# QQ__ = np.concatenate([Qs for i in range(len(MulimSuc))])

Qsv,MulimSucv, decapitatev = np.meshgrid(Qs,MulimSuc, decapitate)
Qsv=Qsv.flatten()
MulimSucv=MulimSucv.flatten()
decapitatev=decapitatev.flatten()
assert len(Qsv) == maxrun,"len(Qsv) != maxrun"

totrun = 0
while totrun < maxrun:
    n_jobs = min((maxrun - totrun),
                     maxcore - 1)
    print("leftToDo:", maxrun - totrun)
    tasks_iterator = (delayed(runSim)
                            (directoryN_ = directoryN,
                             Qmax_ = Qsv[i+totrun],
                             threshold = MulimSucv[i+totrun], 
                             doVTP = False,
                             doDecapitation =decapitatev[i+totrun]) 
                        for i in range(n_jobs))
    # tasks_iterator = (delayed(runSimTest)
                            # (directoryN_ = directoryN,
                             # Qmax_ = Qsv[i+totrun],
                             # threshold = MulimSucv[i+totrun], 
                             # doVTP = False,
                             # doDecapitation =decapitatev[i+totrun]) 
                        # for i in range(n_jobs))
    parallelizer(tasks_iterator)
    totrun += n_jobs