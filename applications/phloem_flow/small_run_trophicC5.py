""" water movement within the root (static soil) """


import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np

from datetime import datetime, timedelta
from small_i_trophicC import runSim
from AuxTrophCDatabase import toTry

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
        try:
            os.remove(results_dir+item)
        except:
            pass

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
totrun = (maxcore - 1)*4
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