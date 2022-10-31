""" water movement within the root (static soil) """


import sys; 
import math
import os
import numpy as np

import time
from masterI import runSim
start_time_ = time.time()

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

    print(results_dir, os.path.exists(results_dir))
    

if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    import shutil
    shutil.rmtree(results_dir)
    os.makedirs(results_dir)
        

isCluster = (os.environ['HOME'] == '/home/m.giraud')
maxcore = 8

PRate =np.array([1+i*2 for i in range(6)])

Mulimaux = np.array([10*(1+i*2) for i in range(6)])
nodeD = np.array([0,3,4,5,6,7,8])

maxrun = len(PRate) * len(Mulimaux) * len(nodeD);maxrun


PRatev,Mulimauxv, nodeDv = np.meshgrid(PRate,Mulimaux, nodeD)
PRatev=PRatev.flatten()
Mulimauxv=Mulimauxv.flatten()
nodeDv=nodeDv.flatten()

assert len(PRatev) == (len(PRate) * len(Mulimaux) * len(nodeD)),"len(Qsv) != maxrun"

totrun = 0
#while totrun < maxrun:
n_jobs = min((maxrun - totrun),
                 maxcore - 1)
temp_time = time.time()
print("leftToDo:", maxrun - totrun,temp_time , start_time_)
print("This sampling run took %5.4f seconds." % (temp_time - start_time_))
#for i in range(n_jobs):
i = 17
runSim(directoryN_ = directoryN, doVTP = False, verbosebase = False,
                         PRate_ = PRatev[i+totrun], thresholdAux = Mulimauxv[i+totrun], 
                         RatiothresholdAux =10,
       Qmax_ = 0, thresholdSuc = 0.8,
       useCWGr = False, UseRatiothresholdAux = False,
                         nodeD = nodeDv[i+totrun], thread = i,
       activeAtThreshold_auxin = True, activeAtThreshold_suc = False,
                         testTime=7, dtBefore = 1/24, dtAfter= 1/(60*24), 
                         start_time = start_time_)
print("did ",i,"runs")
totrun += n_jobs