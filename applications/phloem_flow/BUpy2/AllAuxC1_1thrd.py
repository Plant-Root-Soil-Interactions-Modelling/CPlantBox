import sys; 
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

params = toTry()
Qsv=params['Qsv']
MulimSucv=params['MulimSucv']
nodeDv=params['nodeDv']
thrsholdv= params['Ratiothrshold']
totrun = 0
maxrun = len(Qsv)

n_jobs = min((maxrun - totrun),
                 maxcore - 1)
print("leftToDo:", maxrun - totrun)

def runSim(directoryN_,doVTP, verbosebase,
           PRate_, thresholdAux, RatiothresholdAux,
           Qmax_, thresholdSuc,
           useCWGr, UseRatiothresholdAux,
            nodeD, thread,  
           activeAtThreshold_auxin,activeAtThreshold_suc,
           testTime, dtBefore, dtAfter, start_time, dt_write,dtSIM_write,
           doPrint ,doDict,auxin_D = 0.):
    print(thread ,nodeD,RatiothresholdAux,thresholdSuc,Qmax_)

for i in range(n_jobs):
    runSim(directoryN_ = directoryN, doVTP = False, verbosebase = False,
                         PRate_ = 1., thresholdAux = 0, 
                         RatiothresholdAux =thrsholdv[i + totrun],
       Qmax_ = Qsv[i+totrun], thresholdSuc = MulimSucv[i+totrun], 
       useCWGr = True, UseRatiothresholdAux = True,
                         nodeD = nodeDv[i+totrun], thread = i,
       activeAtThreshold_auxin = True, activeAtThreshold_suc = True,
                         testTime=7, dtBefore = 1/24, dtAfter= 1/24,
                        start_time = start_time_,
                         doPrint = True, doDict = False,
                         dt_write = 0, dtSIM_write = 0, auxin_D=0.)

totrun += n_jobs