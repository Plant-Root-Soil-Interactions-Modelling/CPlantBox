import sys; 
import math
import os
import numpy as np
import time
start_time_ = time.time()

from pathlib import Path
from AllDatabase import toTry
#zip AllAuxC2_getInput AllAuxC2_getInput/*.txt

main_dir=os.environ['PWD']#dir of the file
directoryN = "/"+os.path.basename(__file__)[:-3]+"/"#"/a_threshold/"
results_dir = main_dir +"/results"+directoryN

def getInput(directoryN_,doVTP, verbosebase,
           PRate_, thresholdAux, RatiothresholdAux,
           Qmax_, thresholdSuc,
           useCWGr, UseRatiothresholdAux,
            nodeD, thread,  
           activeAtThreshold_auxin,activeAtThreshold_suc,
           testTime, dtBefore, dtAfter, start_time, dt_write,dtSIM_write,
           doPrint ,doDict,auxin_D = 0.):
    dt_lastWrote = time.time() - dt_write * 2
    dtSIM_lastWrote = -dtSIM_write*2
    temp_time  = time.time()
    outputsDict_array = {}
    outputsDict_float = {}
    #print("num ",thread,temp_time , start_time, "has started. This sampling run took %5.4f seconds." % (temp_time - start_time))
    doDecapitation = (nodeD > 0)
    directoryN = directoryN_
    strPRate_ = str(np.round(PRate_*10))#[:-2]
    strThA = str(np.round(thresholdAux*10))#[:-2]
    strThRatioA = str(np.round(RatiothresholdAux*10))#[:-2]
    
    strQ = str(np.round(Qmax_*1e6))#[:-2]
    strThS = str(np.round(thresholdSuc*10))#[:-2]
    
    strDecap = str(nodeD)#[0]
    
    dir4allResults ="_"+str(thread) #"_"+ strPRate_+ "_"+strThA+ "_"+strThRatioA + "_"+ strQ + "_"+strThS+"_"+strDecap+"_"+str(thread)
    dir4allResults=dir4allResults.replace(".", "o")
    def write_file_array(name, data):
        if doPrint :
            name2 = 'results'+ directoryN+ name+ dir4allResults+ '.txt'
            with open(name2, 'a') as log:
                log.write(','.join([num for num in map(str, data)])  +'\n')
        if doDict:
            if name not in outputsDict_array:
                outputsDict_array[name] = list()
            outputsDict_array[name].append(data)
    
    namef = 'results'+ directoryN+ "running"+ dir4allResults+ '.txt'
    my_file = Path(namef)
    if my_file.is_file():
        print("input", np.array([directoryN_,doVTP, verbosebase,PRate_, thresholdAux, RatiothresholdAux,
       Qmax_, thresholdSuc,useCWGr, UseRatiothresholdAux,nodeD, thread,testTime, dtBefore, dtAfter,
                                   activeAtThreshold_auxin,activeAtThreshold_suc]))
        print(namef)        
        print("runningfile already exists")
        raise Exception
    
    doDict = False
    doPrint = True
    write_file_array("input", np.array(["directoryN_","doVTP", "verbosebase","PRate_", "thresholdAux", "RatiothresholdAux",
           "Qmax_", "thresholdSuc","useCWGr", "UseRatiothresholdAux","nodeD", "thread","testTime", "dtBefore", "dtAfter",
                                       "activeAtThreshold_auxin","activeAtThreshold_suc,auxin_D"])) 
    write_file_array("input", np.array([directoryN_,doVTP, verbosebase,PRate_, thresholdAux, RatiothresholdAux,
           Qmax_, thresholdSuc,useCWGr, UseRatiothresholdAux,nodeD, thread,testTime, dtBefore, dtAfter,
                                       activeAtThreshold_auxin,activeAtThreshold_suc,auxin_D])) 


if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    import shutil
    shutil.rmtree(results_dir)
    os.makedirs(results_dir)

isCluster = (os.environ['HOME'] == '/home/m.giraud')



params = toTry()
Qsv=params['Qsv']
MulimSucv=params['MulimSucv']
nodeDv=params['nodeDv']
thrsholdv= params['Ratiothrshold']
maxcore = 256
totrun = maxcore - 1
maxrun = len(Qsv)
n_jobs = min((maxrun - totrun),
                 maxcore - 1)
print("leftToDo:", maxrun - totrun)
for i in range(n_jobs):
    getInput(directoryN_ = directoryN, doVTP = False, verbosebase = False,
                         PRate_ = 1., thresholdAux = 0, 
                         RatiothresholdAux =thrsholdv[i + totrun],
       Qmax_ = Qsv[i+totrun], thresholdSuc = MulimSucv[i+totrun], 
       useCWGr = True, UseRatiothresholdAux = True,
                         nodeD = nodeDv[i+totrun], thread = i,
       activeAtThreshold_auxin = True, activeAtThreshold_suc = True,
                         testTime=7, dtBefore = 1/24, dtAfter= 1/(60*24),
                        start_time = start_time_,
                         doPrint = True, doDict = False,
                         dt_write = 0, dtSIM_write = 1, auxin_D=0.)
                    
