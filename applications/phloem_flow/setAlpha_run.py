""" water movement within the root (static soil) """


import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np

from datetime import datetime, timedelta
from checkAlpha1601 import runSim

main_dir=os.environ['PWD']#dir of the file
directoryN = "/"+os.path.basename(__file__)[:-3]+"/"#"/a_threshold/"
results_dir = main_dir +"/results"+directoryN

def write_file_array(name, data):
    name2 = 'results'+ directoryN+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(','.join([num for num in map(str, data)])  +'\n')

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
maxcore =256
if maxcore>1:
    parallelizer = Parallel(n_jobs=maxcore - 1)
minAlhpa = 1e-3
maxAlpha = 0.0039249
alphas = np.linspace(minAlhpa,maxAlpha,max(1,maxcore - 1))
maxrun = len(alphas)

totrun = 0
CSVData = open('fig2bRenton2012.txt')
Array2d_result_ = np.loadtxt(CSVData, delimiter="\t")
while totrun < maxrun:
    n_jobs = min((maxrun - totrun),
                    max(1,maxcore - 1))
    print("leftToDo:", maxrun - totrun)
    if maxcore>1:
        tasks_iterator = (delayed(runSim)
                                (directoryN_ = directoryN,
                                 auxin_alpha = alphas[i+totrun],
                                thread = i,
                                verbosebase = False,
                                Array2d_result = Array2d_result_, 
                                 doLogalpha = False) 
                            for i in range(n_jobs))
        result = parallelizer(tasks_iterator)
    else:
        i = 0
        result = list([runSim(directoryN_ = directoryN,
                         auxin_alpha = alphas[i+totrun],
                        thread = i,
                        verbosebase = True,
                        Array2d_result = Array2d_result_,
                              doLogalpha = False)])
    print("resultMSE",result)
    print("alphas",alphas[totrun:(totrun+n_jobs)])
    write_file_array("resultMSE", result)
    write_file_array("alphasMSE", alphas[totrun:(totrun+n_jobs)])
    totrun += n_jobs
    #directoryN_,auxin_alpha =0.075, thread = 0,verbosebase = True, outputSim = np.array([])