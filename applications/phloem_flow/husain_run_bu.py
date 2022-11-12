""" water movement within the root (static soil) """


import sys; 
from joblib import Parallel, delayed
import math
import os
import numpy as np

from datetime import datetime, timedelta
from husain_i import runSim

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
    import shutil
    shutil.rmtree(results_dir)
    os.makedirs(results_dir)

isCluster = (os.environ['HOME'] == '/home/m.giraud')
maxcore = 8
if isCluster:
    maxcore = 256
parallelizer = Parallel(n_jobs=maxcore - 1)

#Qs =np.array([  10., 17.,   28.,   46.,
 #        77.,  120.,  230.,  340.,  450.,  560.,  670.,  780.,  890.,
  #    1000.])*1e-6
#Qs =np.array([ 120.,  230.,  340.,  450.,  560.,  670.,  780.,  890.,
 #     1000.])*1e-6
Qs =np.array([ 100.,  500.,  1000.])*1e-6
#Qs =np.array([ 100.,500.,1000.])*1e-6
#Qs =np.linspace(10,1000,10)*1e-6 #logspace(1,3,num=10)*1e-6
MulimSuc = np.array([0.1*i for i in range(1,9)])
#MulimSuc = np.array([0.1,0.5,0.8])
nodeD = np.array([0,3,4,5,6,7,8])#)[True,False])
kzratio = np.array([0.1,1,10])
rhosratio = np.array([0.1,1,10])
rootGrrate = np.array([0.1,1,10])
maxrun = len(Qs) * len(MulimSuc) * len(nodeD) * len(kzratio) * len(rhosratio) * len(rootGrate)
# QQ_ = np.concatenate((Qs,Qs))
# deca_ = np.concatenate((np.full(True, len(Qs)),np.full(False, len(Qs))))

# QQ__ = np.concatenate([Qs for i in range(len(MulimSuc))])

Qsv,MulimSucv, nodeDv = np.meshgrid(Qs,MulimSuc, nodeD)
Qsv=Qsv.flatten()
MulimSucv=MulimSucv.flatten()
nodeDv=nodeDv.flatten()
assert len(Qsv) == (len(Qs) * len(MulimSuc) * len(nodeD)),"len(Qsv) != maxrun"

totrun = 0
#while totrun < maxrun:
n_jobs = min((maxrun - totrun),
                 maxcore - 1)
print("leftToDo:", maxrun - totrun)
# tasks_iterator = (delayed(runSim)
#                         (directoryN_ = directoryN,
#                          Qmax_ = Qsv[i+totrun],
#                          threshold = MulimSucv[i+totrun], 
#                          doVTP = False,
#                          doDecapitation =decapitatev[i+totrun]) 
#                     for i in range(n_jobs))
tasks_iterator = (delayed(runSim)
                        (directoryN_ = directoryN,
                         Qmax_ = Qsv[i+totrun],
                         threshold = MulimSucv[i+totrun], 
                         doVTP = False,
                         nodeD = nodeDv[i+totrun],
                        thread = i) 
                    for i in range(n_jobs))
parallelizer(tasks_iterator)
totrun += n_jobs