from benchmark_selhausen import run_benchmark
from joblib import Parallel, delayed
import os
import numpy as np
import pickle

maxCore = os.cpu_count() -1   # Change this to the number of available cores on your machine

def run_benchmark_(num):
    fx = run_benchmark()
    return fx

# n_jobs = maxCore
n_jobs = 1
 

runs= 100 #num repetion

parallelizer = Parallel(n_jobs=n_jobs)
tasks_iterator = ( delayed(run_benchmark_)(i) 
                  for i in range(runs))
result = parallelizer( tasks_iterator )
# Merging the output of the jobs
fx = np.vstack(result)


with open('wheat'+'.pkl','wb') as f:
     pickle.dump(fx,f, protocol=pickle.HIGHEST_PROTOCOL)