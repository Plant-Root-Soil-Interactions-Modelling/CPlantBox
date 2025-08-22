"""
    job managment
    
    just put your plan into make_list()
    then run file __main__
    
    use start_jobs() for multiple slurm jobs (for cluster calling sbatch)                     
"""
import sys; sys.path.append("../../../../CPlantBox");  sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../modules"); sys.path.append("../../../build-cmake/cpp/python_binding/");
import numpy as np
import os

import scenario_static

def start_jobs(jobs):
    """ send as individual jobs """

    jobs = np.array(jobs)

    for job in jobs:

        plant, res, soil= job
        job_name = plant + "_" + res + "_resolution_" + soil
        print(job_name)
        job_file = os.path.join("jobs", job_name + ".sh")

        with open(job_file, 'w') as fh:

            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name={:s}\n".format(job_name))
            fh.writelines("#SBATCH --ntasks=1\n")
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --time=48:00:00\n")
            fh.writelines("#SBATCH --mem=200G\n")
            fh.writelines("#SBATCH --partition=cpu256\n")
            fh.writelines("python3 {:s} {:s} {:s}".format(plant, res, soil))

        os.system("sbatch {:s}".format(job_file))


def make_list():
    jobs = []

    # scenario static
    plant = ['maize', 'wheat'] 
    res = ["high", "low"]  # "high resolution = 1cm", "low resolution = 3cm"
    soil = ['hydrus_loam', 'hydrus_clay', 'hydrus_sandyloam']  # 'hydrus_loam', 'hydrus_clay', 'hydrus_sandyloam'

    print("Creating", len(plant) * len(res) * len(soil), "simulations")
    print()

    for p in plant:
        for r in res:
            for s in soil:
                jobs.append([p, r, s])
    return jobs


if __name__ == "__main__":

    # good luck
    jobs = make_list()
    start_jobs(jobs)  