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


def start_jobs(jobs):
    """ send as individual jobs """
    jobs = np.array(jobs)

    for job in jobs:

        plant, size, tube_diam, pdensity = job
        job_name = plant + "_" + size + "_"+ tube_diam +"_" + pdensity
        print(job_name)
        job_file = os.path.join("jobs", job_name + ".sh")

        with open(job_file, 'w') as fh:

            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name={:s}\n".format(job_name))
            fh.writelines("#SBATCH --ntasks=100\n")
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --time=72:00:00\n")
            fh.writelines("#SBATCH --mem=50G\n")
            fh.writelines("#SBATCH --partition=cpu256\n")
            fh.writelines("source /etc/profile.d/modules.sh\n")
            fh.writelines("module load openmpi/4.1.4\n")
            fh.writelines("\n")
            fh.writelines("cd $HOME/Dumux\n")
            fh.writelines("source cpbenv/bin/activate\n")
            fh.writelines("cd $HOME/Dumux/dumux/CPlantBox/tutorial/prld_vrld\n")
            fh.writelines("python3 multiruns_mays_wheat.py {:s} {:s} {:s} {:s}".format(plant, size, tube_diam, pdensity))

        os.system("sbatch {:s}".format(job_file))


def make_list():
    jobs = []

    # scenario static
    plant = ["maize", "wheat"] 
    size = ["small", "large", "complete"] 
    tube_diam = [5.7, 6.4, 7]
    pdensity = ['base', 'alt']

    print("Creating", len(plant) * len(size) * len(tube_diam)* len(pdensity), "simulations")
    print()

    for p in plant:
        for s in size:
            for t in tube_diam:
                for pd in pdensity:
                    jobs.append([p, s, t, pd])

    return jobs


if __name__ == "__main__":

    # good luck
    jobs = make_list()
    start_jobs(jobs)  