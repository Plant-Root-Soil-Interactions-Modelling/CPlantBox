#!/bin/bash
#
#SBATCH --job-name=parallel-job
#SBATCH --ntasks=512
#SBATCH --nodes=2
#SBATCH --partition=cpu256
#SBATCH --time=10:00
#SBATCH --mem=2G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de
#SBATCH --output=./slurmOut/slurm-%j.out
#SBATCH --error=./slurmErr/slurm-%j.err


export OMPI_MCA_btl=^openib
module load openmpi/4.1.4

cd $HOME/CPBLukas/CPlantBox/experimental/wine_fichtl

source $HOME/cpbenv/bin/activate

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
mpirun -n 25 python3 wine_simulation.py B 0 testmpi
