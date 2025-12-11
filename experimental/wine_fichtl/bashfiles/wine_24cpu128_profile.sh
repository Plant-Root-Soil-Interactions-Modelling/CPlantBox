#!/bin/bash
#
#SBATCH --job-name=wine
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=cpu128
#SBATCH --time=20-00:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de
#SBATCH --output=./slurmOut/slurm-%j.out
#SBATCH --error=./slurmErr/slurm-%j.err


export OMPI_MCA_btl=^openib
module load openmpi/4.1.4

cd $HOME/CPBLukas/CPlantBox/experimental/wine_fichtl

source $HOME/cpbenv/bin/activate


export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} py-spy record -o profile.svg -- python3 wine_simulation.py $1 $2 $3
# perf record -g 
# valgrind --tool=callgrind 