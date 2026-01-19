#!/bin/bash
#
#SBATCH --job-name=wine
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=cpu256
#SBATCH --time=20-00:00:00
#SBATCH --mem=9G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de
#SBATCH --output=./slurmOut/slurm-%j.out
#SBATCH --error=./slurmErr/slurm-%j.err


export OMPI_MCA_btl=^openib
module load openmpi/4.1.4

cd $HOME/CPBLukas/CPlantBox/experimental/wine_fichtl

source $HOME/cpbenv/bin/activate

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK python3 wine_simulation.py $1 $2 $3 $4
# OMP_NUM_THREADS=1 python3 wine_simulation.py B 0 dead