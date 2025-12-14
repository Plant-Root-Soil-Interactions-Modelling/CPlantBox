#!/bin/bash
#
#SBATCH --job-name=profile
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=cpu128
#SBATCH --time=20-00:00:00
#SBATCH --mem=60G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END,FAIL,ALL
#SBATCH --mail-user=m.giraud@fz-juelich.de
#SBATCH --output=./slurmOut/slurm-%j.out
#SBATCH --error=./slurmErr/slurm-%j.err


export OMPI_MCA_btl=^openib
module load openmpi/4.1.4

cd $HOME/CPBLukas/CPlantBox/experimental/wine_fichtl

source $HOME/cpbenv/bin/activate


#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} py-spy record --rate 20 --subprocess -o profile2.svg -- python3 wine_simulation.py $1 $2 $3 $4
#OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} perf record -g python3 wine_simulation.py $1 $2 $3 $4
# valgrind --tool=callgrind 

OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK py-spy record \
    --rate 20 \
    --subprocesses \
    -o profile3.svg \
    -- python3 wine_simulation.py $1 $2 $3
#OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK $HOME/valgrind-3.22.0/coregrind/valgrind \
#    --tool=callgrind \
#    python3 wine_simulation.py $1 $2 $3 $4