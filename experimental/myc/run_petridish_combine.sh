#!/bin/bash
#
# Combines the per-seed results of run_petridish.sh, submit with --dependency=afterok:<jobid of the array>
#
#SBATCH --job-name=petricombine
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=01:00:00
#SBATCH --mem=6G
#SBATCH --mail-type=END
#SBATCH --mail-user=an.heck@fz-juelich.de
#SBATCH --output=./slurmOut/slurm-%j.out
#SBATCH --error=./slurmErr/slurm-%j.err

source cpbenv/bin/activate

cd CPlantBox/experimental/myc/
export MPLBACKEND=Agg

python3 petridishcombine.py
