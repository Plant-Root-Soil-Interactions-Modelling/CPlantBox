#!/bin/bash
#
#SBATCH --job-name=petrisimjob
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=10:00:00
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END
#SBATCH --mail-user=an.heck@fz-juelich.de
#SBATCH --output=./slurmOut/slurm-%j.out
#SBATCH --error=./slurmErr/slurm-%j.err

source cpbenv/bin/activate

cd CPlantBox/experimental/myc/

python3 petridishsimauto.py