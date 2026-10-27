#!/bin/bash
#
# Run statistics (petridish_timing.png) and radial profile over all seeds (TipDensOverAllSim*.png) from the per-seed results
# of run_petridish.sh. Submit with --dependency=afterany:<jobid of the array>, so it also runs if single seeds failed or timed out.
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

python3 petridishtiming.py
python3 petridishcombine.py
