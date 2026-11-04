#!/bin/bash
#
# Run statistics (petridish_timing.png), radial and spatial/temporal profiles over all seeds (TipDensOverAllSim*.png,
# SpatialTemporalProfile.png) and wave statistics (petridish_stats.png, petridish_wave_fits.csv) from the per-seed results
# of run_petridish.sh, all in ~/petridish_results. Submit with --dependency=afterany:<jobid of the array>, so it also runs if single seeds failed or timed out.
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
# same output folder as run_petridish.sh: reads densities/ there, writes plots/ and petridish_wave_fits.csv there
export PETRIDISH_OUT="$PWD/petridish_results"

cd CPlantBox/experimental/myc/
export MPLBACKEND=Agg

python3 petridishtiming.py
python3 petridishcombine.py
python3 petridishstats.py
