#!/bin/bash
#
# Run statistics (RunTiming.png), radial and spatial/temporal profiles over all seeds (TipDensOverAllSim*.png,
# SpatialTemporalProfile.png) and wave statistics (WaveStatistics.png, WaveFitsPerRing.csv) from the per-seed results
# of run_petridish.sh, all in ~/petridish_results/plots/<run>. Submitted by submit_petridish.sh with --dependency=afterany:<jobid of the array>,
# so it also runs if single seeds failed or timed out. Run = $PETRIDISH_RUN, without it the newest run in petridish_results/densities.
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
# same output folder as run_petridish.sh: reads densities/<run> there, writes plots/<run> (plots, csv tables, run_info.json)
export PETRIDISH_OUT="$PWD/petridish_results"
echo "run ${PETRIDISH_RUN:-(newest)}"

cd CPlantBox/experimental/myc/
export MPLBACKEND=Agg

python3 petridishtiming.py
python3 petridishcombine.py
python3 petridishstats.py
python3 petridishwave.py
