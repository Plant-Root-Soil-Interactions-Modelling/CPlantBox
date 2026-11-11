#!/bin/bash
#
# Compares several runs (petridishcompare.py): C of f = C rho against R = distTH, from the realized anastomoses and from the tip distances,
# writes ~/petridish_results/plots/compare_<date>/. The names of the runs (folders in petridish_results/densities) are the arguments.
# Submitted by submit_petridish_compare.sh after the combine jobs of the runs, or by hand when the runs are done:
#   sbatch run_petridish_compare.sh <run> <run> ...
#
#SBATCH --job-name=petricompare
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

if [ "$#" -lt 1 ]; then
    echo "Usage: sbatch run_petridish_compare.sh <run> <run> ..." >&2
    exit 1
fi
source cpbenv/bin/activate
# same output folder as run_petridish.sh: reads densities/<run> there, writes plots/compare_<date>
export PETRIDISH_OUT="$PWD/petridish_results"
unset PETRIDISH_RUN # the compare folder gets its own name

cd CPlantBox/experimental/myc/
export MPLBACKEND=Agg

echo "compare runs $*"
python3 petridishcompare.py "$@"
