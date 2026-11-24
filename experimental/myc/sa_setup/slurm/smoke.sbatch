#!/bin/bash
#
# Smoke test of the sensitivity analysis: 0 = nominal, 1 = low (all factors low), 2 = high (all high), in parallel.
# Submitted by submit_smoke.sh from the home directory (like run_petridish.sh), together with the report job.
#
#SBATCH --job-name=sa_smoke
#SBATCH --array=0-2
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=6:00:00           # generous: the smoke test is what measures the real runtime
#SBATCH --mem=8G
#SBATCH --output=./slurmOut/sa_smoke-%A_%a.out
#SBATCH --error=./slurmErr/sa_smoke-%A_%a.err

set -euo pipefail
source "$HOME/cpbenv/bin/activate"            # additionally needs SALib, scikit-learn, networkx
cd "${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
export MPLBACKEND=Agg

POINTS=(nominal low high)
echo "Smoke test ${POINTS[$SLURM_ARRAY_TASK_ID]} ${QUICK:-} on $(hostname), start $(date)"
python3 smoke_test.py "${POINTS[$SLURM_ARRAY_TASK_ID]}" ${QUICK:-}
echo "End $(date)"
