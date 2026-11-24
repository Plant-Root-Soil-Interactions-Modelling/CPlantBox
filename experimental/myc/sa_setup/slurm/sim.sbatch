#!/bin/bash
#
# One simulation run per array task (row of the design table $DESIGN), submitted by submit_sa.sh from the home directory.
#
#SBATCH --job-name=sa_sim
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=4:00:00           # smoke test: max. 0.75 h per run; margin for seeds, factor space and fine dx in Phase 0
#SBATCH --mem=4G                 # smoke test: max. 0.85 GB; margin for dx = 0.0025 in Phase 0. For Morris after Phase 0 maybe 2G
#SBATCH --output=./slurmOut/sa_sim-%A_%a.out
#SBATCH --error=./slurmErr/sa_sim-%A_%a.err

set -euo pipefail
source "$HOME/cpbenv/bin/activate"
cd "${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
export MPLBACKEND=Agg

python3 run_one.py "$DESIGN" "$(( ${OFFSET:-0} + SLURM_ARRAY_TASK_ID ))"   # OFFSET: sub-arrays for > 1000 runs
