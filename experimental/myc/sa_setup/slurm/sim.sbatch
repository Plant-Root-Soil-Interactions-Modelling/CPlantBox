#!/bin/bash
#
# One simulation run per array task (row of the design table $DESIGN), submitted by submit_sa.sh from the home directory.
#
#SBATCH --job-name=sa_sim
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=8:00:00           # 120 h after t*: nominal median 1.6 h; margin for slow seeds, the factor space, dx = 0.005 and dt = 0.25 h
#SBATCH --mem=6G                 # 60 h: max. 0.85 GB; more hyphae over 120 h and fine dx in Phase 0
#SBATCH --output=./slurmOut/sa_sim-%A_%a.out
#SBATCH --error=./slurmErr/sa_sim-%A_%a.err

set -euo pipefail
source "$HOME/cpbenv/bin/activate"
cd "${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
export MPLBACKEND=Agg

python3 run_one.py "$DESIGN" "$(( ${OFFSET:-0} + SLURM_ARRAY_TASK_ID ))"   # OFFSET: sub-arrays for > 1000 runs
