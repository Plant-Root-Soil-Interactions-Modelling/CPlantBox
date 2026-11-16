#!/bin/bash
#
# Auswertung eines Laufs je Array-Task, eingereicht von submit_sa.sh aus dem Home-Verzeichnis.
#
#SBATCH --job-name=sa_eval
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=0:30:00
#SBATCH --mem=4G
#SBATCH --output=./slurmOut/sa_eval-%A_%a.out
#SBATCH --error=./slurmErr/sa_eval-%A_%a.err

set -euo pipefail
source "$HOME/cpbenv/bin/activate"
cd "${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
export MPLBACKEND=Agg

python3 evaluate.py "$DESIGN" "$(( ${OFFSET:-0} + SLURM_ARRAY_TASK_ID ))" ${FORCE:-}   # OFFSET: Teil-Arrays bei > 1000 Läufen
