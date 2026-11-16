#!/bin/bash
#
# Ein Simulationslauf je Array-Task (Zeile der Designtabelle $DESIGN), eingereicht von submit_sa.sh aus dem Home-Verzeichnis.
#
#SBATCH --job-name=sa_sim
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=4:00:00           # Smoke-Test: max. 0,75 h je Lauf; Puffer für Seeds, Faktorraum und feines dx in Phase 0
#SBATCH --mem=4G                 # Smoke-Test: max. 0,85 GB; Puffer für dx = 0,0025 in Phase 0. Für Morris nach Phase 0 ggf. 2G
#SBATCH --output=./slurmOut/sa_sim-%A_%a.out
#SBATCH --error=./slurmErr/sa_sim-%A_%a.err

set -euo pipefail
source "$HOME/cpbenv/bin/activate"
cd "${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
export MPLBACKEND=Agg

python3 run_one.py "$DESIGN" "$(( ${OFFSET:-0} + SLURM_ARRAY_TASK_ID ))"   # OFFSET: Teil-Arrays bei > 1000 Läufen
