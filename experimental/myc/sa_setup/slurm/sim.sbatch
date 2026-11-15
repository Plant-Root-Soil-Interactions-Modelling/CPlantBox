#!/bin/bash
#
# Ein Simulationslauf je Array-Task (Zeile der Designtabelle $DESIGN), eingereicht von submit_sa.sh aus dem Home-Verzeichnis.
#
#SBATCH --job-name=sa_sim
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=4:00:00           # max. Laufzeit EINES Laufs (+ Puffer), Seeds brauchen bis zu 3 h; Phase 0 mit kleinem dx länger
#SBATCH --mem=4G                 # nach dem Smoke-Test an den echten Bedarf anpassen
#SBATCH --output=./slurmOut/sa_sim-%A_%a.out
#SBATCH --error=./slurmErr/sa_sim-%A_%a.err

set -euo pipefail
source "$HOME/cpbenv/bin/activate"
cd "${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
export MPLBACKEND=Agg

python3 run_one.py "$DESIGN" "$SLURM_ARRAY_TASK_ID"
