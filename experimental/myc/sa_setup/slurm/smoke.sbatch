#!/bin/bash
#
# Smoke-Test der Sensitivitätsanalyse: 0 = nominal, 1 = low (alle Faktoren unten), 2 = high (alle oben), parallel.
# Einreichen aus dem Home-Verzeichnis (wie run_petridish.sh), die Datei dorthin kopieren:
#   mkdir -p slurmOut slurmErr
#   sbatch smoke.sbatch                                  # volle Läufe (HOURS_HYPHAE nach t*)
#   sbatch --export=ALL,QUICK=--quick smoke.sbatch       # nur 12 h nach t*
# danach: cd CPlantBox/experimental/myc/sa_setup && python3 smoke_test.py --report
#
#SBATCH --job-name=sa_smoke
#SBATCH --array=0-2
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=6:00:00           # großzügig: der Smoke-Test soll die echte Laufzeit erst messen
#SBATCH --mem=8G
#SBATCH --output=./slurmOut/sa_smoke-%A_%a.out
#SBATCH --error=./slurmErr/sa_smoke-%A_%a.err

set -euo pipefail
source "$HOME/cpbenv/bin/activate"            # braucht zusätzlich SALib, scikit-learn, networkx
cd "${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
export MPLBACKEND=Agg

POINTS=(nominal low high)
echo "Smoke-Test ${POINTS[$SLURM_ARRAY_TASK_ID]} ${QUICK:-} auf $(hostname), Start $(date)"
python3 smoke_test.py "${POINTS[$SLURM_ARRAY_TASK_ID]}" ${QUICK:-}
echo "Ende $(date)"
