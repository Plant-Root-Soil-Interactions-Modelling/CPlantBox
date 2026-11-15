#!/bin/bash
#SBATCH --job-name=sa_smoke
#SBATCH --array=0-2              # 0 = nominal, 1 = low (alle Faktoren unten), 2 = high (alle oben)
#SBATCH --time=06:00:00          # großzügig: der Smoke-Test soll die echte Laufzeit erst messen
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%A_%a.out
# TODO ggf.: #SBATCH --partition=...   #SBATCH --account=...
#
#   mkdir -p logs && sbatch slurm/smoke.sbatch                   # volle Läufe (HOURS_HYPHAE nach t*)
#   mkdir -p logs && sbatch --export=ALL,QUICK=--quick slurm/smoke.sbatch   # nur 12 h nach t*
#   danach: python smoke_test.py --report

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# TODO: Umgebung laden (wie sim.sbatch), z. B.
# module load python/3.11
# source ~/venvs/sim/bin/activate

POINTS=(nominal low high)
python smoke_test.py "${POINTS[$SLURM_ARRAY_TASK_ID]}" ${QUICK:-}
