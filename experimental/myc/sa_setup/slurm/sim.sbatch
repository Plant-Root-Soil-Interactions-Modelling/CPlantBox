#!/bin/bash
#SBATCH --job-name=sa_sim
#SBATCH --time=00:30:00          # max. Laufzeit EINES Laufs (+ Puffer)
#SBATCH --mem=4G                 # TODO: nach Phase 0 an echten Bedarf anpassen
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%A_%a.out
# TODO ggf.: #SBATCH --partition=...   #SBATCH --account=...

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# TODO: Umgebung laden, z. B.
# module load python/3.11
# source ~/venvs/sim/bin/activate

python run_one.py "$DESIGN" "$SLURM_ARRAY_TASK_ID"
