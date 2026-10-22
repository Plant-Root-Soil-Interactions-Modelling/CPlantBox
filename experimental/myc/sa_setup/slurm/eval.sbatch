#!/bin/bash
#SBATCH --job-name=sa_eval
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%A_%a.out

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
# TODO: gleiche Umgebung wie in sim.sbatch laden

python evaluate.py "$DESIGN" "$SLURM_ARRAY_TASK_ID" ${FORCE:-}
