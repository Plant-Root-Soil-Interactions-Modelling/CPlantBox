#!/bin/bash
#SBATCH --job-name=sa_collect
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.out

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
# TODO: gleiche Umgebung wie in sim.sbatch laden

python collect.py "$DESIGN"
