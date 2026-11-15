#!/bin/bash
#
# Sammelt alle outputs.csv einer Designtabelle, eingereicht von submit_sa.sh aus dem Home-Verzeichnis.
#
#SBATCH --job-name=sa_collect
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=0:30:00
#SBATCH --mem=8G
#SBATCH --output=./slurmOut/sa_collect-%j.out
#SBATCH --error=./slurmErr/sa_collect-%j.err

set -euo pipefail
source "$HOME/cpbenv/bin/activate"
cd "${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"

python3 collect.py "$DESIGN"
