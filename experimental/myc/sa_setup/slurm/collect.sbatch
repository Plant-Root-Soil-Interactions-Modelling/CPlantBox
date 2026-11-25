#!/bin/bash
#
# Collects all outputs.csv of a design table and starts the matching analysis (analyze_phase0.py
# or analyze_morris.py), submitted by submit_sa.sh from the home directory.
#
#SBATCH --job-name=sa_collect
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --output=./slurmOut/sa_collect-%j.out
#SBATCH --error=./slurmErr/sa_collect-%j.err

set -euo pipefail
source "$HOME/cpbenv/bin/activate"
cd "${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
export MPLBACKEND=Agg

python3 collect.py "$DESIGN"
case "$(basename "$DESIGN")" in
    design_phase0.csv) python3 analyze_phase0.py ;;
    design_long.csv)   python3 analyze_long.py ;;
    design_morris.csv) python3 evaluate_pooled.py "$DESIGN"; python3 analyze_morris.py ;;
    design_inoculum.csv) python3 evaluate_pooled.py "$DESIGN"; python3 analyze_morris.py inoculum ;;
esac
