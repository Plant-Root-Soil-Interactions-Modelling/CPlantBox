#!/bin/bash
#
# Sammelt alle outputs.csv einer Designtabelle und startet die passende Analyse (analyze_phase0.py
# bzw. analyze_morris.py), eingereicht von submit_sa.sh aus dem Home-Verzeichnis.
#
#SBATCH --job-name=sa_collect
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=1:00:00
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
    design_morris.csv) python3 analyze_morris.py ;;
esac
