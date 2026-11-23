#!/bin/bash
#
# Kalibrierung (calibrate.py) als EIN Job: die Seeds je Auswertung laufen parallel auf den Kernen des Jobs.
# Einreichen aus dem Home-Verzeichnis, Argumente werden an calibrate.py weitergereicht:
#   sbatch calibrate.sbatch
#   sbatch calibrate.sbatch --seeds 16 --max-iter 6
# Abgebrochen (Timeout)? Einfach erneut einreichen: fertige Läufe werden übernommen.
# Ausgabe: slurmOut/sa_cal-<job>.out, Ergebnisse in sa_setup/calibration/
#
#SBATCH --job-name=sa_cal
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16       # = N_SEEDS in calibrate.py, dann eine Runde je Auswertung (~0,75 h)
#SBATCH --partition=cpu256
#SBATCH --time=48:00:00          # ~ 5 + 3 × MAX_ITER Auswertungen im ungünstigsten Fall
#SBATCH --mem=32G                # Smoke-Test: max. 0,85 GB je Lauf
#SBATCH --output=./slurmOut/sa_cal-%j.out
#SBATCH --error=./slurmErr/sa_cal-%j.err

set -euo pipefail
source "$HOME/cpbenv/bin/activate"
cd "${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
export MPLBACKEND=Agg

echo "python3 calibrate.py $* in $PWD auf $(hostname), $(date)"
python3 calibrate.py "$@"
