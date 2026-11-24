#!/bin/bash
#
# Kalibrierung (calibrate.py) als EIN Job: die Seeds je Auswertung laufen parallel auf den Kernen des Jobs.
# Einreichen aus dem Home-Verzeichnis, Argumente werden an calibrate.py weitergereicht:
#   sbatch calibrate.sbatch
#   sbatch calibrate.sbatch --seeds 16 --max-nfev 10
# Abgebrochen (Timeout)? Einfach erneut einreichen: fertige Läufe werden übernommen.
# Ausgabe: slurmOut/sa_cal-<job>.out, Ergebnisse in sa_setup/calibration/
#
#SBATCH --job-name=sa_cal
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128      # Jacobi: 16 Differenzenpunkte × 16 Seeds = 256 Läufe -> 3 Runden à ~0,75 h
#SBATCH --partition=cpu256
#SBATCH --time=72:00:00          # je Iteration ~4 Runden ≈ 3 h; python calibrate.py --check zeigt den Aufwand
#SBATCH --mem=200G               # Smoke-Test: max. 0,85 GB je Lauf, 128 gleichzeitig
#SBATCH --output=./slurmOut/sa_cal-%j.out
#SBATCH --error=./slurmErr/sa_cal-%j.err

set -euo pipefail
source "$HOME/cpbenv/bin/activate"
cd "${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
export MPLBACKEND=Agg

echo "python3 calibrate.py $* in $PWD auf $(hostname), $(date)"
python3 calibrate.py "$@"
