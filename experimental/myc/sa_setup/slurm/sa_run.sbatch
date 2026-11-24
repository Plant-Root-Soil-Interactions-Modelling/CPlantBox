#!/bin/bash
#
# Führt ein Python-Skript aus sa_setup als kurzen Job aus, für alles, was nicht auf dem Login-Knoten laufen darf.
# Einreichen aus dem Home-Verzeichnis, Skript und Argumente hinter der sbatch-Datei:
#   sbatch sa_run.sbatch check_config.py
#   sbatch sa_run.sbatch design.py phase0
#   sbatch sa_run.sbatch smoke_test.py --report
#   sbatch sa_run.sbatch analyze_morris.py
# Ausgabe: slurmOut/sa_run-<job>.out (z. B. cat slurmOut/sa_run-*.out | tail)
#
#SBATCH --job-name=sa_run
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --output=./slurmOut/sa_run-%j.out
#SBATCH --error=./slurmErr/sa_run-%j.err

set -euo pipefail
if [ $# -eq 0 ]; then
    echo "give a script, e.g. sbatch sa_run.sbatch check_config.py" >&2
    exit 1
fi
source "$HOME/cpbenv/bin/activate"
cd "${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
export MPLBACKEND=Agg

echo "python3 $* in $PWD on $(hostname), $(date)"
python3 "$@"
