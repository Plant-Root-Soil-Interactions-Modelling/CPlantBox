#!/bin/bash
# Reicht Simulation, Auswertung und Sammeln als verkettete SLURM-Jobs ein.
# Aus dem Home-Verzeichnis (wie submit_petridish.sh): dieses Skript und die sbatch-Dateien dorthin kopieren
# (z. B. als submit_sa.sh), der Designpfad ist relativ zu sa_setup.
#
#   bash submit_sa.sh results/design_phase0.csv
#   bash submit_sa.sh results/design_morris.csv
#   bash submit_sa.sh results/design_morris.csv --eval-only   # nur neu auswerten (z. B. neue Outputs)
#   SIMTIME=12:00:00 bash submit_sa.sh results/design_long.csv  # Test mit 120 h, längeres Zeitlimit
#
# Abgebrochene Läufe einfach durch erneutes Einreichen nachholen: fertige werden übersprungen.
# sa_setup liegt standardmäßig in ~/CPlantBox/experimental/myc/sa_setup, sonst SA_DIR=... davor setzen.

set -euo pipefail
JOBS="$(cd "$(dirname "$0")" && pwd)"   # hier liegen sim.sbatch, eval.sbatch, collect.sbatch
export SA_DIR="${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
DESIGN="$1"; MODE="${2:-}"
MAXPAR="${MAXPAR:-100}"            # max. gleichzeitig laufende Tasks

N=$(( $(wc -l < "$SA_DIR/$DESIGN") - 1 ))
LAST=$(( N - 1 ))
mkdir -p slurmOut slurmErr         # Ausgaben der Jobs im Verzeichnis, aus dem eingereicht wird
echo "Design: $SA_DIR/$DESIGN  ($N Läufe)"

if [[ "$MODE" == "--eval-only" ]]; then
  EVAL=$(sbatch --parsable --array=0-$LAST%$MAXPAR --export=ALL,DESIGN="$DESIGN",FORCE=--force "$JOBS/eval.sbatch")
else
  # SIMTIME=12:00:00 bash submit_sa.sh ... überschreibt --time von sim.sbatch (z. B. für den Test mit 120 h)
  SIM=$(sbatch --parsable --array=0-$LAST%$MAXPAR ${SIMTIME:+--time=$SIMTIME} --export=ALL,DESIGN="$DESIGN" "$JOBS/sim.sbatch")
  echo "Simulation: Job $SIM"
  # afterany: Auswertung startet, wenn alle Simulationen beendet sind (auch bei Timeouts einzelner
  # Läufe — die werden in evaluate.py übersprungen und in status_*.csv als fehlend gemeldet).
  EVAL=$(sbatch --parsable --array=0-$LAST%$MAXPAR --dependency=afterany:$SIM \
         --export=ALL,DESIGN="$DESIGN" "$JOBS/eval.sbatch")
fi
echo "Auswertung: Job $EVAL"
COL=$(sbatch --parsable --dependency=afterany:$EVAL --export=ALL,DESIGN="$DESIGN" "$JOBS/collect.sbatch")
echo "Sammeln + Analyse: Job $COL  ->  slurmOut/sa_collect-$COL.out"
echo "Status:     squeue -u \$USER   |   danach: $SA_DIR/results/ (status_*.csv, Analysen, plots/)"
