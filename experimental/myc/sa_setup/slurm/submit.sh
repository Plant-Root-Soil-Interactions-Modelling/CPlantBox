#!/bin/bash
# Reicht Simulation, Auswertung und Sammeln als verkettete SLURM-Jobs ein.
#
#   bash slurm/submit.sh results/design_phase0.csv
#   bash slurm/submit.sh results/design_morris.csv
#   bash slurm/submit.sh results/design_morris.csv --eval-only   # nur neu auswerten (z. B. neue Outputs)
#
# Abgebrochene Läufe einfach durch erneutes Einreichen nachholen: fertige werden übersprungen.

set -euo pipefail
cd "$(dirname "$0")/.."
DESIGN="$1"; MODE="${2:-}"
MAXPAR="${MAXPAR:-100}"            # max. gleichzeitig laufende Tasks

N=$(( $(wc -l < "$DESIGN") - 1 ))
LAST=$(( N - 1 ))
mkdir -p logs
echo "Design: $DESIGN  ($N Läufe)"

if [[ "$MODE" == "--eval-only" ]]; then
  EVAL=$(sbatch --parsable --array=0-$LAST%$MAXPAR --export=ALL,DESIGN="$DESIGN",FORCE=--force slurm/eval.sbatch)
else
  SIM=$(sbatch --parsable --array=0-$LAST%$MAXPAR --export=ALL,DESIGN="$DESIGN" slurm/sim.sbatch)
  echo "Simulation: Job $SIM"
  # afterany: Auswertung startet, wenn alle Simulationen beendet sind (auch bei Timeouts einzelner
  # Läufe — die werden in evaluate.py übersprungen und in status_*.csv als fehlend gemeldet).
  EVAL=$(sbatch --parsable --array=0-$LAST%$MAXPAR --dependency=afterany:$SIM \
         --export=ALL,DESIGN="$DESIGN" slurm/eval.sbatch)
fi
echo "Auswertung: Job $EVAL"
COL=$(sbatch --parsable --dependency=afterany:$EVAL --export=ALL,DESIGN="$DESIGN" slurm/collect.sbatch)
echo "Sammeln:    Job $COL"
echo "Status:     squeue -u \$USER   |   danach: results/status_*.csv"
