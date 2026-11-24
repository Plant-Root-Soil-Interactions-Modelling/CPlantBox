#!/bin/bash
# Reicht Simulation, Auswertung und Sammeln als verkettete SLURM-Jobs ein.
# Aus dem Home-Verzeichnis (wie submit_petridish.sh): dieses Skript und die sbatch-Dateien dorthin kopieren
# (z. B. als submit_sa.sh), der Designpfad ist relativ zu sa_setup.
#
#   bash submit_sa.sh results/design_phase0.csv
#   bash submit_sa.sh results/design_morris.csv
#   bash submit_sa.sh results/design_morris.csv --eval-only   # nur neu auswerten (z. B. neue Outputs)
#   SIMTIME=12:00:00 bash submit_sa.sh results/design_long.csv  # Test mit 120 h, längeres Zeitlimit
#   SIMTIME=24:00:00 SIMMEM=16G bash submit_sa.sh results/design_phase0.csv   # fehlende Läufe mit mehr Zeit/Speicher
#
# Abgebrochene Läufe einfach durch erneutes Einreichen nachholen: fertige werden übersprungen.
# sa_setup liegt standardmäßig in ~/CPlantBox/experimental/myc/sa_setup, sonst SA_DIR=... davor setzen.

set -euo pipefail
JOBS="$(cd "$(dirname "$0")" && pwd)"   # hier liegen sim.sbatch, eval.sbatch, collect.sbatch
export SA_DIR="${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
DESIGN="$1"; MODE="${2:-}"
MAXPAR="${MAXPAR:-100}"            # max. gleichzeitig laufende Tasks (insgesamt)
CHUNK="${CHUNK:-1000}"             # max. Tasks je Array: SLURM begrenzt Array-Indizes (MaxArraySize, oft 1001)

N=$(( $(wc -l < "$SA_DIR/$DESIGN") - 1 ))
NCH=$(( (N + CHUNK - 1) / CHUNK ))
PAR=$(( (MAXPAR + NCH - 1) / NCH ))   # je Teil-Array, damit insgesamt höchstens ~MAXPAR laufen
mkdir -p slurmOut slurmErr         # Ausgaben der Jobs im Verzeichnis, aus dem eingereicht wird
echo "Design: $SA_DIR/$DESIGN  ($N runs, $NCH sub-array(s))"

# reicht $1 (sbatch-Datei) in Teil-Arrays ein, Zeile = OFFSET + Array-Index; gibt die Job-IDs als a:b:c aus
submit_arrays() {
  local script="$1"; shift
  local ids=() off n
  for (( off = 0; off < N; off += CHUNK )); do
    n=$(( N - off < CHUNK ? N - off : CHUNK ))
    ids+=( "$(sbatch --parsable --array=0-$(( n - 1 ))%$PAR "$@" --export=ALL,DESIGN="$DESIGN",OFFSET=$off${FORCE_EXPORT:-} "$JOBS/$script")" )
  done
  local IFS=:; echo "${ids[*]}"
}

if [[ "$MODE" == "--eval-only" ]]; then
  EVAL=$(FORCE_EXPORT=",FORCE=--force" submit_arrays eval.sbatch)
else
  # SIMTIME=12:00:00 bash submit_sa.sh ... überschreibt --time von sim.sbatch (z. B. für den Test mit 120 h)
  # SIMMEM=16G überschreibt --mem (z. B. für feines dx in Phase 0)
  SIM=$(submit_arrays sim.sbatch ${SIMTIME:+--time=$SIMTIME} ${SIMMEM:+--mem=$SIMMEM})
  echo "Simulation: job(s) $SIM"
  # afterany: Auswertung startet, wenn alle Simulationen beendet sind (auch bei Timeouts einzelner
  # Läufe — die werden in evaluate.py übersprungen und in status_*.csv als fehlend gemeldet).
  EVAL=$(submit_arrays eval.sbatch --dependency=afterany:$SIM)
fi
echo "Evaluation: job(s) $EVAL"
COL=$(sbatch --parsable --dependency=afterany:$EVAL --export=ALL,DESIGN="$DESIGN" "$JOBS/collect.sbatch")
echo "Collect + analysis: job $COL  ->  slurmOut/sa_collect-$COL.out"
[[ -n "${ID_FILE:-}" ]] && echo "$COL" > "$ID_FILE"   # für submit_all.sh: nächste Stufe hängt am Sammel-Job
echo "Status:     squeue -u \$USER   |   afterwards: $SA_DIR/results/ (status_*.csv, analyses, plots/)"
