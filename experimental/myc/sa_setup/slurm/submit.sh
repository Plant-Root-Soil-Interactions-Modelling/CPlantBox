#!/bin/bash
# Submits simulation, evaluation and collection as chained SLURM jobs.
# From the home directory (like submit_petridish.sh): copy this script and the sbatch files there
# (e.g. as submit_sa.sh), the design path is relative to sa_setup.
#
#   bash submit_sa.sh results/design_phase0.csv
#   bash submit_sa.sh results/design_morris.csv
#   bash submit_sa.sh results/design_morris.csv --eval-only   # only re-evaluate (e.g. new outputs)
#   SIMTIME=24:00:00 bash submit_sa.sh results/design_long.csv  # long test (180 h), longer time limit
#   SIMTIME=24:00:00 SIMMEM=16G bash submit_sa.sh results/design_phase0.csv   # missing runs with more time/memory
#
# Catch up on aborted runs by simply resubmitting: finished ones are skipped.
# sa_setup is in ~/CPlantBox/experimental/myc/sa_setup by default, otherwise prefix SA_DIR=...

set -euo pipefail
JOBS="$(cd "$(dirname "$0")" && pwd)"   # sim.sbatch, eval.sbatch, collect.sbatch are here
export SA_DIR="${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
DESIGN="$1"; MODE="${2:-}"
MAXPAR="${MAXPAR:-100}"            # max. concurrently running tasks (in total)
CHUNK="${CHUNK:-1000}"             # max. tasks per array: SLURM limits array indices (MaxArraySize, often 1001)

N=$(( $(wc -l < "$SA_DIR/$DESIGN") - 1 ))
NCH=$(( (N + CHUNK - 1) / CHUNK ))
PAR=$(( (MAXPAR + NCH - 1) / NCH ))   # per sub-array, so that at most ~MAXPAR run in total
mkdir -p slurmOut slurmErr         # job output in the directory submitted from
echo "Design: $SA_DIR/$DESIGN  ($N runs, $NCH sub-array(s))"

# submits $1 (sbatch file) in sub-arrays, row = OFFSET + array index; prints the job IDs as a:b:c
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
  # SIMTIME=12:00:00 bash submit_sa.sh ... overrides --time of sim.sbatch (e.g. for the long test)
  # SIMMEM=16G overrides --mem (e.g. for fine dx in Phase 0)
  SIM=$(submit_arrays sim.sbatch ${SIMTIME:+--time=$SIMTIME} ${SIMMEM:+--mem=$SIMMEM})
  echo "Simulation: job(s) $SIM"
  # afterany: evaluation starts when all simulations have ended (also on timeouts of single
  # runs — those are skipped in evaluate.py and reported as missing in status_*.csv).
  EVAL=$(submit_arrays eval.sbatch --dependency=afterany:$SIM)
fi
echo "Evaluation: job(s) $EVAL"
COL=$(sbatch --parsable --dependency=afterany:$EVAL --export=ALL,DESIGN="$DESIGN" "$JOBS/collect.sbatch")
echo "Collect + analysis: job $COL  ->  slurmOut/sa_collect-$COL.out"
[[ -n "${ID_FILE:-}" ]] && echo "$COL" > "$ID_FILE"   # for submit_all.sh: the next stage depends on the collect job
echo "Status:     squeue -u \$USER   |   afterwards: $SA_DIR/results/ (status_*.csv, analyses, plots/)"
