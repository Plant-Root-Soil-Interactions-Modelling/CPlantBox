#!/bin/bash
# Submits everything in the right order with SLURM dependencies, directly from the repo (nothing to copy):
#
#   cd ~ && bash CPlantBox/experimental/myc/sa_setup/slurm/submit_all.sh [options]
#
#   1. Smoke test: 3 real runs (nominal, all factors low, all high), then the gate
#      (smoke_test.py --report --gate). The gate fails if the nominal run does not finish cleanly;
#      SLURM then cancels all following jobs (--kill-on-invalid-dep) instead of computing broken runs.
#   2. after the gate, in parallel:
#      - Phase 0 (noise, dx, factor space) -> collect + analyze_phase0 -> Morris -> collect + analyze_morris
#      - 120 h test at the nominal point -> analyze_long
#      - calibration (calibrate.sbatch, one job with many cores)
#
# Options:
#   --quick-smoke      smoke runs only 12 h after t* (faster, checks the wiring)
#   --no-smoke         without smoke test and gate (e.g. if it already passed)
#   --no-morris        stop after Phase 0 (look at Phase 0 first, Morris later with submit.sh)
#   --no-long          without the 120 h test
#   --no-calibration   without calibration
#   --no-sa            only smoke test and calibration
#   --cal-args "..."   arguments for calibrate.py, e.g. --cal-args "--seeds 8 --max-nfev 10"
# Environment: SA_DIR (default ~/CPlantBox/experimental/myc/sa_setup), MAXPAR (concurrent SA runs, 100),
#              SIMTIME_LONG (time limit of the 120 h runs, 12:00:00)
#
# Output in ./slurmOut, overview: squeue -u $USER. Just resubmit anything aborted: finished runs
# are skipped (resubmit with --no-smoke if the smoke test already passed).

set -euo pipefail
JOBS="$(cd "$(dirname "$0")" && pwd)"
export SA_DIR="${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"

SMOKE=1; QUICK=""; MORRIS=1; LONG=1; CAL=1; SA=1; CAL_ARGS=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --quick-smoke)    QUICK="--quick" ;;
        --no-smoke)       SMOKE=0 ;;
        --no-morris)      MORRIS=0 ;;
        --no-long)        LONG=0 ;;
        --no-calibration) CAL=0 ;;
        --no-sa)          SA=0 ;;
        --cal-args)       CAL_ARGS="$2"; shift ;;
        -h|--help)        sed -n '2,27p' "$0"; exit 0 ;;
        *) echo "unknown option $1 (--help)" >&2; exit 1 ;;
    esac
    shift
done

[[ -d "$SA_DIR" ]] || { echo "SA_DIR $SA_DIR missing" >&2; exit 1; }
mkdir -p slurmOut slurmErr
echo "sa_setup: $SA_DIR, sbatch files: $JOBS, output: $PWD/slurmOut"
echo "Commit: $(git -C "$SA_DIR" log -1 --format='%h %s' 2>/dev/null || echo '?')"

AFTER=()   # Abhängigkeit aller folgenden Jobs
if [[ $SMOKE == 1 ]]; then
    S=$(sbatch --parsable --export=ALL,QUICK="$QUICK" "$JOBS/smoke.sbatch")
    G=$(sbatch --parsable --dependency=afterany:$S "$JOBS/sa_run.sbatch" smoke_test.py --report --gate $QUICK)
    echo "1. Smoke runs: job $S, gate: job $G  ->  slurmOut/sa_run-$G.out"
    AFTER=(--dependency=afterok:$G --kill-on-invalid-dep=yes)
fi

if [[ $SA == 1 ]]; then
    NEXT=""; [[ $MORRIS == 1 ]] && NEXT=morris
    P=$(sbatch --parsable "${AFTER[@]}" --export=ALL,STAGE=phase0,NEXT="$NEXT",JOBS="$JOBS",SA_DIR="$SA_DIR" "$JOBS/stage.sbatch")
    echo "2. Phase 0$([[ -n $NEXT ]] && echo ' -> Morris'): stage job $P  ->  slurmOut/sa_stage-$P.out (submits the runs)"
    if [[ $LONG == 1 ]]; then
        L=$(sbatch --parsable "${AFTER[@]}" --export=ALL,STAGE=long,NEXT=,JOBS="$JOBS",SA_DIR="$SA_DIR" "$JOBS/stage.sbatch")
        echo "   120 h test: stage job $L"
    fi
fi

if [[ $CAL == 1 ]]; then
    C=$(sbatch --parsable "${AFTER[@]}" "$JOBS/calibrate.sbatch" $CAL_ARGS)
    echo "3. Calibration: job $C  ->  $SA_DIR/calibration/report.txt"
fi

cat <<EOF

Where to look:
  Smoke gate:     slurmOut/sa_run-*.out          (OK / PROBLEMS per point)
  Calibration:    $SA_DIR/calibration/report.txt
  Phase 0:        slurmOut/sa_collect-*.out, $SA_DIR/results/phase0_*.csv
  Morris:         $SA_DIR/results/morris_indices.csv, results/plots/
  Overview:       squeue -u \$USER (running/pending), $SA_DIR/results/status_*.csv (missing = aborted)
  Aborts:         grep -h -i -E 'TIME LIMIT|oom' slurmErr/*.err | sort | uniq -c
If jobs hang with reason 'DependencyNeverSatisfied', a stage before them failed.
EOF
