#!/bin/bash
# Submits everything in the right order with SLURM dependencies, directly from the repo (nothing to copy):
#
#   cd ~ && bash CPlantBox/experimental/myc/sa_setup/slurm/submit_all.sh [options]
#
#   1. Smoke test: 3 real runs (nominal, all factors low, all high), then the gate
#      (smoke_test.py --report --gate). The gate fails if the nominal run does not finish cleanly;
#      SLURM then cancels all following jobs (--kill-on-invalid-dep) instead of computing broken runs.
#   2. after the smoke gate, in parallel:
#      - Phase 0 (noise, dx, dt, layer, anaMinDist, CRN, factor space) -> collect + analyze_phase0
#      - long test (150 h) at the nominal point -> collect + analyze_long
#   3. Phase 0 gate (gate.sbatch -> phase0_gate.py, results/gate_report.txt) once both are collected:
#      GO   -> Morris (collect + analyze_morris) and the calibration (calibrate.sbatch, with between_cv from Phase 0)
#      STOP -> nothing more; in both cases an e-mail to SA_MAIL with the report
#
# Options:
#   --quick-smoke      smoke runs only 12 h after t* (faster, checks the wiring)
#   --no-smoke         without smoke test and gate (e.g. if it already passed)
#   --no-morris        the Phase 0 gate does not start Morris (Morris later with submit.sh)
#   --no-long          without the long test (150 h); the Phase 0 gate then checks Phase 0 only
#   --no-calibration   the Phase 0 gate does not start the calibration
#   --no-sa            only smoke test and calibration (no Phase 0 and gate, calibration right after the smoke gate)
#   --cal-args "..."   arguments for calibrate.py, e.g. --cal-args "--seeds 8 --max-nfev 10"
#   --mail ADDRESS     e-mail for the Phase 0 gate (or the environment variable SA_MAIL; without: no e-mail)
# Environment: SA_DIR (default ~/CPlantBox/experimental/myc/sa_setup), MAXPAR (concurrent SA runs, 100),
#              SIMTIME_LONG (time limit of the long runs over LONG_HOURS = 150 h, 16:00:00), SA_MAIL
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
        --mail)           export SA_MAIL="$2"; shift ;;
        -h|--help)        sed -n '2,31p' "$0"; exit 0 ;;
        *) echo "unknown option $1 (--help)" >&2; exit 1 ;;
    esac
    shift
done

[[ -d "$SA_DIR" ]] || { echo "SA_DIR $SA_DIR missing" >&2; exit 1; }
mkdir -p slurmOut slurmErr
echo "sa_setup: $SA_DIR, sbatch files: $JOBS, output: $PWD/slurmOut"
echo "Commit: $(git -C "$SA_DIR" log -1 --format='%h %s' 2>/dev/null || echo '?')"

AFTER=()   # dependency of all following jobs
if [[ $SMOKE == 1 ]]; then
    S=$(sbatch --parsable --export=ALL,QUICK="$QUICK" "$JOBS/smoke.sbatch")
    G=$(sbatch --parsable --dependency=afterany:$S "$JOBS/sa_run.sbatch" smoke_test.py --report --gate $QUICK)
    echo "1. Smoke runs: job $S, gate: job $G  ->  slurmOut/sa_run-$G.out"
    AFTER=(--dependency=afterok:$G --kill-on-invalid-dep=yes)
fi

if [[ $SA == 1 ]]; then
    P=$(sbatch --parsable "${AFTER[@]}" --export=ALL,STAGE=phase0,NEXT=,JOBS="$JOBS",SA_DIR="$SA_DIR" "$JOBS/stage.sbatch")
    echo "2. Phase 0: stage job $P  ->  slurmOut/sa_stage-$P.out (submits the runs)"
    L=""
    if [[ $LONG == 1 ]]; then
        L=$(sbatch --parsable "${AFTER[@]}" --export=ALL,STAGE=long,NEXT=,JOBS="$JOBS",SA_DIR="$SA_DIR" "$JOBS/stage.sbatch")
        echo "   long test (150 h): stage job $L"
    fi
    # afterany: also if a stage job fails, the gate then reports STOP by e-mail
    G=$(sbatch --parsable --dependency=afterany:$P${L:+:$L} \
        --export=ALL,MODE=arm,JOBS="$JOBS",SA_DIR="$SA_DIR",P0_STAGE="$P",LONG_STAGE="$L",GO_MORRIS="$MORRIS",GO_CAL="$CAL",CAL_ARGS="$CAL_ARGS" \
        "$JOBS/gate.sbatch")
    echo "3. Phase 0 gate: job $G (waits for the collection$([[ -n $L ]] && echo ' of Phase 0 and long test'))"
    echo "   GO -> $([[ $MORRIS == 1 ]] && echo 'Morris ')$([[ $CAL == 1 ]] && echo 'calibration')  STOP -> nothing;" \
         "e-mail: ${SA_MAIL:-none (--mail or SA_MAIL not set)}; report: $SA_DIR/results/gate_report.txt"
elif [[ $CAL == 1 ]]; then
    C=$(sbatch --parsable "${AFTER[@]}" "$JOBS/calibrate.sbatch" $CAL_ARGS)
    echo "2. Calibration (without Phase 0 gate): job $C  ->  $SA_DIR/calibration/report.txt"
fi

cat <<EOF

Where to look:
  Smoke gate:     slurmOut/sa_run-*.out          (OK / PROBLEMS per point)
  Phase 0 gate:   $SA_DIR/results/gate_report.txt, slurmOut/sa_gate-*.out
  Calibration:    $SA_DIR/calibration/report.txt
  Phase 0:        slurmOut/sa_collect-*.out, $SA_DIR/results/phase0_*.csv
  Morris:         $SA_DIR/results/morris_indices.csv, results/plots/
  Overview:       squeue -u \$USER (running/pending), $SA_DIR/results/status_*.csv (missing = aborted)
  Aborts:         grep -h -i -E 'TIME LIMIT|oom' slurmErr/*.err | sort | uniq -c
If jobs hang with reason 'DependencyNeverSatisfied', a stage before them failed.
EOF
