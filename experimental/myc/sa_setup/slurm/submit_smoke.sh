#!/bin/bash
# Submits the smoke test: the three runs (smoke.sbatch) and then the report (sa_run.sbatch smoke_test.py --report).
# From the home directory, smoke.sbatch and sa_run.sbatch are next to this script:
#   bash submit_smoke.sh             # full runs (HOURS_HYPHAE after t*)
#   bash submit_smoke.sh --quick     # only 12 h after t*
# Report afterwards in slurmOut/sa_run-<job>.out, the runs in slurmOut/sa_smoke-<job>_<0..2>.out

set -euo pipefail
JOBS="$(cd "$(dirname "$0")" && pwd)"
QUICK="${1:-}"
mkdir -p slurmOut slurmErr
SMOKE=$(sbatch --parsable --export=ALL,QUICK="$QUICK" "$JOBS/smoke.sbatch")
echo "Smoke runs: job $SMOKE"
# afterany: the report also comes if a run fails or hits the time limit
REPORT=$(sbatch --parsable --dependency=afterany:$SMOKE "$JOBS/sa_run.sbatch" smoke_test.py --report)
echo "Report:     job $REPORT  ->  slurmOut/sa_run-$REPORT.out"
