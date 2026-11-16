#!/bin/bash
# Reicht den Smoke-Test ein: die drei Läufe (smoke.sbatch) und danach den Bericht (sa_run.sbatch smoke_test.py --report).
# Aus dem Home-Verzeichnis, smoke.sbatch und sa_run.sbatch liegen neben diesem Skript:
#   bash submit_smoke.sh             # volle Läufe (HOURS_HYPHAE nach t*)
#   bash submit_smoke.sh --quick     # nur 12 h nach t*
# Bericht danach in slurmOut/sa_run-<job>.out, die Läufe in slurmOut/sa_smoke-<job>_<0..2>.out

set -euo pipefail
JOBS="$(cd "$(dirname "$0")" && pwd)"
QUICK="${1:-}"
mkdir -p slurmOut slurmErr
SMOKE=$(sbatch --parsable --export=ALL,QUICK="$QUICK" "$JOBS/smoke.sbatch")
echo "Smoke-Läufe: Job $SMOKE"
# afterany: der Bericht kommt auch, wenn ein Lauf fehlschlägt oder ins Zeitlimit läuft
REPORT=$(sbatch --parsable --dependency=afterany:$SMOKE "$JOBS/sa_run.sbatch" smoke_test.py --report)
echo "Bericht:     Job $REPORT  ->  slurmOut/sa_run-$REPORT.out"
