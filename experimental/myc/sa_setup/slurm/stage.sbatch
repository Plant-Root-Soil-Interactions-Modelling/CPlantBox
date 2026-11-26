#!/bin/bash
#
# One stage of submit_all.sh: creates the design (design.py $STAGE) and then submits simulation,
# evaluation and collection with submit.sh. After phase0, if NEXT is set, it submits the next stage, which depends
# on the collect job of Phase 0. As a job because Python should not run on the login node and the design only
# exists afterwards. Submitted by submit_all.sh with --export=ALL,STAGE=...,JOBS=...,SA_DIR=...[,NEXT=...]
#
#SBATCH --job-name=sa_stage
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --output=./slurmOut/sa_stage-%j.out
#SBATCH --error=./slurmErr/sa_stage-%j.err

set -euo pipefail
source "$HOME/cpbenv/bin/activate"
export MPLBACKEND=Agg
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"   # slurmOut/ and slurmErr/ are there
echo "Stage $STAGE on $(hostname), $(date)"

(cd "$SA_DIR" && python3 design.py "$STAGE")

case "$STAGE" in
    phase0) DESIGN=results/design_phase0.csv ;;
    long)   DESIGN=results/design_long.csv; export SIMTIME="${SIMTIME_LONG:-24:00:00}" ;;
    morris) DESIGN=results/design_morris.csv ;;
    *) echo "unknown stage $STAGE" >&2; exit 1 ;;
esac

cd "$SUBMIT_DIR"
ID_FILE="$(mktemp)"
ID_FILE="$ID_FILE" bash "$JOBS/submit.sh" "$DESIGN"
COL="$(cat "$ID_FILE")"; rm -f "$ID_FILE"

if [[ -n "${NEXT:-}" ]]; then
    # next stage only once collection and analysis of this stage succeeded
    NXT=$(sbatch --parsable --dependency=afterok:$COL --kill-on-invalid-dep=yes \
          --export=ALL,STAGE="$NEXT",NEXT=,JOBS="$JOBS",SA_DIR="$SA_DIR" "$JOBS/stage.sbatch")
    echo "Next stage $NEXT: job $NXT (after collect job $COL)"
fi
