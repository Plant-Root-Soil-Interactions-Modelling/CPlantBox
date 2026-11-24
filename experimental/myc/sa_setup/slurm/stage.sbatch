#!/bin/bash
#
# Eine Stufe von submit_all.sh: erzeugt das Design (design.py $STAGE) und reicht dann mit submit.sh Simulation,
# Auswertung und Sammeln ein. Nach phase0 reicht sie, falls NEXT gesetzt ist, die nächste Stufe ein, die am
# Sammel-Job von Phase 0 hängt. Als Job, weil Python nicht auf dem Login-Knoten laufen soll und das Design erst
# danach existiert. Eingereicht von submit_all.sh mit --export=ALL,STAGE=...,JOBS=...,SA_DIR=...[,NEXT=...]
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
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"   # dort liegen slurmOut/ und slurmErr/
echo "Stufe $STAGE auf $(hostname), $(date)"

(cd "$SA_DIR" && python3 design.py "$STAGE")

case "$STAGE" in
    phase0) DESIGN=results/design_phase0.csv ;;
    long)   DESIGN=results/design_long.csv; export SIMTIME="${SIMTIME_LONG:-12:00:00}" ;;
    morris) DESIGN=results/design_morris.csv ;;
    *) echo "unbekannte Stufe $STAGE" >&2; exit 1 ;;
esac

cd "$SUBMIT_DIR"
ID_FILE="$(mktemp)"
ID_FILE="$ID_FILE" bash "$JOBS/submit.sh" "$DESIGN"
COL="$(cat "$ID_FILE")"; rm -f "$ID_FILE"

if [[ -n "${NEXT:-}" ]]; then
    # nächste Stufe erst, wenn Sammeln und Analyse dieser Stufe erfolgreich waren
    NXT=$(sbatch --parsable --dependency=afterok:$COL --kill-on-invalid-dep=yes \
          --export=ALL,STAGE="$NEXT",NEXT=,JOBS="$JOBS",SA_DIR="$SA_DIR" "$JOBS/stage.sbatch")
    echo "Nächste Stufe $NEXT: Job $NXT (nach Sammel-Job $COL)"
fi
