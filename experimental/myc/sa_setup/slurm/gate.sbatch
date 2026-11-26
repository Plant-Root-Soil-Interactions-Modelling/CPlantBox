#!/bin/bash
#
# Gate after Phase 0 and the long test (submitted by submit_all.sh), two steps:
#   MODE=arm    after the stage jobs: reads the ids of their collect jobs (slurmOut/collect_<stage job>.id) and submits
#               itself with MODE=check after them (afterany: also when a collection or analysis fails, then the gate says STOP)
#   MODE=check  python phase0_gate.py -> results/gate_report.txt. GO: submits the Morris stage (GO_MORRIS=1) and the
#               calibration (GO_CAL=1, CAL_ARGS). STOP: submits nothing. In both cases an e-mail to SA_MAIL with the report.
# Exports: MODE, JOBS, SA_DIR, P0_STAGE, LONG_STAGE (empty = without long test), GO_MORRIS, GO_CAL, CAL_ARGS, SA_MAIL
#
#SBATCH --job-name=sa_gate
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=0:30:00
#SBATCH --mem=4G
#SBATCH --output=./slurmOut/sa_gate-%j.out
#SBATCH --error=./slurmErr/sa_gate-%j.err

set -uo pipefail
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$SUBMIT_DIR"

# e-mail with the file $2 as body; `mail`/`sendmail` if the node has them, otherwise a tiny job whose start SLURM
# announces by mail (subject only, the report stays in results/gate_report.txt)
send_mail() {
    local subject="$1" body="$2"
    if [[ -z "${SA_MAIL:-}" ]]; then
        echo "SA_MAIL not set: no e-mail ($subject)"; return
    fi
    if command -v mail >/dev/null 2>&1; then
        mail -s "$subject" "$SA_MAIL" < "$body" && { echo "e-mail sent (mail): $subject"; return; }
    fi
    if command -v sendmail >/dev/null 2>&1; then
        { echo "Subject: $subject"; echo; cat "$body"; } | sendmail "$SA_MAIL" && { echo "e-mail sent (sendmail): $subject"; return; }
    fi
    local jn; jn="$(echo "$subject" | tr -c 'A-Za-z0-9_\n' '_' | cut -c1-60)"
    sbatch --parsable --job-name="$jn" --mail-type=BEGIN --mail-user="$SA_MAIL" --time=0:01:00 --mem=100M \
        --partition=cpu256 --output=/dev/null --error=/dev/null --wrap="true" >/dev/null \
        && echo "e-mail via SLURM (subject only): $subject"
}

if [[ "$MODE" == "arm" ]]; then
    DEPS=()
    for st in $P0_STAGE ${LONG_STAGE:-}; do
        f="slurmOut/collect_$st.id"
        if [[ ! -s "$f" ]]; then
            msg="$(mktemp)"; echo "Stage job $st did not submit its runs ($f missing), see slurmOut/sa_stage-$st.out" > "$msg"
            send_mail "SA gate STOP: stage $st failed" "$msg"; rm -f "$msg"
            exit 1
        fi
        DEPS+=("$(cat "$f")")
    done
    D=$(IFS=:; echo "${DEPS[*]}")
    G=$(sbatch --parsable --dependency=afterany:$D --export=ALL,MODE=check "$JOBS/gate.sbatch")
    echo "gate check: job $G (after collect job(s) $D)"
    exit 0
fi

# ---- MODE=check ----
source "$HOME/cpbenv/bin/activate"
export MPLBACKEND=Agg
LONGFLAG=""; [[ -z "${LONG_STAGE:-}" ]] && LONGFLAG="--no-long"
(cd "$SA_DIR" && python3 phase0_gate.py $LONGFLAG)
RC=$?
REPORT="${SA_RESULTS_DIR:-$SA_DIR/results}/gate_report.txt"   # SA_RESULTS_DIR as in config.py (tests)
[[ -s "$REPORT" ]] || { REPORT="$(mktemp)"; echo "phase0_gate.py failed (exit $RC), see slurmOut/sa_gate-$SLURM_JOB_ID.out" > "$REPORT"; RC=2; }

if [[ $RC -eq 0 ]]; then
    STARTED=""
    if [[ "${GO_MORRIS:-0}" == 1 ]]; then
        M=$(sbatch --parsable --export=ALL,STAGE=morris,NEXT=,JOBS="$JOBS",SA_DIR="$SA_DIR" "$JOBS/stage.sbatch")
        STARTED+="Morris: stage job $M (slurmOut/sa_stage-$M.out). "
    fi
    if [[ "${GO_CAL:-0}" == 1 ]]; then
        C=$(sbatch --parsable "$JOBS/calibrate.sbatch" ${CAL_ARGS:-})
        STARTED+="Calibration: job $C ($SA_DIR/calibration/report.txt). "
    fi
    echo "$STARTED"
    BODY="$(mktemp)"; { echo "Started: ${STARTED:-nothing (--no-morris --no-calibration)}"; echo; cat "$REPORT"; } > "$BODY"
    send_mail "SA gate GO $(head -1 "$REPORT" | sed 's/^Phase 0 gate: [A-Z]* *//'): Morris/calibration started" "$BODY"; rm -f "$BODY"
else
    BODY="$(mktemp)"; { echo "Nothing started. Fix, then resubmit (finished runs are reused)."; echo; cat "$REPORT"; } > "$BODY"
    send_mail "SA gate STOP $(head -1 "$REPORT" | sed 's/^Phase 0 gate: [A-Z]* *//'): nothing started" "$BODY"; rm -f "$BODY"
    exit 1
fi
