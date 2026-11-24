#!/bin/bash
# Reicht alles in der richtigen Reihenfolge mit SLURM-Abhängigkeiten ein, direkt aus dem Repo (nichts kopieren):
#
#   cd ~ && bash CPlantBox/experimental/myc/sa_setup/slurm/submit_all.sh [Optionen]
#
#   1. Smoke-Test: 3 echte Läufe (nominal, alle Faktoren unten, alle oben), danach das Gate
#      (smoke_test.py --report --gate). Das Gate schlägt fehl, wenn der Nominal-Lauf nicht sauber durchläuft;
#      dann bricht SLURM alle folgenden Jobs ab (--kill-on-invalid-dep), statt kaputte Läufe zu rechnen.
#   2. nach dem Gate, parallel:
#      - Phase 0 (Rauschen, dx, Faktorraum) -> Sammeln + analyze_phase0 -> Morris -> Sammeln + analyze_morris
#      - 120-h-Test am Nominalpunkt -> analyze_long
#      - Kalibrierung (calibrate.sbatch, ein Job mit vielen Kernen)
#
# Optionen:
#   --quick-smoke      Smoke-Läufe nur 12 h nach t* (schneller, prüft die Verkabelung)
#   --no-smoke         ohne Smoke-Test und Gate (z. B. wenn er schon erfolgreich gelaufen ist)
#   --no-morris        nach Phase 0 aufhören (Phase 0 erst ansehen, Morris später mit submit.sh)
#   --no-long          ohne 120-h-Test
#   --no-calibration   ohne Kalibrierung
#   --no-sa            nur Smoke-Test und Kalibrierung
#   --cal-args "..."   Argumente für calibrate.py, z. B. --cal-args "--seeds 8 --max-nfev 10"
# Umgebung: SA_DIR (Standard ~/CPlantBox/experimental/myc/sa_setup), MAXPAR (gleichzeitige SA-Läufe, 100),
#           SIMTIME_LONG (Zeitlimit der 120-h-Läufe, 12:00:00)
#
# Ausgaben in ./slurmOut, Übersicht: squeue -u $USER. Abgebrochenes einfach erneut einreichen: fertige Läufe
# werden übersprungen (bei erneutem Einreichen mit --no-smoke, wenn der Smoke-Test schon OK war).

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
        *) echo "unbekannte Option $1 (--help)" >&2; exit 1 ;;
    esac
    shift
done

[[ -d "$SA_DIR" ]] || { echo "SA_DIR $SA_DIR fehlt" >&2; exit 1; }
mkdir -p slurmOut slurmErr
echo "sa_setup: $SA_DIR, sbatch-Dateien: $JOBS, Ausgaben: $PWD/slurmOut"
echo "Commit: $(git -C "$SA_DIR" log -1 --format='%h %s' 2>/dev/null || echo '?')"

AFTER=()   # Abhängigkeit aller folgenden Jobs
if [[ $SMOKE == 1 ]]; then
    S=$(sbatch --parsable --export=ALL,QUICK="$QUICK" "$JOBS/smoke.sbatch")
    G=$(sbatch --parsable --dependency=afterany:$S "$JOBS/sa_run.sbatch" smoke_test.py --report --gate $QUICK)
    echo "1. Smoke-Läufe: Job $S, Gate: Job $G  ->  slurmOut/sa_run-$G.out"
    AFTER=(--dependency=afterok:$G --kill-on-invalid-dep=yes)
fi

if [[ $SA == 1 ]]; then
    NEXT=""; [[ $MORRIS == 1 ]] && NEXT=morris
    P=$(sbatch --parsable "${AFTER[@]}" --export=ALL,STAGE=phase0,NEXT="$NEXT",JOBS="$JOBS",SA_DIR="$SA_DIR" "$JOBS/stage.sbatch")
    echo "2. Phase 0$([[ -n $NEXT ]] && echo ' -> Morris'): Stufen-Job $P  ->  slurmOut/sa_stage-$P.out (reicht die Läufe ein)"
    if [[ $LONG == 1 ]]; then
        L=$(sbatch --parsable "${AFTER[@]}" --export=ALL,STAGE=long,NEXT=,JOBS="$JOBS",SA_DIR="$SA_DIR" "$JOBS/stage.sbatch")
        echo "   120-h-Test: Stufen-Job $L"
    fi
fi

if [[ $CAL == 1 ]]; then
    C=$(sbatch --parsable "${AFTER[@]}" "$JOBS/calibrate.sbatch" $CAL_ARGS)
    echo "3. Kalibrierung: Job $C  ->  $SA_DIR/calibration/report.txt"
fi

cat <<EOF

Morgen ansehen:
  Smoke-Gate:     slurmOut/sa_run-*.out          (OK / PROBLEME je Punkt)
  Kalibrierung:   $SA_DIR/calibration/report.txt
  Phase 0:        slurmOut/sa_collect-*.out, $SA_DIR/results/phase0_*.csv
  Morris:         $SA_DIR/results/morris_indices.csv, results/plots/
  Übersicht:      squeue -u \$USER (laufend/wartend), $SA_DIR/results/status_*.csv (fehlend = abgebrochen)
  Abbrüche:       grep -h -i -E 'TIME LIMIT|oom' slurmErr/*.err | sort | uniq -c
Hängen Jobs mit Reason 'DependencyNeverSatisfied', ist eine Stufe davor fehlgeschlagen.
EOF
