#!/bin/bash
#
# Submits several petri dish runs with different hyphal parameters as ONE job array (run_petridish.sh) over all runs and seeds,
# then one combine job per run (run_petridish_combine.sh) and a compare job (run_petridish_compare.sh) after all combine jobs.
# Each argument (except the options) is the parameter set of one run: name=value pairs separated by spaces (quoted), or xml for the
# values of the xml. Build CPlantBox by hand before, or with --build as a job before the array.
#
# Usage, from the home directory (where the other job scripts are):
#   bash submit_petridish_compare.sh xml ana=0 distTH=0.005 distTH=0.02 distTH=0.05
#   bash submit_petridish_compare.sh --seeds=10-59 --throttle=20 ana=0 "distTH=0.02 ana=0.5"
# Options:
#   --build          build job (build_plantbox.sh) before, the array waits for it
#   --seeds=A-B      seeds A to B for every run (default 10-109)
#   --throttle=N     at most N tasks of the array at the same time, over all runs (default 50)
#   --no-compare     no compare job (used by submit_petridish.sh for a single run)
# The runs are named <date and time>_<parameters> (xml: only the date and time); the task list and the run names are written to
# petridish_results/tasks/<date and time>.txt and petridish_results/tasks/<date and time>_runs.txt
# Note: the number of tasks (runs x seeds) must be below MaxArraySize of the cluster (scontrol show config | grep MaxArraySize)
#
set -e
build=0
compare=1
first=10
last=109
throttle=50
sets=()
for a in "$@"; do
    case "$a" in
        --build) build=1 ;;
        --no-build) build=0 ;; # still accepted, the default
        --no-compare) compare=0 ;;
        --seeds=*-*) first="${a#--seeds=}"; first="${first%-*}"; last="${a##*-}" ;;
        --throttle=*) throttle="${a#--throttle=}" ;;
        --*) echo "unknown option $a" >&2; exit 1 ;;
        *) sets+=("$a") ;;
    esac
done
if [ "${#sets[@]}" -lt 1 ]; then
    echo "Usage: bash submit_petridish_compare.sh [--build] [--seeds=A-B] [--throttle=N] <parameter set> ...  (parameter set: xml or \"name=value ...\")" >&2
    exit 1
fi
mkdir -p slurmOut slurmErr petridish_results/tasks

# run names and task list: one line "<run> <seed> [name=value ...]" per run and seed
stamp=$(date +%Y-%m-%d_%H-%M-%S)
tasks="$PWD/petridish_results/tasks/$stamp.txt"
: > "$tasks"
runs=()
for s in "${sets[@]}"; do
    params=()
    [ "$s" == "xml" ] || read -r -a params <<< "$s"
    run="$stamp"
    for p in "${params[@]}"; do
        if [[ "$p" != *=* ]]; then
            echo "parameter $p is not name=value" >&2; exit 1
        fi
        run="${run}_${p/=/-}"
    done
    for r in "${runs[@]}"; do
        if [ "$r" == "$run" ]; then
            echo "parameter set $s given twice" >&2; exit 1
        fi
    done
    runs+=("$run")
    for seed in $(seq "$first" "$last"); do
        echo "$run $seed ${params[*]}" >> "$tasks"
    done
done
printf "%s\n" "${runs[@]}" > "${tasks%.txt}_runs.txt"
n=$(wc -l < "$tasks")

dependency=""
if [ "$build" == 1 ]; then
    buildid=$(sbatch --parsable build_plantbox.sh)
    dependency="--dependency=afterok:$buildid"
    echo "build job $buildid"
fi
arrayid=$(sbatch --parsable --array=0-$((n - 1))%$throttle "--export=ALL,PETRIDISH_TASKS=$tasks" $dependency run_petridish.sh)
echo "array job $arrayid: $n tasks (${#runs[@]} runs x seeds $first-$last), at most $throttle at the same time, task list $tasks"

combineids=()
for run in "${runs[@]}"; do
    cid=$(sbatch --parsable --job-name="pcomb${run#$stamp}" "--export=ALL,PETRIDISH_RUN=$run" --dependency=afterany:$arrayid run_petridish_combine.sh)
    combineids+=("$cid")
    echo "run $run: combine job $cid"
done

if [ "$compare" == 1 ] && [ "${#runs[@]}" -gt 1 ]; then
    compareid=$(sbatch --parsable --dependency="afterany:$(IFS=:; echo "${combineids[*]}")" run_petridish_compare.sh "${runs[@]}")
    echo "compare job $compareid after the combine jobs"
fi
