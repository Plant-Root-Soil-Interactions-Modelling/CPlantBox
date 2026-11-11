#!/bin/bash
#
# Submits several petri dish runs with different hyphal parameters (each with submit_petridish.sh: seed array and combine job)
# and a compare job (run_petridish_compare.sh) that runs when all combine jobs are done. Build CPlantBox by hand before,
# or with --build once as a job before all runs (the seed arrays wait for it).
# Each further argument is the parameter set of one run: name=value pairs separated by spaces (quoted), or xml for the values of the xml.
#
# Usage, from the home directory (where the other job scripts are):
#   bash submit_petridish_compare.sh xml ana=0 distTH=0.005 distTH=0.02 distTH=0.05
#   bash submit_petridish_compare.sh --build ana=0 "distTH=0.02 ana=0.5"
# The run names are printed and written to petridish_results/compare_<date>_runs.txt
#
set -e
build=0
sets=()
for a in "$@"; do
    if [ "$a" == "--build" ]; then
        build=1
    elif [ "$a" == "--no-build" ]; then # still accepted, the default now
        build=0
    else
        sets+=("$a")
    fi
done
if [ "${#sets[@]}" -lt 1 ]; then
    echo "Usage: bash submit_petridish_compare.sh [--build] <parameter set> <parameter set> ...  (parameter set: xml or \"name=value ...\")" >&2
    exit 1
fi
mkdir -p slurmOut slurmErr petridish_results

after=()
if [ "$build" == 1 ]; then
    buildid=$(sbatch --parsable build_plantbox.sh)
    after=("--after=$buildid")
    echo "build job $buildid"
fi

runs=()
combineids=()
for s in "${sets[@]}"; do
    params=()
    [ "$s" == "xml" ] || read -r -a params <<< "$s"
    read -r run cid <<< "$(bash submit_petridish.sh "${after[@]}" --print-ids "${params[@]}")"
    echo "run $run (${s}): combine job $cid"
    runs+=("$run")
    combineids+=("$cid")
done

dependency="afterany:$(IFS=:; echo "${combineids[*]}")"
compareid=$(sbatch --parsable --dependency=$dependency run_petridish_compare.sh "${runs[@]}")
list="petridish_results/compare_$(date +%Y-%m-%d_%H-%M-%S)_runs.txt"
printf "%s\n" "${runs[@]}" > "$list"
echo "compare job $compareid after the combine jobs, runs in $list"
