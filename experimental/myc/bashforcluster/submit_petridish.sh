#!/bin/bash
#
# Submits one complete petri dish run: build of CPlantBox (skip with --no-build), the seed array (run_petridish.sh)
# and the combine job (run_petridish_combine.sh). All jobs get the same run name, date and time of the submission,
# so all outputs go into petridish_results/densities/<run>, plots/<run> and animation/<run> (each with run_info.json).
# Further arguments name=value set hyphal parameters for all seeds (petridishsimauto.py), they are appended to the run name.
#
# Usage, from the home directory (where the other job scripts are):
#   bash submit_petridish.sh                        # build first
#   bash submit_petridish.sh --no-build             # without build
#   bash submit_petridish.sh --no-build ana=0       # network without anastomoses (tip distances for any R, petridishcompare.py)
#   bash submit_petridish.sh --no-build distTH=0.02 # other anastomosis distance [cm]
# Several runs and their comparison in one go: submit_petridish_compare.sh
# Options used by submit_petridish_compare.sh:
#   --after=<jobid>  the seed array waits for this job (e.g. a build job), implies --no-build
#   --print-ids      prints only "<run> <combine job id>"
#
set -e
build=1
after=""
print_ids=0
args=()
for a in "$@"; do
    if [ "$a" == "--no-build" ]; then
        build=0
    elif [[ "$a" == --after=* ]]; then
        after="${a#--after=}"
        build=0
    elif [ "$a" == "--print-ids" ]; then
        print_ids=1
    elif [[ "$a" == *=* ]]; then
        args+=("$a")
    else
        echo "unknown argument $a (expected --no-build, --after=<jobid>, --print-ids or name=value)" >&2
        exit 1
    fi
done
mkdir -p slurmOut slurmErr
export PETRIDISH_ARGS="${args[*]}"
suffix=""
for a in "${args[@]}"; do
    suffix="${suffix}_${a/=/-}"
done
export PETRIDISH_RUN=$(date +%Y-%m-%d_%H-%M-%S)$suffix
# --export passes the variables to the jobs, also on clusters where SBATCH_EXPORT=NONE is the default
EXPORT="--export=ALL,PETRIDISH_RUN=$PETRIDISH_RUN,PETRIDISH_ARGS=$PETRIDISH_ARGS"

dependency=""
if [ "$build" == 1 ]; then
    buildid=$(sbatch --parsable build_plantbox.sh)
    dependency="--dependency=afterok:$buildid"
    [ "$print_ids" == 1 ] || echo "build job $buildid"
elif [ -n "$after" ]; then
    dependency="--dependency=afterok:$after"
fi
jobid=$(sbatch --parsable "$EXPORT" $dependency run_petridish.sh)
combineid=$(sbatch --parsable "$EXPORT" --dependency=afterany:$jobid run_petridish_combine.sh)
if [ "$print_ids" == 1 ]; then
    echo "$PETRIDISH_RUN $combineid"
else
    echo "run $PETRIDISH_RUN: seed array $jobid, combine job $combineid"
fi
# the run names differ by at least one second
sleep 1
