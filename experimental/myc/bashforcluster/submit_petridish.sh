#!/bin/bash
#
# Submits one petri dish run: the seed array (run_petridish.sh) and the combine job (run_petridish_combine.sh), with --build
# a build job (build_plantbox.sh) before, else build by hand before. Arguments name=value set hyphal parameters for all seeds
# (petridishsimauto.py), they are appended to the run name. Outputs go into petridish_results/densities/<run>, plots/<run>
# and animation/<run> (each with run_info.json). Same as submit_petridish_compare.sh with one parameter set and no compare job.
#
# Usage, from the home directory (where the other job scripts are):
#   bash submit_petridish.sh                 # values of the xml
#   bash submit_petridish.sh --build         # build job first, the array waits for it
#   bash submit_petridish.sh ana=0           # network without anastomoses (tip distances for any R, petridishcompare.py)
#   bash submit_petridish.sh distTH=0.02     # other anastomosis distance [cm]
#   bash submit_petridish.sh --seeds=10-59 --throttle=20 distTH=0.02
# Several runs in one array and their comparison: submit_petridish_compare.sh
#
set -e
options=()
params=()
for a in "$@"; do
    case "$a" in
        --*) options+=("$a") ;;
        *=*) params+=("$a") ;;
        *) echo "unknown argument $a (expected options --build, --seeds=A-B, --throttle=N or name=value)" >&2; exit 1 ;;
    esac
done
set_="xml"
[ "${#params[@]}" -eq 0 ] || set_="${params[*]}"
bash "$(dirname "$0")/submit_petridish_compare.sh" --no-compare "${options[@]}" "$set_"
