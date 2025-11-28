#!/bin/bash

# Array to store all job IDs
all_jobs=()

# Submit B jobs
for i in {0..7}; do
    jid=$(sbatch --parsable wine_24cpu128.sh B "$i" "$1")
    all_jobs+=($jid)
done

# Submit D jobs
for i in {0..7}; do
    jid=$(sbatch --parsable wine_24cpu128.sh D "$i" "$1")
    all_jobs+=($jid)
done

# Submit E jobs
for i in {0..7}; do
    jid=$(sbatch --parsable wine_24cpu128.sh E "$i" "$1")
    all_jobs+=($jid)
done

# Submit the final job that depends on all previous jobs
dep=$(IFS=:; echo "${all_jobs[*]}")  # join job IDs with colon
echo "Dependency string: $dep"
sbatch --dependency=afterok:$dep gatherAndPlotOutputs.sh $1