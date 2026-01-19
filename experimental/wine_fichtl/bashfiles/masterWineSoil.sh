#!/bin/bash

# Array to store all job IDs
all_jobs=()

# Submit B jobs
for i in {1,2,3,5,8,13,20,32,50}; do
    jid=$(sbatch --nodelist=node04 --parsable wineSoil_24cpu256.sh B 0 "$i" 0 "$1")
    all_jobs+=($jid)
done

# Submit D jobs
for i in {1,2,3,5,8,13,20,32,50}; do
    jid=$(sbatch --nodelist=node04 --parsable wineSoil_24cpu256.sh D 0 "$i" 0 "$1")
    all_jobs+=($jid)
done

# Submit E jobs
for i in {1,2,3,5,8,13,20,32,50}; do
    jid=$(sbatch --nodelist=node04 --parsable wineSoil_24cpu256.sh E 0 "$i" 0 "$1")
    all_jobs+=($jid)
done

for i in {1,2,3,5,8,13,20,32,50}; do
    jid=$(sbatch --nodelist=node04 --parsable wineSoil_24cpu256.sh B 0 "$i" 1 "$1")
    all_jobs+=($jid)
done

# Submit D jobs
for i in {1,2,3,5,8,13,20,32,50}; do
    jid=$(sbatch --nodelist=node04 --parsable wineSoil_24cpu256.sh D 0 "$i" 1 "$1")
    all_jobs+=($jid)
done

# Submit E jobs
for i in {1,2,3,5,8,13,20,32,50}; do
    jid=$(sbatch --nodelist=node05 --parsable wineSoil_24cpu256.sh E 0 "$i" 1 "$1")
    all_jobs+=($jid)
done

# Submit the final job that depends on all previous jobs
dep=$(IFS=:; echo "${all_jobs[*]}")  # join job IDs with colon
echo "Dependency string: $dep"
sbatch --nodelist=node16 --dependency=afterany:$dep gatherAndPlotOutputsSoil.sh $1