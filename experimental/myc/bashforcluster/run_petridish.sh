#!/bin/bash
#
# Job array: SLURM starts this script once per entry of SEEDS, at most 5 copies at the same time (%5).
# Each copy is a normal 1-CPU job like before, $SLURM_ARRAY_TASK_ID (0, 1, 2, ...) picks the seed from SEEDS.
# --array must be 0-(number of seeds - 1), e.g. 100 seeds -> 0-99
#
# Submit from the home directory (same place as the old single job), after C++ changes with the build first:
#   mkdir -p slurmOut slurmErr
#   buildid=$(sbatch --parsable build_plantbox.sh)
#   jobid=$(sbatch --parsable --dependency=afterok:$buildid run_petridish.sh)
#   sbatch --dependency=afterany:$jobid run_petridish_combine.sh
# Seeds that do not reach t* within max_wait_days (petridishsimauto.py) end normally with status "not_crossed".
#
#SBATCH --job-name=petrisimjob
#SBATCH --array=0-99%5
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=10:00:00
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_50,END
#SBATCH --mail-user=an.heck@fz-juelich.de
#SBATCH --output=./slurmOut/slurm-%A_%a.out
#SBATCH --error=./slurmErr/slurm-%A_%a.err

# fixed list of seeds that are simulated
SEEDS=($(seq 10 109))
# SEEDS=(10 25 45 50 70 105 115 125)

if [ "$SLURM_ARRAY_TASK_ID" -ge "${#SEEDS[@]}" ]; then
    echo "Array index $SLURM_ARRAY_TASK_ID out of range, only ${#SEEDS[@]} seeds" >&2
    exit 1
fi
SEED=${SEEDS[$SLURM_ARRAY_TASK_ID]}

source cpbenv/bin/activate

cd CPlantBox/experimental/myc/
mkdir -p plots animation petridish_results/densities
export MPLBACKEND=Agg

echo "Task $SLURM_ARRAY_TASK_ID: seed $SEED on $(hostname), start $(date)"
# runtime per phase and peak memory are stored by petridishsimauto.py itself (see petridishtiming.py),
# since slurm accounting is not available on the cluster (and /usr/bin/time is not installed on the nodes)
python3 petridishsimauto.py "$SEED"
echo "Seed $SEED done $(date)"
