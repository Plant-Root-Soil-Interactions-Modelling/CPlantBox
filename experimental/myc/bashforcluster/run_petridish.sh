#!/bin/bash
#
# Job array: SLURM starts this script once per entry of SEEDS, at most 5 copies at the same time (%5).
# Each copy is a normal 1-CPU job like before, $SLURM_ARRAY_TASK_ID (0, 1, 2, ...) picks the seed from SEEDS.
# --array must be 0-(number of seeds - 1), e.g. 100 seeds -> 0-99
#
# Submit from the home directory with submit_petridish.sh (build, this array and the combine job, all with the same run name):
#   bash submit_petridish.sh [--no-build]
# Outputs go into petridish_results/densities/<run>, plots/<run>, animation/<run>, run = date and time of the submission ($PETRIDISH_RUN)
# Seeds that do not reach t* within max_wait_days (petridishsimauto.py) end normally with status "not_crossed".
#
#SBATCH --job-name=petrisimjob
#SBATCH --array=0-99%50
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu256
#SBATCH --time=5:00:00
#SBATCH --mem=3G
#SBATCH --mail-type=BEGIN,TIME_LIMIT_80,END
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
# all outputs (npz, plots, vtp) go there, not into the repository (amfio.petridishOutputDirs), absolute path since we cd below
export PETRIDISH_OUT="$PWD/petridish_results"
# all tasks of the array write into the same run folders; submitted without submit_petridish.sh, the array job id is the run name
export PETRIDISH_RUN="${PETRIDISH_RUN:-job_$SLURM_ARRAY_JOB_ID}"
cd CPlantBox/experimental/myc/
export MPLBACKEND=Agg

echo "Task $SLURM_ARRAY_TASK_ID: seed $SEED, run $PETRIDISH_RUN, parameters ${PETRIDISH_ARGS:-(xml)} on $(hostname), start $(date)"
# runtime per phase and peak memory are stored by petridishsimauto.py itself (see petridishtiming.py),
# since slurm accounting is not available on the cluster (and /usr/bin/time is not installed on the nodes)
# PETRIDISH_ARGS: hyphal parameters name=value of submit_petridish.sh (empty: values of the xml), not quoted, so each is one argument
python3 petridishsimauto.py "$SEED" $PETRIDISH_ARGS
echo "Seed $SEED done $(date)"
 