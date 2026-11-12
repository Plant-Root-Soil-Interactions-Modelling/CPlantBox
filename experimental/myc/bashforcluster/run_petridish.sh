#!/bin/bash
#
# Job array over a task list: one line per array task, "<run> <seed> [name=value ...]" (run name, seed, hyphal parameters).
# The task list ($PETRIDISH_TASKS) and --array=0-(number of lines - 1)%<throttle> are set by submit_petridish_compare.sh
# (or submit_petridish.sh for a single run), all runs of a submission are in this one array.
# Each task is a normal 1-CPU job, $SLURM_ARRAY_TASK_ID (0, 1, 2, ...) picks the line of the task list.
#
# Submit from the home directory:
#   bash submit_petridish.sh [--build] [name=value ...]                          # one run
#   bash submit_petridish_compare.sh [--build] xml ana=0 distTH=0.02 ...        # several runs and their comparison
# Outputs go into petridish_results/densities/<run>, plots/<run>, animation/<run>
# Seeds that do not reach t* within max_wait_days (petridishsimauto.py) end normally with status "not_crossed".
#
#SBATCH --job-name=petrisim
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

if [ -z "$PETRIDISH_TASKS" ] || [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "no task list or no array task: submit with submit_petridish.sh or submit_petridish_compare.sh" >&2
    exit 1
fi
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$PETRIDISH_TASKS")
if [ -z "$line" ]; then
    echo "Array index $SLURM_ARRAY_TASK_ID out of range of the task list $PETRIDISH_TASKS" >&2
    exit 1
fi
read -r RUN SEED ARGS <<< "$line"

source cpbenv/bin/activate
# all outputs (npz, plots, vtp) go there, not into the repository (amfio.petridishOutputDirs), absolute path since we cd below
export PETRIDISH_OUT="$PWD/petridish_results"
export PETRIDISH_RUN="$RUN"
cd CPlantBox/experimental/myc/
export MPLBACKEND=Agg

echo "Task $SLURM_ARRAY_TASK_ID: seed $SEED, run $RUN, parameters ${ARGS:-(xml)} on $(hostname), start $(date)"
# runtime per phase and peak memory are stored by petridishsimauto.py itself (see petridishtiming.py),
# since slurm accounting is not available on the cluster (and /usr/bin/time is not installed on the nodes)
# ARGS: hyphal parameters name=value (empty: values of the xml), not quoted, so each is one argument
python3 petridishsimauto.py "$SEED" $ARGS
echo "Seed $SEED done $(date)"
