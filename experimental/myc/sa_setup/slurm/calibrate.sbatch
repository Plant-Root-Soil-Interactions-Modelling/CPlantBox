#!/bin/bash
#
# Calibration (calibrate.py) as ONE job: the seeds of each evaluation run in parallel on the cores of the job.
# Submit from the home directory, arguments are passed on to calibrate.py:
#   sbatch calibrate.sbatch
#   sbatch calibrate.sbatch --seeds 16 --max-nfev 10
# Aborted (timeout)? Just resubmit: finished runs are reused.
# Output: slurmOut/sa_cal-<job>.out, results in sa_setup/calibration/
#
#SBATCH --job-name=sa_cal
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128      # Jacobian: 16 difference points × 16 seeds = 256 runs -> 3 rounds of ~0.75 h
#SBATCH --partition=cpu256
#SBATCH --time=72:00:00          # 120 h after t*: ~1.6 h per run, per iteration ~3 rounds ≈ 5 h; resubmit to continue
#SBATCH --mem=200G               # smoke test: max. 0.85 GB per run, 128 at once
#SBATCH --output=./slurmOut/sa_cal-%j.out
#SBATCH --error=./slurmErr/sa_cal-%j.err

set -euo pipefail
source "$HOME/cpbenv/bin/activate"
cd "${SA_DIR:-$HOME/CPlantBox/experimental/myc/sa_setup}"
export MPLBACKEND=Agg

echo "python3 calibrate.py $* in $PWD on $(hostname), $(date)"
python3 calibrate.py "$@"
