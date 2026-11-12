#!/bin/bash
# Slurm array for sample_genus_10.jl on hydra, from julia/:
#   sbatch --array=0-113 studies/g2/curves/sample_genus_10.slurm
# Task t of T runs the experiments BASE + t + i T, i = 0, 1, ..., in Julia
# processes of at most BATCH experiments each until the time limit; a process
# ends early once its peak RSS passes 6000 MiB (polymake memory is not
# returned to the system), and the experiments it leaves out are not run.
# Fresh processes peaked below 9 GB, hence 10G per task. Results go to OUTPUT, by default
# results/curves/genus_10_sampling; the progress bars are dropped from the logs.
#SBATCH --job-name=genus10-sampling
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=14:00:00
#SBATCH --output=work/genus_10_sampling/logs/task-%a.out

set -uo pipefail
export PATH="$HOME/.juliaup/bin:$PATH"
cd "$SLURM_SUBMIT_DIR"
BASE=${BASE:-100000}
BATCH=${BATCH:-10}
OUTPUT=${OUTPUT:-results/curves/genus_10_sampling}
tasks=$SLURM_ARRAY_TASK_COUNT
task=$SLURM_ARRAY_TASK_ID
failures=0
for ((batch = 0; ; batch++)); do
    first=$((BASE + task + batch * BATCH * tasks))
    julia --project=. -t 2 --heap-size-hint=2500M studies/g2/curves/sample_genus_10.jl --first=$first \
        --experiments=$BATCH --stride=$tasks --task=$task --output=$OUTPUT 2>&1 |
        stdbuf -oL tr '\r' '\n' | grep --line-buffered -v '^Progress:'
    status=${PIPESTATUS[0]}
    echo "batch $batch from experiment $first: exit status $status"
    if ((status != 0)); then
        ((++failures >= 3)) && { echo "three failed batches, stopping"; exit 1; }
    else
        failures=0
    fi
done
