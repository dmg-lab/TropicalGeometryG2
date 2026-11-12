#!/bin/bash
# Variants of benchmark_oscar_prevariety.jl in parallel on one hydra node, one
# single-threaded process each, from julia/:
#   sbatch studies/g2/supplemental/benchmark_oscar_prevariety.slurm
# The cpu partition was down on 2026-11-06, so this asks for the gpu partition
# without a GPU. Logs: work/oscar_prevariety/<job>/<variant>.log
#SBATCH --job-name=oscar-prevariety
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=47:00:00
#SBATCH --output=work/oscar_prevariety/slurm-%j.out

set -uo pipefail
export PATH="$HOME/.juliaup/bin:$PATH"
cd "$SLURM_SUBMIT_DIR"
logs=work/oscar_prevariety/$SLURM_JOB_ID
mkdir -p "$logs"
lscpu | grep -E "Model name|^CPU\(s\)|MHz"
free -g | head -2
git rev-parse --short HEAD 2>/dev/null
# Precompile once, before the parallel processes need it.
julia --project=. -e 'using Oscar' 2>&1 | tail -3

variants=(
    "--case=flag78 --order=input --mode=oscar"
    "--case=flag78 --order=input --mode=steps"
    "--case=flag78 --order=greedy --mode=steps"
    "--case=flag72 --order=input --mode=steps"
    "--case=flag72 --order=greedy --mode=steps"
    "--case=dressian42 --order=input --mode=oscar --hull=ppl"
    "--case=dressian42 --order=input --mode=oscar --hull=cdd"
    "--case=dressian42 --order=input --mode=oscar --hull=lrs"
    "--case=dressian42 --order=input --mode=oscar --hull=beneath_beyond"
    "--case=dressian42 --order=greedy --mode=steps"
)
for variant in "${variants[@]}"; do
    name=$(echo "$variant" | sed 's/--[a-z]*=//g; s/ /_/g')
    # The script prints its own times and peak memory; the nodes have no /usr/bin/time.
    julia --project=. -t 1 studies/g2/supplemental/benchmark_oscar_prevariety.jl $variant \
        > "$logs/$name.log" 2>&1 &
done
wait
grep -h "RESULT\|MISMATCH\|ERROR" "$logs"/*.log
