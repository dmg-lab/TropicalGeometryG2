#!/bin/bash
#SBATCH --job-name=mukai-curves
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=work/curves/%x-%j.out

# Submit from mukai+project/julia with the arguments of run_experiments.jl, e.g.
# 10 genus-10 experiments with translations from the shift basis, then 10 with
# uniformly random translations:
#   mkdir -p work/curves
#   sbatch mukai_curves/run_experiments.slurm 10 10 certified results/curves/genus_10/basis_shifts --shift-basis
#   sbatch mukai_curves/run_experiments.slurm 10 10 certified results/curves/genus_10/uniform_shifts
# The threads of Julia run the feasibility filter, one linear-space cell each.
# Use the solver certified rather than hybrid: the SoPlex fallback of hybrid
# crashed with a segmentation fault in both batches of the first runs here.
set -euo pipefail

export PATH="$HOME/.juliaup/bin:$PATH"
cd "$SLURM_SUBMIT_DIR"
julia --project=. --threads=$SLURM_CPUS_PER_TASK mukai_curves/run_experiments.jl "$@"
