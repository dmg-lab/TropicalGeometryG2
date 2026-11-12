#!/bin/sh
###############################################################################
# The table of the three cubics of the tropical basis of the G2-Grassmannian,
# with the simple G2-tree of a ray over GF(2) under each monomial: runs
# cubic_tree_table.jl, which certifies the trees and draws the table.
#
#   sh julia/to_tex/cubic_tree_table.sh [--edge-length=LENGTH] [--standalone] [--output=FILE]
###############################################################################
set -eu

to_tex=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec julia --project="$to_tex/.." "$to_tex/cubic_tree_table.jl" "$@"
