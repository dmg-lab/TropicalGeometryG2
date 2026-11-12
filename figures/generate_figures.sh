#!/usr/bin/env bash
###############################################################################
# Writes every figure and table of main.tex that a script draws, with the
# commands of the % COMPUTATION comments of main.tex.
#
#   bash figures/generate_figures.sh [DIR]
#
# The files go to DIR/figures/ and DIR/tables/, by default those of the
# repository. Not drawn by a script: the prime of gamma_1' in
# figures/square_cell.tex, added by hand, and the hand-drawn figures
# (taxon_hexagon, root_lattice_V2, algebraic_hexagon, triangle_incidence*,
# triangle_projection*). The curve portraits also write their previews next to
# the curves, in julia/curves_g10/overnight_20261028/.
###############################################################################
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$ROOT}"
mkdir -p "$OUT/figures" "$OUT/tables"
OUT="$(cd "$OUT" && pwd)"
F="$OUT/figures"
T="$OUT/tables"
cd "$ROOT"

J=(julia --project=julia)
TREE=("${J[@]}" julia/to_tex/render_maximal_cell_tree.jl --color --scale=0.5 --taxon-scale=0.6)
CURVES=julia/curves_g10/overnight_20261028

# Section 2: the G2-trees of a ray, an edge, a cell and the lineality space.
"${TREE[@]}" --ray=8 --bounded-length=0.8 --output="$F/ray_tree.tex"
"${TREE[@]}" --edge=7,41 --bounded-length=0.8 --output="$F/edge_tree.tex"
"${TREE[@]}" --rays=7,18,41 --bounded-length=0.8 --output="$F/triangle_tree.tex"
"${TREE[@]}" --lineality=1 --bounded-length=0.8 --output="$F/lineality_star_1.tex"
"${TREE[@]}" --lineality=2 --bounded-length=0.8 --output="$F/lineality_star_2.tex"
"${TREE[@]}" --lineality=3 --bounded-length=0.8 --output="$F/lineality_star_3.tex"
"${TREE[@]}" --ray=53 --field=GF2 --bounded-length=1 --output="$F/ray_GF2.tex"

# Section 3: the cubics of the tropical basis and the trees of their monomials.
sh julia/to_tex/cubic_tree_table.sh --output="$T/proof_int_tropical_basis.tex"

# Section 5: the G2-forest.
"${J[@]}" julia/to_tex/build_styled_ray_orbit_tables.jl --color --label=false --scale=0.95 --taxon-scale=0.75 \
    --bounded-length=1.1 --output="$F/ray_orbits_styled_catalog.tex"
"${J[@]}" julia/to_tex/toric_degeneration_tree_table.jl --color --scale=0.88 --taxon-scale=0.75 --bounded-length=1 \
    --output="$F/toric_degeneration_trees.tex"
"${J[@]}" julia/to_tex/toric_degeneration_tree_table.jl --non-toric --color --scale=0.88 --taxon-scale=0.75 \
    --bounded-length=1 --output="$F/non_toric_cell_trees.tex"
"${J[@]}" julia/to_tex/square_cell_figure.jl --output="$F/square_cell.tex"
"${J[@]}" julia/to_tex/cell_ray_orbit_table.jl --output="$T/cell_ray_orbits.tex"

# Section 6: the Newton-Okounkov bodies.
"${J[@]}" julia/to_tex/toric_degeneration_table.jl --output="$T/toric_degeneration_table.tex"

# Section 7: the G2-orchard. The prime edge trees check their panels against
# figures/toric_degeneration_trees.tex of the repository.
"${J[@]}" julia/to_tex/toric_degeneration_prime_edge_trees.jl --fragment \
    --output="$F/toric_degeneration_prime_edge_trees.tex"
"${J[@]}" julia/to_tex/decorated_trees_figure.jl --scale=0.7 --taxon-scale=0.55 --output="$F/decorated_trees.tex"
"${J[@]}" julia/to_tex/triangle_cell_figure.jl --scale=0.5 --taxon-scale=0.5 --ray-scale=0.8 --name-below \
    --no-cell-letters --cell-overhang --height=4.3 --output="$F/triangle_cell.tex"
"${J[@]}" julia/to_tex/triangle_cell_figure.jl --orbit=5 --scale=0.5 --taxon-scale=0.5 --ray-scale=0.8 --name-below \
    --no-cell-letters --cell-overhang --height=4.3 --output="$F/triangle_cell_5.tex"
"${J[@]}" julia/to_tex/local_flag_complex_table.jl --output="$T/local_flag_complexes.tex"

# Section 9: the tropical hexagon.
"${J[@]}" julia/studies/hexagon/figures/hexagon_tree_ring.jl --output="$F/hexagon_tree_ring.tex"

# Section 10: the two curves of genus 10.
"${J[@]}" julia/studies/g2/curves/curve_portraits.jl "$CURVES/wpart_2722/decorated_graph_1_reduced.csv" \
    --output-prefix="$CURVES/wpart_2722/wpart_2722_row" --outer-face=1 --seed=0 \
    --side-by-side --frames --trivalent-nodes --fragments="$F/genus_10_curve"
D="$CURVES/experiment_100128"
"${J[@]}" julia/studies/g2/curves/curve_portraits.jl "$D/curve_reduced.csv" \
    --output-prefix="$D/experiment_100128" --points="$D/point_places.csv" --side-by-side --chain-power=0.6 \
    --crossing=2 --outer-face=1 --seed=0 --trivalent-nodes --connected-regions --region-drawing="$D/region_drawing.csv" \
    --tree-scale=4 --overlapping-trees --fragments="$F/genus_10_curve_100128"
