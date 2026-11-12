# Historical fine Gröbner census (noncanonical)

This directory holds the tree renderers and TeX atlases of the 251-cell fine
G2 fan, the fan `julia/data/g2/grassmannian/01_raw_fan.mrdi` of the census:
`balanced_phylogenetic_trees.tex`, the trees of its 51 rays, and
`balanced_phylogenetic_triangles.tex`, the trees of its 251 maximal triangles.

These files are reference material only. They are not imported by the current
census, they do not define the symmetry group of the coarsest fan, and they are
never stage-6 input. The maintained stage 6 instead computes one exact
canonical-interior weight for each of the 31 coarsest-cell orbits.

The renderers are:

- `find_better_gens.jl`, the exact tree/gauge helper (`TreeMetric`, tree
  reconstruction, Coxeter balancing);
- `build_balanced_fan_and_atlas.jl`, which writes
  `balanced_phylogenetic_trees.tex`;
- `build_balanced_triangle_atlas.jl`, the trees of the triangles, which the
  atlas scripts of `../../atlas/` include; it wrote
  `balanced_phylogenetic_triangles.tex`.

`build_balanced_fan_and_atlas.jl` reads the fan of the census and writes it
with the balanced ray generators to
`julia/data/g2/reference/fine_groebner/normalized_fine_{fan,rays}.mrdi`, which
the atlas scripts read; `--no-pdf` skips compiling its atlas. These scripts
are deliberately kept here, outside the maintained census driver.
