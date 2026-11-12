# Historical fine Gröbner census (noncanonical)

This directory preserves the scientifically useful report, CSV classification,
and TeX atlases from the 251-cell fine G2 fan. The exact fan, rays, weights,
incidence, initial-ideal census, and converted Newton–Okounkov body data are centralized under
`../../../../data/g2/reference/fine_groebner/`.

These files are reference material only. They are not imported by the current
census, they do not define the symmetry group of the coarsest fan, and they are
never stage-6 input. The maintained stage 6 instead computes one exact
canonical-interior weight for each of the 31 coarsest-cell orbits.

The restored historical renderers are:

- `find_better_gens.jl`, the exact tree/gauge helper used by the old atlas;
- `build_balanced_fan_and_atlas.jl`, which writes
  `balanced_phylogenetic_trees.tex`;
- `build_balanced_triangle_atlas.jl`, which writes
  `balanced_phylogenetic_triangles.tex`.

Their machine inputs and derived matrices are resolved through the central
`data/g2/reference/fine_groebner/` directory. These scripts are deliberately
kept here, outside the maintained census driver.
