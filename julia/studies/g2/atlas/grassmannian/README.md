# G2 Grassmannian tree and orbit atlases

This directory contains the QQ/GF(2) orbit atlas and the W(G2) orbit CSVs
that it, the paper ray catalog, and the [flag atlas](../flag/) read.
The shared G2-tree renderer is in [`to_tex/`](../../../../to_tex/).

Run Julia commands from `julia/` with `--project=.`.
Inputs come from `data/g2/grassmannian/`,
`data/g2/grassmannian_gf2/`, and the fine fan
`data/g2/reference/fine_groebner/normalized_fine_{fan,rays}.mrdi`,
written by
[`build_balanced_fan_and_atlas.jl --no-pdf`](../../reference/fine_groebner/build_balanced_fan_and_atlas.jl).
TeX sources stay in this directory; PDF outputs go to `../output/pdf/`. The
field atlas uses `work/tex/` for compilation files.

## Scripts

| Script | Purpose |
|---|---|
| [`build_field_comparison_atlas.jl`](build_field_comparison_atlas.jl) | Assemble the three-page QQ/GF(2) orbit atlas. |
| [`export_field_atlas_assets.jl`](export_field_atlas_assets.jl) | Export its exact tree fragments, memberships, metrics, and source hashes. |
| [`to_tex/render_maximal_cell_tree.jl`](../../../../to_tex/render_maximal_cell_tree.jl) | Render the tree of a cone over QQ or GF(2), or a supplied G2 tree through the Julia API. |
| [`to_tex/build_styled_ray_orbit_tables.jl`](../../../../to_tex/build_styled_ray_orbit_tables.jl) | Build the paper ray catalog with the shared renderer. |
| [`export_maximal_cell_tree_metrics.jl`](export_maximal_cell_tree_metrics.jl) | Write the distances `d_ij` of the trees of the 35 maximal-cell orbits of both fields, numbered as in the paper tables, to `data/g2/maximal_cell_tree_metrics.csv`. It includes only [`../flag/cell_orbit_trees.jl`](../flag/cell_orbit_trees.jl), no renderer of `to_tex`. |
| [`compute_label_orbits.jl`](compute_label_orbits.jl) | Compute the W(G2) label action and write the orbit CSVs. |
| [`atlas_common.jl`](atlas_common.jl) | Shared definitions only: orbit CSV readers, rooted tree signatures, the 248 coarsest support cells, and the censuses of both fields in atlas rays. |
| [`ray_orbit_labels.csv`](ray_orbit_labels.csv) | Kept by hand: the label of the paper and the catalog position of each ray orbit, read by the paper ray catalog and the stage zenodo of the census. |

### Code layout

The include graph of the atlas scripts is:

```text
reference/fine_groebner/build_balanced_triangle_atlas.jl   fine-fan trees (via find_better_gens.jl)
└── atlas_common.jl
    ├── compute_label_orbits.jl
    ├── ../flag/cell_orbit_trees.jl  (+ ../flag/cell_orbit_labels.jl)
    │   ├── ../../../../to_tex/render_maximal_cell_tree.jl
    │   │   ├── ../../../../to_tex/build_styled_ray_orbit_tables.jl
    │   │   └── export_field_atlas_assets.jl  <- build_field_comparison_atlas.jl
    │   └── export_maximal_cell_tree_metrics.jl
    └── ../flag/local_flag_geometry.jl  <- build_local_flag_cell_atlas.jl,
                                           with ../flag/local_flag_degenerations.jl
```


[`atlas_common.jl`](atlas_common.jl) is the only atlas file that includes the
historical [fine Groebner reference](../../reference/fine_groebner/) scripts.
Those scripts provide the Coxeter balancing, the balanced fine-fan rays, and
the fine-triangle classification; `TreeMetric` and the tree reconstruction
are those of the package, in [`src/trees.jl`](../../../../src/trees.jl).


## QQ/GF(2) orbit atlas

The main output is `../output/pdf/g2_grassmannian_qq_gf2_orbit_atlas.pdf`.
Its three A3 landscape pages contain:

1. Twelve ray classes in a 3×4 table, with GF(2)-specific `R12` highlighted.
2. All 28 common maximal-cell orbits on one page, with paired census orbit
   IDs, sizes, ray-class multisets, and exact tree representatives.
3. A colored tetrahedron comparison of the three QQ-only boundary-face
   orbits and four GF(2)-only center-to-edge triangle orbits.

| Fan | Rays | Ray orbits | Maximal cells | Maximal-cell orbits |
|---|---:|---:|---:|---:|
| QQ | 51 | 11 | 248 | 31 |
| GF(2) | 54 | 12 | 254 | 32 |
| Common | 51 | 11 | 236 | 28 |

To rebuild the document using the saved, verified tree assets:

```sh
julia --project=. studies/g2/atlas/grassmannian/build_field_comparison_atlas.jl --reuse-assets
```

Omit `--reuse-assets` for a full build, which re-exports the assets only when
their inputs, the exporter, or the renderer changed. Use `--force-assets`
(without `--reuse-assets`) to regenerate every tree fragment. `--no-pdf`
writes the ray page and standalone wrapper without compiling. PDF
compilation requires `pdflatex`.

The standalone document it writes is `field_comparison_atlas.tex`.
It inputs `field_ray_classes.tex`, `field_common_cells.tex`, and
`field_tetrahedron_comparison.tex`, using `assets/tree_atlas_style.tex`.

The exporter reads stage-3 normalized fans and rays, stage-4 normalized cell
trees, and stage-5 symmetry data. It checks the field correspondence
against [`09_qq_comparison.json`](../../../../data/g2/grassmannian_gf2/09_qq_comparison.json).
Its output directory, [`field_comparison_assets/`](field_comparison_assets/),
contains `rays/`, `common/`, `exceptional/`, and `tetra/` tree fragments, plus
[`atlas_data.toml`](field_comparison_assets/atlas_data.toml) with exact metrics,
memberships, index correspondences, and input hashes. The exporter itself
accepts `--force` to regenerate all assets.

`R1`–`R11` follow the QQ census ray-orbit order. `R12` consists of GF(2) rays
`3,36,39`, the centers of the three Dressian tetrahedra; the tree renderer
numbers them `52,53,54`. In each tetrahedron,
the four QQ boundary triangles are replaced by six GF(2) triangles joining
the center to the six edges.

Ray representatives use primitive integral Coxeter-balanced metrics with
minimum pendant length zero. Cell representatives use the equal positive
sum of their normalized ray generators. Drawn edge lengths are schematic;
numeric edge labels give the exact metric.

## Render one tree

The shared tree renderer is
[`to_tex/render_maximal_cell_tree.jl`](../../../../to_tex/render_maximal_cell_tree.jl).
It draws the tree of a generator of the lineality space, of a ray, of an
edge, or of a maximal cell, over QQ or over GF(2). Its options, its numbering of rays and
cells, and its Julia API are described in
[`to_tex/README.md`](../../../../to_tex/README.md).

## Paper ray catalog

The catalog contains the twelve ray orbits, the eleven over QQ and the one
that occurs only over GF(2), in the article order and Greek labels of
[`ray_orbit_labels.csv`](ray_orbit_labels.csv), described in
[`ray_orbit_article_labels.md`](ray_orbit_article_labels.md). These labels
differ from the field atlas's `R1`–`R12` notation. Its builder is
[`to_tex/build_styled_ray_orbit_tables.jl`](../../../../to_tex/build_styled_ray_orbit_tables.jl),
described in [`to_tex/README.md`](../../../../to_tex/README.md).

## W(G2) orbit CSVs

[`compute_label_orbits.jl`](compute_label_orbits.jl) tests every permutation
of taxa `1,...,6` on the exact balanced ray trees of the historical
251-triangle fan and writes
[`label_permutation_group.csv`](label_permutation_group.csv),
[`ray_label_orbits.csv`](ray_label_orbits.csv), and
[`triangle_label_orbits.csv`](triangle_label_orbits.csv). The paper catalog
and the flag atlas read these files. The maintained census symmetry stage
does not read these CSVs.
