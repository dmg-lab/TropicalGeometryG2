# Canonical data

This folder is written by the scripts below, or unpacked from the Zenodo tarball
([10.5281/zenodo.24104003](https://doi.org/10.5281/zenodo.24104003)), which is a copy of it as the census writes
it, with curated aliases; both give the same tree. Git keeps the inputs, which no script below rewrites: this
README, the notes [`g2/SUMMARY.md`](g2/SUMMARY.md) (the seven `08_summary.toml` in one table), the tree metrics
`g2/maximal_cell_tree_metrics.csv`, `chevalley/`, the shift basis `mukai_curves/genus_10_shift_basis.mrdi`, and
`fan_uniform_linear_space_16_21.mrdi`.


- `make_paper_data.jl`: the single driver, which runs the census in order and writes every other file below.
  The stage `zenodo` names the census orbits as the paper does from the label tables of `studies/g2/atlas/`,
  which git also keeps.
- `studies/g2/run_census.jl`: the targets `dressian`, `flag_dressian`, `base` (`g2/grassmannian/`) and `flag`
  over QQ, and the stage `zenodo` of the target `zenodo`.
- `studies/g2/run_census_gf2.jl`: the targets `base` and `flag` over GF(2).
- `studies/g2/run_census_slong.jl`: `g2/grassmannian_Slong/` and `g2/grassmannian_Slong/short_long/`.
- The generators named below, for the files outside the numbered stages.

Study scripts resolve these paths with `data_path` of the package in `src/`; data are not duplicated into study
directories.

## `g2/<target>/`

One directory per census target: the G2-Dressian `dressian` and the G2 flag Dressian `flag_dressian`,
`grassmannian` (G2/P2) and `flag` (G2/B) over QQ, `grassmannian_gf2` and `flag_gf2` over GF(2), and
`grassmannian_Slong`, G2/P2 in the coordinates of S_long. The Dressians are the prevarieties of the generators,
which have coefficients 1 and -1, so they are the same over every field and stored only here. A file name starts
with the number of the stage that writes it, and `provenance/<stage>.toml` records each stage.

| Files | Stage | Content |
|---|---|---|
| `00_*` | `tropicalize` | The ideal; for the flag targets also its purity certificate and the verified prime initial ideal. No Dressian. |
| `01_*` | `tropicalize` | The tropical variety `01_raw_fan.mrdi` (with the complex `01_variety_complex.mrdi` for `grassmannian`), or the prevariety `01_prevariety_fan.mrdi` of a Dressian. Over GF(2) also the tropical basis and its certificates. |
| `02_*` | `coarsen` | The coarsest fan, its cells, and the report `02_coarsest_report.md`. |
| `03_*` | `normalize` | The coarsest fan, with each ray normalized modulo the lineality space, and the matrix of these rays. |
| `04_*` | `trees` | The canonical weight and the tree of each maximal cell. Not for `grassmannian_Slong`. |
| `05_*` | `symmetry` | The action of the Weyl group D12 on the rays and the cells, and the orbits. |
| `06_*` | `degenerations` | One row per orbit of maximal cells: the initial ideal at a canonical interior weight and its invariants. In `06_canonical_interior_degenerations_objects/` the ideals and, for the toric orbits, the Newton–Okounkov bodies (`.poly`, Oscar mrdi). No Dressian. |
| `07_*` | `projections` | Flag targets: exact rational witnesses of the projection onto the maximal cells of the Grassmannian, one per orbit (270 over QQ), and the target cell and Weyl-group transport of each maximal cell (2714 over QQ). `grassmannian_Slong`: `07_projections.{csv,toml}`. The checks `studies/g2/supplemental/verify_gf2_census.jl` of the census over GF(2) write `grassmannian_gf2/07_verification.toml`, which is not part of the data. |
| `08_summary.toml` | `summarize` | The statements of the paper about the target: dimension, f-vector and orbit f-vector modulo lineality, the combinatorial types of the maximal cells, the reduced integral homology, the size of the stage-1 fan and the counts of stage 6. Over GF(2) also the comparison of the maximal cells and toric degenerations with those over QQ. |

Over GF(2), `census.md` is rewritten after every stage.

## Other files

- The stage `zenodo` (`studies/g2/census/zenodo_census.jl`) writes the files of the tarball that no other stage
  writes: `g2/pointspace{,_gf2}/00_ideal.mrdi`, `g2/hexagon/tropicalized_hexagon.mrdi`, the orbit correspondences
  `g2/{grassmannian,flag}_orbits.csv`, `g2/dressian/dressian_orbits.csv` and
  `g2/flag_dressian/flag_dressian_orbits.csv`, and the two curves of genus ten
  `mukai_curves/genus_10_method_{1,2}.mrdi`: that of Method 1 with its pieces `genus_10_method_1_pieces.mrdi`,
  computed from its shift vector `curves_g10/overnight_20261028/wpart_2722/translation_1.mrdi`, and that of
  Method 2 from its run `curves_g10/overnight_20261028/experiment_100128/curve.mrdi`; and the tropical K3 surface
  of Method 1, `mukai_curves/genus_10_method_1_k3.mrdi`, polygon by polygon (`studies/g2/curves/k3_surface_polygons.jl`,
  7 min), whose bounded part `paper_results.jl` checks to collapse onto a 2-sphere (`studies/g2/curves/k3_fins.jl`).
  Provenance in `g2/provenance/zenodo.toml`.
- `g2/maximal_cell_tree_metrics.csv`: the 35 orbits of maximal cells of `g2/grassmannian/` and
  `g2/grassmannian_gf2/` together, numbered as in the two tables of maximal trees of the paper, with their census
  orbits, the rays of the cell drawn there, and the 21 distances `d_ij` of its tree. Kept in git, and written again
  with the label tables by the scripts of `studies/g2/atlas/`, last
  `studies/g2/atlas/grassmannian/export_maximal_cell_tree_metrics.jl`. Those scripts also write
  `g2/reference/fine_groebner/normalized_fine_{fan,rays}.mrdi`, the fine fan of the Grassmannian whose rays name the
  51 rays (`studies/g2/reference/fine_groebner/build_balanced_fan_and_atlas.jl --no-pdf`).
- `chevalley/`: Chevalley polytopes of Spacek and Wang, as native polymake files, written by
  `generate_chevalley_polytopes.jl` and, for the 18 of G2, `g2_chevalley_study.jl`, both in
  `studies/g2/degenerations_grassmannian/Spacek_Wang/`.
- `mukai_curves/`: `genus_10_shift_basis.mrdi`, written by `studies/g2/curves/shift_basis.jl`, the two
  curves of genus ten and the K3 surface of Method 1.
- `fan_uniform_linear_space_16_21.mrdi`: the Bergman fan B_4 of the uniform matroid of rank 17 on 21 elements (16 unit
  vectors per maximal cone, min convention) of the curves of genus ten, written by
  `studies/g2/curves/uniform_bergman_fan.jl`.

