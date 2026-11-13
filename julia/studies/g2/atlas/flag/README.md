# Local G2 flag cell atlas

`build_local_flag_cell_atlas.jl` draws the local cells of the coarsest G2/B
flag fans over QQ and over GF(2). It writes `g2_tree_atlas.tex`,
`local_flag_cell_atlas.csv`, and `local_flag_intersection_certificate.tsv`.

| File | Contents |
|---|---|
| `local_flag_geometry.jl` | Exact geometry of a flag fan over its support cells and the ray labels. It includes the shared Grassmannian definitions in [`../grassmannian/atlas_common.jl`](../grassmannian/atlas_common.jl). |
| `local_flag_fields.jl` | The second field: the GF(2) support cells in the ray numbering of the atlas, matched to the QQ support orbits. |
| `local_flag_degenerations.jl` | The stage-6 censuses of both fields, the page order, the polytope notation, and the polymake drawings. |

The tables of main.tex number the support orbits by the labels of
`local_flag_cell_atlas.csv`. Two files hold this numbering for the scripts of
`to_tex/` and for `../grassmannian/export_maximal_cell_tree_metrics.jl`;
the build does not use them.

| File | Contents |
|---|---|
| `cell_orbit_labels.jl` | The stage-6 census rows of both fields, and the label, the census orbits, the orbit size and the support cell of every support orbit, in atlas order. It reads CSV files only. |
| `cell_orbit_trees.jl` | The trees of the lineality generators, the rays, the edges and the maximal cells of the coarsest Grassmannian fans of both fields, which `to_tex/render_maximal_cell_tree.jl` draws, and the tree of the support cell of every numbered orbit. |

## The two fields

Rays `R_1`, ..., `R_51` are those of the QQ census. The three rays that occur
only over GF(2) are `R_52`, `R_53`, `R_54`; they form the ray orbit `O12` and
are the centres of the three prevariety tetrahedra. A support cell present in
both fields has the same rays in both, and the build checks that it carries
the same flag cells in both. There are 28 support orbits over both fields, 3
only over QQ, and 4 only over GF(2).

| Field | Flag rays | Maximal flag cells | Cell orbits | Prime, all toric |
|---|---:|---:|---:|---:|
| QQ | 141 | 2714 | 270 | 81 |
| GF(2) | 156 | 2792 | 281 | 83 |

Run the census stages through `../../run_census.jl` (QQ) and
`../../run_census_gf2.jl` (GF(2)). The renderer reads the coarsest fans
`02_coarsest_fan.mrdi` of `data/g2/flag/` and `data/g2/flag_gf2/` and the
Grassmannian data of both fields through `TropicalGeometryG2.data_path`.

## Pages

Page 2 shows the twelve ray orbits, the one over GF(2) only last. Page 3
draws the three- and four-dimensional smallest pieces with polymake's TikZ
backend, the latter as Schlegel diagrams, and defines the polytope notation.

The support orbits are labelled 1 to 35 in page order. Those with a toric
degeneration come first: over both fields (the Chevalley simplex, the other
simplex, then by increasing smallest piece and description), then only over
QQ, then only over GF(2). The others follow: present over both fields, only
over QQ, only over GF(2), the prime ones first. Each page says over which
fields its orbit occurs.

Each support page states the status of its initial ideal and, when it is prime
toric, the type, f-vector, and lattice-point count of its Newton--Okounkov
polytope `P_O`. The polytope of a prime flag cell is a Minkowski sum: the type
of the summand from the point coordinates, plus the summand from the Pluecker
coordinates. A status or a summand line names its field only when the two
fields of the page disagree.

## Certificates of the build

- The polytope types are recomputed with `RecursivePolytopeTypes` and must
  equal the census strings of both fields.
- The Pluecker summand of every prime flag cell is unimodularly equivalent to
  `P_O` of its support orbit, in each field.
- A support orbit that is toric over both fields has unimodularly equivalent
  polytopes in both.
- The polytope that the census marks as the Chevalley polytope of `G2/P2` is
  unimodularly equivalent to the lecture hall simplex of `(1, 3, 2, 3, 1)`.
- Every drawn smallest piece is combinatorially isomorphic to the census piece.

## Output

`local_flag_cell_atlas.csv` has one row per drawn flag cell. `label` is the
atlas label, `fields` the fields of the support orbit, `qq_orbit` and
`gf2_orbit` its census numbers, and `qq_flag_rays` and `gf2_flag_rays` the rays
of the cell in each flag fan.

Run the renderer from `julia/`:

```sh
julia --project=. studies/g2/atlas/flag/build_local_flag_cell_atlas.jl --no-pdf
```

It takes about two minutes and under three gigabytes of memory. When built
without `--no-pdf`, the PDF is written to `../output/pdf/g2_tree_atlas.pdf`.
