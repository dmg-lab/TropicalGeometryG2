# G2 census

Every target keeps what it produces in one directory `../../data/g2/<target>/`:
the numbered stage artifacts, the coarsening report `02_coarsest_report.md`,
and the stage provenance `provenance/<stage>.toml`. The directories are
`grassmannian`, `flag`, `grassmannian_gf2`, `flag_gf2`, and `dressian` and
`flag_dressian` for [the two Dressians](#the-two-dressians), which are the
same over every field and are stored once.

Run the commands of this file from `julia/`.

The independent characteristic-two census uses
[`run_census_gf2.jl`](run_census_gf2.jl) and stores all canonical artifacts,
certificates, reports, and provenance in
[`../../data/g2/grassmannian_gf2/`](../../data/g2/grassmannian_gf2/).
The polynomial coefficients are in GF(2), with the trivial valuation and the
min convention. Its tropicalization starts from the G2-Dressian in
`../../data/g2/dressian/` after verifying that reduction modulo 2 preserves
every defining polynomial support. It certifies whole cones with GF(2) Groebner bases and
torus saturation, then runs the usual census stages with its own cell and
orbit counts.

```bash
# Status only.
julia --project=. studies/g2/run_census_gf2.jl
# Complete census; use --resume after an interrupted computation.
julia --project=. studies/g2/run_census_gf2.jl --target=base --from=tropicalize --through=degenerations --bits=0
# Independently verify the saved certificates and census artifacts.
julia --project=. studies/g2/supplemental/verify_gf2_census.jl
```

The same script runs the G2/B flag over GF(2) with `--target=flag`; its
artifacts are in [`../../data/g2/flag_gf2/`](../../data/g2/flag_gf2/).

```bash
julia --project=. studies/g2/run_census_gf2.jl --target=flag --stage=tropicalize --workers=2 --threads=8 --bits=0
julia --project=. studies/g2/run_census_gf2.jl --target=flag --from=coarsen --through=projections --resume
```

As over QQ, the flag is tropicalized by balanced support. An exact prime
binomial initial ideal over GF(2), at the weight used over QQ, proves purity
in dimension 8. The candidate fan is the prevariety of the 78 defining
equations together with the certified tropical basis of the Grassmannian over
GF(2), whose 12 further polynomials lie in the flag ideal. It has dimension 8
and no oversized cell, and exact interior membership keeps all of its cells:
156 rays and 2792 maximal cells in 281 orbits, already coarsest. Stage 6
classifies one initial ideal per orbit over GF(2). The projections stage
certifies that forgetting the point maps every maximal cell onto a maximal cell
of the Grassmannian fan over GF(2) and that all 254 of them are reached.

`flag_gf2/01_prevariety_fan.mrdi` is the prevariety of these 90 polynomials,
not the flag Dressian, which is in `../../data/g2/flag_dressian/`.

`run_census.jl` is the maintained entry point for the G2/P2 Grassmannian
(`--target=base`) and G2/B flag (`--target=flag`). Both targets participate in
stage 6, which classifies the initial ideal at one canonical interior weight
per maximal-cell orbit.

The flag tropicalizer uses all **28 original coordinates and 78 equations**.
It checks the higher-dimensional cells first, proves whole-relative-interior
exclusions, and then tests one representative per orbit of the verified
order-12 signed G2 action.
An exact prime-initial certificate proves purity in dimension 8 before
balanced support selection.
The actual four-dimensional lineality is retained throughout. Its checkpoint
directory is `../../work/g2/flag/tropicalize_symmetry_full28/`.
The output is a coarsened tropical support fan; stage 6 computes initial
ideals at its stored canonical interior weights.

```bash
JULIA_DEPOT_PATH=/tmp/mukai-flag-julia-depot:$HOME/.julia \
  julia --startup-file=no --project=. studies/g2/run_census.jl \
  --target=flag --stage=tropicalize --workers=2 --threads=8 --bits=0 --force
```

For support membership, an all-binomial weight-initial Groebner basis proves
monomial-freeness directly. Other homogeneous initial ideals use Oscar's
cyclic variable-saturation test, with a fresh degree-reverse-lexicographic
ring; its exactness is [Bohm--Keicher--Ren, Proposition 3.1 and Algorithm
3.3](https://arxiv.org/pdf/1603.09241). General product saturation remains the
fallback for nonhomogeneous ideals or unavailable cyclic support.

Exact support comparisons and benchmark hashes are recorded in
[`full28_benchmark.toml`](supplemental/full28_benchmark.toml).

Canonical data are in `../../data/g2/`. The stage-2 and stage-3 base fans have
the same quotient support and incidence, but remain separate artifacts because
stage 3 records the tree gauge representatives.

Stage 4 stores two vectors for every maximal cell, in every target.
`04_canonical_cell_weights.mrdi` is the equal positive sum of its normalized
rays, the weight at which stage 6 computes. `04_normalized_cell_trees.mrdi` is
that sum normalized as the rays are: shifted until its shortest pendant edge
has length 0 and divided by the gcd of its edge lengths. The `cell` rows of
`04_trees.csv` and the atlases show the normalized tree.

The order-12 taxon group fixes taxon 7. The base coarsest fan has 248 maximal
cells in 31 group orbits. Stage 6 computes at the exact equal-ray-sum weight of
one representative per orbit. The flag has 2714 maximal cells in 270 orbits;
its stage 6 records initial-ideal invariants without a degree-one body because
the flag coordinate ring is multigraded. For the flag, monomial-freeness follows
from the certified tropical fan, so stage 6 tests primality and binomiality
without recomputing torus saturation or the degree for every orbit. It saves
each orbit as it finishes; use `--resume` after an interrupted run.

Over `CC`, the 17 Grassmannian and 81 flag toric orbit counts remain the same.
The Grassmannian's one prime nonbinomial orbit over `QQ` splits into two
nonbinomial components over `CC`. See
[`supplemental/geometric_toric_report.md`](supplemental/geometric_toric_report.md)
for the check and its command.

The exact dihedral representation of the order-12 taxon group is analyzed in
[`fan_symmetries/`](fan_symmetries/).

## Multiplicities of the Grassmannian

Oscar's `tropical_variety` does not compute multiplicities for an ideal that is
not principal, binomial or affine linear: it sets every multiplicity to 1
(`compute_weights_and_construct_tropical_variety` in
`src/TropicalGeometry/variety_prime.jl`, Oscar 1.8.2, "TODO: compute
multiplicities"). The documentation does not mention this. The multiplicities
of the Grassmannian fan are therefore computed separately in
[`weights/weights.jl`](weights/weights.jl), following Maclagan--Sturmfels,
Definition 3.4.3: for `w` in the relative interior of a maximal cone, sum the
lengths `deg(Q)/deg(P)` of the primary components `Q` of the torus saturation
of `in_w(I)` at its minimal primes `P`. Over QQbar each QQ-prime contributes
once per conjugate geometric prime (`absolute_primary_decomposition`).

```bash
julia --project=. studies/g2/weights/weights.jl --absolute
```

It treats the 251 maximal cones of `../../data/g2/grassmannian/01_raw_fan.mrdi`
at the sum of their primitive rays, appends one row per cone to
`weights/weights.csv` (an interrupted run resumes), and prints a summary.
Without `--absolute` only the multiplicities over QQ are computed.

Result (about a minute in total): over QQ all 251 maximal cones have
multiplicity 1, and the torus saturation of every initial ideal is prime
(137 initial ideals are already saturated). Over QQbar, 249 cones have
multiplicity 1 and the two cones 158 and 159 have multiplicity 2: these are the
cells of type 21, whose prime splits into two conjugate components, as in the
example in `main.tex`.

[`weights/balancing.jl`](weights/balancing.jl) checks the balancing condition
(MS15, Definition 3.3.1) at all 244 walls with exact lattice normal vectors.
With the QQbar multiplicities the fan is balanced at every wall. With the QQ
multiplicities (all 1) it fails at exactly the six walls of cones 158 and 159,
where four maximal cones meet. So the multiplicities over the algebraically
closed field are the right ones, as in MS15.

```bash
julia --project=. studies/g2/weights/balancing.jl
```

[`weights/saturation_primality.jl`](weights/saturation_primality.jl) checks,
for every face of the raw fan (not only the maximal cones), whether the torus
saturation `in_w(I) : (x_1 ⋯ x_n)^∞` is prime over QQ (`is_prime`) and over
QQbar (`absolute_primary_decomposition`: one component, reduced, one conjugate).
It also checks the 11 initial forms of the quadric q.

```bash
julia --project=. studies/g2/weights/saturation_primality.jl
```

Result, in `weights/saturation_primality.csv` (one row per face):

| faces | number | prime over QQ | prime over QQbar |
|---|---|---|---|
| lineality | 1 | 1 | 1 |
| rays | 51 | 51 | 51 |
| walls | 244 | 244 | 244 |
| maximal cones | 251 | 251 | 249 |

The two exceptions are faces 454 and 455, i.e. maximal cones 158 and 159
(face = cone + 296), of type 21: there the saturation splits over QQbar into
two conjugate components. All 11 initial forms of q are prime over QQ and QQbar.

## The two Dressians

The G2-Dressian and the G2 flag Dressian are the tropical prevarieties of the
defining generators of the Grassmannian and of the flag variety. All their
coefficients are 1 and -1, so the two fans are the same over every field. They
are computed once, as the targets `dressian` and `flag_dressian` of
`run_census.jl`, and stored in `../../data/g2/dressian/` and
`../../data/g2/flag_dressian/`. The censuses of the tropical varieties, over QQ
and over GF(2), keep no copy of them.

```bash
julia --project=. studies/g2/run_census.jl --target=dressian --from=tropicalize --through=symmetry
for stage in tropicalize coarsen normalize trees symmetry; do
  julia --project=. studies/g2/run_census.jl --target=flag_dressian --stage=$stage --threads=8
done
```

Stage 1 writes `01_prevariety_fan.mrdi` after checking the coefficients. For
the Grassmannian it calls Oscar. For the flag it is
[`census/compute_flag_prevariety.jl`](census/compute_flag_prevariety.jl),
which can also be run alone: it uses the census's Gfan intersection and
coordinate-lifting routines, in the min convention with all 28 flag
coordinates. It checks that eliminating the six binomial equations preserves
the generator prevariety and verifies every maximal cone against all 78
original defining equations. The input, native Gfan output, and log are
retained in `../../work/g2/flag/prevariety/`; matching completed checkpoints
are reused.

Neither fan is pure, so stage 2 uses `coarsen_prevariety_fan`.
`02_coarsest_cells.csv` lists the stage-1 cones in each cell. Stages 3 and 4
use the tree gauge and the tree census of the corresponding variety. Stage 5
finds the order-12 taxon group; its additional output
`05_symmetry_face_orbits.csv` counts the faces and their orbits in each
dimension modulo lineality.

| | G2-Dressian | G2 flag Dressian |
|---|---|---|
| stage-1 fan: rays, maximal cones | 51, 239 | 147, 2726 |
| coarsest fan: rays, maximal cones | 51, 239 | 141, 2627 |
| maximal cones by affine dimension | 236 of dimension 6, 3 of dimension 7 | 2594 of dimension 8, 33 of dimension 9 |
| lineality dimension | 3 | 4 |
| f-vector modulo lineality | `(51, 241, 248, 3)` | `(141, 1352, 3584, 2729, 33)` |
| orbit f-vector | `(11, 33, 31, 1)` | `(27, 173, 388, 274, 6)` |
| ray orbits, maximal-cell orbits | 11, 29 | 27, 259 |

The stage-1 fan of the Grassmannian is already the coarsest one. The Gfan fan
of the flag is not: 99 pairs of maximal cones merge and six rays disappear. Its
f-vector modulo lineality, including the zero cone, is
`(1, 147, 1424, 3755, 2837, 36)`.

The rays of the coarsest flag Dressian are those of the tropical flag variety.
Its 2714 maximal cells are all cells of this fan: 2594 are maximal, and 120
are facets of the 33 nine-dimensional simplicial cones. Eleven of these cones
lie over each of the three tetrahedra of the G2-Dressian, and their only
further facets are 15 walls, each in three of them.

Stages 6 and 7 are not defined for these targets. The initial ideals in the
cones above the dimension of the variety contain monomials, and no projection
certificate from the flag Dressian onto the G2-Dressian is computed.

## Maximal-cell projection certificates

```bash
julia --project=. studies/g2/run_census.jl --target=flag --stage=projections
```

This seventh stage uses the normalized coarsest fans and symmetry data of both
targets. It certifies that forgetting the seven point coordinates, retaining
Pluecker coordinates `8:28`, maps every maximal flag cell **onto** a maximal
Grassmannian cell. It checks one representative in each of the 270 flag
Weyl-group orbits and records transports for all 2714 cells; their images cover
all 248 Grassmannian maximal cells. Cached certificates are verified again.

The outputs in `../../data/g2/flag/` are:

- [`07_projection_certificate.toml`](../../data/g2/flag/07_projection_certificate.toml):
  exact rational ray and lineality witnesses, orbit certificates, and transports.
- [`07_projection_orbits.csv`](../../data/g2/flag/07_projection_orbits.csv):
  one flag representative, its rays, and its target maximal cell per orbit.
- [`07_projection_cells.csv`](../../data/g2/flag/07_projection_cells.csv):
  the target cell for every flag maximal cell, together with the flag orbit and
  group-element index into the certificate's `group_permutations` list.

All cell and ray indices are one-based and refer to the stored stage-3 fans.
For each flag ray `F_i`, a witness gives
`F_i[8:28] = scale * R_j + coefficients * L`, where `R_j` is a normalized
Grassmannian ray and the rows of `L` are its stored fan's lineality basis.
`scale` is strictly positive for a nonvertical ray; a vertical ray uses target
index 0 and scale 0. The seven vertical rays lie in target lineality.

The certificate also expresses each projected flag lineality basis vector in
`L` and checks that these images span `L`. For each representative, its
nonvertical ray images are exactly the extreme-ray set of the recorded target
maximal cell. Positive scales prove both cone inclusions, and surjectivity on
lineality proves equality before taking the quotient. Both groups' ray and
cell actions are verified against the exact fans, and coordinate projection
commutes with those actions. The recorded group elements therefore carry each
representative equality to every member of its orbit. This stage makes no
claim about lower-dimensional flag cells.

`../../data/g2/flag/provenance/projections.toml` records SHA-256 hashes of all
inputs and outputs.


The local flag cell atlas lives under `atlas/flag/`.
The Grassmannian atlases live under `atlas/grassmannian/` and
`reference/fine_groebner/`.


## The Grassmannian in the coordinates of S_long

[`run_census_slong.jl`](run_census_slong.jl) tropicalizes G2/P2 over QQ in the
coordinates of the long tripod of Section 8: the 18 Pluecker coordinates other
than p14, p25, p36, and u = p25 - p36, v = p36 - p14, w = p14 - p25 in their
slots ([`ideals/G2GrassmannianSlongIdeal.jl`](ideals/G2GrassmannianSlongIdeal.jl)).
The target `slong` writes to `../../data/g2/grassmannian_Slong/`. The target
`short_long` tropicalizes in all 24 coordinates, the Pluecker coordinates and
u, v, w, and writes to `../../data/g2/grassmannian_Slong/short_long/`.

```bash
for stage in tropicalize coarsen normalize symmetry degenerations summarize; do
  julia --project=. studies/g2/run_census_slong.jl --target=slong --stage=$stage --workers=1
done
for stage in tropicalize coarsen normalize symmetry; do
  julia --project=. studies/g2/run_census_slong.jl --target=short_long --stage=$stage --workers=1
done
julia --project=. studies/g2/run_census_slong.jl --target=slong --stage=projections
```

Both use `run_tropicalization` with Gfan, with the Weyl group as signed
variable permutations. A reflection that exchanges p25 and p36 with sign c maps
u to -c u. Points of trop_S_long are not tree metrics, so rays are normalized
by orthogonal projection and there is no tree stage. The stage `projections`
projects one cell per orbit of the 24-coordinate fan to trop(G2/P2) and to
trop_S_long, and checks that the pieces over a subdivided cell cover it
(`07_projections.csv`, `07_projections.toml`).

## Summaries for the paper

```bash
julia --project=. studies/g2/run_census.jl --target=all --stage=summarize
julia --project=. studies/g2/run_census_gf2.jl --target=all --stage=summarize
```

The last stage writes `08_summary.toml` in each of the six target directories
with what the paper states about that target. It examines one maximal cell per
orbit of the taxon group and multiplies by the orbit sizes; no cone is
enumerated globally. The faces of a representative come from the facets of its
own cone, and a face orbit is identified by its smallest image under the ray
actions. From these the stage computes the f-vector and orbit f-vector modulo
lineality, the combinatorial types of the maximal cells and the polygons among
the 2-faces. The reduced integral homology is computed from the simplicial
complex whose facets are the vertex sets of all maximal cells, which has the
same nerve. The stage checks it against the f-vector by comparing Euler
characteristics. The file also counts the rays and maximal cones of the stage-1
fan and, when there is a stage 6, its toric, prime and binomial orbits and
cells. Over GF(2), each representative is looked up among the maximal cells
over QQ modulo lineality, and the stage counts the toric degenerations of both
fields and of only one.

## Files for Zenodo

```bash
julia --project=. studies/g2/run_census.jl --target=zenodo --stage=zenodo
bash ../production/make_zenodo_tarball.sh
```

The target `zenodo` of `run_census.jl` has the single stage `zenodo`, which
runs after both censuses, over QQ and over GF(2). It writes in
`../../data/g2/` the files of the Zenodo tarball that no other stage writes,
from the artifacts of the six targets ([`census/zenodo_census.jl`](census/zenodo_census.jl)):

- `pointspace{,_gf2}/00_ideal.mrdi`: the quadric Q, the one generator of the
  flag ideal in x1, ..., x7 alone.
- `hexagon/tropicalized_hexagon.mrdi`: the hexagon of the Example in Section 9,
  completed from its flags 1 and 4 by `tropical_complete_hexagon` (70 s,
  2.5 GB) and checked against the table of `main.tex`.
- `grassmannian_orbits.csv`, `flag_orbits.csv`, `dressian/dressian_orbits.csv`
  and `flag_dressian/flag_dressian_orbits.csv`: the labels of `main.tex` of
  the ray and cell orbits. They also read `atlas/grassmannian/ray_label_orbits.csv`,
  `atlas/grassmannian/ray_orbit_labels.csv` and
  `atlas/flag/local_flag_cell_atlas.csv` here, and
  `../../data/g2/maximal_cell_tree_metrics.csv`.
- `../mukai_curves/genus_10_method_{1,2}.mrdi`: the two genus-10 curves of
  Section 10 as polyhedral complexes. That of Method 1 is cut from the raw fan
  by B + q, from its shift q in
  `../../curves_g10/overnight_20261028/wpart_2722/translation_1.mrdi`, by
  `curves/method_1_curve.jl` (about 1 min, cell by cell by
  `curves/bergman_curve.jl`), which also writes its 454 pieces to
  `genus_10_method_1_pieces.mrdi`; that of Method 2 is the stored curve of
  `../../curves_g10/overnight_20261028/experiment_100128/`. Both become Oscar
  complexes by `curves/curve_complex.jl`.

A part whose files all exist is skipped unless `--force`, even when its inputs
are missing, and its stored files are checked again. The provenance is `../../data/g2/provenance/zenodo.toml`.
