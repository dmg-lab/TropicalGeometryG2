# Tropical Geometry of G2: code

Code for the paper *Tropical Geometry of G2* by Michael Joswig, Kevin Kuehn, Bernd Sturmfels and Mathieu Vallée,
[arXiv:2611.24456](https://arxiv.org/abs/2611.24456). It computes the data of the paper, which are also on Zenodo,
[10.5281/zenodo.24104003](https://doi.org/10.5281/zenodo.24104003), and checks every computational statement of the
paper against them. Section and theorem numbers below are those of arXiv:2611.24456v1; comments in the code that
mention `main.tex` refer to its source.

Run everything from `julia/`: `cd julia` after cloning; the commands and paths below assume it.

## Setup

Julia 1.13 (the manifest was resolved with 1.13.3: Oscar 1.8.2, Polymake.jl 0.13.14, HiGHS 1.26.1) and Gfan 0.8, whose
`gfan_tropicalprevariety` must be on the PATH:

```sh
git clone https://github.com/dmg-lab/TropicalGeometryG2.git && cd TropicalGeometryG2/julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
bash scripts/install_gfan.sh   # links gfan_tropicalprevariety into ~/bin, which must be on the PATH
```

The scripts write scratch files to `work/` and results to `results/`; git ignores both.

## The data: download or compute

The scripts read the data from `data/`. Fill it in one of two ways; both give the same files.

1. **Download** `TropicalGeometryG2_data.tar.gz` from Zenodo into the repository root (`..`) and unpack its folder
   `data/` into `julia/`:

   ```sh
   tar -xzf ../TropicalGeometryG2_data.tar.gz --strip-components=1 TropicalGeometryG2_data/data
   ```

2. **Compute** it: the census, from the ideals, the shift vector of the first curve of genus ten and the run of the
   second, about 1 h 40 min and at most 6 GB on a laptop. The tables that name the census orbits as the paper does are
   kept in git (`studies/g2/atlas/`).

   ```sh
   julia --project=. make_paper_data.jl --list   # the steps, their time, what is stored
   julia --project=. make_paper_data.jl          # the missing steps, one Julia process each
   bash ../production/make_zenodo_tarball.sh     # optional: the Zenodo tarball from data/, at the repository root
   ```

## Check the paper

```sh
julia --project=. paper_results.jl                 # every section, about 7 minutes and 4.5 GB
julia --project=. paper_results.jl --quick         # without the slowest checks, about 5 minutes
julia --project=. paper_results.jl --sections=3,9
```

It goes through the paper in order, computes each statement from `data/` and the certificates stored next to
the scripts, compares it with the paper, and prints the tables that the paper describes in words. The exit status is 1
if a statement fails. Every check runs in both ways of getting `data/` (computed, or unpacked from Zenodo); only
`--quick` leaves some out. With the Zenodo tarball at the repository root, Appendix A also compares its files with
`data/`.

## Quick examples

In `julia --project=.`, from `julia/`; times on a laptop (Intel i7-12700H).

```julia
using Oscar, TropicalGeometryG2                           # 15 s

# The ideal of G2/P2 in P^20: 7 linear forms and 35 quadrics in the Plücker coordinates p12, ..., p67.
include("studies/g2/ideals/G2GrassmannianIdeal.jl"); using .G2GrassmannianIdeal
G = build_g2_grassmannian_ideal()                         # (ring, ideal, generators)
dim(G.ideal), degree(G.ideal)                             # (6, 18)

# The G2-Dressian, the tropical prevariety of the 42 generators (Theorem 2.4): Gfan, 1 s.
D = gfan_prevariety(G.generators; threads=8)
f_vector(D)[lineality_dim(D)+1:end]                       # [51, 241, 248, 3] before coarsening
# Oscar alone: tropical_prevariety(G.generators, tropical_semiring_map(QQ, min)), 13 min.
# The flag Dressian: include("studies/g2/ideals/FlagG2Ideal.jl"), then
# gfan_prevariety(FlagG2Ideal.build_full_flag_g2_ideal().generators; threads=8), 10 s.

# trop(G2/P2) (Theorem 3.1): Oscar, 3 min and 2.5 GB.
T = tropical_variety(G.ideal, tropical_semiring_map(QQ, min); skip_decomposition=true)

# The stored fans have the coarsest structure; their f-vector modulo lineality and the homology of the link.
F = load("data/g2/grassmannian/03_normalized_coarsest_fan.mrdi")
f_vector(F)[lineality_dim(F)+1:end]                       # [51, 241, 248]
H = TropicalGeometryG2.reduced_integral_homology(TropicalGeometryG2.fan_cells(F))
TropicalGeometryG2.describe_reduced_homology(H)           # "H_2 = Z^57"

# The action of the Weyl group D12: the orbits of rays and maximal cells.
S = TropicalGeometryG2.load_symmetry_artifacts("data/g2/grassmannian/05_symmetry")
length(S.ray_orbits), length(S.cell_orbits)               # (11, 31)

# Everything the census states about a fan, from its stored files (8 s).
census_summary("data/g2/grassmannian"; stage1_fan="data/g2/grassmannian/01_raw_fan.mrdi")

# A Newton-Okounkov polytope of a toric degeneration (Theorem 6.3).
P = load("data/g2/grassmannian/06_canonical_interior_degenerations_objects/orbit_001_newton_okounkov_body.poly")
dim(P), normalized_volume(P), has_idp(P)                  # (5, 18, true)

# The tropicalized hexagon of Section 9, and the curve of genus ten of Method 2 (Section 10).
h = load("data/g2/hexagon/tropicalized_hexagon.mrdi")     # (points, lines), min-convention
C = load("data/mukai_curves/genus_10_method_2.mrdi")      # (complex, multiplicities)
n_maximal_polyhedra(C.complex)                            # 447
```

<details>
<summary>Layout of the repository</summary>

```text
julia/
  Project.toml, Manifest.toml   the environment, which is also the package TropicalGeometryG2
  src/                          TropicalGeometryG2: tropicalization, coarsening, trees, symmetry, degenerations
  make_paper_data.jl            writes data/
  paper_results.jl              checks the paper, with paper_results/section_*.jl
  data/                         the data (see data/README.md); git keeps only the inputs
  studies/g2/                   G2/P2 and G2/B
    run_census.jl               the census over QQ; run_census_gf2.jl over GF(2); run_census_slong.jl Section 8
    census/                     the code the drivers include; zenodo_census.jl writes the rest of the Zenodo tarball
    ideals/                     the ideals of G2/P2 and G2/B over QQ and GF(2)
    char-p/                     the certificates of Section 4 (other characteristics)
    weights/                    multiplicities and balancing
    atlas/                      the names of the ray and cell types (Section 5): tables kept in git, their scripts
    degenerations_grassmannian/ comparisons with Spacek-Wang, Christophersen-Ilten, Feigin-Fourier-Littelmann
    homogeneous_spaces/         Section 8
    curves/                     the curves of genus ten and the Mukai model (Section 10)
    supplemental/, reports/     further checks and their reports
  studies/hexagon/              Section 9: completion of hexagons and the certificate of Theorem 9.1
  curves_g10/overnight_20261028/       the runs of the two curves of the paper (shift vector q of Method 1)
  to_tex/                       drawing code that the atlas, curve and check scripts load
production/make_zenodo_tarball.sh   packs julia/data/ into the Zenodo tarball
```

`data/`, once filled:

```text
g2/<target>/          one folder per fan: dressian, flag_dressian, grassmannian (G2/P2 over QQ), grassmannian_gf2,
                      flag (G2/B over QQ), flag_gf2, grassmannian_Slong (Section 8), with the files of the stages
  00_ideal.mrdi                         the ideal
  01_raw_fan.mrdi / 01_prevariety_fan.mrdi   the fan as computed
  02_coarsest_*                         the coarsest fan structure and its report
  03_normalized_coarsest_fan.mrdi       the fan of the paper, rays normalized; 03_normalized_rays.mrdi
  04_trees.csv, 04_*                    the trees of the cells
  05_symmetry_*                         the action of D12 on rays and cells, and the orbits
  06_canonical_interior_degenerations*  initial ideals, toric degenerations, Newton-Okounkov polytopes
  07_projection_*                       the projection of trop(G2/B) onto trop(G2/P2)
  08_summary.toml                       the statements of the paper about the fan
  provenance/                           versions, options and hashes of each stage
g2/{grassmannian,flag}_orbits.csv       the orbit types of the paper, by census orbit
g2/pointspace{,_gf2}/, g2/hexagon/      the quadric of P^6, the tropicalized hexagon
mukai_curves/genus_10_method_*.mrdi     the two curves of genus ten and the pieces of the first
```

</details>

<details>
<summary>The census</summary>

`studies/g2/run_census.jl` (over QQ), `run_census_gf2.jl` (over GF(2)) and `run_census_slong.jl` (Section 8) compute
their targets in stages: `tropicalize`, `coarsen`, `normalize`, `trees`, `symmetry`, `degenerations`, `projections`
and `summarize`. Without `--stage` or `--from` they only print what is stored. A stage whose files are stored is
skipped unless `--force` is given. For example (`make_paper_data.jl --dry-run --force` prints all the commands it
runs, in order):

```sh
julia --project=. studies/g2/run_census.jl          # status
julia --project=. studies/g2/run_census.jl --target=base --from=tropicalize --through=summarize
julia --project=. studies/g2/run_census.jl --target=flag --stage=tropicalize --workers=1
julia --project=. studies/g2/run_census_gf2.jl --target=base --from=tropicalize --through=summarize --bits=0
julia --project=. studies/g2/run_census.jl --target=zenodo --stage=zenodo
```

The targets of `run_census.jl` are `base` (G2/P2), `dressian`, `flag` (G2/B), `flag_dressian` and `zenodo`. The other
options are `--engine=auto|oscar|parallel`, `--workers`, `--threads`, `--bits`, `--resume` and `--no-pdf`. G2/P2 uses
Oscar's `tropical_variety`; G2/B uses an exact, resumable parallel engine with Gfan and Oscar, whose checkpoints go to
`work/`. See `src/README.md` and `studies/g2/README.md`.

</details>

<details>
<summary>Section by section: the scripts behind each statement</summary>

`paper_results.jl` runs or reads all of these. The certificates and reports they wrote are stored next to them.

- **Section 1.** The model of G2: `studies/g2/supplemental/g2_split_symmetry.jl`,
  `quadric_computation.jl`; the stabilizer of Φ and the representation E: `studies/g2/char-p/adjoint_char.jl`;
  dimension 6 and degree 18 over QQ, GF(2), GF(3), GF(5): `studies/g2/char-p/supplementary/tetrahedra/g2check.jl`;
  the 78 generators of the flag ideal: `studies/g2/char-p/supplementary/flag/gens.jl`.
- **Sections 2 and 3.** The f-vectors, orbits and homology: the census, `08_summary.toml`. The action of D12:
  `studies/g2/fan_symmetries/analyze_dihedral_action.jl`. Multiplicities and balancing:
  `studies/g2/weights/`. The three cubics of the tropical basis: `sh to_tex/cubic_tree_table.sh`.
  Cross-checks of trop(F): `studies/g2/char-p/supplementary/flag/`.
- **Section 4.** Other characteristics: the certificates of `studies/g2/char-p/` (see its README) and the scan
  `supplementary/char2/scan.jl`; characteristic two: `run_census_gf2.jl`, checked by
  `studies/g2/supplemental/verify_gf2_census.jl`.
- **Section 5.** The ray and cell types and their trees: `studies/g2/atlas/`.
- **Section 6.** The toric degenerations: the stage `degenerations` of the census;
  `studies/g2/supplemental/toric_body_equivalence.jl`, `verify_geometric_toric.jl`, `g2_prime_cone.jl`; the
  comparisons with the literature: `studies/g2/degenerations_grassmannian/`.
- **Section 7.** The fiber complexes: `to_tex/local_flag_complex_table.jl`; the 807 prime cones: the census of
  G2/B and the atlas `studies/g2/atlas/flag/`, counted by `paper_results.jl`.
- **Section 8.** `studies/g2/homogeneous_spaces/` (see its README) and `run_census_slong.jl`.
- **Section 9.** Theorem 9.1: the certificate `studies/hexagon/certificate/hexagon_certificate.mrdi`, checked by
  `verify_certificate.jl` and computed by `compute_certificate.jl` there (hours; see the README there). Lemma 9.2:
  `studies/hexagon/tests/`. The hexagon of the example: `to_tex/hexagon_tree_ring.jl --complete`.
- **Section 10.** In `studies/g2/curves/` (see its README): the Mukai model, `mukai_g10_degenerations.jl` and
  `mukai_models.jl`. Method 1: the search `core_trivalence.jl` with `point_filter.jl`; the curve from its shift q:
  `method_1_curve.jl`, cell by cell by `bergman_curve.jl`. Method 2: the search `sample_genus_10.jl`, the curve by
  `hyperplane_intersection_census.jl` and `curve_details.jl`. Both curves as complexes: `curve_complex.jl`. Rays and
  degree: `curve_rays.jl`; the K3 surface of Method 1: `k3_surface_polygons.jl`, `k3_fins.jl`. The runs:
  `curves_g10/overnight_20261028/`.
- **Appendix A.** `../production/make_zenodo_tarball.sh` and the stage `zenodo` of the census.

</details>
