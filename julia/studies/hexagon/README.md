# Completing hexagons

The code of Section `sec:hexagons` of main.tex (Theorem `thm:completing_hexagons`, Lemma `lem:complete_hexagon`), with
the tropically faithful representatives of the proof of Theorem `thm:completing_hexagons`.
Two flags `(x1, p12)`, `(x4, p45)` with opposite points extend uniquely to a hexagon with points `x1, ..., x6` and lines
`p12, p23, ..., p61`; each half is given by Lemma `lem:complete_hexagon` applied to `(p12, x4)` and `(p45, x1)`. Min
convention throughout.

## Layout

| Path | Content |
|---|---|
| `polynomial_expressions.jl` | Definitions: the 56 coordinates f of Lemma `lem:complete_hexagon` (`hexagon_expressions`, `evaluate_expressions`), the ideal I of linespace × pointspace (`line_point_ideal`), the ring `EXPRESSION_RING` (p_12, …, p_67, x_1, …, x_7) |
| `complete_hexagon.jl` | Classical completion over any field, including characteristic 2: `complete_hexagon`, `random_flag` |
| `tropically_faithful_representatives.jl` | The algorithm for tropically faithful representatives, for any homogeneous ideal: `faithful_setup` (steps 1-3), `tropically_faithful_representative` (the loop) |
| `tropical_complete_hexagon.jl` | Definitions: the cones of Trop = trop(linespace) × trop(pointspace) (`raw_fan_faces`, `quadric_cone_point`); the lookup of f_σ in the certificate (`certificate_representative`); and the tropical completion at a point: `tropical_complete_hexagon`, `random_tropical_flag`, `tropicalize` |
| `certificate/` | **The proof of Theorem `thm:completing_hexagons`**: `hexagon_certificate.mrdi`, `compute_certificate.jl`, `verify_certificate.jl`; see `certificate/README.md` |
| `tests/test_tropical_hexagon.jl` | Membership tests (`is_tropical_hexagon`, `check_tropical_incidences`) and tests of the above: `run_all_tests` |
| `tests/hexagon_expressions.jl` | Tests of the incidences of Lemma `lem:complete_hexagon` for `hexagon_expressions`, modulo I and at random points over several fields |
| `../../to_tex/hexagon_tree_ring.jl` | Writes `../figures/hexagon_tree_ring.tex` of the repository (`fig:tropical_hexagon`) |
| `figures/hexagon_example_figure.jl` | Writes `figures/hexagon_example.tex` (six decorated G2-trees per hexagon) and `figures/hexagon_examples.md` (the coordinates) |
| `flag_projection_rank.jl` | Rank of the x- and p-parts of the span of each maximal cone of trop(F): bad cones (x-rank < 6), see below |
| `flag_residue_test.jl`, `flag_residue_test/` | The certificate's representatives on the residues that flags actually reach, on the bad cones; results, see below |
| `flag_lift_check.jl`, `flag_lift_check/` | The same, independently, by explicit flags over QQ(t); pilot of representatives modulo the ideal of linespace × flagspace |
| `flag_example_check.jl` | The residue test at the two flags of Example `ex:tropical_hexagon` |
| `line_flag_certificate.jl`, `line_flag_certificate/` | Pilot of a certificate modulo the ideal of linespace × flagspace on sampled cone pairs, for an estimate of the full run |

Include chain: `tropical_complete_hexagon.jl` → `complete_hexagon.jl` → `polynomial_expressions.jl`, and
`tropical_complete_hexagon.jl` → `tropically_faithful_representatives.jl`. Everything else includes
`tropical_complete_hexagon.jl` (or `complete_hexagon.jl`).

## The proof

For each coordinate f, the certificate gives closed cells σ, each a full-dimensional subcone of a maximal cone C of Trop,
and a representative f_σ = x^(-u) g ∈ f + I on each. `certificate/verify_certificate.jl` checks, from the definitions
above only:

- (a) the cells cover Trop,
- (b) f_σ ∈ f + I,
- (c) f_σ is tropically faithful on relint σ (σ lies in a Gröbner cone and in one normal cone of Newt(f_σ)),
- (d) J = in_w(I) : (x_1 ⋯ x_28)^∞ is prime on all 6017 cones of Trop, also the non-maximal ones.

The fan of trop(linespace), `../../data/g2/grassmannian/01_raw_fan.mrdi`, is taken as given.

## Open: the residues of flags (2026-11-14)

The certificate proves faithfulness for generic residues of the mixed pair (p12, x4) in G_d × Q_a. In the theorem the
residue of x4 comes from a flag (x4, p45): it is the x-part of a point of F_w = V(in_w I_F) ∩ T, w = (a, d45).

- `flag_projection_rank.jl` (seconds): of the 2714 maximal cones of trop(F), all have p-rank 6, but only 1311 have
  x-rank 6; 1200 have x-rank 5 and 203 have x-rank 4 (the *bad* cones). On a bad cone the components of F_w are orbits
  of the torus of span σ, so the x-parts of residues lie in a set O of dimension 5 or 4 inside Q_a (dimension 6).
- `flag_residue_test.jl` (about 10 minutes with 4 processes; results in `flag_residue_test/`, format in the file):
  for each bad cone, a random point w, the points of F_w on a slice (one component each), and two random maximal cones
  of trop(G) with a random point d each: whether the initial form of the certificate's representative (or of any other
  representative of that coordinate) vanishes on G_d × O, which is then a residue of every lift. Result: the 140 good
  control cones pass all 280 tests; of the 1383 tested bad cones (20 have points of degree 2, not tested), 611 fail in
  every test, 460 in some, 312 in none. A failing test has on average 3.6 failing coordinates, mostly of p34.
  Whether a test fails depends also on the line cone of d, of which only two per bad cone were sampled.
  The test is at one point x̄ of O. With `--homogeneity` (`flag_residue_test/homogeneity.csv`): in 6090 of the 6091
  failing (test, coordinate) pairs in_v(f~) is homogeneous for the x-part of span σ, so it vanishes on all of G_d × O,
  and the valuation is larger than the certificate's value for every lift of the flags (the exception: cone 1975,
  coordinate 40). All slices have determinant 1 and one point, so F_w is a single orbit on the tested cones.
- `flag_lift_check.jl` (about 15 minutes for 80 trials), independent of the residues: explicit special flags (x4, p45)
  over QQ(t) and random flags (x1, p12), the valuations of the 56 coordinates f(p12, x4) against the certificate's
  values. Result (seed 1, 80 trials; `flag_lift_check/trials.csv`): 59 flags in good cones, all agree; 8 in bad cones, 4
  of them with valuations above the certificate's value (trials 17, 43, 61, 74; e.g. cone 1119, x-rank 4: 11
  coordinates), never below; 13 on the boundary of cones. For the 8 flags in bad cones, at the same tropical point:
  the residue of the lift lies in F_w and G_d, and the coordinates with a larger valuation are exactly those with
  in_v(f~) = 0 at the residue of the lift, and exactly those failing the residue test.

So on bad cones the certificate's values are in general not the valuations of the hexagon, and the proof of Theorem
`thm:completing_hexagons` does not hold there. A failure needs the mixed pair v = (d, a) on a wall of the certificate's
cells: in the relative interior of a cell, in_v(f~) is homogeneous for the span of the cell and V(in_v I) ∩ T is a union
of orbits of its torus, so in_v(f~) vanishes nowhere, and the certificate's value is the valuation for every lift.

Pilot of a repair (`flag_lift_check.jl --pilot=2`, `flag_lift_check/pilot.csv`): tropically faithful representatives
modulo the ideal of linespace × flagspace, I_G(p12) + I_F(x4, p45) (49 variables), at W = (d12, a4, d45), for trials 17
and 43. All 56 coordinates have one, found after 0 steps of the algorithm where the certificate is right and after 1
step (one multiplication by p45_12) where it is wrong; its value is the valuation at the explicit lift. Setup 5-8 s per
point, then at most 0.1 s per coordinate (one 1.6 s), at most 183 terms.

`flag_example_check.jl` (about 2 minutes): Example `ex:tropical_hexagon` is not affected. Its flags (a1, d12) and (a4,
d45) lie in the relative interiors of the good cones 797 and 648, and the residue test passes for all 56 coordinates of
both halves.

## Running

`tropical_complete_hexagon` looks up the faithful representatives in `certificate/hexagon_certificate.mrdi`
(`certificate_representative`): loading the certificate takes about a minute once per session, a hexagon then a fraction
of a second. With `method=:algorithm` they are computed by the algorithm instead; its steps 1-3 only depend on the Gröbner
cone of the weight, so they are done once per cone in a session, which takes a few seconds. The tropical code reads the
fans in `../../data/g2/flag/` and `../../data/g2/grassmannian/`, and is over QQ. From the `julia` folder:

```sh
julia --project=. studies/hexagon/tests/test_tropical_hexagon.jl
julia --project=. studies/hexagon/tests/hexagon_expressions.jl
julia --project=. to_tex/hexagon_tree_ring.jl [--complete]
julia --project=. studies/hexagon/figures/hexagon_example_figure.jl
bash studies/hexagon/certificate/run_compute_overnight.sh          # many hours
bash studies/hexagon/certificate/run_verification_overnight.sh     # about 15 h
```
