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
