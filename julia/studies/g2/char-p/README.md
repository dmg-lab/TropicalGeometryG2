# Characteristic p: certificates for the tropical G2-Grassmannian

Code and data for the proof of Theorem `thm:char2` of `main.tex` (Section "Characteristic
Two"): trop_p(G) = trop(G) for every odd prime p, and the fan trop_2(G); plus cross-checks for
Theorem `thm:flag`.

## Running

Oscar master (here commit 43cd571f68, 2026-09-30, "Oscar 1.9.0-DEV"), Julia 1.13:

```
env JULIA_DEPOT_PATH=/cache/oscar-dev-depot: julia-1.13 --project=/cache/oscar-system/Oscar.jl SCRIPT.jl ARGS
```

All scripts locate the census data in `julia/data/g2` relative to their own directory.
Conventions: all weights and fans are in the MIN convention, as in `julia/data/g2`;
`main.tex` uses the MAX convention (a weight w here is -w there).  The ideal I over ZZ is
given by the 42 generators of `julia/data/g2/grassmannian/00_ideal.mrdi` (7 linear forms,
35 Pluecker quadrics, all coefficients +-1).

## First inclusion, odd p: trop_p(G) in trop(G)

The proof uses the tropical basis of Corollary `cor:tropicalbasis`: the 42 generators have
coefficients +-1 and the three cubics of `eq:threecubics` have leading coefficient 2.
`char2.jl` (below) checks the three cubics of the paper directly.  `witness.jl` is the
earlier version of the same check with Weyl images of one cubic:

| statement in the proof | file |
|---|---|
| generators have coefficients +-1, so their prevariety is the G2-Dressian over every K | `witness.jl` -> `results/witness.txt` |
| the 12 taxon permutations act by signed permutations p_ij -> s_i s_j p_s(i)s(j) preserving I over ZZ | `witness.jl` -> `results/witness.txt`, `certificates/witness.mrdi` (`weyl_sigma`, `weyl_signs`) |
| f = 2 p17 p27 p37 + p36 p47 p56 - p14 p45 p67 - p26 p45 p57 - p35 p46 p67 lies in I over ZZ | `certificates/witness.mrdi` (`f`, `f_cofactors`) |
| for each tetrahedron a Weyl image of f whose term +-2 p17 p27 p37 is strictly minimal on the open cell | `certificates/witness.mrdi` (`tetrahedron_witnesses`, `tetrahedron_witness_cofactors`), checked in `witness.jl` |

Runtime: about 30 s.

## Second inclusion: the 251 cones of trop(G) lie in trop_p(G)

Certificates of the proof, one per maximal cone of the raw fan over QQ.

`compute.jl CHUNK NCHUNKS` computes, for each of the 251 maximal cones C of
`julia/data/g2/grassmannian/01_raw_fan.mrdi`, a certificate file
`certificates/cones/cone_NNN.mrdi` (a `Dict{String,Any}`, load with Oscar's `load`):

- `w`: interior point of C, the sum of its primitive rays;
- `G`: strong Groebner basis of I over ZZ for the weight order of w refined by degrevlex
  (`weight_ordering(c .- w, degrevlex(R))`, c = max(w) + 1);
- `cofactors`: G[k] = sum_j cofactors[k][j] * gens_j over ZZ;
- `pm1_point`: a common zero in {+-1}^21 of all in_w(g) (237 cones), or empty;
- for the other 14 cones (raw cones 71, 72, 101, 103, 121, 133, 134, 135, 136, 141, 158,
  159, 162, 163; 158 and 159 are the cones of type 21, of multiplicity 2 in
  `../weights/weights.csv`, and the other twelve are the boundary triangles of the three
  tetrahedra):
  `slice_S`, `slice_gens` (in_w(G), p_s - 1 for s in S, y prod_{t not in S} p_t - 1) and
  `slice_G`, a strong Groebner basis of the slice ideal over ZZ (degrevlex).

`verify.jl CHUNK NCHUNKS` re-checks every certificate WITHOUT computing a Groebner basis,
writing `results/verify-CHUNK.csv`:

| column | check | used in the proof for |
|---|---|---|
| `in_I` | G[k] = sum_j cofactors[k][j] gens_j (polynomial identity over ZZ) | G and the generators generate the same ideal over ZZ |
| `I_in_G` | every generator of I is a ZZ-combination of G (division by G, identity checked) | G and the generators generate the same ideal over ZZ |
| `spairs_checked`, `spairs_failed` | every S-polynomial over ZZ of two elements of G with non-coprime leading monomials equals sum q_g g with lm(q_g) lm(g) < lcm of the leading monomials (identity and order checked in Julia) | S-polynomials have standard representations over ZZ (coprime pairs: product criterion); G reduces to a Groebner basis modulo p |
| `lcs`, `lcs_powers_of_2` | leading coefficients of G | G reduces to a Groebner basis modulo every p not dividing them |
| `lin_hom` | every g is homogeneous for the lineality grading | in_v(g) = in_w(g) on relint C |
| `H`, `I` | the terms of in_w(g) tie at every ray of C and are <= all other terms there | in_v(g) = in_w(g) on relint C |
| `torus`, `torus_kind`, `slice_lcs` | +-1 point: all in_w(g) vanish there; slice: `slice_gens` recomputed from in_w(G), `slice_gens` in <`slice_G`>, S-polynomials of `slice_G` with standard representations, no constant, leading coefficients powers of 2 | the in_w(g) have a common zero in the torus over the algebraic closure |

Result (2026-10-02): all 251 cones pass every check; 40271 S-polynomials checked, none
failed; leading coefficients of G in {1} (246 cones) or {1,2} (5 cones); torus points:
237 by +-1 points, 14 by slice certificates with leading coefficients in {1}, {1,2},
{1,8}, {1,64}.

Note.  The verifier never uses that G is a strong Groebner basis.  It checks directly that
G and I generate the same ideal over ZZ and that all S-polynomials have standard
representations over ZZ, which is what the proof uses (Buchberger's criterion modulo p).
For the slice ideals only `slice_gens` in <`slice_G`> is needed (if 1 were in the slice
ideal over K, it would lie in <`slice_G`> over K, impossible since that is a Groebner basis
without constant leading monomial).  The 12 slice cones with a leading coefficient 2 (in G
or in `slice_G`) are exactly the facets of the three tetrahedra; the type-21 cones 158, 159
have all leading coefficients 1, so their certificates hold modulo 2 as well.

Runtime: `compute.jl` and `verify.jl` about 2 minutes each with 10 parallel
chunks (`for c in 1..10: julia compute.jl c 10`).  Logs in `logs/` (not tracked).

## Characteristic 2: `char2.jl`

`char2.jl` -> `results/char2.txt`, `results/char2.csv`, `certificates/gf2/cone_NNN.mrdi`
(NNN = index in `julia/data/g2/grassmannian_gf2/01_raw_fan.mrdi`), about 1 minute, 2.4 GB:

- the 257 maximal cones of the raw fan over GF(2) against the 251 over QQ (same lineality,
  rays compared modulo lineality): 239 coincide, the 12 QQ cones missing over GF(2) are the
  facets of the three Dressian tetrahedra (71, 72, 101, 103, 121, 133, 134, 135, 136, 141,
  162, 163), and the 18 new cones are spanned by the new ray epsilon of a tetrahedron (in its
  relative interior, not the barycenter) and one of its six edges;
- the three cubics of `eq:threecubics`: in I over ZZ; the 2-term is strictly dominant on the
  open tetrahedron (cubics 1, 2, 3 for tetrahedra 3, 2, 1); the four other terms tie at
  epsilon; at each vertex one of them is strictly non-minimal and the other three tie with the
  2-term; each term is the non-minimal one at exactly one vertex (MIN convention);
- for each of the 18 new cones: `G`, a Groebner basis of I over GF(2) for the weight order of
  `w` (sum of the primitive rays, refined by degrevlex), the checks `lin_hom`, `H`, `I` of
  `verify.jl`, and the torus point: the all-ones vector is a common zero of all in_w(g) for
  all 18 cones (`ones_point`); the slice fallback was not needed.

## Supplementary

- `supplementary/g2lib.jl`: Bayer's monomial-freeness test (`in_trop`, MIN convention),
  used by the scripts below.
- `supplementary/char2/`: not used in the proof.  `scan.jl FIELD` tests the 248 canonical cell
  weights (`04_canonical_cell_weights.mrdi`) and an interior point of each tetrahedron;
  `scan-FIELD.csv` for QQ, GF2, GF3, GF5, GF7, GF11, GF13.  Over GF2 exactly the 12
  boundary triangles of the tetrahedra (cells 70, 71, 100, 102, 120, 132, 133, 134, 135,
  140, 159, 160) fail; over QQ and every odd p all 248 cells are in trop.
- `supplementary/tetrahedra/`: Bernd's witness and points inside the tetrahedra over QQ,
  GF2, GF3, GF5 (`g2check.jl`, `g2tet.jl`), results in `RESULTS-2026-10-01.md`.
- `supplementary/flag/`: checks for Theorem thm:flag.
  - `projection.jl`, `rays.jl`, `linrays.jl` (outputs `*.txt`): maximal cells and rays of
    the coarsest flag fan over the coarsest Grassmannian fan (decorated G2-trees).
  - `census_types.jl` -> `census_types.txt`: cell types of the census flag fan
    (2588 tetrahedra, 108 square pyramids, 18 triangular prisms).
  - `direct.jl` -> `direct-*.csv`, `direct_fvector.jl` -> `direct_fvector.txt`: the flag
    fan constructed directly from the Grassmannian via the fiber description; every cell
    equals a census cell; f-vector (141, 1352, 3584, 2714), 27 ray and 270 cell orbits.
  - `fstep2.jl` -> `fstep2-*.csv`: the method of thm:char2 for the 3879 raw flag cones
    (leading coefficients in {1,2,4}).
  - `fiber.jl` -> `fiber-*.csv`: Bayer tests of the fiber description over QQ, GF3, GF5, GF2.
  - `gens.jl`: the flag ideal has 7 linear + 35 Pluecker + 35 incidence + 1 quadric generators.
  Run the chunked scripts with at most 10 parallel jobs (16 jobs ran out of memory).
