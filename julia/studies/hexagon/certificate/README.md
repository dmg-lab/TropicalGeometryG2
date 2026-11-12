# Certificate for Theorem `thm:completing_hexagons`

For each of the 56 coordinate functions f of Lemma `lem:complete_hexagon` (the coordinates of x^(2), x^(3), p^(23),
p^(34), as polynomials in K[p_ij, x_k]), a subdivision of Trop = trop(linespace) × trop(pointspace) into closed cells σ
with a tropically faithful representative f_σ ∈ f + I on each, where I is the ideal of linespace × pointspace (Remark
`rem:hexagon_reduction` of tropicalizing morphisms). Min convention throughout.

## Files

| File | Content |
|---|---|
| `hexagon_certificate.mrdi` | the certificate (Oscar, `load`), format below |
| `GUIDE.md` | a short guide: loading the certificate in Julia, reading a tree, checking a cell by hand |
| `tree_example.tex` | the tree of p23_13 on the maximal cone (297, [1, 4]), drawn (compile with pdflatex) |
| `compute_certificate.jl` | computes it from the definitions (stages `cones`, `setups`, `greedy`, `merge`, `cover`, `assemble`) |
| `verify_certificate.jl` | checks (a)-(d) below from the definitions only |
| `run_compute_overnight.sh`, `run_verification_overnight.sh` | run these on buckow |
| `verification/` | the results of `verify_certificate.jl` |
| `work/` | intermediate files of `compute_certificate.jl` (not part of the certificate) |

The definitions, in `..`: `polynomial_expressions.jl` (the 56 f: `hexagon_expressions()`; I: `line_point_ideal()`),
`tropical_complete_hexagon.jl` (the faces of the raw fan of trop(linespace), `../../../data/g2/grassmannian/01_raw_fan.mrdi`,
and the cones of the tropical quadric), `tropically_faithful_representatives.jl` (Algorithm `alg:faithful_representative`,
used only by `compute_certificate.jl`).

## Format of `hexagon_certificate.mrdi`

A Dict with

- `"coordinate"`: the names `x2_1, …, x2_7, x3_1, …, x3_7, p23_1, …, p23_21, p34_1, …, p34_21` (the indices of p are the
  pairs i < j in lexicographic order), and `"f"`: the 56 polynomials f, in the variables p_12, …, p_67, x_1, …, x_7.
- `"g"`, `"u"`: for each coordinate the list of its representatives h_j = x^(-u_j) g_j.
- `"face"`, `"quadric"`: the 1506 maximal cones C of Trop. C is the product of the face `face` of the raw fan of
  trop(linespace) (numbered as in `raw_fan_faces()`) and the cone `quadric` of the tropical quadric (the terms of
  q = x1x4 + x2x5 + x3x6 − x7² of minimal weight, 1 = x1x4, 2 = x2x5, 3 = x3x6, 4 = x7²).
- `"tree"`: for each maximal cone and each coordinate (in the order above) the tree of the subdivision of C:

      node  := "."  |  j "(" child ("," child)* ")"
      child := t ":" node

  The root is C. A node `j(...)` is a cell σ that is cut by the normal fan of h_j: its children are the cells σ ∩ N_v,
  one for each vertex v of the normal fan of Newt(h_j) on σ (i.e. each v whose normal cone N_v meets σ in a
  full-dimensional cone). The vertex is given by `t`: the position of its exponent vector in the lexicographically sorted
  list of the exponent vectors of g_j, counted from 1 (the smallest one, if several exponents have the same pairing with
  the span of σ). A child `.` is a cell σ of the cover, with f_σ = h_j.

  Example: `1(4:.,9:2(3:.,7:.))` cuts C by h_1 into two cells; on the first h_1 is faithful, the second is cut by h_2 into
  two cells, on which h_2 is faithful.

- `"description"`: a short version of this.

The cells are not stored (with their rays they take 2.2 GB); they are recomputed from the trees.

## What is verified

- (a) ∪ σ = Trop: every cone of Trop lies in a maximal one, there is a tree for every maximal cone and coordinate, and the
  children of each node are exactly the vertices of the normal fan on it.
- (b) f_σ ∈ f + I: x^(u) f − g ∈ I; and the f of the certificate are those of `polynomial_expressions.jl`.
- (c) f_σ is tropically faithful on relint σ: in_w(g) ∉ J_C at a point w of relint σ, where J_C = in_w(I) : (x_1 ⋯ x_28)^∞
  on relint C. σ lies in a Gröbner cone: C does, as in_w(G) is constant on relint C for the reduced Gröbner basis G at a
  point of C. σ lies in one normal cone of Newt(f_σ): it is a child σ' ∩ N_v; so in_w(g) is constant on relint σ.
- (d) J is prime over QQ on every cone of Trop (all 6017, also the non-maximal ones), as needed by the continuity lemma
  `lem:continuity_faithful` for the faces of the cells. Primality over QQ suffices: the components of V(in_w I) are
  Galois conjugate and in_w(g) has rational coefficients.

From the `julia` folder:

    julia --project=. studies/hexagon/certificate/verify_certificate.jl --part=1 --parts=3     (and parts 2, 3)
    julia --project=. studies/hexagon/certificate/verify_certificate.jl --summary

The results are written to `verification/`; the summary ends with `VERIFIED` if all checks pass. With 3 processes on
buckow this takes about 15 h.

## History

The certificate was computed on 2026-11-05/06 (about 11 h on buckow) and verified on 2026-11-07;
`compute_certificate.jl` merges the scripts of that computation with the same logic.


On 2026-11-08/09 `run_compute_overnight.sh` recomputed the certificate on buckow from scratch (about 12 h); the result is
identical to `hexagon_certificate.mrdi`. `run_verification_overnight.sh` then verified it (about 15 h): `verification/`,
summary `6017 cones (1506 maximal), 1632641 cells, 403 representatives`, `VERIFIED`.
