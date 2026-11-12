# A short guide to `hexagon_certificate.mrdi`

What the certificate proves, and its exact format: `README.md`. Min convention throughout.

For each of the 56 coordinates f of Lemma `lem:complete_hexagon` and each maximal cone C of
Trop = trop(linespace) × trop(pointspace), the certificate subdivides C into cells σ and gives a tropically faithful
representative f_σ ∈ f + I on each. The representatives of f are h_1, h_2, ..., with h_j = x^(-u_j) g_j: a monomial in
all 28 variables p_12, ..., p_67, x_1, ..., x_7 times a polynomial.

## Loading it

From the `julia` folder, `julia --project=.`, then

```julia
include("studies/hexagon/certificate/verify_certificate.jl")   # definitions and helpers; runs nothing
cert = hexagon_certificate();                                  # about a minute
```

| | |
|---|---|
| `cert.names[i]` | the name of coordinate i = 1..56: `x2_1..7, x3_1..7, p23_1..21, p34_1..21` |
| `cert.reps[i][j]` | h_j of coordinate i: `.g` = g_j, `.u` = u_j, `.exponents` = the exponent vectors of g_j in lexicographic order (this numbers the terms of g_j) |
| `cert.cones[n]` | the maximal cone C number n = 1..1506, as `(k, S)`: k is a face of the raw fan of trop(linespace), S the terms of q of minimal weight (1 = x1x4, 2 = x2x5, 3 = x3x6, 4 = x7²) |
| `cert.trees[n][i]` | the subdivision of cone n for coordinate i |
| `trop_cone(k, S)` | the cone C = {w : E w = 0, A w ≤ 0}, with fields `.cone`, `.E`, `.A`, `.point` (in relint C), `.dim` |
| `certificate_representative(i, w)` | f_σ = h_j at a point w of Trop: `.j`, `.g`, `.u`, `.value` = trop(f_σ)(w) |

(`cert.parsed` is only a cache of parsed trees.)

Example: `cert.cones[1]` is `(k = 297, S = [1, 2])`. For w = (p, x) ∈ R^21 × R^7 this is the cone of all w with p in
face 297 of the raw fan (a maximal face, of dimension 6) and x1 + x4 = x2 + x5 ≤ x3 + x6, 2 x7 (dimension 6). So C has
dimension 12.

## Why these numbers

- **56 coordinates**: x^(2), x^(3) have 7 coordinates each, p^(23), p^(34) have 21 Plücker coordinates each.
- **6017 cones of Trop**: 547 faces of the raw fan of trop(linespace) times 11 cones of the tropical quadric (the sets S
  of at least two of the 4 terms of q).
- **1506 maximal cones** (dimension 12): 251 maximal faces times 6 maximal quadric cones (the pairs S). Only these are
  subdivided; the other cones are their faces. J is prime on all 6017 cones (check (d) of `README.md`).
- **1,632,641 cells** over all coordinates and maximal cones (about 29,000 per coordinate), stored as 1506 × 56 trees.
- **403 representatives**, 1 to 25 per coordinate.

## Reading a tree

`cert.trees[n][i]` is a string:

- `j(...)`: this piece of C is cut along the normal fan of Newt(g_j), h_j = `cert.reps[i][j]`, i.e. into the parts on
  which one and the same term of g_j has minimal weight;
- `t:` before a child: the part on which the t-th term of g_j has minimal weight;
- `.`: a cell σ; f_σ is the h_j of the node just above it.

The subdivision is lazy: a piece is cut by h_1; where h_1 is faithful it stops, the other parts are cut by h_2, and so on.

`cert.trees[1][15]`, p23_1 on `cert.cones[1]`:

    1(14:.,16:.,25:.,36:.,37:.,39:.,51:.,52:.,69:.,77:.,82:.,85:.,86:.,92:.,95:.,97:.,98:.,103:.,106:.,109:.)

h_1 cuts C into 20 cells, one for each term of g_1 that has minimal weight somewhere on C, and is faithful on all of
them. So on C the valuation of f at a generic lift of w is trop(h_1)(w), linear on each cell.

`cert.trees[3][27]`, p23_13 on `cert.cones[3]` = (297, [1, 4]), drawn in `tree_example.tex`:

    1(17:.,29:.,33:.,50:2(49:3(45:4(48:.,61:.,63:.,92:.)),
                         52:3(39:4(35:.,48:.,61:.,63:.),45:4(35:.,48:.,61:.,63:.)),
                         53:3(39:4(37:.,48:.,61:.,63:.),45:4(37:.,48:.,61:.,63:.))),
      55:.,57:.,68:.,71:.,74:.,89:.,97:.,99:.)

h_1 cuts C into 12 parts and is faithful on 11 of them. The part where the 50th term of g_1 is minimal is cut by h_2,
then by h_3 (on the part `49:` only one term of g_3 is ever minimal, so it is not cut), then by h_4, which is faithful on
all 20 cells it produces. So f_σ = h_1 on 11 cells and f_σ = h_4 on 20 cells, e.g. on the cell where the 50th term of
g_1, the 52nd of g_2, the 39th of g_3 and the 35th of g_4 have minimal weight. h_2 and h_3 only cut.

## Checking a cell by hand

This builds that cell, takes a point w in its relative interior, and checks what the verification checks there:

```julia
n, i = 3, 27                                    # p23_13 on cert.cones[3]
c = trop_cone(cert.cones[n].k, cert.cones[n].S)
I = line_point_ideal()
f = verification_coordinates()[i][2]
J = check_trop_cone(I, c).J                     # J on C (also returns whether J is prime)
# the part where the t-th term of g_j has minimal weight: <e_t - e, w> ≤ 0 for all other exponent vectors e
minimal(j, t) = (e = cert.reps[i][j].exponents; matrix(QQ, e[t:t, :] .- e[setdiff(axes(e, 1), t), :]))
σ = cone_from_inequalities(vcat(c.A, minimal(1, 50), minimal(2, 52), minimal(3, 39), minimal(4, 35)), c.E)
dim(σ)                                          # 12
w = sum(integer_rays(σ)[1])                     # a point of relint σ
r = certificate_representative(i, w)            # r.j == 4: f_σ = h_4
v = gens(base_ring(I))
iszero(normal_form(prod(v .^ r.u) * f - r.g, I))                          # f_σ ∈ f + I
!(initial_form(r.g, w) in J)                                              # f_σ is faithful at w
initial_form(cert.reps[i][1].g, w) in J                                   # h_1 is not
r.value == tropically_faithful_representative(f, cone_setup(w); w).value  # the algorithm agrees
```

All of these return `true`. Any other path of a tree works the same way.

## Use

`certificate_representative(i, w)` finds a maximal cone containing w and walks down its tree, always to the child whose
term of g_j has minimal weight at w. Only minima of linear forms are computed, no Gröbner bases. On a wall between cells
any of them may be taken: their representatives are all faithful there. `tropical_complete_hexagon` uses exactly this.
