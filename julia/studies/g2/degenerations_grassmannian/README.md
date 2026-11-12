# Toric degenerations of the G2-Grassmannian: comparisons and observations

This folder compares the Newton–Okounkov polytopes of the toric degenerations
of the census of trop(G2/P2) (and, for FFL, of trop(G2/B)) with the
degenerations in the literature. The census data are in
[`../../../data/g2/`](../../../data/g2/): `grassmannian/` and `flag/` over QQ,
`grassmannian_gf2/` and `flag_gf2/` over GF(2). Stage 6 of the census,
`06_canonical_interior_degenerations.csv`, gives one initial ideal per orbit of
maximal cells and its Newton–Okounkov body when the initial ideal is prime.

Run every script from the repository root with `julia --project=julia <path>`.

## Orbits and types

The scripts print census orbit numbers, but the paper uses the types \celltype{n}.
Types 1–20 are the toric types of
[`tables/toric_degeneration_table.tex`](../../../../tables/toric_degeneration_table.tex);
types 21–35 are the other atlas pages
([`../atlas/flag/local_flag_cell_atlas.csv`](../atlas/flag/local_flag_cell_atlas.csv)).

The tropical G2-Grassmannian over QQ has 248 maximal cells in 31 orbits of
the order-12 group. 18 orbits are prime: 17 are toric (135 cells), and orbit
31 (type 21) is prime but not binomial. The bodies are 5-dimensional with
normalized volume 18 (the degree), except type 21, where it is 9.

| Type | QQ orbit | Cells | Newton–Okounkov polytope | f-vector | Lattice points | Fiber polytope over the hexagon |
|---:|---:|---:|---|---|---:|---|
| 1 | 29 | 6 | 5-simplex (Chevalley, lecture hall) | 6 15 20 15 6 | 14 | 3-cube |
| 2 | 27 | 6 | 5-simplex | 6 15 20 15 6 | 16 | 3-simplex |
| 3 | 9 | 12 | Pyramid(Pyramid(Pyramid(4-gon))) | 7 19 26 19 7 | 14 | 3-cube |
| 4 | 20 | 12 | Pyramid(Pyramid(Pyramid(4-gon))) | 7 19 26 19 7 | 14 | Δ1 × Δ2 |
| 5 | 7 | 6 | Pyramid(Wedge(Pyramid(4-gon); segment)) | 9 26 35 24 8 | 14 | 3-cube |
| 6 | 23 | 3 | Pyramid(Wedge(Pyramid(4-gon); segment)) | 9 26 35 24 8 | 15 | 3-cube |
| 7 | 22 | 3 | Wedge(Wedge(Pyramid(4-gon); segment); 3-gon) | 13 38 48 30 9 | 14 | 3-cube |
| 8 | 2 | 12 | Pyramid(Pyramid(M(3,5))) | 8 23 32 23 8 | 15 | 7 vertices, 6 facets |
| 9 | 10 | 12 | Pyramid(Pyramid(M(3,5))) | 8 23 32 23 8 | 16 | 10 vertices, 8 facets |
| 10 | 14 | 6 | Pyramid(Wedge(M(3,5); segment)) | 11 33 44 29 9 | 14 | 5-gonal prism |
| 11 | 21 | 6 | Pyramid(Wedge(M(3,5); segment)) | 11 33 44 29 9 | 14 | 10 vertices, 7 facets |
| 12 | 19 | 6 | Wedge(Pyramid(M(3,5)); 3-gon) | 11 33 44 29 9 | 14 | 3-cube |
| 13 | 25 | 6 | Pyramid(Pyramid(4-gonal antiprism)) | 10 33 50 37 12 | 15 | 14 vertices, 11 facets |
| 14 | 26 | 3 | Pyramid(Pyramid(4-gonal antiprism)) | 10 33 50 37 12 | 15 | 4-gonal antiprism |
| 15 | 4 | 12 | Pyramid(M(4,7)) | 9 27 38 27 9 | 14 | 15 vertices, 11 facets |
| 16 | 1 | 12 | Pyramid(4D, 8 vertices, 8 facets) | 9 27 38 27 9 | 14 | 14 vertices, 10 facets |
| 17 | 8 | 12 | Pyramid(4D, 9 vertices, 8 facets) | 10 30 41 28 9 | 14 | 11 vertices, 8 facets |
| 21 | 31 | 2 | 5-simplex (not toric, volume 9) | 6 15 20 15 6 | 13 | 3-simplex |

M(d,n) is the multiplex; "segment" and "3-gon" are the faces of the wedges.
The last column is from OBSERVATIONS.md, with type 14 updated (see the
caveats below).
Types 18–20 are toric over GF(2) only: GF(2) orbits 9, 20 and 1. The
non-prime QQ orbits are 3, 5, 6, 11, 12, 13, 15, 16, 17, 18, 24, 28 and 30,
which are types 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33 and 34.

## Top-level files: polytope types and the hexagon

[`OBSERVATIONS.md`](OBSERVATIONS.md) (2026-10-06):

- The lineality space of trop(G2/P2) is 2-dimensional. Each of its 21
  columns is 0 or a root of G2, with α short and β long
  ([`roots.png`](roots.png)). Their convex hull is a hexagon H whose
  G2-lattice points are exactly the roots and 0.
- For every prime cell C, the projection π45 of P(C) onto the two lineality
  coordinates is onto H. The table above gives the Billera–Sturmfels fiber
  polytope of π45 (3-dimensional).

How the tables were made:

- [`check_types.jl`](check_types.jl) classifies every stored body with
  `RecursivePolytopeTypes`.
- [`check_fiber_type.jl`](check_fiber_type.jl) does the same for the fiber
  polytope of the projection onto coordinates 4 and 5.
- Both scripts save the parts they cannot name to [`unrecognized/`](unrecognized/)
  through `save_unrecognized` in [`utils.jl`](utils.jl). `utils.jl` also has
  `fiber_polytope` and a catalog test for 3-polytopes (Platonic, Archimedean,
  Catalan, pyramids, prisms, bipyramids, Johnson).
- They use the relative path `julia/data/...`, so run them from the
  repository root.

Caveats, checked 2026-10-20:

- OBSERVATIONS.md was made when the census stored a polytope for all 31
  orbits. Its rows for the 13 non-prime orbits are value polytopes, not
  Newton–Okounkov bodies. Since 2026-10-09 only the 18 prime orbits have a
  `.poly` file.
- The stored bodies now use intrinsic lattice coordinates. There,
  coordinates 4 and 5 no longer project onto the hexagon: the images are
  triangles and quadrilaterals.
- Rerunning the fiber classification on these files gives the same types as
  the table, except orbit 26 (type 14). OBSERVATIONS.md left it unrecognized
  (8 vertices, 10 facets); the current recognizer finds the 4-gonal
  antiprism.
- `unrecognized/` (2026-10-07) holds the bases of the unrecognized bodies of
  orbits 1, 3, 5, 8, 11 and 13. It also holds the unrecognized fibers of
  orbits 1, 2, 4, 5, 8, 10, 11, 21 and 25.

## Spacek_Wang: Chevalley polytopes

This folder implements the Chevalley polytopes of Spacek and Wang
(arXiv:2411.10276) from scratch. The implementation is
[`chevalley_polytopes.jl`](Spacek_Wang/chevalley_polytopes.jl), which
`TropicalCensus` includes. It computes the root system, the minimal coset
representative, V_ϖ, the torus expansion of u_-·v_ϖ, the minimal terms and
the hull.

- [`appendix_chevalley.jl`](Spacek_Wang/appendix_chevalley.jl) checks every
  numerical claim of their Appendix A; the output is in
  [`test_summary.txt`](Spacek_Wang/test_summary.txt).
  - **Results:** 32 checks pass, 0 fail. The checks are degrees
    (Borel–Hirzebruch), volumes, IDP, vertex sets, and the Minkowski
    decompositions of A.1.
  - **Two published volumes are not reproduced:** F4/P4 with s_2 gives 48
    (the paper says 63), and SG(2,8) with s_2 gives 112 (the paper says 105).
  - **Agreement with the paper:** the cases F4 s_2 and C4 s_1, s_2 have
    volume below the degree, so they are not Newton–Okounkov bodies, as the
    paper states.
- **G2/P2 (A.4):** with the word 21212, the Chevalley polytope is the
  lecture hall simplex P_(1,3,2,3,1). It has 6 vertices, 14 lattice points
  and volume 18, and it is IDP. It is census type 1 (QQ orbit 29, column
  `chevalley_equivalent`).
- [`generate_chevalley_polytopes.jl`](Spacek_Wang/generate_chevalley_polytopes.jl)
  writes one polymake file per commutation class for the five spaces of
  Appendix A to `julia/data/chevalley/`. Its run log is
  [`appendix_summary.txt`](Spacek_Wang/appendix_summary.txt).
- [`test_chevalley_triangulations.jl`](Spacek_Wang/test_chevalley_triangulations.jl)
  needs UniTriSat, which is not a project dependency. It checks that the G2
  lecture hall simplex has 418 unimodular triangulations, all regular, 180 of
  them flag (equivalently quadratic). The 4-dimensional facet x_1 = 1 has the
  same counts.
- [`g2_chevalley_study.jl`](Spacek_Wang/g2_chevalley_study.jl) writes the G2
  polytopes for twelve weights to `julia/data/chevalley/`. Its printed table
  was not saved. From the stored files (2026-10-20):
  - Every polytope has volume equal to the degree and exactly dim V_ϖ
    lattice points.
  - G2/P1 and G2/P2 give 5-simplices; G2/B gives polytopes with f-vector
    22 76 120 104 50 12, for both words 121212 and 212121.
  - At ϖ = ω1 + ω2: 64 lattice points and volume 720.
- [`compare_flag.jl`](Spacek_Wang/compare_flag.jl) writes
  [`compare_flag.txt`](Spacek_Wang/compare_flag.txt). It compares the G2/B
  polytopes at ω1 + ω2 with the flag census, as for FFL, through
  [`flag_census.jl`](flag_census.jl).
  - **Summands:** the valuation of a word is defined on all of C[G2/B]. Its
    values on the sections of ω1 (7) and of ω2 (14) give the two summands,
    and their sums are its 64 values on ω1 + ω2. For both words, the summands
    are lattice equivalent to the Chevalley polytopes of G2/P1 (volume 2) and
    G2/P2 (the lecture hall simplex).
  - **Match over QQ:** of the 81 toric orbits, 4 have the f-vector
    (136, 137, 149, 152). For both words, only orbit 149 (12 cells) is lattice
    equivalent, and the pair of summands matches under one map. So the two
    words give the same toric variety in P^6 × P^13.
  - **Over GF(2):** orbit 165 matches; the other orbits with the f-vector
    are 155, 156 and 168.
  - **Where orbit 149 lies:** over type 1 (Grassmannian orbit 29, 6 cells, two
    flag cells over each). Over the atlas cell 121, the tree is the
    caterpillar 16 – 2 – 7 – 5 – 34. Orbit 149 puts the apple on the pendant
    edge of leaf 3 or of leaf 6 (atlas `(:ray, 3)` and `(:ray, 6)`), that is,
    on a leaf of one of the two cherries.

## Christophersen_Ilten: Stanley–Reisner degenerations

Christophersen and Ilten (Math. Z. 278 (2014), Theorem 3.1, g = 10) give a
7×7 skew matrix whose 35 Pfaffians cut out the G2-Grassmannian (dimension 6,
degree 18).

- [`compare.jl`](Christophersen_Ilten/compare.jl) writes
  [`compare.txt`](Christophersen_Ilten/compare.txt), and
  [`verify.jl`](Christophersen_Ilten/verify.jl) writes
  [`verify.txt`](Christophersen_Ilten/verify.txt).
- **Translation to our coordinates:** the matrix is our Plücker matrix on the
  linear space, with u = p17, v = p47, w = p14, r = p25, x00 = p27, x01 = p37,
  x10 = p67, x11 = p13, y00 = p24, y01 = p34, y10 = p16, y11 = p15, z0 = p26,
  z1 = p35 and −w−r = p36.
- **Eight initial ideals:** their term-order conditions leave three
  comparisons open, which gives 8 squarefree quadratic initial ideals with 28
  generators each. Each is the Stanley–Reisner ideal of K ∗ Δ2, where
  Δ2 = {u, v, r} and K is a flag 2-sphere with 11 vertices and 27 edges.
- **T11:** K is T11 for exactly two of the eight. T11 has valencies
  4 4 5 5 5 5 5 5 5 5 6.
- **The starring rule:** read literally, CI's rule fails at the step T5 → T6
  (it gives valencies 3 3 4 4 5 5, not the octahedron). Starting from the
  octahedron, the rule works as stated.
- **Cells below each order:** a cell is below an order when its canonical
  weight lies in the closed Gröbner cone. 12 of the 248 cells are below at
  least one order:
  - type 7: cells 68, 69;
  - type 10: cells 77, 125;
  - type 11: cells 43, 181;
  - type 12: cells 35, 177;
  - type 17: cells 39, 203;
  - type 28 (non-toric QQ orbit 15): cells 101, 131.
- **Below every order:** cells 69, 77, 101, 125 and 131. Each toric type 7,
  10, 11, 12 and 17 also lies below a T11 order, which gives the normality
  statement and the T11 ∗ Δ2 triangulations of Section 6.
- **Consistency checks:**
  - The Gröbner test and the polyhedral test (the regular subdivision is a
    unimodular triangulation with the facets of K ∗ Δ2, Sturmfels Theorem
    8.3) agree on all 135 toric cells.
  - The tests in 14 and in 21 variables agree on all cells.
- **Numbering of the orders:** the two output files number the eight orders
  differently. `verify.txt` uses sign vectors; `compare.txt` uses M1–M8,
  sorted by diagonals. The cell sets agree up to this relabelling; for
  instance the T11 orders are 3 and 6 in one file and M6 and M3 in the other.

## Feigin_Fourier_Littelmann: the polytopes of Definition 11.5

[`verify.jl`](Feigin_Fourier_Littelmann/verify.jl) writes
[`verify.txt`](Feigin_Fourier_Littelmann/verify.txt). Its flag comparison is
in [`flag_census.jl`](flag_census.jl), shared with
[`Spacek_Wang/compare_flag.jl`](Spacek_Wang/compare_flag.jl).

**Conventions** (FFL, Transform. Groups 22 (2017), Section 11):

- Roots β1, …, β6 = 3a1+2a2, 3a1+a2, 2a1+a2, a1+a2, a2, a1, with a1 short.
- The order on PBW monomials is homogeneous reverse lexicographic.
- The essential exponents are computed from the Lie bracket in Oscar. With
  homogeneous lexicographic order they do not give Definition 11.5.

**Grassmannian, ϖ = ω2:**

- **Lattice points:** the essential exponents of V(ω2) (14) and V(2ω2) (77)
  are the lattice points of Definition 11.5, and S(2ω2) = S(ω2) + S(ω2).
- **The polytope:** P(ω2) lies on p6 = 0. It has dimension 5, f-vector
  13 38 48 30 9, 14 lattice points and volume 18.
- **Match:** among types 1–20, it is lattice equivalent only to type 7 (QQ
  orbit 22).

**Flag, ϖ = ω1 + ω2:**

- **Lattice points:** the essential exponents of V(ω1) (7) and of
  V(ω1+ω2) ⊂ V(ω1) ⊗ V(ω2) (64) are the lattice points of Definition 11.5,
  and S(ω1+ω2) = S(ω1) + S(ω2).
- **The polytope:** P(ω1+ω2) has dimension 6, f-vector
  35 123 185 146 62 13, 64 lattice points, volume 720 and Ehrhart
  polynomial (t+1)^6.
- **Match over QQ:** of the 81 toric flag orbits, only orbit 64 (6 cells) is
  lattice equivalent to it. Its summands are Pyramid(Pyramid(Pyramid(4-gon)))
  from the vector coordinates and the type 7 body from the Plücker
  coordinates.
- **The pair of summands also matches:** P(ω1) and P(ω2) map to the two
  summands of orbit 64 under one lattice map. This is tested on the Cayley
  polytopes, keeping each summand in its layer. So the two toric varieties in
  P^6 × P^13 are isomorphic, not only the polytopes.
- **The near miss:** orbit 67 has the same f-vector, the same Ehrhart
  polynomial and the same summand types, but it is not combinatorially
  isomorphic to P(ω1+ω2): the numbers of lattice points on the facets differ.
- **Over GF(2):** of the 83 toric orbits, orbit 88 matches; orbit 91 is the
  near miss.
- **Where orbit 64 lies:** over the type 7 cells (Grassmannian orbit 22, 3
  cells, two flag cells over each). The type 7 tree is the caterpillar
  12 – 3 – 7 – 6 – 45 (cell 59), with inner edges 1, 3, 3, 1. Orbit 64
  puts the apple on a cherry edge, 12|34567 or 45|12367 (atlas
  `(:edge, 8, 9)` and `(:edge, 11, 12)`). The apple position was read from
  the canonical weight through d_k8 = −val(x_k); `verify.jl` does not compute
  it. Orbit 67 puts the apple on the pendant edges of leaves 3 and 6.
