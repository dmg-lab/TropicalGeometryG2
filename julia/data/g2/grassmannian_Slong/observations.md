# Tropical G2-Grassmannian in the coordinates of S_long

Coordinates: the 18 Pluecker coordinates other than p14, p25, p36, and the
zero-weight forms u = p25 - p36, v = p36 - p14, w = p14 - p25 of the long tripod
S_long (Section 8), in the slots of p14, p25, p36. Characteristic 0. The
directory `short_long/` holds the tropicalization in all 24 coordinates
(S_short and S_long together), of which both are coordinate projections.

## The fans

| | trop(L), Pluecker | trop_S_long(L) | both, 24 coordinates |
|---|---|---|---|
| f-vector | (51, 241, 248) | (42, 211, 275) | 57 rays, 383 maximal cells |
| orbit f-vector | (11, 33, 31) | (11, 31, 35) | 13 ray orbits, 48 cell orbits |
| maximal cells | 245 triangles, 3 quadrilaterals | 263 triangles, 12 quadrilaterals | |
| homology | Z^57 in dim 2 | Z^105 in dim 2 | |
| toric degenerations | 135 in 17 orbits | 21 in 5 orbits | |

The symmetry group of trop_S_long(L) is the Weyl group, of order 12. Its 35
orbits of initial ideals are monomial-free of degree 18. 23 cells are prime
and 27 binomial, and none of the 5 toric bodies is the Chevalley polytope. The
2 cells coming from type 21 stay prime and not binomial.

## trop(L) against trop_S_long(L)

- 36 rays (9 orbits) are rays of both. The 3 rays of type alpha_4 of trop(L)
  (splits 14|, 25|, 36|) go to the lineality space in S_long, and the 12 rays
  of type beta_1 land inside maximal cells of trop_S_long. Conversely, the 3
  rays of u, v, w go to the lineality space of trop(L), and 3 rays of
  trop_S_long land inside the cells of type 33.
- Over trop(L), 245 cells lift to one cell of the 24-coordinate fan and the 3
  cells of type 33 are cut into 3 each. 129 = 3 x 43 cells are vertical: those
  through the rays of u, v, w.
- Over trop_S_long(L), 230 cells lift to one cell and 45 cells (5 orbits) are
  cut into 2 each. 63 = 3 x 21 cells are vertical: the 21 cells of trop(L)
  through the alpha_4 ray of p36 of Remark rem:weight-basis, once for each of
  p14, p25, p36.

## A PL map trop(L) -> trop_S_long(L)

The 54 rays of the 24-coordinate fan other than those of u, v, w project to
distinct rays of trop(L): its 51 rays and the 3 rays inside the cells of type
33. So the cells without a ray of u, v, w form the graph of a continuous PL
section over trop(L), and projecting gives a PL map to trop_S_long. With
d = -val, it fixes the 18 other coordinates and sets
d_u = d_v = d_w = max(d14, d25, d36). The exception is the cells of type 33,
where one of u, v, w cancels to leading order at the interior ray. The map
collapses the 63 cells through the alpha_4 rays, and its image is the 146 cells
of trop_S_long without a ray of u, v, w; it misses the other 129 (43 per ray).
The map back is the same with the roles exchanged.

No PL homeomorphism exists, since H_2 is Z^57 against Z^105. In slot
coordinates the alpha_4 ray of p14 and the ray of u are the same vector, but no
cell contains two alpha_4 rays or two rays of u, v, w. The identity maps 30 of
the 51 rays of trop(L) outside trop_S_long.

## Toric degenerations

No toric degeneration of trop(L) stays toric in S_long. Types 6, 9, 11, 13, 14,
15 (42 cells) collapse. Types 1-5, 7, 8, 10, 12, 16, 17 (93 cells) map onto or
into cells of trop_S_long that are not toric. The 21 toric cells of
trop_S_long lie over walls of trop(L) (18 cells) or inside the cells of type
33 (3 cells). In the weight basis, by comparison, 90 of the 135 stay toric.

## Computation

`julia/studies/g2/run_census_slong.jl`, with Gfan through
`run_tropicalization` and the Weyl group as signed variable permutations: a
reflection exchanging p25 and p36 with sign c maps u to -c u. The comparison is
`07_projections.csv` (one cell per orbit of the 24-coordinate fan) and
`07_projections.toml` (per target orbit). Commands are in
`julia/studies/g2/README.md`.
