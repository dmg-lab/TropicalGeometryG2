# The tropical K3 surface of shift wpart 2722

X18 = trop(G2/P2) ∩ (U(18,21) + t), the tropical K3 surface that contains the
genus-10 curve X17 = trop(G2/P2) ∩ (U(17,21) + t) of this folder.

- trop(G2/P2): the 248 maximal cells of `julia/data/g2/grassmannian/02_coarsest_fan.mrdi`.
- U(18,21): the Bergman fan of the uniform matroid, {x : min x attained at least 4 times},
  with cones B_S = {x_S constant, x_j ≥ x_S off S}, |S| ≥ 4.
- t: trial 2722 of the `wpart` night search (row 2722 of
  `julia/results/curves/genus_10/core_trivalence/core_trivalence_wpart.csv`).
  `translation_1.mrdi` here is 2t, from the certified run. Both sets are fans, so 2t gives the same combinatorics.

All counts are modulo the diagonal (1,...,1). X18 has no lineality beyond the diagonal.

## f-vectors

|                | vertices | edges | polygons |
|----------------|---------:|------:|---------:|
| X18            |      454 |  2603 |     3153 |
| bounded part   |      454 |  1165 |      713 |
| unbounded      |          |  1438 |     2440 |

## Bounded part

- Euler characteristic 454 − 1165 + 713 = 2.
- Rational homology: b0 = 1, b1 = 0, b2 = 1, so it is a homology 2-sphere.
- It is not a surface. Bounded edges by the number of bounded polygons containing them:
  1^422 2^369 3^170 4^129 5^75. It is a 2-sphere with fins, see below.
- Every vertex lies on a bounded edge.
- Bounded polygons by number of sides: 3^337 4^334 5^39 6^2 7^1.

Bounded edges by the number of polygons of X18 containing them: 3^941 4^14 5^210.
- 3 or 4 polygons: the 955 edges on a wall of trop(G2/P2) (3 or 4 cells around the wall).
- 5 polygons: the 210 edges inside a cell (|S| = 5).

## The 2-sphere without the fins (2026-11-10)

H2 of the bounded part over QQ is spanned by one cycle z, with coefficients ±1 on 27 polygons. These 27 polygons
form a 2-sphere: 33 vertices, 58 edges, each edge in exactly 2 of them, every vertex link one circle, Euler
characteristic 2. The other 686 bounded polygons are fins. The 40 edges of X18 on the core of the curve all lie
on the sphere, which they cut into 11 discs, the faces of the planar core. The other 170 bounded edges of the curve
lie on fins only. The core (18 trivalent vertices, 27 edges) is simple and 3-connected, so by Steinitz's theorem
it is the edge graph of a simple 3-polytope, as the last paragraph of Section 10 of main.tex says (fig:planar_curve).

    julia --project=julia julia/studies/g2/curves/k3_fins.jl julia/curves_g10/overnight_20261028/wpart_2722

reads the output of k3_surface_polygons.jl (below) and checks all of this.

## Method (2026-10-28)

The polygons are the pieces C ∩ (B_S + t) with C a maximal cell and |S| = 4:
1. A floating-point pass selects the pairs (C, S) that may meet, with a margin of 1e-6
   (t scaled to max |t| = 1). This gives 3153 candidates.
2. Of the 495636 pairs of deficient rank, none meets in exact arithmetic.
3. Each candidate is computed exactly in ℚ in its plane x_S − t_S = 0 (vertices, rays, edges).
   All 3153 are 2-dimensional.
4. The polygons are glued along their exact vertices, and the homology is computed from the
   oriented boundary matrices over ℚ.

These numbers agree with the earlier count of X18, which used max-slack LPs over the pairs
(face of the coarsest fan, S) (min slack 1e-4, no undecided pairs). The two methods are independent.

A run of the exact step on 4 Julia threads lost 2 unbounded polygons (3151). The cause was not
found. The single-thread run gives the counts above.

## The region atlas of Figure 15 against X18 (2026-10-30)

Figure 15 of main.tex (fig:planar_curve, `curve_portraits.jl --outer-face=1 --seed=0 --side-by-side --frames
--trivalent-nodes`) draws the planar core of the curve with a coloured region for each of the 5 maximal cells it runs
through: 35 and 177 (type 12), 43 and 181 (type 11), 68 (type 7). The regions come from the drawing alone: every face
is cut among the cells of its border. Two checks against the exact polygons of X18.

1. The regions as a picture of the surface: right along the core only.
   - Every polygon holding a core edge lies in the cell of that edge (ring 0 below). One, two and three polygons
     further, 65%, 28% and 0% of the area lies in the cells of the core.
   - The bounded part meets 83 maximal cells; the 5 of the core hold 47.5% of its area, cells 215 and 45 (type 29)
     and 146 (type 30) 29%.
   - Along every core edge 5 polygons meet, 4 or 5 of them bounded. Cut along the core, the bounded part falls into
     8 pieces: a small disc in each of the faces 3 and 6 to 11 of the drawing (1 to 3 polygons, 0.9% of the area
     together, in the cells of the core on their borders), and one piece of 702 polygons, 99.1% of the area, along
     all 40 core edges. The outer face and the faces 2, 4 and 5, the largest, have no piece of their own.
2. The regions as a map of the cells of the core: right, but not to scale.
   - The part of the bounded part in each cell of the core is connected and holds all the core edges of that cell,
     as each region does.
   - The cells of the core that meet along walls of X18 are exactly those whose regions meet: 35–68, 68–177,
     177–181, 181–43 and 43–35, a cycle. The other 5 pairs meet in neither.
   - The areas and the lengths of the borders are not those of X18: cell 35 holds 37% of the area of the parts but
     9% of the atlas, the cells 43 and 181 4% and 2% against 15% and 14%.

Areas and lengths are Euclidean in the slice sum(x) = 0. The faces are numbered as `portrait_faces` returns them, the
outer one 1.

| piece of the cut | polygons | area | core edges on its border | faces holding them all | area in the cells of the core |
|---|---:|---:|---:|---|---:|
| 1 | 702 | 99.1% | 40 | none | 47.0% |
| 2 | 1 | 0.3% | 7 | 10 | 100% |
| 3 | 3 | 0.2% | 10 | 7 | 100% |
| 4 | 2 | 0.1% | 6 | 8 | 100% |
| 5 | 1 | 0.1% | 4 | 6 | 100% |
| 6 | 1 | 0.1% | 4 | 3 | 100% |
| 7 | 2 | 0.0% | 5 | 9 | 100% |
| 8 | 1 | 0.0% | 4 | 11 | 100% |

| ring from the core | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|---:|
| polygons | 136 | 176 | 168 | 136 | 78 | 17 | 2 |
| area | 18.5% | 30.4% | 33.5% | 12.3% | 4.7% | 0.5% | 0.2% |
| in the cells of the core | 100% | 64.7% | 27.8% | 0% | 0% | 0% | 0% |

| cell | type | bounded polygons | area on X18 | area of the region |
|---:|---:|---:|---:|---:|
| 35 | 12 | 54 | 37.0% | 9.1% |
| 43 | 11 | 25 | 3.6% | 15.0% |
| 68 | 7 | 45 | 16.9% | 22.1% |
| 177 | 12 | 91 | 40.4% | 40.0% |
| 181 | 11 | 10 | 2.1% | 13.8% |

| cells | wall on X18 | border in the atlas |
|---|---:|---:|
| 68, 177 | 38.5% | 23.1% |
| 35, 68 | 36.2% | 24.2% |
| 35, 43 | 10.0% | 8.4% |
| 177, 181 | 7.7% | 23.4% |
| 43, 181 | 7.6% | 20.9% |

The polygons of X18 are found again from the curve alone, crossing the edges of the polygons from the sheets along
the curve: the same 3153 polygons, 713 bounded, with the f-vector above.

    julia --project=julia julia/studies/g2/curves/k3_surface_polygons.jl julia/curves_g10/overnight_20261028/wpart_2722
    julia --project=julia julia/studies/g2/curves/k3_core_regions.jl julia/curves_g10/overnight_20261028/wpart_2722

The first takes 16 minutes on one thread and writes julia/results/curves/genus_10/k3/wpart_2722_x18.mrdi; the second
prints the checks.

## The unbounded rays of X17 (2026-10-30)

main.tex says that the curves of fig:planar_curve and fig:nonplanar_curve, of degree 18 in R^21/1, have up to
21 × 18 = 378 unbounded rays. The bound is 15 × 18 = 270, counted with multiplicity: trop(G2/P2) lies in x12 = x67,
x13 = x57, x17 = x56, x23 = x47, x27 = x46 and x37 = x45 (checked on its rays and lineality, and no other pair), so a
ray points along e12 + e67, never along e12 alone. X17 reaches the bound.

| curve | rays | with multiplicity | multiplicities |
|---|---:|---:|---|
| X17, wpart 2722 (fig:planar_curve) | 239 | 270 | 1^208 2^31 |
| experiment 100128 (fig:nonplanar_curve) | 210 | 240 | 1^185 2^23 4^1 5^1 |

- X17: 15 directions, each of total multiplicity 18: e_i for the 9 coordinates 14, 15, 16, 24, 25, 26, 34, 35, 36,
  and e_ab + e_cd for the 6 pairs.
- 100128 (julia/curves_g10/overnight_20261028/experiment_100128): the 9 single directions of total multiplicity 18,
  the 6 pairs 12 each, and one more direction, the sum of the 12 paired coordinates, of total multiplicity 6:
  9 · 18 + 6 · 12 + 6 = 240.
- Both have degree 18 in each of the 21 coordinates (the sum over the rays of m (u_i − min u), u primitive) and are
  balanced at infinity.

X17 is recomputed exactly from translation_1.mrdi (2t), as in core_trivalence.jl, every piece of C ∩ (B_S + t) of
multiplicity m_C [Z^21 : L_C + L_S], with L_C = span(C) ∩ Z^21 and L_S = {x_S constant}. All 449 pieces are
1-dimensional, of multiplicities 1^408 2^41, the curve is balanced with them, and its 239 rays have the same vertices
and cells as in decorated_graph_1.mrdi. The certified graph there has 215 bounded segments against the 210 recomputed:
5 of its vertices have valence 2, with collinear edges, on the diagonal of the quadrilateral 146 of the coarsest fan,
a wall of the raw fan where the certified run cut the curve (checked 2026-11-10). Since then the stage zenodo merges
them, and julia/data/mukai_curves/genus_10_method_1.mrdi has 449 cells and 201 vertices. For 100128 the multiplicities are those that
hyperplane_intersection_census.jl stored in curve.mrdi.

    julia --project=julia -t 8 julia/studies/g2/curves/curve_rays.jl julia/curves_g10/overnight_20261028/wpart_2722
    julia --project=julia -t 8 julia/studies/g2/curves/curve_rays.jl julia/curves_g10/overnight_20261028/experiment_100128

The first takes 94 s, the second 54 s, both under 2.7 GB.
