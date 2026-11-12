# G2 census: what the paper states

Computed 2026-10-15 by the census stage `summarize` from the stored census
artifacts. Each row links to the `08_summary.toml` of its target, which has
all the numbers below and more: the polygons among the 2-faces, the counts of
prime and binomial orbits, and the size of the stage-1 fan.

All counts are modulo lineality, where the fan is a polyhedral complex: a
vertex is a ray and a face of dimension k is a cone of dimension k + 1. The
f-vector starts with the vertices. Orbits are those of the Weyl group D12 of
order 12. Homology is reduced and integral. It has no torsion in any case.

| case | dim | f-vector | orbit f-vector | maximal cells | homology | toric degenerations |
|---|---:|---|---|---|---|---|
| [G2-Dressian](dressian/08_summary.toml) | 3 | (51, 241, 248, 3) | (11, 33, 31, 1) | 233 triangles, 3 quadrilaterals, 3 tetrahedra | ℤ⁵⁴ in dim 2 | – |
| [G2 flag Dressian](flag_dressian/08_summary.toml) | 4 | (141, 1352, 3584, 2729, 33) | (27, 173, 388, 274, 6) | 2468 tetrahedra, 108 square pyramids, 18 triangular prisms, 33 4-simplices | ℤ³²⁴ in dim 3 | – |
| [tropical G2-Grassmannian, char 0](grassmannian/08_summary.toml) | 2 | (51, 241, 248) | (11, 33, 31) | 245 triangles, 3 quadrilaterals | ℤ⁵⁷ in dim 2 | 135 in 17 orbits |
| [tropical G2-Grassmannian in S_long, char 0](grassmannian_Slong/08_summary.toml) | 2 | (42, 211, 275) | (11, 31, 35) | 263 triangles, 12 quadrilaterals | ℤ¹⁰⁵ in dim 2 | 21 in 5 orbits |
| [tropical G2 flag variety, char 0](flag/08_summary.toml) | 3 | (141, 1352, 3584, 2714) | (27, 173, 388, 270) | 2588 tetrahedra, 108 square pyramids, 18 triangular prisms | ℤ³⁴² in dim 3 | 807 in 81 orbits |
| [tropical G2-Grassmannian, char 2](grassmannian_gf2/08_summary.toml) | 2 | (54, 253, 254) | (12, 36, 32) | 251 triangles, 3 quadrilaterals | ℤ⁵⁴ in dim 2 | 147 in 20 orbits |
| [tropical G2 flag variety, char 2](flag_gf2/08_summary.toml) | 3 | (156, 1445, 3758, 2792) | (31, 192, 416, 281) | 2630 tetrahedra, 144 square pyramids, 18 triangular prisms | ℤ³²⁴ in dim 3 | 825 in 83 orbits |

The maximal cells of the two Dressians are not all of top dimension. In the
flag Dressian, the 2729 faces of dimension 3 are the 2594 maximal ones listed
above and the 135 facets of the 33 4-simplices.

## Characteristic 2 compared with characteristic 0

Each maximal cell over GF(2) is looked up among the maximal cells over QQ,
modulo lineality.

| | maximal cells in both | toric in both | toric only in char 2 | toric only in char 0 |
|---|---|---|---|---|
| Grassmannian | 236 in 28 orbits | 135 in 17 orbits | 12 in 3 orbits | none |
| flag variety | 2594 in 253 orbits | 807 in 81 orbits | 18 in 2 orbits | none |

## How it is computed

```bash
julia --project=julia julia/studies/g2/run_census.jl --target=all --stage=summarize
julia --project=julia julia/studies/g2/run_census_gf2.jl --target=all --stage=summarize
```

The stage reads the normalized coarsest fan, the symmetry artifacts and the
stage-6 census of each target. It examines one maximal cell per orbit and
multiplies by the orbit sizes. Faces come from the facets of each
representative's cone. Homology comes from the simplicial complex on the
vertex sets of the maximal cells, which has the same nerve, and is checked
against the f-vector through the Euler characteristic. The code is in
`julia/src/summary.jl`; see `julia/studies/g2/README.md`.
