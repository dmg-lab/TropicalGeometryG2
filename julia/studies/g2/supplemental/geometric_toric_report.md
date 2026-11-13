# Geometric primality of the canonical G2 degenerations

The stage-6 `prime` and `binomial` columns are computed over `QQ`. For the
canonical interior weights, the toric classification is unchanged over `CC`:

| Variety | Maximal-cell orbits | Prime over `QQ` | Prime over `CC` | Toric over `CC` |
| --- | ---: | ---: | ---: | ---: |
| `G2/P2` Grassmannian | 31 | 18 | 17 | 17 |
| `G2/B` flag | 270 | 81 | 81 | 81 |

All 98 `QQ`-prime binomial ideals have saturated exponent-difference lattices:
the nonzero Smith normal form entries are all 1. Since these ideals avoid the
coordinate hyperplanes, they remain prime binomial ideals over `CC`.

The only `QQ`-prime nonbinomial orbit is Grassmannian orbit 31 (representative
cell 155, cell orbit `{155, 156}`). OSCAR's
`absolute_primary_decomposition` gives one conjugacy class of extension degree
2 over `QQ(sqrt(-3))`. Thus the ideal has **two** conjugate prime components
over `CC`. Each component has degree 9 and is nonbinomial in the Pluecker
coordinates. The original ideal has degree 18. Its reduced Groebner basis has
three-term polynomials over `QQ`; this remains true after extension to `CC`.
The older fine Groebner atlas numbers the corresponding two triangles 158 and
159; those are the two `QQ`-prime nonbinomial entries in that atlas.

Reproduce the geometric check from julia/ with:

```sh
julia --project=. studies/g2/supplemental/verify_geometric_toric.jl
```

The Smith normal form criterion is the lattice-ideal primality criterion of
[Eisenbud and Sturmfels, *Binomial Ideals*](https://arxiv.org/abs/alg-geom/9401001)
(1996). The script verifies every
stage-6 toric orbit and computes the absolute primary decomposition of the
exceptional orbit directly.
