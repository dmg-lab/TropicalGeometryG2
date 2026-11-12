# Functions in `TropicalCensus` to use instead of Oscar

This package was written for the tropical census of G2 and B3 homogeneous
spaces, but several of its functions are general. Each one below either
replaces an Oscar function that is too slow or too restrictive, or fills a gap
that Oscar 1.8 leaves. All arithmetic is exact over `QQ`.

```julia
# from the repository root
julia --project=julia
julia> using TropicalCensus, Oscar
```

Tested with Oscar 1.8.2 and gfan 0.8beta, on a 20-core machine with 16 GB of memory.

| Use | Instead of (Oscar) | Gain |
|---|---|---|
| [`gfan_prevariety`](#gfan_prevariety) | `tropical_prevariety` | Same fan in 0.8 s, where Oscar did not finish in 5 min |
| [`run_tropicalization`](#run_tropicalization) | `tropical_variety` | 2.3× faster on G2/P2; runs in parallel and can resume after a crash |
| [`coarsen_fan`, `coarsen_prevariety_fan`](#coarsen_fan-and-coarsen_prevariety_fan) | none: Oscar only returns refined structures | Coarsest polyhedral structure, with a certificate |
| [`canonical_initial_ideal`](#canonical_initial_ideal) | `initial` + `is_prime` + `saturation` + `degree` | One call, and detects ideals prime over `QQ` that split over `C` |
| [`tree_from_min_pluecker`](#tree_from_min_pluecker) | `phylogenetic_tree(M, taxa)` | Any tree metric, not only ultrametric ones |
| [`degree_of_embedding`](#degree_of_embedding) | `degree(I)` of the Plücker ideal of `G/P` | Closed formula, no Gröbner basis |
| [`determine_polytope`](#determine_polytope) | `Polymake.polytope.isomorphic(P, Q)` against candidates | Names the combinatorial type, recursively |

## `gfan_prevariety`

```julia
F = gfan_prevariety(polynomials; threads=Sys.CPU_THREADS, bits=64, symmetry=nothing)
```

Computes the tropical prevariety of polynomials over `QQ` with trivial
valuation and returns a `PolyhedralFan` in the min convention, i.e. the
convention of `tropical_semiring_map(QQ, min)`. It needs
`gfan_tropicalprevariety` on the `PATH`, or its location passed as `gfan=...`.

Oscar's `tropical_prevariety` intersects the tropical hypersurfaces one at a
time, using a single-threaded polymake call on every pair of cells. Gfan
instead intersects all the hypersurfaces at once, on several threads, with
machine integers.

Example: the 10 quadrics in 14 variables of one edge of the G2 hexagon.


| Polynomials | 2 | 3 | 4 | 5 | 10 |
|---|---|---|---|---|---|
| `tropical_prevariety` | 0.5 s | 3.4 s | 20 s | 60 s | not finished after 5 min |
| `gfan_prevariety` | | | | | 0.8 s |

On the first 3 polynomials both functions return the same 164 maximal cones.

- Gfan uses the max convention. `gfan_prevariety` negates the rays, so its
  result can be compared directly with Oscar's.
- Oscar returns a `PolyhedralComplex`; this returns a `PolyhedralFan`. With
  trivial valuation the cells are the same.
- With `bits=64`, gfan uses machine integers. If they overflow, the computation
  is retried automatically with arbitrary precision.
- `symmetry` takes a `PermGroup`, or a list of generators given as
  `PermGroupElem`s or permutation vectors of `1:nvars`, of a group permuting
  the variables. Gfan then computes only one maximal cone per orbit
  (`--symmetry`), and the orbits are expanded so that the full fan is still
  returned. The group must permute the tropical hypersurfaces of the
  polynomials; this is not checked.
- By default (`non_redundant=true`) polymake takes gfan's rays and maximal
  cones as they are, instead of spending minutes recomputing them on a large
  fan (15 minutes for 37,986 cones). This requires that the rays and the
  incidence matrix really describe a fan: distinct rays modulo lineality, each
  cone listed by exactly its extreme rays, no cone inside another, and cones
  meeting in common faces. Gfan's output satisfies this, but polymake does not
  check it; pass `non_redundant=false` to let polymake reduce the input itself.
- Gfan has no memory limit and can fill all the memory on large inputs. To
  keep it from freezing the machine, run Julia under a cap, for example
  `systemd-run --user --scope -p MemoryMax=8G -p MemorySwapMax=0 julia ...`.

## `run_tropicalization`

```julia
problem = TropicalProblem(
    name="G2/P2", slug="g2_p2",
    provider_file=study_path(:g2, "ideals", "G2GrassmannianIdeal.jl"),
    provider_module=:G2GrassmannianIdeal,
    builder=:build_g2_grassmannian_ideal,
    output_prefix="TGrass_G2",
    default_output="work/g2_p2",
)
run_tropicalization(problem, TropicalizationConfig(output="work/g2_p2", workers=2, threads=8))
fan = load("work/g2_p2/TGrass_G2_Fan.mrdi")
```

Computes the tropical variety of a homogeneous ideal over `QQ` with trivial
valuation. It returns `nothing` and writes its results to the `output`
directory. `<output_prefix>_Fan.mrdi` holds the fan, with cones merged into
Gröbner cones.

The ideal is read from a provider file: a module with a builder function that
returns a named tuple with fields `ring`, `ideal` and `generators`. Copy
`julia/studies/g2/ideals/G2GrassmannianIdeal.jl` as a template. The variable
names are read from `ring`. Only a provider that works in a reduced torus
chart, such as `build_flag_g2_ideal` in
`julia/studies/g2/ideals/FlagG2Ideal.jl`, also returns
`original_variable_names` and `lift_indices`: the fan is lifted back to the
original coordinates, the `i`-th taking the value of ring coordinate
`lift_indices[i]`. A provider may also return the Krull `dimension` of its
ideal when computing it is too expensive.
With the default certification it only decides which cells are probed first.

How it works, round by round:

1. Gfan computes the prevariety of the current list of polynomials.
2. Worker processes check every cone in parallel, using the initial ideal at an
   interior weight. A cone is either confirmed, or a polynomial that cuts it
   away is added to the list.
3. The loop ends when every cone is confirmed. The list is then a tropical
   basis.

Every round is saved to `output`, so after a crash `resume=true` continues
where the computation stopped.

With a symmetry group, the list is closed under the group and only one cone
per orbit is checked. Give it as gfan does, by generators permuting the
variables. A generator `p` sends variable `k` to variable `p[k]`. If the
variables also change sign, write it `p => signs`: it sends variable `k` to
`signs[k]` times variable `p[k]`, like gfan's `--symsigns`.

```julia
problem = TropicalProblem(
    name="G2/P2 in S_long", slug="g2_p2_slong",
    provider_file=study_path(:g2, "ideals", "G2GrassmannianSlongIdeal.jl"),
    provider_module=:G2GrassmannianSlongIdeal,
    builder=:build_g2_grassmannian_slong_ideal,
    output_prefix="TGrass_G2_Slong",
    default_output="work/g2_p2_slong",
    symmetry_candidates=:SLONG_WEYL_GENERATORS,   # rho and tau, which generate the Weyl group of G2
    symmetry_action=:variables,
)
```


- The workers close the generators to the whole group and check exactly that
  every element preserves the ideal.
- `gfan_symmetry=true` also gives the group to gfan (`--symmetry`). Gfan then
  computes one cone per orbit, and the orbits are expanded before the checks.
  On the hexagon path A-B-C this cut the prevariety from 71 to 15 minutes.
- `coarsen_prevariety=false` checks the orbits of gfan's cones without merging
  them first. With the default certification the merge only orders the checks,
  and for the 69,876 cones of that path it needed more than 7 GB.
- The default `symmetry_action=:taxa` instead takes the whole
  group of permutations of the taxa `1:n`, acting on Plücker coordinates
  `p_ij`, with signed images. The G2 and flag censuses use it.

Example: G2/P2 (21 variables, 42 generators), 2 workers, 8 gfan threads.

| | Time | Peak memory | Maximal cones |
|---|---|---|---|
| `tropical_variety(I, nu)` | 173 s | 2.7 GB | 251 |
| `run_tropicalization` | 76 s, including worker start-up | 6.2 GB | 251 |

Each worker is a full Julia and Oscar process, so count about 2 GB of memory
per worker.

## `coarsen_fan` and `coarsen_prevariety_fan`

```julia
result = coarsen_fan(fan)               # pure fans, e.g. a tropical variety
result = coarsen_prevariety_fan(fan)    # fans that may not be pure, e.g. a prevariety
result.fan, result.coarsest_certified, result.original_cells
```

`tropical_variety` returns the cells of a Gröbner fan, and
`tropical_prevariety` returns a common refinement. Both usually cut cells that
belong together, and Oscar has no function to merge them again. These
functions merge neighbouring cones whenever their union is convex, without
changing the support or the lineality space.

- Every merged cone is proved convex in exact rational arithmetic.
- `original_cells` lists, for each output cone, the input cones it is made of.
- `coarsest_certified == true` proves that no two output cones can be merged.
  `false` means only that this proof was not found.
- `open_equal_span_walls(fan)` returns the number of walls that might still be
  removable; `0` proves the fan is already coarsest.
- Otherwise `open_wall_certificate(fan)` decides the walls one by one, and the
  coarseners use it for `coarsest_certified`. Cones merged into a convex cone
  are connected through open walls, so it checks every set of at least two
  cones connected through them; when no such union is convex, no coarsening
  exists. It gives up (inconclusive) past `max_unions=100_000` sets in a
  component. On the prevariety of the 84 Pfaffians in `Equation_ab_surfaces`,
  the 216 open walls of the 5,640 coarsened cones form 72 triangles, and the
  288 unions checked certify it coarsest.
- `coarsen_prevariety_fan(fan; symmetry=generators)` takes coordinate
  permutations `p`, acting by `v -> v[p]`, of a group preserving the fan. The
  cones with a common linear span are coarsened together; with a symmetry, only
  one such set per orbit is coarsened, and its merges are carried to the rest
  of the orbit. Each generator is checked to permute the rays and the maximal
  cones; `transported_span_partitions` counts the sets that were carried.

Example: the prevariety of the hexagon edge above goes from 912 to 852 maximal
cones in 50 s, and the result is certified coarsest.

## `canonical_initial_ideal`

```julia
r = canonical_initial_ideal(I, w)       # w: weight vector, min convention
r.ideal, r.prime, r.geometrically_prime, r.toric, r.dimension, r.degree
```

Computes the initial ideal of `I` at `w`, then its saturation by the product of
the variables, primality, dimension and degree, all in one call. It also
answers a question Oscar does not: `is_prime` works over `QQ`, so it accepts
ideals that split over `C`. For a binomial initial ideal, `canonical_initial_ideal`
decides primality over `C` from the index of the exponent lattice, using a
Smith normal form.

```julia
R, (x, y) = polynomial_ring(QQ, [:x, :y])
r = canonical_initial_ideal(ideal(R, [x^2 - 2y^2]), [0, 0])
r.prime                    # true: prime over QQ, as Oscar's is_prime says
r.geometrically_prime      # false: x^2 - 2y^2 = (x - √2 y)(x + √2 y)
r.binomial_lattice_index   # 2
```

- `geometrically_prime` is `missing` when a prime initial ideal is not binomial.
- `toric` is `true` exactly when the initial ideal is binomial, prime over `C`
  and contains no monomial.
- `compute_degree=false` skips the degree.
- `compute_saturation=false` skips the saturation, but only together with
  `certified_tropical=true`, i.e. when `w` is already known to lie in the
  tropical variety.

## `tree_from_min_pluecker`

```julia
T = tree_from_min_pluecker(w)    # w: point of Trop Gr(2,n), min convention
internal_splits(T), pendant_lengths(T), tree_edge_lengths(T), tree_distance(T, i, j)
```

Reconstructs the phylogenetic tree whose distances are `-w`, with exact
rational edge lengths. The coordinates of `w` are ordered
`(1,2), (1,3), …, (n-1,n)`. The function first checks the tropical Plücker
relations, then rebuilds the tree by exact neighbour joining.

Oscar's `phylogenetic_tree(M, taxa)` goes through polymake, which rejects
every distance matrix that is not ultrametric ("The dissimilarity is not
ultrametric"). Most points of `Trop Gr(2,n)` are not ultrametric.

```julia
pendant = [1, 2, 3, 4, 5]        # caterpillar tree: splits 12|345 (length 1), 123|45 (length 2)
d(i, j) = pendant[i] + pendant[j] + ((i <= 2) != (j <= 2)) + 2((i <= 3) != (j <= 3))
T = tree_from_min_pluecker([-d(i, j) for i in 1:4 for j in i+1:5])
internal_splits(T)               # [12|345 of length 1, 45|123 of length 2]
```

## `degree_of_embedding`

```julia
degree_of_embedding(root_system(:G, 2), [0, 1])    # (5, 18)
```

Returns the dimension and degree of `G/P` in `P(V_ϖ)` for the dominant weight
`ϖ` (in fundamental-weight coordinates), using the Borel–Hirzebruch formula.
Oscar's `degree(I)` needs the ideal of `G/P` and a Gröbner basis of it, which
took 9 s for G2/P2 in a fresh session. Checked values: 2 for G2/P1 (a quadric),
18 for G2/P2, 2 for Gr(2,4).

The same module computes the Chevalley polytopes of Spacek–Wang
(arXiv:2411.10276), which have no Oscar counterpart:
`chevalley_polytope(chevalley_data(:G, 2, [0, 1]))`. Its source is
`studies/g2/degenerations_grassmannian/Spacek_Wang/chevalley_polytopes.jl`.

## `determine_polytope`

This module is not loaded by `using TropicalCensus`; include it directly:

```julia
include("julia/src/recursive_polytope_types.jl")
using .RecursivePolytopeTypes
describe_type(determine_polytope(cyclic_polytope(4, 7)))   # "cyclic polytope C(4, 7)"
describe_type(determine_polytope(pyramid(cube(2))))        # "Pyramid(4-gon)"
describe_type(determine_polytope(simplex(2) * simplex(2))) # "Δ2 × Δ2"
```

Oscar 1.8 has no combinatorial isomorphism test for polytopes, and polymake's
`isomorphic(P, Q)` only compares two polytopes you already have. This function
names the combinatorial type of a bounded polytope. It recognizes pyramids,
prisms, wedges over a face, duals and a catalog of named types, recursively.
Every decomposition is proved by rebuilding the full vertex–facet incidence.
`unrecognized_parts(P, r)` returns the pieces it could not name.
