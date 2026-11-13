# Timings of the tropicalization of the G2 ideals outside the census

The prevarieties are computed with Gfan's threaded intersection
`gfan_tropicalprevariety`, and the tropical G2 Grassmannian with Oscar's
`tropical_variety`. Oscar's `tropical_prevariety` is timed for comparison.
All timings in this table were measured on a laptop. Oscar's
`tropical_prevariety` was run to the end on a cluster node, see the last
section.

| Computation | Method | Time per run | Result |
| --- | --- | --- | --- |
| G2-Dressian | Gfan, 1 thread | 1.64 s, 1.63 s, 1.63 s | 51 rays, 239 maximal cones |
| G2-Dressian | Gfan, 8 threads | 0.68 s, 0.68 s, 0.70 s | 51 rays, 239 maximal cones |
| G2 flag Dressian | Gfan, 1 thread | 34.4 s, 34.3 s | 147 rays, 2726 maximal cones |
| G2 flag Dressian | Gfan, 8 threads | 9.8 s, 9.8 s | 147 rays, 2726 maximal cones |
| G2 Grassmannian | `tropical_variety`, default keywords | 160 s, 142 s | 51 rays, 251 maximal cones |
| G2 Grassmannian | `tropical_variety` with `skip_decomposition=true`, the census call | 162 s, 143 s | 51 rays, 251 maximal cones |
| G2-Dressian | Oscar's `tropical_prevariety` | 777 s, 770 s, 770 s | 51 rays, 239 maximal cones |
| G2 flag Dressian | Oscar's `tropical_prevariety` | stopped after 89 min, not finished; finished on a cluster node, see below | |

The G2-Dressian and the G2 Grassmannian use the 42 generators of the
Grassmannian ideal in 21 variables. The G2 flag Dressian uses the 78
generators of the flag ideal in all 28 variables. The G2 flag variety is not
timed: `tropical_variety` takes too long on it.

## Setup

- Gfan 0.8beta with `--bits 64`, Oscar 1.8.2 and Julia 1.13.1 with one thread,
  on a laptop with a 12th Gen Intel Core i7-12700H and 16 GB of memory.
- The ideals are over `QQ` with the trivial valuation. The results are in the
  min convention.
- The Gfan times are those of `gfan_tropicalprevariety` alone. The Oscar times
  are those of the Oscar call alone.
- A small example, the Pluecker ideal of Gr(2,5), runs first, so that
  compilation is not counted. Each run starts from a newly built ideal, so
  that no Groebner basis is reused.
- The numbers of rays and of maximal cones are those of the census fans
  `../../../data/g2/dressian/01_prevariety_fan.mrdi`,
  `../../../data/g2/flag_dressian/01_prevariety_fan.mrdi` and
  `../../../data/g2/grassmannian/01_raw_fan.mrdi`.

## Reproduction

From `julia/`, one case per process. The arguments are the case,
the number of runs and the number of Gfan threads:

```sh
julia --project=. studies/g2/supplemental/benchmark_tropicalization.jl dressian_gfan 3 1
julia --project=. studies/g2/supplemental/benchmark_tropicalization.jl dressian_gfan 3 8
julia --project=. studies/g2/supplemental/benchmark_tropicalization.jl flag_dressian_gfan 2 1
julia --project=. studies/g2/supplemental/benchmark_tropicalization.jl flag_dressian_gfan 2 8
julia --project=. studies/g2/supplemental/benchmark_tropicalization.jl grassmannian 2
julia --project=. studies/g2/supplemental/benchmark_tropicalization.jl grassmannian_skip_decomposition 2
julia --project=. studies/g2/supplemental/benchmark_tropicalization.jl dressian 3
julia --project=. studies/g2/supplemental/benchmark_tropicalization.jl flag_dressian 1
```

The script prints the timings and writes nothing to `data/`.

## Oscar's `tropical_prevariety` on a cluster node

Hardware: a cluster node with 2 x Intel Xeon Gold 6226R at 2.90 GHz (32 cores,
64 threads), 754 GB of memory. Oscar 1.8.2, Julia 1.13.1, one thread per
process. Ten single-threaded processes ran at the same time (12 CPUs, 96 GB),
each with the warm-up on Gr(2,5) first.

`tropical_prevariety` takes no keywords (`skip_saturation` and
`skip_decomposition` belong to `tropical_variety`): it is
`reduce(common_refinement, ...)` over the tropical hypersurfaces, in the order
given. Polymake's `common_refinement` intersects every pair of maximal cells,
each intersection a polymake cone with its own convex hull, 20 ms to 55 ms
per pair here. What can be varied from outside is the order of the
generators, the coordinates and polymake's convex hull backend.

| Computation | Variant | Time | Result |
| --- | --- | --- | --- |
| G2-Dressian | input order, ppl (polymake's default) | 1262 s | 51 rays, 239 maximal cones |
| G2-Dressian | input order, lrs | 1424 s | 51 rays, 239 maximal cones |
| G2-Dressian | input order, cdd | 1490 s | 51 rays, 239 maximal cones |
| G2-Dressian | greedy order | 886 s | 51 rays, 239 maximal cones |
| G2 flag Dressian, 71 generators, 22 coordinates | input order | 17,578 s (4 h 53 min) | 147 rays, 2726 maximal cones |
| G2 flag Dressian, 71 generators, 22 coordinates | greedy order | 15,634 s (4 h 21 min) | 147 rays, 2726 maximal cones |

- Input order: linear forms, Pluecker relations, incidence relations, quadric
  q. Greedy order: the linear forms, then always the generator with the
  fewest variables not met so far, the two variables of a binomial linear form
  counting as one. The 22 coordinates are those left by the six binomial
  linear forms, as in the Gfan computation of the census; two Pluecker
  relations become equal there, leaving 71 distinct generators.
- Peak memory was 2.2 GiB to 2.8 GiB in every run.

Reproduction, from `julia/` on a Slurm cluster (logs in
`work/oscar_prevariety/<job>/`):

```sh
mkdir -p work/oscar_prevariety
sbatch studies/g2/supplemental/benchmark_oscar_prevariety.slurm
julia --project=. studies/g2/supplemental/benchmark_oscar_prevariety.jl --case=flag72 --order=greedy --mode=steps
```
