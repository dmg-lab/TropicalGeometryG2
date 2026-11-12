# Timings of the tropicalization of the G2 ideals outside the census

The prevarieties are computed with Gfan's threaded intersection
`gfan_tropicalprevariety`, and the tropical G2 Grassmannian with Oscar's
`tropical_variety`. Oscar's `tropical_prevariety` is timed for comparison.
All timings in this table were measured on 2026-10-12 on a laptop. Oscar's
`tropical_prevariety` was run to the end on a hydra node on 2026-11-07, see
the last section.

| Computation | Method | Time per run | Result |
| --- | --- | --- | --- |
| G2-Dressian | Gfan, 1 thread | 1.64 s, 1.63 s, 1.63 s | 51 rays, 239 maximal cones |
| G2-Dressian | Gfan, 8 threads | 0.68 s, 0.68 s, 0.70 s | 51 rays, 239 maximal cones |
| G2 flag Dressian | Gfan, 1 thread | 34.4 s, 34.3 s | 147 rays, 2726 maximal cones |
| G2 flag Dressian | Gfan, 8 threads | 9.8 s, 9.8 s | 147 rays, 2726 maximal cones |
| G2 Grassmannian | `tropical_variety`, default keywords | 160 s, 142 s | 51 rays, 251 maximal cones |
| G2 Grassmannian | `tropical_variety` with `skip_decomposition=true`, the census call | 162 s, 143 s | 51 rays, 251 maximal cones |
| G2-Dressian | Oscar's `tropical_prevariety` | 777 s, 770 s, 770 s | 51 rays, 239 maximal cones |
| G2 flag Dressian | Oscar's `tropical_prevariety` | stopped after 89 min, not finished; finished on hydra, see below | |

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

## Observations

- On one thread Gfan is about 470 times faster than Oscar's
  `tropical_prevariety` on the G2-Dressian. On the G2 flag Dressian the Oscar
  call was stopped after 89 minutes, with 2.7 GB in use; Gfan takes 34 s.
- `TropicalCensus.gfan_prevariety` runs Gfan, reads its output and returns an
  Oscar fan. The whole call took 2.0 s and 1.0 s on the G2-Dressian, with 1
  and 8 threads, and 34 s to 35 s and 10.0 s on the G2 flag Dressian. It
  checks once per second whether Gfan has finished.
- The primary decomposition makes no measurable difference on the
  Grassmannian: the two variants agree within the variation between runs.
- The peak memory of the Julia process was about 2.5 GiB in every finished
  run. Most of it is Oscar itself.
- `../../../data/g2/flag/provenance/full28_benchmark.toml` records 173 s as
  `initial_intersection_seconds` for the same 78 generators with 8 threads and
  64 bits. Gfan alone takes 9.8 s here. The two were not measured under the
  same conditions, and what the recorded value includes has not been
  determined.

## Reproduction

From the repository root, one case per process. The arguments are the case,
the number of runs and the number of Gfan threads:

```sh
julia --project=julia julia/studies/g2/supplemental/benchmark_tropicalization.jl dressian_gfan 3 1
julia --project=julia julia/studies/g2/supplemental/benchmark_tropicalization.jl dressian_gfan 3 8
julia --project=julia julia/studies/g2/supplemental/benchmark_tropicalization.jl flag_dressian_gfan 2 1
julia --project=julia julia/studies/g2/supplemental/benchmark_tropicalization.jl flag_dressian_gfan 2 8
julia --project=julia julia/studies/g2/supplemental/benchmark_tropicalization.jl grassmannian 2
julia --project=julia julia/studies/g2/supplemental/benchmark_tropicalization.jl grassmannian_skip_decomposition 2
julia --project=julia julia/studies/g2/supplemental/benchmark_tropicalization.jl dressian 3
julia --project=julia julia/studies/g2/supplemental/benchmark_tropicalization.jl flag_dressian 1
```

The script prints the timings and writes nothing to `data/`.

## Oscar's `tropical_prevariety` on a hydra node (2026-11-07)

Hardware: hydra node gpu04 (partition `gpu`, no GPU used; the `cpu` partition
was down), 2 x Intel Xeon Gold 6226R at 2.90 GHz (32 cores, 64 threads),
754 GB of memory. Oscar 1.8.2, Julia 1.13.1, one thread per process. Ten
single-threaded processes ran at the same time (Slurm job 378795, 12 CPUs,
96 GB), each with the warm-up on Gr(2,5) first.

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
| G2-Dressian | input order, beneath_beyond | 1910 s | 57 rays, 245 maximal cones: wrong |
| G2-Dressian | greedy order | 886 s | 51 rays, 239 maximal cones |
| G2 flag Dressian, 78 generators, 28 coordinates | input order, the plain call | running | |
| G2 flag Dressian, 78 generators, 28 coordinates | input order, step by step | running: step 66 of 78 after 17,436 s | |
| G2 flag Dressian, 78 generators, 28 coordinates | greedy order | running: step 70 of 78 after 16,880 s | |
| G2 flag Dressian, 71 generators, 22 coordinates | input order | 17,578 s (4 h 53 min) | 147 rays, 2726 maximal cones |
| G2 flag Dressian, 71 generators, 22 coordinates | greedy order | 15,634 s (4 h 21 min) | 147 rays, 2726 maximal cones |

- Input order: linear forms, Pluecker relations, incidence relations, quadric
  q. Greedy order: the linear forms, then always the generator with the
  fewest variables not met so far, the two variables of a binomial linear form
  counting as one. The 22 coordinates are those left by the six binomial
  linear forms, as in the Gfan computation of the census; two Pluecker
  relations become equal there, leaving 71 distinct generators.
- Peak memory was 2.2 GiB to 2.8 GiB in every run.
- The incidence relations take the time. With 22 coordinates and the greedy
  order: Pluecker relations 691 s, the 35 incidence relations 13,281 s, the
  quadric 1193 s. A step costs about the number of cells times the number of
  terms in intersections; incidence steps grow from 23 s to 520 s as the
  cells grow from 717 to 2627.
- Greedy is 30% faster than the input order on the G2-Dressian and 11%
  faster on the flag Dressian in 22 coordinates. Ordering by leaves 1..7 (largest
  last) was worse in a local test: in 22 coordinates its first relations share
  few variables and the cells tripled at every step.
- ppl, polymake's default, is the fastest backend. beneath_beyond gives a
  different fan from Gfan and the other three backends (not investigated).
- With all 77 other generators, the fold has 2627 maximal cells; the quadric
  q then splits 99 of them, giving 2726. Checked with Gfan: the prevariety of
  the 77 generators other than q is exactly the coarsest flag Dressian
  `../../../data/g2/flag_dressian/02_coarsest_fan.mrdi` (same lineality, rays
  and cells). So q does not change the support, and the 99 pairs merged in the
  coarsening come from q.
- Gfan takes 34 s on one thread on the laptop: about 450 times less than the
  fastest Oscar run.
- The G2-Dressian takes 1262 s here and 770 s on the laptop: the node is older
  and ran ten processes at once.

Reproduction, from `julia/` on hydra (logs in `work/oscar_prevariety/<job>/`):

```sh
sbatch --exclude=gpu02 studies/g2/supplemental/benchmark_oscar_prevariety.slurm
julia --project=. studies/g2/supplemental/benchmark_oscar_prevariety.jl --case=flag72 --order=greedy --mode=steps
```

gpu02 was excluded because its `/tmp` was full and reading the home
directory hung there on 2026-11-07.
