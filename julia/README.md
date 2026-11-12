# TropicalCensus Julia project

This directory is both the Julia environment and the package `TropicalCensus`. The stored data live below `data/`;
the studies below `studies/` read and write them through the path helpers of the package (`data_path`,
`study_path`), so scripts can be run from any working directory. The commands of the repository README are run
from the repository root with `julia --project=julia`.

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia julia/make_paper_data.jl --list
```

`make_paper_data.jl` fills `data/` by running the census below, one Julia process per
step; the Zenodo tarball holds the same files (see the repository README).

## The census

`studies/g2/run_census.jl` (over QQ) and `studies/g2/run_census_gf2.jl` (over GF(2)) compute the targets `dressian`,
`flag_dressian`, `base` (G2/P2) and `flag` (G2/B) in stages: `tropicalize`, `coarsen`, `normalize`, `trees`,
`symmetry`, `degenerations`, `projections` (flags only) and `summarize`. Each target has one directory
`data/g2/<target>/` with its numbered stage files, its coarsening report `02_coarsest_report.md`, its summary
`08_summary.toml` and the provenance of each stage, `provenance/<stage>.toml`; see `data/README.md`. The stage
`zenodo` of the target `zenodo` writes the remaining files of the Zenodo tarball.

```sh
julia --project=julia julia/studies/g2/run_census.jl --target=base --stage=normalize,trees,symmetry
julia --project=julia julia/studies/g2/run_census.jl --target=flag --from=coarsen --through=degenerations
julia --project=julia julia/studies/g2/run_census.jl --target=zenodo --stage=zenodo
```

Other options are `--engine=auto|oscar|parallel`, `--resume`, `--force`, `--workers`, `--threads`, `--bits` and
`--no-pdf`. Oscar's `tropical_variety` is the default engine for G2/P2; G2/B uses an exact, resumable parallel
engine with Gfan and Oscar, whose checkpoints go to `work/`.

The stages 3–7 use the coarsest fan structure of stage 2. For each maximal cell C, the degeneration is computed at
the sum of the normalized rays of C, once per orbit of the symmetry group, after checking that the group acts on the
ideal and on the weights modulo lineality.

## Layout

- `src/`: the package: tropicalization, exact coarsening, trees, symmetry, degenerations, Chevalley polytopes,
  provenance and atlas support. `src/README.md` lists the functions that replace Oscar ones.
- `studies/g2/`: G2/P2 and G2/B: the census drivers, the ideals, the characteristic-p certificates (`char-p/`),
  the multiplicities (`weights/`), the atlases of trees (`atlas/`), the toric degenerations
  (`degenerations_grassmannian/`), Section 8 (`homogeneous_spaces/`), the curves of genus ten (`curves/`) and
  supplemental checks. See `studies/g2/README.md`.
- `studies/hexagon/`: the completion of hexagons and the certificate of Theorem 9.1.
- `to_tex/`: the scripts that write the figures and tables of the paper from the stored data; each prints its TeX
  or writes it to `--output=FILE`.
- `data/`: the data, written by `make_paper_data.jl` or unpacked from the Zenodo tarball; git keeps only the inputs,
  listed in `data/README.md`.
- `mukai_curves/`, `intersection_lib.jl`: the intersection of trop(G2/P2) with translated tropical linear spaces;
  `curves_g10/`: the runs that gave the two curves of the paper; the stage `zenodo` computes the curve of Method 1
  again from its shift vector q and takes that of Method 2 from its run.
- `scripts/`: the installer of Gfan and the generator of the uniform tropical linear space.
