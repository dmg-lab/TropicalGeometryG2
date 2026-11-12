# Tropicalizing homogeneous spaces (Section 8)

The computations of Section 8 of main.tex. Run them from the repository root.

| File | What it checks | Output |
|---|---|---|
| [`tropicalizing_reductive_groups_tits.jl`](tropicalizing_reductive_groups_tits.jl) | The Tits lifts: four signed lifts in G2 of the rotation and of the reflection each; with the signed diagonal elements they generate N(T)(ZZ), of order 48. | [`tropicalizing_reductive_groups_tits.txt`](tropicalizing_reductive_groups_tits.txt) |
| [`tropicalizing_reductive_groups.jl`](tropicalizing_reductive_groups.jl) | Forgetting p_36 keeps 227 of the 248 maximal cells of trop(G) and contracts the others; section 6 computes the initial ideals in the weight basis on the 227 cells. Exact over QQ. | [`tropicalizing_reductive_groups.txt`](tropicalizing_reductive_groups.txt) |
| [`long_root_check.jl`](long_root_check.jl) | No Plücker coordinate in any basis of K^7 restricts to a long root line on G: the Pfaffian ideal is the unit ideal for p14 − p25, and for p14 its only zero is a = 0. | printed |

```sh
julia --project=julia julia/studies/g2/homogeneous_spaces/tropicalizing_reductive_groups_tits.jl
julia --project=julia julia/studies/g2/homogeneous_spaces/tropicalizing_reductive_groups.jl julia/data/g2/grassmannian
julia --project=julia julia/studies/g2/homogeneous_spaces/long_root_check.jl
```

`tropicalizing_reductive_groups.jl` reads `00_ideal.mrdi`, `02_coarsest_cells.csv`, `04_trees.csv` and
`06_canonical_interior_degenerations.csv` of `julia/data/g2/grassmannian/`, and the labels of the paper from
`../atlas/flag/local_flag_cell_atlas.csv` and `../atlas/grassmannian/ray_label_orbits.csv`; it takes about 20 s. The
Tits lifts need only the standard library and take a few seconds, the long-root check a few seconds after loading
Oscar; it stops with an error if either ideal is not the expected one.

The tropicalizations in the bases S_long and S_long ∪ S_short are targets of
[`../run_census_slong.jl`](../run_census_slong.jl), written to `julia/data/g2/grassmannian_Slong/`.
