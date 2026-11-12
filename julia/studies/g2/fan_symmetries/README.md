# Fan symmetries

This study tests whether the order-12 taxon symmetry of the coarsest
`Trop(G2/P2)` fan is dihedral and whether its abstract reflections are literal
linear reflections of the ambient Pluecker space `R^21`.

Run

```sh
julia --project=julia julia/studies/g2/fan_symmetries/analyze_dihedral_action.jl
```

The script uses exact rational arithmetic, verifies the action on the fan,
isolates the canonical `G2` plane outside lineality, detects its finite
`D4^(3)` (twisted affine G2) root window from cone incidence, and regenerates
[`REPORT.md`](REPORT.md) and [`projected_rays.csv`](projected_rays.csv).
