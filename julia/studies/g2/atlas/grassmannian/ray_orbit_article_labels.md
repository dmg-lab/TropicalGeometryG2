# Grassmannian ray-orbit labels for the article

The data orbit numbers 1 to 11 are the `orbit` values in [`ray_label_orbits.csv`](ray_label_orbits.csv). Orbit 12 holds the rays 52, 53, 54, which occur only over GF(2). The Greek labels are assigned in catalog order. They do not renumber the data or the ray generators.

| Article label | Data orbit | Size | Representative ray | Tree type |
|:--------------|-----------:|-----:|-------------------:|:----------|
| $\alpha_1$ | 3 | 6 | 3 | Split |
| $\alpha_2$ | 5 | 3 | 8 | Split |
| $\alpha_3$ | 4 | 6 | 6 | Split |
| $\alpha_4$ | 6 | 3 | 9 | Split |
| $\alpha_5$ | 9 | 2 | 15 | Split |
| $\beta_1$ | 1 | 12 | 1 | Caterpillar, taxon 7 at an end |
| $\beta_2$ | 2 | 6 | 2 | Caterpillar, taxon 7 at an end |
| $\gamma_1$ | 7 | 6 | 13 | Caterpillar, taxon 7 in the middle |
| $\gamma_2$ | 8 | 3 | 14 | Caterpillar, taxon 7 in the middle |
| $\gamma_3$ | 10 | 3 | 18 | Caterpillar, taxon 7 in the middle |
| $\delta$ | 11 | 1 | 37 | Snowflake |
| $\epsilon$ | 12 | 3 | 52 | Snowflake, taxon 7 on an arm; characteristic 2 only |

The catalog order is **$\alpha_1,\ldots,\alpha_5$, $\beta_1,\beta_2$, $\gamma_1,\gamma_2,\gamma_3$, $\delta$, $\epsilon$**. The canonical order and labels used by the catalog builder are `PAPER_RAY_CATALOG_ORDER` and `PAPER_RAY_ARTICLE_LABELS` in [`build_styled_ray_orbit_tables.jl`](../../../../to_tex/build_styled_ray_orbit_tables.jl).
