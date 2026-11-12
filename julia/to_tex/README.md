# TeX fragments for the paper

These scripts write figures and tables of `main.tex` from the stored data.
Run them from the repository root with `--project=julia`. Each prints its
TeX, or writes it to `--output=FILE`.

The tables have the booktabs form, and `main.tex` loads booktabs for them.
`\toprule` and `\bottomrule` frame a table, and a single `\midrule`
separates its heading from its body; a table without heading has no
`\midrule`. Rows end in `\tabularnewline`, and `\arraystretch` is `0.9`.
A table of fixed width is a `tabular*` of `\linewidth` with
`@{\extracolsep{\fill}}`, whose column gaps stretch to the width, and the
heading of a right-aligned column is centred by `\multicolumn{1}{c}{...}`.

| Script | Output |
|---|---|
| [`render_maximal_cell_tree.jl`](render_maximal_cell_tree.jl) | The G2-tree of a cone of the tropical G2-Grassmannian over QQ or GF(2). |
| [`render_simple_tree.jl`](render_simple_tree.jl) | The same trees in a plain style: taxa as numbers, an edge of length `n` as `n` parallel lines. |
| [`build_styled_ray_orbit_tables.jl`](build_styled_ray_orbit_tables.jl) | The catalog of the twelve ray orbits of the tropical G2-Grassmannian, the last of which occurs only in characteristic 2. |
| [`cubic_tree_table.jl`](cubic_tree_table.jl) | The three cubics of the tropical basis, with the simple tree of a ray over GF(2) under each monomial; [`cubic_tree_table.sh`](cubic_tree_table.sh) runs it. |
| [`toric_degeneration_table.jl`](toric_degeneration_table.jl) | The table of Newton--Okounkov polytopes of the toric degenerations of the G2-Grassmannian in characteristics 0 and 2. |
| [`toric_degeneration_tree_table.jl`](toric_degeneration_tree_table.jl) | The G2-trees of the maximal cells that yield these toric degenerations, with the labels and in the order of that table; with `--non-toric`, those of the other maximal cells. |
| [`cell_ray_orbit_table.jl`](cell_ray_orbit_table.jl) | The ray orbits of each of the 35 orbits of maximal cells, numbered as in the two tables of trees. |
| [`square_cell_figure.jl`](square_cell_figure.jl) | A quadrilateral cell over QQ: the G2-trees of its four rays in circles at the corners of a square, and the G2-tree of the cell in its middle. |
| [`triangle_cell_figure.jl`](triangle_cell_figure.jl) | A triangular cell, by default of orbit 21: the G2-trees of its three rays in circles at the corners of a triangle and the G2-tree of the cell in its middle, without lengths, with the inner vertices of the rays lettered in all four trees. |
| [`lettered_square_figure.jl`](lettered_square_figure.jl) | The quadrilateral cell as the triangle figure draws a triangle: every ray its own tree with lettered inner vertices, the cell in the middle, and the hook for drawing over it. |
| [`decorated_trees_figure.jl`](decorated_trees_figure.jl) | Three decorated G2-trees: the star tree of the lineality space with an apple on its pendant edge to taxon 1 and, of the triangle of orbit 21, the ray α3 with an apple at its inner vertex A and the tree of the cell with an apple on its edge AB, letters as in the triangle figure but not drawn. |
| [`toric_degeneration_apple_trees.jl`](toric_degeneration_apple_trees.jl) | A standalone page of the trees of `toric_degeneration_tree_table.jl`, in color and without edge lengths, with an apple on every edge whose decorated tree is a prime cone of the tropical G2 flag variety; `--pdf` compiles it to `output/pdf/g2_toric_degeneration_apple_trees.pdf`. |
| [`toric_degeneration_prime_edge_trees.jl`](toric_degeneration_prime_edge_trees.jl) | The figure of `figures/toric_degeneration_trees.tex` in black, with the edges that have a prime cone of the tropical G2 flag variety over them, and their lengths, in green, the others in red, and the number of prime cones and a green apple at the top right of every tree; it checks that the rest is that figure. With `--fragment` it is the figure of main.tex (`figures/toric_degeneration_prime_edge_trees.tex`), with `--standalone` that figure as the document of `standalone_figure.jl`, whose PDF main.tex includes, otherwise a standalone page with its caption (`figures/toric_degeneration_prime_edge_trees_standalone.tex`). |
| [`standalone_figure.jl`](standalone_figure.jl) | A figure of main.tex as a standalone document: the preamble of `main.tex` in the class `standalone`, cropped to the figure at its size in main.tex. The `--standalone` of the other scripts, and the document of a hand-drawn figure. |
| [`local_flag_complex_table.jl`](local_flag_complex_table.jl) | The isomorphism types of the local flag complexes above the 35 orbits of maximal cells of both fields, with their maximal cells, f-vectors and orbits. |

## G2-trees

`render_maximal_cell_tree.jl` draws the tree of one cone of the coarsest fan,
or of one generator of its lineality space. Select it with exactly one of:

| Option | Cone |
|---|---|
| `--lineality=K` | Generator `K` of the lineality space, `K` in `1:4`: a star tree. |
| `--ray=K` | Ray `K`. |
| `--edge=I,J` | The two-dimensional cone spanned by rays `I` and `J`. |
| `--cell=INDEX` | The maximal cell with this index in the census of the field. |
| `--rays=I,J,K[,L]` | The maximal cell with these rays: a triangle or a square. |

`--field=QQ` (the default) or `--field=GF2` selects the fan.

| Fan | Rays | Edges | Maximal cells |
|---|---:|---:|---:|
| QQ | 51 | 241 | 248 |
| GF(2) | 54 | 253 | 254 |

```sh
julia --project=julia julia/to_tex/render_maximal_cell_tree.jl \
  --lineality=1 --color --label=true --output=/tmp/g2_lineality_1.tex

julia --project=julia julia/to_tex/render_maximal_cell_tree.jl \
  --ray=52 --field=GF2 --color --label=true --output=/tmp/g2_ray_52.tex

julia --project=julia julia/to_tex/render_maximal_cell_tree.jl \
  --edge=13,14 --color --label=true --output=/tmp/g2_edge.tex

julia --project=julia julia/to_tex/render_maximal_cell_tree.jl \
  --rays=13,14,17,18 --color --label=true --scale=1.4 \
  --bounded-length=1.0 --edge-labels=aligned --output=/tmp/g2_square.tex
```

### Generators of the lineality space

The lineality space is the same over both fields. Its points are the star
trees, with no bounded edge, and a generator is drawn as the star tree
labelled by its seven pendant lengths, which may be negative.

| `K` | Pendant lengths of taxa 1, ..., 7 |
|---|---|
| 1 | `1, -1, 0, -1, 1, 0, 0` |
| 2 | `-1, 0, 1, 1, 0, -1, 0` |
| 3 | `0, 1, -1, 0, -1, 1, 0` |
| 4 | `1/2` for every taxon; the Pluecker vector is the all-ones vector |

Generators 1, 2 and 3 are the three short coroots of the maximal torus of G2
in `main.tex`, and 3 is minus the sum of 1 and 2. The Pluecker vectors of any
two of them and of generator 4 are a lattice basis of the integer points of
the lineality space. The rows are `LINEALITY_SPACE_STAR` of
[`find_better_gens.jl`](../studies/g2/reference/fine_groebner/find_better_gens.jl),
in the order `G2_TREE_LINEALITY_STARS` of `render_maximal_cell_tree.jl`.

### Numbering

Both fields share one ray numbering, that of the
[G2-tree atlas](../studies/g2/atlas/flag/). Rays `1:51` are the rays of the QQ
census, in its order; they are also rays over GF(2). Rays `52`, `53`, `54`
occur only over GF(2): they are rays `3`, `36`, `39` of the GF(2) census, the
centres of the three Dressian tetrahedra. A cone that lies in both fans has
the same rays and the same tree in both fields.

A cell index is the index of the cell in the census of its field, as in the
columns `representative_cell` and `cells` of
`06_canonical_interior_degenerations.csv`. Before the renderer moved here,
its QQ cell index was the historical atlas order, with the three squares
last; fragments that record such an index also record the rays of the cell.

The three listings take `--field` and no other option:

| Option | Lines |
|---|---|
| `--list-rays` | Every ray with its index in the census of the field. |
| `--list-edges` | The two rays of every edge. Two opposite rays of a square do not span an edge. |
| `--list-cells` | Every maximal cell: its index, `triangle` or `square`, and its rays. |

### Metrics

A ray tree has the primitive integral Coxeter-balanced metric with minimum
pendant length zero. The tree of an edge or of a maximal cell is that of the
sum of the balanced generators of its rays. Drawn edge lengths are schematic;
the edge labels give the exact metric.

### Drawing options

| Option | Default | Effect |
|---|---|---|
| `--color` / `--grayscale` | Grayscale | Select the taxon palette. |
| `--label=true\|false` | `false` | Show or hide taxon numbers. |
| `--scale=FACTOR` | `1.0` | Scale tree coordinates, keeping symbol sizes, fonts, and strokes fixed. |
| `--taxon-scale=FACTOR` | `1.0` | Scale all seven taxon symbols and their numbers, keeping tree coordinates and metric fonts fixed. |
| `--bounded-length=LENGTH` | `0.70` | Set the drawn bounded-edge length before coordinate scaling, keeping exact metrics fixed. |
| `--edge-labels=disk\|aligned\|none` | `disk` | Show metric lengths in disks, sloped above edges, or hide them. |
| `--standalone` | Off | Supply the document of `standalone_figure.jl` instead of a fragment: its PDF is the figure at its size in main.tex. |
| `--output=FILE` | Standard output | Write the generated TeX to a file. |

Scale factors and bounded lengths must be positive and finite.
`--taxon-scale=0.7` gives 70% of the default symbol size; `1.5` gives 150%.
Pendant metric-label placement accounts for symbol size. A length on a disk
that would come within a point of a taxon symbol, a vertex dot, or another
disk moves along its edge until it is clear, and beside its edge only when
the edge is too short. Bounded lengths
that push internal vertices outside the taxon heptagon are rejected.

Fragments can be included with `\input{/tmp/g2_square.tex}` in a document
loading TikZ and xcolor. `--standalone` supplies the preamble of `main.tex`.

### Julia API

Include the renderer and pass the same options as Julia keywords:

```julia
include("julia/to_tex/render_maximal_cell_tree.jl")
source = load_g2_maximal_cell_source(; field=:gf2)
tex = g2_maximal_cell_tree_tex(1; source, palette=:color, label=true,
    taxon_scale=0.7, bounded_length=1.0, edge_labels=:aligned)
write_g2_edge_tree_tex("/tmp/g2_edge.tex", (6, 52); source, taxon_scale=0.7)
```

The fragments are `g2_lineality_tree_tex`, `g2_ray_tree_tex`,
`g2_edge_tree_tex`, and `g2_maximal_cell_tree_tex`, each with a writer
`write_...(path, ...)`. `g2_lineality_tree_tex` takes the number of the
generator. The others take `field=:qq` or `field=:gf2`, or a loaded
`source`, which then fixes the field. The same names without `_tex` return
the tree itself, as a `TreeMetric` with its description. `load_g2_ray_source(; field)` loads the
rays only; `load_g2_maximal_cell_source(; field)` also loads the cells and the
edges. Reuse a loaded source when drawing many trees.
`g2_tree_tikz(tree; kwargs...)` draws a supplied seven-taxon `TreeMetric`.
`g2_tree_commands(io, tree, colors; kwargs...)` writes the same drawing
without its picture, its heptagon centred at the origin, for a picture that
holds several trees; `colors` comes from `g2_minimal_tree_styles`.
Its `vertex_labels`, a dictionary from inner vertices to text, writes each
text in white in a larger dot, which the lengths on disks keep clear of.
`g2_heptagon_layout(tree, order)` gives the `layout_coordinates` that put the
taxa on the heptagon in a cyclic `order`, clockwise from the top, such as
the `topology_compatible_leaf_order` of another tree; every split of `tree`
must be contiguous in it.

The renderer includes the shared Grassmannian definitions in
[`atlas_common.jl`](../studies/g2/atlas/grassmannian/atlas_common.jl). The
[Grassmannian atlases](../studies/g2/atlas/grassmannian/) include the renderer.


## Simple trees

`render_simple_tree.jl` draws the same trees in black, without the taxon
symbols and the metric labels:

- A node with taxa is the sequence of their numbers. An edge of length zero
  is contracted, so taxa at distance zero share a node: `13` stands for taxa
  1 and 3 at distance zero. Such a node can be an inner node.
- An inner node without taxa is not marked.
- An edge of length `n` is `n` parallel lines. With `--max-lines=N`, an
  edge longer than `N` is one line with `n` written beside it, in
  `\scriptsize`, on the side where it is farther from the rest of the tree.
- All edges are drawn with the same length, measured between the labels or
  unmarked nodes at their ends. The position of a taxon has no meaning.
- The edges around a node are equally spaced, except that an edge to an
  outer node with `k` taxa takes `k` shares of the angle. A star with outer
  nodes `1`, `2`, `3`, `4` is a cross; with `12`, `3`, `4`, the edge to `12`
  makes an angle of 135 degrees with each of the other two.

Among the orders of the edges around the nodes, the script takes the one
that keeps the parts of the drawing farthest apart, and it turns the drawing
so that its longest direction is horizontal.

```sh
julia --project=julia julia/to_tex/render_simple_tree.jl \
  --rays=7,18,41 --edge-length=0.8 --output=/tmp/g2_simple_triangle.tex
```

The tree is selected as above, with `--lineality=K`, `--ray=K`, `--edge=I,J`,
`--cell=INDEX`, or `--rays=I,J,K[,L]`, and `--field`; the three listings are
also available. The other options are:

| Option | Default | Effect |
|---|---|---|
| `--edge-length=LENGTH` | `1.0` | The drawn length of every edge between its two ends, in centimetres. Line widths, the distance between parallel lines, and the labels keep their size. |
| `--max-lines=N` | None | Draw an edge longer than `N` as one line with its length beside it. |
| `--standalone` | Off | Supply the document of `standalone_figure.jl` instead of a fragment: its PDF is the figure at its size in main.tex. |
| `--output=FILE` | Standard output | Write the generated TeX to a file. |

The trees of the cells have edges of length up to 10; a tree with long edges
needs a larger `--edge-length`, or `--max-lines`, to stay legible. Edge lengths must be
integers. A lineality star may have edges of negative length, whose lines
are dashed; generator 4 has pendant lengths `1/2` and is rejected.

The fragment places the labels in em, so it follows the font size of the
document that inputs it; the order of the edges around the nodes is chosen
for a 10pt font.

From Julia, `g2_simple_tree_tikz(tree; edge_length, max_lines)` draws a
`TreeMetric`, and `g2_simple_tree_tex(g2_ray_tree(8); edge_length, max_lines,
standalone)` a selected tree with its description. To put several trees in
one picture, `g2_simple_tree_drawing` lays a tree out,
`g2_simple_tree_commands(io, drawing; prefix)` writes its nodes and lines
with node names starting with `prefix`, `g2_simple_tree_mirror(drawing)`
reflects it in a vertical line, `g2_simple_tree_styles()` gives the
styles they use, and `g2_simple_tree_outline(drawing, em_cm)` gives points
around which a circle encloses the drawing.

## Ray catalog

`build_styled_ray_orbit_tables.jl` writes a `tabular` of four rows of three
G2-trees, one for each orbit of rays under the Weyl group: the eleven orbits
of the rays `1:51`, from
[`ray_label_orbits.csv`](../studies/g2/atlas/grassmannian/ray_label_orbits.csv),
and the orbit of the rays `52`, `53`, `54`, which occur only over GF(2) and
whose panel is marked `char 2` at its bottom left and follows a vertical line. Each panel has the label of the orbit,
followed by its size in parentheses, and the tree of its first ray.
The label is a large Greek letter,
thickened by an outline of `pdfrender`, which `main.tex` loads, and the size
is smaller and not bold; the tree rises beside them; with
`--color` both, and the rules and lines of the table, are black, as the
labels of the maximal trees below. As in every tree, the edges, the inner
vertices, and the lengths are black. The labels and the order are
recorded in
[`ray_orbit_article_labels.md`](../studies/g2/atlas/grassmannian/ray_orbit_article_labels.md).
`main.tex` supplies the float and the caption.

```sh
julia --project=julia julia/to_tex/build_styled_ray_orbit_tables.jl \
  --color --label=false --scale=0.95 --taxon-scale=0.75 --bounded-length=1.1 \
  --output=figures/ray_orbits_styled_catalog.tex
```

This is the figure of `main.tex`, which includes it as the PDF of `--standalone`. The drawing options are those
of the G2-trees, except that the bounded length is `0.9` by default, that
`--edge-labels=none` is not available, and that `--scale` also scales the
columns, the row heights, and the gaps of the table. `--standalone` supplies
the preamble of `main.tex`. `--pdf` compiles a preview with that preamble to
`output/pdf/g2_ray_orbits_styled_catalog.pdf` instead of printing the TeX; it
needs `latexmk`.

## Cubics of the tropical basis

`cubic_tree_table.jl` (or `cubic_tree_table.sh`, which runs it) writes a `tabular` with two rows for each of the three
cubics of the tropical basis: its five monomials and, under each monomial,
the simple tree of a ray over GF(2). The table has no heading: `\toprule`
and `\bottomrule` frame it, and the space `CUBIC_TREE_TABLE_TREE_SKIP` of the script, not a
rule, separates the cubics.

Each cubic belongs to one tetrahedron of the G2-Dressian. Its first monomial,
with the coefficient 2, belongs to the ray that exists only over GF(2), in
the interior of the tetrahedron: there it alone attains the minimum. Each
other monomial belongs to a vertex of the tetrahedron: there it alone does
not attain the minimum. The script lists the ray of every monomial and
certifies both statements before it draws the trees, in one Julia session of
about half a minute.

```sh
sh julia/to_tex/cubic_tree_table.sh --edge-length=0.4 --output=FILE
```

`--edge-length` is that of the simple trees, `0.4` by default; five trees
then fit the text width of `main.tex`.

## Toric degeneration table

`toric_degeneration_table.jl` reads the stage-6 censuses of
`data/g2/grassmannian/` and `data/g2/grassmannian_gf2/`. The G2-tree atlas
supplies the orbit labels, their order, and the matching of the orbits of
both fields, through `local_flag_cell_atlas.csv`. A row is numbered by the
atlas label of its orbit, in bold black old-style figures, followed
by the orbit size in parentheses as in the tables of trees. Its last
column lists the characteristics in which the orbit occurs, "all" if both. The script writes a
`tabular*` of the text width only; `main.tex` supplies the float and the
caption. It runs in about two seconds.

```sh
julia --project=julia julia/to_tex/toric_degeneration_table.jl --output=FILE
```

## Trees of the toric degenerations

`toric_degeneration_tree_table.jl` writes a `tabular` of five rows of four
G2-trees, one for each row of the toric degeneration table: its panels are
the rows of that table, in the same order and with the same labels, so both
tables always agree. Each panel has the number of the orbit, set as the
labels of the ray catalog in old-style figures and followed by its size in
parentheses, and the
tree of the cell that the G2-tree atlas draws for the orbit, which is the
`support` column of `local_flag_cell_atlas.csv`. The script certifies that
this cell lies in the census orbit of the row in every field in which the
orbit is toric. A panel marked `char 2` at its bottom left, on a line below
its tree, belongs to an orbit that is toric
only over GF(2). `main.tex` supplies the float and the caption.

```sh
julia --project=julia julia/to_tex/toric_degeneration_tree_table.jl \
  --color --scale=0.88 --taxon-scale=0.75 --bounded-length=1 \
  --output=figures/toric_degeneration_trees.tex
```

This is the fragment in `figures/`. The options are those of the ray catalog,
whose panels the script shares; the trees have the automatic layout of
`render_maximal_cell_tree.jl`, and the bounded length is `0.70` by default.

Unlike those of the ray catalog, the panels are fitted to their trees. A
panel is as high as its heading and its tree, and the rows follow each other
at a quarter of the scale, in centimetres. A panel is as wide as the widest
tree that the options allow, and at least `3.9` cm, which the heading
`18 char 2 6` needed at the 12pt of `main.tex` when it still read `size`: below scale `0.99`
the table is therefore `15.6` cm wide plus its three gaps. The scale is
limited by the height of the page: with these options the five rows and a
caption of two lines fit a page of `main.tex` up to scale `0.90`, and a
caption of three lines up to `0.87`.

With `--non-toric` the script writes the same table for the other orbits of
maximal cells, those that are not toric in the census of some field. Their
numbers, `21` to `35` from the atlas labels, continue those
of the toric orbits, so both tables together have every orbit of maximal
cells of both fields once. Orbit
`21` is among them: at its canonical weight the initial ideal is prime in
the census, but not binomial. A panel marked `char ≠ 2` belongs to an orbit
of maximal cells over QQ only, and a panel marked `char 2` to one over GF(2)
only. A square at the top right marks a quadrilateral cell: orbit `30`, with
the rays `13, 14, 17, 18`. All other maximal cells are triangles.

In both tables a vertical line separates two panels of a row whose orbits
belong to different fields: the orbits of both fields come first, then those
of one field. The line stands between `17` and `18`, between `31` and `32`, and
between `34` and `35`. Where such a change falls on the end
of a row, there is no line.

```sh
julia --project=julia julia/to_tex/toric_degeneration_tree_table.jl --non-toric \
  --color --scale=0.88 --taxon-scale=0.75 --bounded-length=1 \
  --output=figures/non_toric_cell_trees.tex
```

This table has five rows of three panels: its widest heading,
`33 char ≠ 2 3`, needed `4.7` cm with `size`, and four panels of this width would
exceed the text width of `main.tex`. Its preview is
`output/pdf/g2_non_toric_cell_trees.pdf`.
The preview of `--pdf` is `output/pdf/g2_toric_degeneration_trees.pdf`.

## Ray orbits of the maximal cells

`cell_ray_orbit_table.jl` lists, for each of the 35 orbits of maximal cells of
both fields, the ray orbits of its rays, with the numbers of the two tables of
trees and the labels of the ray catalog. It reads the support cells and their
`polygon_ray_orbits` from the atlas `local_flag_cell_atlas.csv`, and translates
them with `PAPER_RAY_ARTICLE_LABELS`. The table has seven rows and five
columns, filled column by column. The numbers and the labels are black, as
the cell and ray labels of `main.tex`. The labels of a triangle are sorted in the
order of the catalog; those of the quadrilateral `30` follow its cyclic order.
The table is a `tabular*` of the text width without heading, between
`\toprule` and `\bottomrule`.

```sh
julia --project=julia julia/to_tex/cell_ray_orbit_table.jl \
  --output=tables/cell_ray_orbits.tex
```

## The quadrilateral cell

`square_cell_figure.jl` draws a square cell of the fan over QQ, by default
the one with the rays `13, 14, 17, 18` of orbit `30`. Its corners are circles
that hold the G2-trees of its rays `R_1`, ..., `R_4`, clockwise from the top
left in their cyclic order around the cell, starting from the smallest ray;
the cyclic order comes from the edges of the fan. Each ray is named outside
its circle, on the diagonal of the square, with the label of its orbit in
the ray catalog in parentheses, such as `R_1 (gamma_1)`. The label of an
orbit has the font of the tables of trees, bold and thickened by the outline
of `styled_tree_label` (pdfrender, which `main.tex` loads), and the number of the cell is in old-style figures. The G2-tree of the
cell, the sum of the balanced generators of its rays, stands in the middle;
it is named `S (30)` in green inside the bottom side of the square, with the
number of the orbit of the squares in the tables of maximal trees. All
squares of both fields form this one orbit, and the script reads its number
from the only orbit of squares of the atlas `local_flag_cell_atlas.csv`.

The trees are drawn by `g2_tree_commands` with the color palette and the
lengths on disks, in black, and so are the names of the rays and of the cell,
as the labels of `main.tex`. Every tree has its taxa at the places
of the heptagon that the tree of the cell gives them, by
`g2_heptagon_layout`: the splits of the rays are splits of the cell, so the
trees of the rays stay planar. Every ray is moreover drawn on the vertices
and edges of the tree of the cell: a bounded edge has the length of its
split in the ray, and `0` for a split of the cell that the ray lacks, whose
ends the ray merges; the pendant lengths are those of the ray. The script
certifies that each tree so drawn has the metric of its ray. The inner
vertices of every tree are lettered `A`, `B`, ... from top to bottom and
from left to right, in white in larger dots; as every tree has the vertices
of the cell at the same places, a letter names the same vertex in all five. The sides of the square and the circles are black, and the inside of
the square is a light tint of the blue of the taxa, `g2TaxonBlue!10`,
which is also the color of the disks under the lengths of the tree of the cell.

The taxa of every tree lie on the same heptagon, so the circles have one
radius, that of the heptagon plus the reach of a symbol. The square is the
smallest that leaves a side of `1` cm between two circles and keeps the tree
of the cell clear of the circles, the sides, and its name.

```sh
julia --project=julia julia/to_tex/square_cell_figure.jl \
  --output=figures/square_cell.tex
```

| Option | Default | Effect |
|---|---|---|
| `--rays=I,J,K,L` | `13,14,17,18` | The rays of the square. |
| `--scale=FACTOR` | `1.0` | That of the G2-trees: the size of the trees and of their edges. |
| `--taxon-scale=FACTOR` | `0.75` | That of the G2-trees. |
| `--bounded-length=LENGTH` | `1.1` | That of the G2-trees: the drawn length of the bounded edges, before the scale. |
| `--ray-scale=FACTOR` | `0.8` | The trees of the rays are this factor of the size of the tree of the cell: it multiplies their scale and their taxon scale. |
| `--ray-bounded-length=LENGTH` | `1.1` | The bounded length of the trees of the rays instead. |
| `--simple` | Off | Draw simple trees instead, with the two options below. |
| `--edge-length=LENGTH` | `0.8` | With `--simple`: that of the simple trees. |
| `--max-lines=N` | `3` | With `--simple`: an edge longer than `N` is one line with its length. |
| `--standalone` | Off | Supply the document of `standalone_figure.jl` instead of a fragment: its PDF is the figure at its size in main.tex. |
| `--output=FILE` | Standard output | Write the generated TeX to a file. |

The lettered dots need longer edges than the plain ones: hence the scale
`1.0` and the bounded length `1.1`. A bounded length of `1.3` pushes the
outer inner vertices of the tree of the cell, a chain of four bounded edges,
against the taxa, and leaves no room for the lengths of their pendant edges.
The trees of the rays have the four bounded edges of the cell, and the same
bounded length: they are copies of the tree of the cell at `0.8` of its
size, whose lettered dots, lengths, and fonts keep their full size, so that
the circles shrink but stay legible. Their bounded length `1.3` also pushes
`D` and `E` against the taxa.

With `--simple`, a tree whose taxon `1` lies right of its taxon `4` is
mirrored, so that every tree has `1` on the left of `4`, and each circle is
the smallest around its tree; the circles take the largest radius. The
inner vertices are not marked and have no letters. The fragment then sets
its font to the `12pt` of `main.tex`, for which these sizes are computed,
since the simple trees place their labels in em.

With the defaults the figure is about `13` cm wide and `11.5` cm high; it
takes half a minute.

## Local flag complexes

Above a maximal cell `S` of the tropical G2-Grassmannian stands the local
flag complex: the maximal cells of the tropical G2 flag variety that project
onto `S`, with all their faces. Its support is the preimage of `S`.
`local_flag_complex_table.jl` computes it for one cell of every orbit of
maximal cells, in each field in which the orbit occurs, and writes a
`tabularx` of the text width with one row per f-vector modulo lineality:
the f-vector, the maximal cells, and the list of the isomorphism types with
this f-vector, each label followed by the orbits of that type, numbered as in
the tables of trees. `main.tex` loads `tabularx` for it. The only rule inside
the table is the `\midrule` under its heading.

The script reads, for each field, the normalized coarsest flag fan, the
certified projection census `07_projection_cells.csv` (the Grassmannian cell
under every flag cell), and the stage-6 Grassmannian census for the orbits;
the atlas `local_flag_cell_atlas.csv` gives the orbit numbers. The faces of
every flag cell come from the facets of its own cone. Two complexes are
isomorphic when the incidences of their vertices with their faces are, by
polymake's `isomorphic`; the faces are cones over their rays, so these
incidences determine the face lattices.

The types are ordered by f-vector, then by their least orbit. Types with the
same f-vector share a letter and are numbered: `\mathfrak B_1`, ...,
`\mathfrak B_6`; a letter with a single type has no index (`\mathfrak D`,
`\mathfrak F`). In the data the f-vector also fixes the maximal cells; the
script checks it. The list of types wraps in its column, only between two
types, which are a `\quad` apart. An orbit whose complex had different types
in the two fields would carry its characteristic; none does. The script
checks that the complexes of all orbits, counted with the orbit sizes, hold
every maximal flag cell once: 2714 over QQ, 2792 over GF(2).

```sh
julia --project=julia julia/to_tex/local_flag_complex_table.jl \
  --output=tables/local_flag_complexes.tex
```

Its TeX comments list, for each type, the census orbit of each field, the
representative Grassmannian cell and the number of flag cells over it. It
takes about 40 seconds and 2.3 GB of memory.

## A triangular cell

`triangle_cell_figure.jl` draws a triangular cell as `square_cell_figure.jl`
draws the square, by default the atlas cell of orbit `21` over QQ, with the
rays `7, 12, 40`. The corners of an equilateral triangle are circles holding
the G2-trees of the rays, clockwise from the top in increasing order, each
named outside its circle by the label of its orbit in the ray catalog; the
G2-tree of the cell stands in the middle, named by the number of its orbit
inside the bottom side. No edge length is shown.

Every ray keeps its own tree: it is not drawn on the vertices and edges of
the cell, and has no edge of length zero. Its taxa are where the tree of the
cell has them. The inner vertices of the rays are lettered `A`, `B`, ..., ray
after ray, and from top to bottom and left to right in each. A letter is not
in its dot but beside it, at the end of each bounded edge of its vertex. A
split of a ray is a split of the cell, and in the tree of the cell the
letter stands at the end of the edge of that split on the same side. A
vertex of the cell then has several letters, each beside its own edge.
Rays that share a split share its edge, and their letters at one end are
written together, such as `B,E`.

The two letters of an edge stand on the same side of it, each near its end,
in the angle that the edge makes on that side with the next edge at that
end, at most `30` degrees from the edge; an angle holds at most one letter.
Among all choices of the sides, which the at most four bounded edges of a
tree allow to try, the script takes one that puts the letters of as many
edges as possible on one side, then the one whose narrowest angle with a
letter is the widest. The tree of the cell chooses first, and the trees of
the rays then put the letters of each split on the side that the cell chose
for it. The TeX comments of the fragment list the side of each letter, as
the taxa beyond its vertex.

To leave room for two letters on one side, the bounded edges are long: `1.6`
in the cell and `2.0` in the rays before scaling, that is `1.6` cm in all four
trees with the defaults, against `1.1` in the square figure. A tree whose
inner vertices would then come farther from the centre than `0.8` of the
radius of the heptagon of the taxa gets shorter edges; an explicit
`--bounded-length` or `--ray-bounded-length` is taken as it is. The TeX
comments record the lengths used. In orbits whose rays share a split, such
as `3` and `18`, the joined letters at the two ends of a shortened edge can
touch.

```sh
julia --project=julia julia/to_tex/triangle_cell_figure.jl \
  --scale=0.5 --taxon-scale=0.5 --ray-scale=0.8 --name-below \
  --no-cell-letters --cell-overhang --height=4.3 --output=figures/triangle_cell.tex
```

`--orbit=N` draws the atlas cell of another orbit of triangles, over QQ if it
occurs there, and `--rays=I,J,K` any triangle; `--field=GF2` selects the
other field. The drawing options are those of the square, without the edge
lengths and the simple trees: `--scale`, `--taxon-scale`, `--bounded-length`,
`--ray-scale`, `--ray-bounded-length`, `--standalone`, and `--output`.
`--name-below` puts the orbit number of the cell under the bottom side,
between the circles, instead of inside it under the tree, which lets the
triangle shrink around the same trees. With
the defaults the figure is `14.9` cm wide and `13.4` cm high, and with
`--scale=0.9` `14.4` cm by `13.1` cm; the paper uses `--scale=0.7
--taxon-scale=0.55 --name-below`, a triangle of circumradius `4.1` cm
instead of `4.96` without `--name-below`. The triangle is the
smallest whose sides keep the tree of the cell `0.3` cm away. With
`--cell-overhang` the tree of the cell may reach past the sides and only keeps
clear of the circles and of the label, and `--no-cell-letters` leaves it
without letters; the letters of the rays then simply stand at their
vertices, each once, as in `lettered_square_figure.jl`, and the rules
above for the ends and the sides of the edges no longer apply. `--height=CM` gives the triangle
that height, if the trees fit in it, instead of the smallest. The paper adds
these three with `--scale=0.5 --taxon-scale=0.5`, a triangle `4.3` cm high,
the height of the square of `lettered_square_figure.jl`; these are the largest
trees that fit in it. It takes half a minute.

The picture ends with `\ifdefined\gtwocelloverlay\gtwocelloverlay\fi`:
a document that defines this macro draws over the figure, in its coordinates
with the centre of the triangle at the origin. The paper's
`figures/triangle_projection.tex` defines it to put the hand-drawn incidence
graph of the tetrahedra of the local flag complex and their shared triangles
(`figures/triangle_incidence.tex`, a scope) on
the left of the triangle and lower than its picture, with a straight
two-headed arrow to the left side that stops just outside it, and draws the
whole at `0.95`. The positions are set by hand for the triangle of the paper's
options.
`figures/triangle_projection_5.tex` does the same for orbit `5`, the triangle
`figures/triangle_cell_5.tex` drawn with `--orbit=5` and the same options, and
the incidence graph `figures/triangle_incidence_5.tex` of its local flag
complex, 10 tetrahedra and a triangular prism, on the left at the height of
the arrow and lower than the triangle. The two graphs share their shapes,
their edge lengths, and their scale `0.64`.

## Decorated trees

A decorated G2-tree carries the free taxon `8` at a point of the tree, with a
pendant edge of length `0`: at an inner vertex it is a vertex of the tropical
G2 flag variety, inside a bounded edge of a maximal tree a maximal cell.
`decorated_trees_figure.jl` draws the free taxon as an apple. It takes the
trees of `triangle_cell_figure.jl` for a triangle (`g2_triangle_trees`, which
both scripts share), the ray at the size of the cell, and draws them side by
side without letters, under their orbit names: in the middle the ray with the
inner vertex lettered `--vertex` (default `A`), with the apple at it; on the
right the tree of the cell with the apple in the middle of the bounded edge
whose ends are lettered `--edge` (default `AB`). The letters are those of the
triangle figure. On the left, under the name lineality, is the star tree, its
taxa where the other two have them, with the apple in the middle of its
pendant edge to the taxon `--pendant` (default `1`, the dark square). The apple hangs straight down from its point, or turns in
steps of 5 degrees until its body clears the edges and the taxa; its leaf
takes the side where the edges through the point leave more room. For orbit
`21` all three apples hang straight down: from the pendant edge of the star
to taxon `1`, from vertex `8` of ray `7` (the side `{3,4}`) and from the edge
between vertices `10` and `11` of the cell, the tetrahedron `ABDE` of the
atlas.

```sh
julia --project=julia julia/to_tex/decorated_trees_figure.jl \
  --scale=0.7 --taxon-scale=0.55 --output=figures/decorated_trees.tex
```

`--orbit=N`, `--scale`, `--taxon-scale` and `--standalone` are as in the
triangle figure; the apple follows `--taxon-scale`. `--apple-style` draws
the apple `flat` (the default, as the taxa: no outline, no highlight),
`outlined` (its body and leaf outlined in dark green) or `glossy` (outlined,
with a white highlight on its body, the former style); the scripts
`toric_degeneration_prime_edge_trees.jl` and
`toric_degeneration_apple_trees.jl` take the same option. The figure is about
`11.5` cm wide with these options, `16` cm with the defaults; it takes about a
minute.

## A lettered quadrilateral

`lettered_square_figure.jl` draws the quadrilateral of
`square_cell_figure.jl` the way `triangle_cell_figure.jl` draws a triangle:
the four rays clockwise from the top left, in their cyclic order, each its own
tree in a circle at a corner, with its inner vertices lettered `A`, ..., `L`
in that order; the tree of the cell in the middle, without letters unless
`--cell-letters`, and then each letter of a ray once, beside its vertex in
the widest free angle; the orbit names; and the hook `\gtwocelloverlay`, the
centre of the square at the origin. The trees and their letters come from
`g2_triangle_trees`, which takes a square as well as a triangle. The
placement of the letters avoids overlaps, with each other and with the taxa,
before it puts the letters of an edge on one side. Without letters in the
cell, `g2_triangle_trees` puts those of the rays at their vertices, for the
square as for a triangle with `--no-cell-letters`.

```sh
julia --project=julia julia/to_tex/lettered_square_figure.jl \
  --scale=0.56 --taxon-scale=0.5 --name-below --ray-bounded-length=1.3 \
  --output=figures/square_lettered.tex
```

