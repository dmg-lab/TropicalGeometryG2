###############################################################################
# The G2-trees of the cones of the coarsest fans of the tropical
# G2-Grassmannian over QQ and over GF(2), in the rays of the atlas: the
# generators of the lineality space, the rays, the edges and the maximal
# cells; and the tree of the support cell of each orbit of maximal cells
# that cell_orbit_labels.jl numbers.
#
# Included by to_tex/render_maximal_cell_tree.jl, which draws these trees,
# and by export_maximal_cell_tree_metrics.jl. The public entry points are
# load_g2_ray_source, load_g2_maximal_cell_source, the trees
# g2_lineality_tree, g2_ray_tree, g2_edge_tree and g2_maximal_cell_tree, and
# toric_tree_table_rows.
#
# Both fields share one ray numbering: rays 1:51 are the rays over QQ, and
# rays 52:54 occur only over GF(2). A cell index is the index of the cell in
# the census of its field.
###############################################################################

include(joinpath(@__DIR__, "..", "grassmannian", "atlas_common.jl"))
# The censuses, the atlas labels and their numbers, unless a script that
# includes this file has included them already.
isdefined(@__MODULE__, :toric_table_support_rows) ||
    include(joinpath(@__DIR__, "cell_orbit_labels.jl"))

const G2_TREE_FAN_SIZES = (qq=(rays=51, cells=248), gf2=(rays=54, cells=254))

# Pendant lengths of the star trees that generate the lineality space: the
# three short coroots of the maximal torus of G2 in main.tex, the third being
# minus the sum of the first two, which span it modulo the all-ones vector,
# and the star of the all-ones Pluecker vector.
const G2_TREE_LINEALITY_STARS = (
    QQ.(LINEALITY_SPACE_STAR[2, :]),
    QQ.(LINEALITY_SPACE_STAR[3, :]),
    QQ.(LINEALITY_SPACE_STAR[4, :]),
    QQ.(LINEALITY_SPACE_STAR[1, :]) .// 2,
)

g2_tree_field_name(field::Symbol) = field == :qq ? "QQ" : "GF(2)"

# Fragments over QQ keep the descriptions they had before the fields.
g2_tree_field_suffix(field::Symbol) = field == :qq ? "" : " over GF(2)"

"""
    load_g2_ray_source(; field=:qq)

Load the balanced ray generators of the G2 Grassmannian over `field`, `:qq`
or `:gf2`, without loading its fan: 51 rays over QQ, and three more over
GF(2). `ray_index` maps the rays of the census of `field` to these rays.
"""
function load_g2_ray_source(; field::Symbol=:qq)
    balanced_rays = load(DEFAULT_BALANCED_RAYS_PATH)
    nrows(balanced_rays) == 51 || error("expected 51 balanced rays")
    census = census_atlas_rays(field, balanced_rays)
    nrows(census.rays) == G2_TREE_FAN_SIZES[field].rays ||
        error("unexpected number of rays over $(g2_tree_field_name(field))")
    field == :gf2 && (balanced_rays = matrix(QQ, census.rays))
    return (field=field, balanced_rays=balanced_rays, ray_index=census.ray_index)
end

"""
The edges of the fan with maximal cells `cells`: the pairs of rays that span
one of its two-dimensional cones. All pairs of a triangle do. Of a square
only the four sides do; in the linear relation among its rays modulo
lineality, a side pairs two coefficients of opposite sign.
"""
function g2_fan_edges(cells, rays)
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, balanced_lineality_basis()); side=:right))
    projected = Matrix{QQFieldElem}(rays) * quotient
    edges = Set{Tuple{Int,Int}}()
    for cell in cells
        positive = falses(length(cell))
        if length(cell) == 4
            relation = kernel(matrix(QQ, projected[collect(cell), :]); side=:left)
            nrows(relation) == 1 || error("square $cell has no single relation among its rays")
            positive = [relation[1, position] > 0 for position in 1:4]
            count(positive) == 2 && all(position -> !iszero(relation[1, position]), 1:4) ||
                error("the relation among the rays of square $cell does not alternate")
        end
        for first in 1:length(cell), second in (first + 1):length(cell)
            length(cell) == 4 && positive[first] == positive[second] && continue
            push!(edges, (cell[first], cell[second]))
        end
    end
    return sort!(collect(edges))
end

"""
    load_g2_maximal_cell_source(; field=:qq)

Load the coarsest fan of the G2 Grassmannian over `field`: its rays as in
`load_g2_ray_source`, its maximal cells in the order of the census of
`field`, each the sorted tuple of its rays, and its edges. Over QQ these are
248 cells, over GF(2) 254. Reuse the returned source when drawing many trees.
"""
function load_g2_maximal_cell_source(; field::Symbol=:qq)
    source = load_g2_ray_source(; field=field)
    cells = census_atlas_cells(field, source.ray_index)
    length(cells) == G2_TREE_FAN_SIZES[field].cells ||
        error("unexpected number of maximal cells over $(g2_tree_field_name(field))")
    all(cell -> length(cell) in (3, 4), cells) ||
        error("a maximal cell is neither a triangle nor a square")
    lookup = Dict{Tuple{Vararg{Int}},Int}()
    for (index, cell) in enumerate(cells)
        haskey(lookup, cell) && error("duplicate maximal cell $cell")
        lookup[cell] = index
    end
    edges = g2_fan_edges(cells, source.balanced_rays)
    return merge(source, (cells=cells, lookup=lookup, edges=edges))
end

"""Resolve a 1-based cell index or a set of 3/4 ray indices."""
function g2_resolve_maximal_cell(source, cell)
    if cell isa Integer && !(cell isa Bool)
        1 <= cell <= length(source.cells) ||
            throw(ArgumentError("cell index must lie in 1:$(length(source.cells))"))
        return (index=Int(cell), rays=source.cells[Int(cell)])
    end
    cell isa Tuple || cell isa AbstractVector ||
        throw(ArgumentError("cell must be a 1-based index or ray-index collection"))
    all(value -> value isa Integer && !(value isa Bool), cell) ||
        throw(ArgumentError("ray indices must be integers"))
    rays = Tuple(sort!(Int.(collect(cell))))
    length(rays) in (3, 4) && length(unique(rays)) == length(rays) ||
        throw(ArgumentError("a maximal cell needs 3 or 4 distinct ray indices"))
    index = get(source.lookup, rays, 0)
    index > 0 || throw(ArgumentError(
        "$rays is not a maximal cell of the coarsest fan over $(g2_tree_field_name(source.field))",
    ))
    return (index=index, rays=source.cells[index])
end

"""Resolve the two ray indices of an edge of the fan of `source`."""
function g2_resolve_edge(source, edge)
    edge isa Tuple || edge isa AbstractVector ||
        throw(ArgumentError("edge must be a collection of two ray indices"))
    all(value -> value isa Integer && !(value isa Bool), edge) ||
        throw(ArgumentError("ray indices must be integers"))
    rays = Tuple(sort!(Int.(collect(edge))))
    length(rays) == 2 && rays[1] != rays[2] ||
        throw(ArgumentError("an edge needs 2 distinct ray indices"))
    rays in source.edges || throw(ArgumentError(
        "$rays is not an edge of the coarsest fan over $(g2_tree_field_name(source.field))",
    ))
    return rays
end

"""The star tree with the seven pendant `lengths`, which may be negative."""
function g2_star_tree(lengths)
    length(lengths) == 7 || throw(ArgumentError("a seven-taxon star has 7 pendant lengths"))
    centre = 8
    adjacency = Dict{Int,Vector{Tuple{Int,QQFieldElem}}}(
        taxon => [(centre, QQ(lengths[taxon]))] for taxon in 1:7
    )
    adjacency[centre] = [(taxon, QQ(lengths[taxon])) for taxon in 1:7]
    return TreeMetric(7, adjacency)
end

###############################################################################
# The tree of a selection, with the description that heads its fragment
###############################################################################

"""
    g2_lineality_tree(generator)

The star tree of a generator of the lineality space, which is the same over
both fields, with its pendant lengths. Generators 1, 2 and 3 are the short
coroots `(1,-1,0,-1,1,0,0)`, `(-1,0,1,1,0,-1,0)` and `(0,1,-1,0,-1,1,0)`,
which sum to zero; generator 4 has all pendant lengths `1/2`, and its
Pluecker vector is the all-ones vector. The Pluecker vectors are primitive
and span the lineality space.
"""
function g2_lineality_tree(generator::Integer)
    !(generator isa Bool) && generator in eachindex(G2_TREE_LINEALITY_STARS) ||
        throw(ArgumentError(
            "lineality generator must lie in 1:$(length(G2_TREE_LINEALITY_STARS))",
        ))
    lengths = G2_TREE_LINEALITY_STARS[generator]
    tree = g2_star_tree(lengths)
    pluecker = star_pluecker(lengths)
    lineality = matrix(QQ, balanced_lineality_basis())
    realizes(tree, pluecker) && primitive_integer(pluecker) == pluecker &&
        rank(vcat(lineality, matrix(QQ, permutedims(pluecker)))) == rank(lineality) ||
        error("star $generator is not a primitive vector of the lineality space")
    description = "Lineality space of the G2 Grassmannian: generator $generator, " *
                  "pendant lengths $(join(rational_key.(lengths), ","))"
    return (tree=tree, description=description)
end

"""
    g2_ray_tree(ray; field=:qq, source=nothing)

The tree of a ray index in `1:51` over QQ or in `1:54` over GF(2), with the
exact balanced ray metric. A preloaded `source` fixes the field.
"""
function g2_ray_tree(ray::Integer; field::Symbol=:qq, source=nothing)
    source === nothing && (source = load_g2_ray_source(; field=field))
    !(ray isa Bool) && 1 <= ray <= nrows(source.balanced_rays) ||
        throw(ArgumentError("ray index must lie in 1:$(nrows(source.balanced_rays))"))
    record = balance_outer_edges(row_as_vector(source.balanced_rays, Int(ray)))
    description = "G2 Grassmannian ray $ray" * g2_tree_field_suffix(source.field)
    return (tree=record.tree, description=description)
end

"""
    g2_edge_tree(edge; field=:qq, source=nothing)

The tree of an edge of the coarsest fan over `field`: `edge` holds its two ray
indices, and the tree is that of the sum of the two balanced ray generators.
Pass a `source` of `load_g2_maximal_cell_source` to resolve many edges
without reloading the fan.
"""
function g2_edge_tree(edge; field::Symbol=:qq, source=nothing)
    source === nothing && (source = load_g2_maximal_cell_source(; field=field))
    rays = g2_resolve_edge(source, edge)
    record = normalized_prevariety_cone_record(rays, source.balanced_rays)
    description = "G2 Grassmannian edge" * g2_tree_field_suffix(source.field) *
                  ": rays $(join(rays, ","))"
    return (tree=record.tree, description=description)
end

"""
    g2_maximal_cell_tree(cell; field=:qq, source=nothing)

The tree of a maximal cell of the coarsest fan over `field`: `cell` is its
1-based index in the census of `field` or its 3/4 ray indices, and the tree is
that of the sum of the balanced generators of its rays. Pass a preloaded
`source`, which then fixes the field, to resolve many cells without reloading
the fan.
"""
function g2_maximal_cell_tree(cell; field::Symbol=:qq, source=nothing)
    source === nothing && (source = load_g2_maximal_cell_source(; field=field))
    selected = g2_resolve_maximal_cell(source, cell)
    record = normalized_prevariety_cone_record(selected.rays,
                                               source.balanced_rays)
    description = "Coarsest G2 Grassmannian cell $(selected.index)" *
                  g2_tree_field_suffix(source.field) *
                  ": rays $(join(selected.rays, ","))"
    return (tree=record.tree, description=description)
end

"""
The rows of `toric_table_support_rows`, each with the `tree` of its support
cell. The cell is certified to lie in the census orbit of the row in each of
its fields; it has the same tree in these fields.
"""
function toric_tree_table_rows(; toric::Bool=true)
    sources = Dict(field => load_g2_maximal_cell_source(; field=field)
                   for field in TORIC_TABLE_FIELDS)
    census = toric_table_census()
    return map(toric_table_support_rows(census; toric=toric)) do row
        for (field, orbit) in row.orbits
            cell = g2_resolve_maximal_cell(sources[field], row.support)
            cell.index in parse.(Int, split(census[field][orbit]["cells"], ';')) || error(
                "the support of atlas label $(row.label) is not in " *
                "$(TORIC_TABLE_FIELD_NAMES[field]) orbit $orbit",
            )
        end
        source = sources[first(first(row.orbits))]
        merge(row, (tree=g2_maximal_cell_tree(row.support; source=source).tree,))
    end
end
