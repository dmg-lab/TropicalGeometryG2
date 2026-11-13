###############################################################################
# A triangular cell of the tropical G2-Grassmannian: a triangle whose corners
# are circles holding the G2-trees of its three rays, with the G2-tree of the
# cell in its middle, as square_cell_figure.jl draws a quadrilateral. No edge
# length is shown. Every ray keeps its own tree, and its inner vertices are
# lettered A, B, C, ...: a letter stands next to its vertex, once, at the end
# of one bounded edge of that vertex. The splits of a ray are splits of the
# cell, so a bounded edge of a ray is a bounded edge of the cell, and in the
# tree of the cell the letters stand at the same ends of that edge: a vertex
# of the cell can so receive several letters, each on the side of its edge.
# The two letters of an edge stand on one side of it, the same in the tree of
# the ray and in that of the cell, and the bounded edges are long enough for
# that. With --no-cell-letters the cell has no letters, and those of a ray
# simply stand at their vertices, as in lettered_square_figure.jl.
#
#   julia --project=. to_tex/triangle_cell_figure.jl \
#     [--orbit=N | --rays=I,J,K] [--field=QQ|GF2] [--scale=FACTOR] \
#     [--taxon-scale=FACTOR] [--bounded-length=LENGTH] [--ray-scale=FACTOR] \
#     [--ray-bounded-length=LENGTH] [--name-below] [--standalone] [--output=FILE]
#
# The orbit numbers are those of the tables of maximal trees, from the atlas.
###############################################################################

include(joinpath(@__DIR__, "square_cell_figure.jl"))

# The orbit of maximal cells drawn by default, whose rays have one split each.
const G2_TRIANGLE_ORBIT = 21
# The letters of the inner vertices, in points.
const G2_TRIANGLE_LETTER_FONT_PT = 8
# A letter keeps this far from the edges around its vertex, in centimetres.
const G2_TRIANGLE_LETTER_CLEARANCE_CM = 0.04
# At most this angle from its edge, in degrees.
const G2_TRIANGLE_LETTER_ANGLE = 30
# The drawn lengths of the bounded edges of the cell and of the rays before
# scaling, long enough for the two letters of an edge on one side of it.
const G2_TRIANGLE_BOUNDED_LENGTH = 1.6
const G2_TRIANGLE_RAY_BOUNDED_LENGTH = 2.0
# Shorter if an inner vertex would come farther from the centre than this
# fraction of the radius of the heptagon of the taxa.
const G2_TRIANGLE_INNER_REACH = 0.8

"""
The census `field` of an orbit `number` of the tables of maximal trees, QQ
if it occurs there, and the rays of the triangle that the atlas draws for it.
"""
function g2_triangle_orbit_support(number::Int; field=nothing)
    rows = [row for row in toric_table_csv_rows(TORIC_TABLE_ATLAS)
            if all(f -> f in ("qq", "gf2"), split(row["fields"], '+')) &&
               toric_table_label_number(row["label"]) == number]
    isempty(rows) && throw(ArgumentError("the atlas has no orbit $number of maximal cells"))
    row = first(rows)
    row["kind"] == "triangle" || throw(ArgumentError("orbit $number is a $(row["kind"]), not a triangle"))
    fields = Symbol.(split(row["fields"], '+'))
    field === nothing && (field = :qq in fields ? :qq : :gf2)
    field in fields || throw(ArgumentError("orbit $number does not occur over $(g2_tree_field_name(field))"))
    return field, Tuple(parse.(Int, split(row["support"], ';')))
end

"""The number in the tables of maximal trees of the orbit of the cell `index` of the census of `field`."""
function g2_triangle_orbit_number(field::Symbol, index::Int)
    census = toric_table_csv_rows(TORIC_TABLE_CENSUS[field])
    orbit = only(row["orbit"] for row in census if string(index) in split(row["cells"], ';'))
    numbers = unique(toric_table_label_number(row["label"]) for row in toric_table_csv_rows(TORIC_TABLE_ATLAS)
                     if row["$(field)_orbit"] == orbit)
    return only(numbers)
end

"""The bounded edges of `tree`, as pairs of inner vertices, each in both directions."""
g2_triangle_bounded_ends(tree) =
    [(vertex, neighbor) for (vertex, neighbors) in tree.adjacency if vertex > tree.ntaxa
     for (neighbor, _) in neighbors if neighbor > tree.ntaxa]

"""
The places of the `letters` of a tree drawn at `coordinates`: `letters` maps
an end `(v, w)` of a bounded edge, at `v`, to its text; the ends of a vertex
with the same text are one letter, which stands at one of them. A letter
stands near its end, in the angle that its edge makes on one side with the
next edge at that end, and an angle holds at most one letter. Among the
choices of the ends and the sides, as few letters as possible overlap, as
disks around their texts, each other or the disks `obstacles`, given as
centres and radii, then the letters of as many edges as possible stand on one
side of their edge, then as many letters as possible stand on the side that
`preferred` gives them, and then the narrowest angle that holds a letter is
the widest: all choices are tried, since a tree has at most four bounded
edges. A letter is at most `G2_TRIANGLE_LETTER_ANGLE` from its edge and keeps
clear of both edges of its angle. Returns the places, by the ends that hold
the letters, and the sides: an end `(v, w)` is named by the taxa on the side
of `v`, and its side is `1` counterclockwise from `v -> w` and `-1`
clockwise, as `preferred` names them, so that the trees of a ray and of its
cell, whose taxa are at the same places, can share them.
"""
function g2_triangle_letter_places(tree, coordinates, letters; preferred=Dict{Vector{Int},Int}(),
                                   obstacles=Tuple{Tuple{Float64,Float64},Float64}[])
    isempty(letters) && return Dict{Tuple{Int,Int},Tuple{Float64,Float64}}(), Dict{Vector{Int},Int}()
    angle(v, w) = atan(coordinates[w][2] - coordinates[v][2], coordinates[w][1] - coordinates[v][1])
    font_cm = G2_TRIANGLE_LETTER_FONT_PT * 2.54 / 72.27
    wanted = deg2rad(G2_TRIANGLE_LETTER_ANGLE)
    # The angle at v from the edge to w to the next edge, counterclockwise for
    # side 1 and clockwise for side -1.
    function sector(v, w, side)
        turn(u) = mod(side * (angle(v, u) - angle(v, w)), 2pi)
        next = argmin(turn, [u for (u, _) in tree.adjacency[v] if u != w])
        return (gap=turn(next), key=(v, minmax(w, next)...))
    end
    edges = sort!(unique(minmax(end_...) for end_ in keys(letters)))
    length(edges) <= 6 || error("too many lettered edges to try every choice of sides")
    # The letters in the order of the edges, each with its ends and their
    # sides, as sides at the end: 1 counterclockwise from a -> b. An edge
    # (v, w), v < w, tries its left side first, which is 1 at v and -1 at w.
    slots = Vector{Vector{Tuple{Tuple{Int,Int},Int}}}()
    slot_of = Dict{Tuple{Int,String},Int}()
    for (v, w) in edges, (a, b) in ((v, w), (w, v))
        haskey(letters, (a, b)) || continue
        options = [((a, b), a == v ? side : -side) for side in (1, -1)]
        index = get!(slot_of, (a, letters[(a, b)])) do
            push!(slots, eltype(slots)())
            length(slots)
        end
        append!(slots[index], options)
    end
    # Half the diagonal of a text, as a disk around its centre.
    reach(text) = hypot(.6font_cm * length(text), .7font_cm) / 2 + G2_TRIANGLE_LETTER_CLEARANCE_CM
    # The centre of the letter of the end (a, b) on its side.
    function place(a, b, side)
        offset = min(sector(a, b, side).gap / 2, wanted)
        distance = reach(letters[(a, b)]) / sin(offset)
        direction = angle(a, b) + side * offset
        return (coordinates[a][1] + distance * cos(direction), coordinates[a][2] + distance * sin(direction))
    end
    best = nothing
    for choice in Iterators.product(slots...)
        sectors = [sector(a, b, side) for ((a, b), side) in choice]
        allunique(s.key for s in sectors) || continue
        sides = Dict(end_ => side for (end_, side) in choice)
        # An edge with letters at both ends has them on one side of it when
        # they turn opposite ways.
        one_side = count(haskey(sides, (v, w)) && haskey(sides, (w, v)) && sides[(v, w)] == -sides[(w, v)]
                         for (v, w) in edges)
        agree = count(get(preferred, leaves_on_side(tree, a, b), 0) == side for ((a, b), side) in choice)
        disks = [(place(a, b, side), reach(letters[(a, b)])) for ((a, b), side) in choice]
        overlaps = count(g2_simple_distance(p, q) < r + t for (i, (p, r)) in enumerate(disks)
                         for (q, t) in disks[i + 1:end]) +
                   count(g2_simple_distance(p, q) < r + t for (p, r) in disks for (q, t) in obstacles)
        score = (-overlaps, one_side, agree, round(minimum(s.gap for s in sectors); digits=3))
        (best === nothing || score > best.score) && (best = (; score, choice))
    end
    best === nothing && error("the letters do not fit the angles of their vertices")
    places = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    chosen = Dict{Vector{Int},Int}()
    for ((a, b), side) in best.choice
        chosen[leaves_on_side(tree, a, b)] = side
        places[(a, b)] = place(a, b, side)
    end
    return places, chosen
end

"""
The places of the `letters` of a tree drawn at `coordinates`, each letter at
its vertex rather than at an end: `letters` maps an end `(v, w)` of a bounded
edge, at `v`, to its text, which stands once for `v`. A letter stands on the
bisector of an angle between consecutive edges at its vertex, pendant edges
included, as close as it can while it keeps clear of both edges and of the
dot. Among the choices of the angles, as few letters as possible overlap, as
disks around their texts, each other, the other edges or the disks
`obstacles`, given as centres and radii, and then the narrowest chosen angle
is the widest. Returns the places by one end of each vertex and no sides.
"""
function g2_vertex_letter_places(tree, coordinates, letters;
                                 obstacles=Tuple{Tuple{Float64,Float64},Float64}[])
    places = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    isempty(letters) && return places, Dict{Vector{Int},Int}()
    font_cm = G2_TRIANGLE_LETTER_FONT_PT * 2.54 / 72.27
    dot_cm = G2_TREE_DOT_RADIUS_PT * 2.54 / 72.27
    reach(text) = hypot(.6font_cm * length(text), .7font_cm) / 2 + G2_TRIANGLE_LETTER_CLEARANCE_CM
    ends = Dict{Int,Tuple{Int,Int}}()
    for end_ in sort!(collect(keys(letters)))
        haskey(ends, end_[1]) || (ends[end_[1]] = end_)
        letters[end_] == letters[ends[end_[1]]] || error("vertex $(end_[1]) has two letters")
    end
    vertices = sort!(collect(keys(ends)))
    segments = [(v, w) for v in keys(tree.adjacency) for (w, _) in tree.adjacency[v] if v < w]
    distance_to(p, (v, w)) = g2_segment_distance(p, coordinates[v], coordinates[w])
    # The angles at a vertex, each as its place and its width.
    function options(v)
        r = reach(letters[ends[v]])
        directions = sort!([atan(coordinates[u][2] - coordinates[v][2], coordinates[u][1] - coordinates[v][1])
                            for (u, _) in tree.adjacency[v]])
        map(enumerate(directions)) do (i, a)
            gap = i < length(directions) ? directions[i + 1] - a : directions[1] + 2pi - a
            bisector = a + gap / 2
            distance = max(r / sin(min(gap / 2, pi / 2)), r + dot_cm)
            ((coordinates[v][1] + distance * cos(bisector), coordinates[v][2] + distance * sin(bisector)), gap)
        end
    end
    choices = [options(v) for v in vertices]
    best = nothing
    for choice in Iterators.product(choices...)
        disks = [(p, reach(letters[ends[v]])) for (v, (p, _)) in zip(vertices, choice)]
        overlaps = count(g2_simple_distance(p, q) < r + t for (i, (p, r)) in enumerate(disks)
                         for (q, t) in disks[i + 1:end]) +
                   count(g2_simple_distance(p, q) < r + t for (p, r) in disks for (q, t) in obstacles) +
                   count(distance_to(p, segment) < r for (v, (p, r)) in zip(vertices, disks)
                         for segment in segments if v ∉ segment)
        score = (-overlaps, round(minimum(gap for (_, gap) in choice); digits=3))
        (best === nothing || score > best.score) && (best = (; score, choice))
    end
    for (v, (p, _)) in zip(vertices, best.choice)
        places[ends[v]] = p
    end
    return places, Dict{Vector{Int},Int}()
end

"""
A G2-tree for the triangle: drawn by `g2_tree_commands` without edge
lengths, with the taxa in the cyclic `order` of `g2_heptagon_layout`, and the `letters` of
`g2_triangle_letter_places` with their `preferred` sides, which returns their
`sides`, or with `at_vertices` of `g2_vertex_letter_places`. The writer takes an `io`; `points` are points of
the drawing, around its centre at the origin.
"""
function g2_triangle_tree(tree; order, scale, taxon_scale, bounded_length, letters,
                          preferred=Dict{Vector{Int},Int}(), at_vertices::Bool=false)
    layout = g2_heptagon_layout(tree, order)
    coordinates = g2_minimal_tree_coordinates(tree; scale=scale, bounded_length=bounded_length,
                                               layout_coordinates=layout)
    reach = G2_SQUARE_SYMBOL_REACH_MM / 10 * G2_TREE_SYMBOL_SCALE * taxon_scale
    taxa = [(coordinates[taxon], reach) for taxon in 1:tree.ntaxa]
    places, sides = at_vertices ? g2_vertex_letter_places(tree, coordinates, letters; obstacles=taxa) :
                    g2_triangle_letter_places(tree, coordinates, letters; preferred=preferred, obstacles=taxa)
    function write(io)
        colors = g2_minimal_tree_styles(io, :color)
        g2_tree_commands(io, tree, colors; scale=scale, taxon_scale=taxon_scale,
                         bounded_length=bounded_length, layout_coordinates=layout,
                         edge_labels=:none)
        for (vw, (x, y)) in sort!(collect(places); by=first)
            println(io, "\\node[g2triangle letter] at ($(g2_simple_number(x)),$(g2_simple_number(y))) {$(letters[vw])};")
        end
        println(io, raw"\endgroup")
    end
    points = [coordinates[taxon] .+ reach .* (cos(a), sin(a)) for taxon in 1:tree.ntaxa
              for a in range(0, 2pi; length=17)]
    append!(points, values(places))
    radius = G2_TREE_UNIT_CM * G2_TREE_LEAF_RADIUS * scale + reach
    return (write=write, coordinates=coordinates, points=points, radius=radius, sides=sides)
end

"""
The letters of the trees of the `rays` and of the `cell`, drawn in the
taxon `order`: the inner vertices of the rays, in their order and from top
to bottom and left to right in each, take A, B, C, ...; each end of a bounded
edge at a lettered vertex has its letter, and so has the end of the same
split of the cell on the same side. Returns the letters of each ray, those of
the cell, and lines that record them.
"""
function g2_triangle_letters(rays, cell, order, ray_coordinates)
    ray_letters = Dict{Tuple{Int,Int},String}[]
    cell_letters = Dict{Tuple{Int,Int},String}()
    cell_ends = Dict(leaves_on_side(cell, v, w) => (v, w) for (v, w) in g2_triangle_bounded_ends(cell))
    record = String[]
    count = 0
    for (tree, coordinates) in zip(rays, ray_coordinates)
        inner = sort!([v for v in keys(tree.adjacency) if v > tree.ntaxa];
                      by=v -> (-round(coordinates[v][2]; digits=6), round(coordinates[v][1]; digits=6)))
        letters = Dict{Tuple{Int,Int},String}()
        for v in inner
            count += 1
            count <= 26 || error("more inner vertices than letters")
            letter = string(Char(Int('A') + count - 1))
            sides = String[]
            for (vv, w) in g2_triangle_bounded_ends(tree)
                vv == v || continue
                letters[(v, w)] = letter
                side = leaves_on_side(tree, v, w)
                haskey(cell_ends, side) || error("a split of a ray is not a split of the cell")
                end_ = cell_ends[side]
                cell_letters[end_] = haskey(cell_letters, end_) ? cell_letters[end_] * "," * letter : letter
                push!(sides, "{" * join(side, ",") * "}")
            end
            isempty(sides) && error("an inner vertex of a ray has no bounded edge")
            push!(record, "$letter: the side $(join(sides, " and ")) of its edges")
        end
        push!(ray_letters, letters)
    end
    return ray_letters, cell_letters, record
end

"""
The drawn length of the bounded edges of `trees` before scaling, in the taxon
`order`: `default`, or less if an inner vertex of one of them would come
farther from the centre than `G2_TRIANGLE_INNER_REACH` of the radius of the
heptagon.
"""
function g2_triangle_bounded_length(trees, order, default)
    # The inner vertices are sums of bounded edges, from their mean: their
    # distances grow with the length in proportion.
    reach = maximum(trees) do tree
        coordinates = g2_minimal_tree_coordinates(tree; bounded_length=0.5,
                                                   layout_coordinates=g2_heptagon_layout(tree, order))
        2maximum(hypot(coordinates[v]...) for v in keys(coordinates) if v > tree.ntaxa)
    end
    return min(default, G2_TRIANGLE_INNER_REACH * G2_TREE_LEAF_RADIUS / reach)
end

"""
The trees of the triangular or quadrilateral `cell` of the fan over `field`,
given by its index or its rays, as `g2_triangle_cell_tex` draws them: its
`rays`, those of a triangle in increasing order and those of a square in their
cyclic order from the smallest, the `number` of its orbit, the `selected` tree of the cell,
the taxon `order`, the tree of the cell in the `middle`, the `ray_trees` and their drawings at
the `corners`, from `g2_triangle_tree`, with the `letters` of the cell and `ray_letters` of
the rays from `g2_triangle_letters` and the lines that `record` them, and the
drawn lengths of the bounded edges. Without `cell_letters`, the tree of the
cell has no letters and those of a ray stand at their vertices, by
`g2_vertex_letter_places`.
"""
function g2_triangle_trees(cell; field::Symbol=:qq,
                           scale::Real=1.0, taxon_scale::Real=0.75, bounded_length=nothing,
                           ray_scale::Real=0.8, ray_bounded_length=nothing, cell_letters::Bool=true,
                           source=nothing)
    source === nothing && (source = load_g2_maximal_cell_source(; field=field))
    resolved = g2_resolve_maximal_cell(source, cell)
    rays = collect(resolved.rays)
    length(rays) in (3, 4) || throw(ArgumentError("$(Tuple(rays)) is not a triangle or a square"))
    length(rays) == 4 && (rays = g2_square_cycle(source, cell))
    selected = g2_maximal_cell_tree(rays; source=source)
    number = g2_triangle_orbit_number(source.field, resolved.index)
    # The order of the taxa of square_cell_figure.jl.
    order = topology_compatible_leaf_order(selected.tree)
    neighbour(taxon) = first(only(selected.tree.adjacency[taxon]))
    neighbour(4) == neighbour(5) && (order = replace(order, 4 => 5, 5 => 4))

    ray_trees = [g2_ray_tree(ray; source=source).tree for ray in rays]
    bounded_length === nothing &&
        (bounded_length = g2_triangle_bounded_length([selected.tree], order, G2_TRIANGLE_BOUNDED_LENGTH))
    ray_bounded_length === nothing &&
        (ray_bounded_length = g2_triangle_bounded_length(ray_trees, order, G2_TRIANGLE_RAY_BOUNDED_LENGTH))
    ray_options = (scale=ray_scale * scale, taxon_scale=ray_scale * taxon_scale,
                   bounded_length=ray_bounded_length)
    ray_coordinates = [g2_minimal_tree_coordinates(tree; scale=ray_options.scale,
                           bounded_length=ray_options.bounded_length,
                           layout_coordinates=g2_heptagon_layout(tree, order))
                       for tree in ray_trees]
    ray_letters, letters, record = g2_triangle_letters(ray_trees, selected.tree, order, ray_coordinates)
    if cell_letters
        middle = g2_triangle_tree(selected.tree; order=order, scale=scale,
                                  taxon_scale=taxon_scale, bounded_length=bounded_length,
                                  letters=letters)
        # The rays put their letters on the sides of the edges that the cell chose.
        corners = [g2_triangle_tree(tree; order=order, letters=letters,
                                    preferred=middle.sides, ray_options...)
                   for (tree, letters) in zip(ray_trees, ray_letters)]
    else
        middle = g2_triangle_tree(selected.tree; order=order, scale=scale,
                                  taxon_scale=taxon_scale, bounded_length=bounded_length,
                                  letters=Dict{Tuple{Int,Int},String}())
        corners = [g2_triangle_tree(tree; order=order, letters=letters, at_vertices=true, ray_options...)
                   for (tree, letters) in zip(ray_trees, ray_letters)]
    end
    return (; rays, number, selected, order, middle, ray_trees, corners, letters,
            ray_letters, record, bounded_length, ray_bounded_length)
end

"""
    g2_triangle_cell_tex(cell; field=:qq, scale=1.0, taxon_scale=0.75,
                         bounded_length=nothing, ray_scale=0.8, ray_bounded_length=nothing,
                         name_below=false, standalone=false)

The figure of the triangular `cell` of the fan over `field`, given by its
index or its three rays. Its rays go clockwise from the top corner, in
increasing order; each is labeled outside its circle by its orbit in the ray
catalog, and the cell by its orbit in the tables of maximal trees, inside
the bottom side, or below it, between the circles, when `name_below`. The trees are G2-trees with the options `scale`,
`taxon_scale`, and `bounded_length` of `g2_tree_tikz`, without edge lengths,
all with their taxa where the tree of the cell has them. The trees of the
rays are smaller by the factor `ray_scale`, which multiplies their `scale`
and `taxon_scale`, and their bounded edges have the length
`ray_bounded_length`. Either length, when it is `nothing`, is
`G2_TRIANGLE_BOUNDED_LENGTH` or `G2_TRIANGLE_RAY_BOUNDED_LENGTH`, shortened
by `g2_triangle_bounded_length` for long trees. The letters are those of `g2_triangle_letters`. Every
corner has a circle of the same radius, and the triangle is the smallest
equilateral one whose sides between the circles are `G2_SQUARE_SIDE_MIN_CM`
long and that holds the tree of the cell clear of the circles, the sides,
and the label of the cell inside: with the label below, the triangle can be
smaller. Without `cell_inside`, the tree of the cell may reach past the
sides, and only keeps clear of the circles and of the label. Without
`cell_letters`, the tree of the cell has no letters, and those of a ray
stand at their vertices.
With a `height` in centimetres, the triangle has this height, if the trees
fit in it, instead of the smallest.
"""
function g2_triangle_cell_tex(cell; field::Symbol=:qq,
                              scale::Real=1.0, taxon_scale::Real=0.75, bounded_length=nothing,
                              ray_scale::Real=0.8, ray_bounded_length=nothing,
                              name_below::Bool=false, cell_letters::Bool=true, cell_inside::Bool=true,
                              height=nothing, standalone::Bool=false, source=nothing)
    (; rays, number, selected, order, middle, corners, record, bounded_length, ray_bounded_length) =
        g2_triangle_trees(cell; field, scale, taxon_scale, bounded_length, ray_scale,
                          ray_bounded_length, cell_letters, source)
    length(rays) == 3 || throw(ArgumentError("$(Tuple(rays)) is a square, not a triangle"))

    radius = maximum(corner -> corner.radius, corners) + G2_SQUARE_CIRCLE_PAD_CM
    # The corners at the distance `circumradius` from the centre, the top one first.
    directions = [(cos(a), sin(a)) for a in (pi / 2, -pi / 6, 7pi / 6)]
    name_width, name_height = G2_SQUARE_CELL_NAME_CM
    function clear(circumradius)
        sqrt(3) * circumradius - 2radius >= G2_SQUARE_SIDE_MIN_CM || return false
        bottom = -circumradius / 2
        # The top of the label of the cell, inside the bottom side or below it.
        name_top = name_below ? bottom - G2_SQUARE_NAME_GAP_CM :
                   bottom + G2_SQUARE_NAME_GAP_CM + name_height
        all(middle.points) do point
            # The sides are at the distance circumradius / 2 from the centre,
            # opposite the corners.
            (!cell_inside || all(-(point[1] * d[1] + point[2] * d[2]) <= circumradius / 2 - G2_SQUARE_CENTRE_GAP_CM
                                 for d in directions)) &&
                all(g2_simple_distance(point, circumradius .* d) >= radius + G2_SQUARE_CENTRE_GAP_CM
                    for d in directions) &&
                (name_below && cell_inside || abs(point[1]) > name_width / 2 + G2_SQUARE_CENTRE_GAP_CM ||
                 point[2] >= name_top + G2_SQUARE_CENTRE_GAP_CM)
        end
    end
    if height === nothing
        circumradius = ceil(2radius * 100) / 100
        while !clear(circumradius)
            circumradius += 0.01
        end
    else
        # The height of an equilateral triangle is 3/2 of its circumradius.
        circumradius = 2height / 3
        clear(circumradius) || throw(ArgumentError("the trees do not fit in a triangle of height $height cm"))
    end
    places = [circumradius .* d for d in directions]
    point((x, y)) = "($(g2_simple_number(x)),$(g2_simple_number(y)))"

    tikz = sprint() do io
        println(io, "% Bounded edges of length $(g2_simple_number(bounded_length)) in the cell and " *
                    "$(g2_simple_number(ray_bounded_length)) in the rays, before scaling.")
        println(io, "% Letters of the inner vertices of the rays, in the order of the rays:")
        foreach(line -> println(io, "%   ", line), record)
        println(io, raw"\begingroup")
        println(io, "\\definecolor{g2RayLabel}{HTML}{$G2_SQUARE_RAY_COLOR}")
        println(io, "\\definecolor{g2CellLabel}{HTML}{$G2_SQUARE_CELL_COLOR}")
        println(io, G2_SQUARE_FILL)
        println(io, raw"""\tikzset{
  g2square side/.style={draw=g2RayLabel,line width=.8pt},
  g2square circle/.style={draw=g2RayLabel,fill=white,line width=.8pt},
  g2square cell/.style={fill=g2SquareFill},
  g2square name/.style={text=g2RayLabel,inner sep=0pt,font=\boldmath},
  g2square cell name/.style={g2square name,text=g2CellLabel},""")
        println(io, "  g2triangle letter/.style={text=g2RayLabel,inner sep=0pt," *
                    "font={\\sffamily\\bfseries\\fontsize{$G2_TRIANGLE_LETTER_FONT_PT}" *
                    "{$G2_TRIANGLE_LETTER_FONT_PT}\\selectfont}},")
        println(io, "}")
        println(io, raw"\begin{tikzpicture}[x=1cm,y=1cm,baseline=(current bounding box.center)]")
        triangle = join(point.(places), "--") * "--cycle;"
        println(io, "\\fill[g2square cell] ", triangle)
        println(io, "\\draw[g2square side] ", triangle)
        for (ray, place, direction, corner) in zip(rays, places, directions, corners)
            println(io, "% Ray $ray")
            println(io, "\\draw[g2square circle] $(point(place)) " *
                        "circle[radius=$(g2_simple_number(radius))];")
            # The label touches the circle on the line from the centre, away from the triangle.
            anchor = direction[2] > 0.9 ? "south" : direction[1] > 0 ? "north west" : "north east"
            println(io, "\\node[g2square name,anchor=$anchor] at " *
                        "$(point(place .+ (radius + G2_SQUARE_NAME_GAP_CM) .* direction)) {" *
                        g2_square_table_font("\\($(g2_square_ray_type(ray))\\)") * "};")
            println(io, "\\begin{scope}[shift={$(point(place))}]")
            corner.write(io)
            println(io, raw"\end{scope}")
        end
        println(io, "% The cell")
        name_y = name_below ? -circumradius / 2 - G2_SQUARE_NAME_GAP_CM : -circumradius / 2 + G2_SQUARE_NAME_GAP_CM
        println(io, "\\node[g2square cell name,anchor=$(name_below ? "north" : "south")] at " *
                    "$(point((0.0, name_y))) {" *
                    g2_square_table_font("\\oldstylenums{$number}") * "};")
        middle.write(io)
        # A document draws over the figure, in its coordinates with the centre
        # of the triangle at the origin, by defining this macro.
        println(io, raw"\ifdefined\gtwocelloverlay\gtwocelloverlay\fi")
        println(io, raw"\end{tikzpicture}")
        println(io, raw"\endgroup")
    end
    description = selected.description * ", orbit $number, with its rays clockwise from the top " *
                  "corner; G2-trees with lettered inner vertices, generated by " *
                  "julia/to_tex/triangle_cell_figure.jl"
    return g2_wrap_tree_tikz(tikz, description; standalone=standalone)
end

function g2_triangle_cell_main(args=ARGS)
    script = "julia --project=. to_tex/triangle_cell_figure.jl"
    if "--help" in args
        println("Usage: $script [--orbit=N | --rays=I,J,K] [--field=QQ|GF2] [--scale=FACTOR] " *
                "[--taxon-scale=FACTOR] [--bounded-length=LENGTH] [--ray-scale=FACTOR] " *
                "[--ray-bounded-length=LENGTH] [--name-below] [--no-cell-letters] [--cell-overhang] " *
                "[--height=CM] " *
                "[--standalone] [--output=FILE]")
        return nothing
    end
    valued = ("orbit", "rays", "field", "scale", "taxon-scale", "bounded-length", "ray-scale",
              "ray-bounded-length", "height", "output")
    for arg in args
        arg in ("--standalone", "--name-below", "--no-cell-letters", "--cell-overhang") ||
            any(name -> startswith(arg, "--$name="), valued) ||
            throw(ArgumentError("unknown option $arg; use --help"))
    end
    g2_tree_option(args, "output") == "" && throw(ArgumentError("output path must not be empty"))
    orbit, rays = g2_tree_option(args, "orbit"), g2_tree_option(args, "rays")
    orbit !== nothing && rays !== nothing && throw(ArgumentError("pass --orbit or --rays, not both"))
    field = g2_tree_option(args, "field")
    field === nothing || (field = g2_tree_field(field))
    if rays === nothing
        field, cell = g2_triangle_orbit_support(parse(Int, something(orbit, string(G2_TRIANGLE_ORBIT)));
                                                field=field)
    else
        field, cell = something(field, :qq), Tuple(parse.(Int, split(rays, ',')))
    end
    number(name, default) = parse(Float64, something(g2_tree_option(args, name), default))
    length_option(name) = (value = g2_tree_option(args, name); value === nothing ? nothing : parse(Float64, value))
    tex = g2_triangle_cell_tex(
        cell;
        field=field,
        scale=number("scale", "1.0"),
        taxon_scale=number("taxon-scale", "0.75"),
        bounded_length=length_option("bounded-length"),
        ray_scale=number("ray-scale", "0.8"),
        ray_bounded_length=length_option("ray-bounded-length"),
        name_below="--name-below" in args,
        cell_letters=!("--no-cell-letters" in args),
        cell_inside=!("--cell-overhang" in args),
        height=length_option("height"),
        standalone="--standalone" in args,
    )
    return g2_tree_emit(args, tex)
end

if abspath(PROGRAM_FILE) == @__FILE__
    g2_triangle_cell_main()
end
