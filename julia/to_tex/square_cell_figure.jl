###############################################################################
# A quadrilateral cell of the tropical G2-Grassmannian over QQ: a square whose
# corners are circles holding the G2-trees of its four rays, in their cyclic
# order around the cell, with the G2-tree of the cell in its middle. The trees
# are those of render_maximal_cell_tree.jl, with the taxon shapes, the names
# of the rays and of the cell in black, as the labels of main.tex, and every
# tree has its taxa where the tree of the cell has them; with
# --simple they are the simple trees of render_simple_tree.jl, which this file
# includes, and an edge longer than max_lines is one line with its length
# beside it. The rays are labeled by their orbits in the ray catalog of
# build_styled_ray_orbit_tables.jl, which this file includes as well, and the
# cell by the number of the orbit of the squares in the tables of maximal
# trees, read from the atlas by the helpers of toric_degeneration_table.jl.
#
#   julia --project=julia julia/to_tex/square_cell_figure.jl \
#     [--rays=I,J,K,L] [--scale=FACTOR] [--taxon-scale=FACTOR] \
#     [--bounded-length=LENGTH] [--ray-scale=FACTOR] [--ray-bounded-length=LENGTH] \
#     [--weight-font=PT] [--standalone] [--output=FILE]
#   julia --project=julia julia/to_tex/square_cell_figure.jl --simple \
#     [--rays=I,J,K,L] [--edge-length=LENGTH] [--max-lines=N] \
#     [--standalone] [--output=FILE]
#
# The simple trees place their labels in em: their fragment sets its font to
# G2_SQUARE_FONT_PT, the size of main.tex, for which the circles and the
# square are fitted to them.
###############################################################################

include(joinpath(@__DIR__, "build_styled_ray_orbit_tables.jl"))
include(joinpath(@__DIR__, "render_simple_tree.jl"))
include(joinpath(@__DIR__, "toric_degeneration_table.jl"))

const G2_SQUARE_FONT_PT = 12
const G2_SQUARE_EM_CM = G2_SQUARE_FONT_PT / G2_SIMPLE_PT_PER_CM
# The colors of the names of the rays and of the maximal cells, and of the
# circles and the sides of the polygon: black, as the ray labels (\raytype)
# and the cell labels (\celltype) of main.tex.
const G2_SQUARE_RAY_COLOR = "000000"
const G2_SQUARE_CELL_COLOR = "000000"
# The fill of the polygon: a light tint of the blue of the taxa (g2TaxonBlue).
const G2_SQUARE_FILL = raw"""\definecolor{g2TaxonBlue}{HTML}{1676A3}
\colorlet{g2SquareFill}{g2TaxonBlue!10}"""
# How far a taxon symbol reaches beyond its place on the heptagon, in
# millimetres at symbol scale 1: the corners of the triangles.
const G2_SQUARE_SYMBOL_REACH_MM = 3.8
# Between a circle and the label of its ray, in centimetres.
const G2_SQUARE_NAME_GAP_CM = 0.08
# The label of the cell inside the bottom side of the square: the width and
# the height that it does not exceed, in centimetres at 12pt.
const G2_SQUARE_CELL_NAME_CM = (0.6, 0.5)
# Between a tree and its circle, in centimetres.
const G2_SQUARE_CIRCLE_PAD_CM = 0.15
# The shortest visible side of the square, between two circles.
const G2_SQUARE_SIDE_MIN_CM = 1.0
# Between the tree in the middle and the circles and the sides.
const G2_SQUARE_CENTRE_GAP_CM = 0.3
const G2_SQUARE_RAYS = (13, 14, 17, 18)

"""The smallest circle containing `points`, as its centre and its radius."""
function g2_enclosing_circle(points)
    inside(circle, point) = g2_simple_distance(circle[1], point) <= circle[2] + 1e-9
    diametral(p, q) = (((p[1] + q[1]) / 2, (p[2] + q[2]) / 2), g2_simple_distance(p, q) / 2)
    function through(p, q, r)
        determinant = 2 * (p[1] * (q[2] - r[2]) + q[1] * (r[2] - p[2]) + r[1] * (p[2] - q[2]))
        # Three points on a line: the circle on the two farthest apart.
        abs(determinant) < 1e-12 &&
            return argmax(circle -> circle[2], (diametral(p, q), diametral(q, r), diametral(p, r)))
        square(point) = point[1]^2 + point[2]^2
        centre = ((square(p) * (q[2] - r[2]) + square(q) * (r[2] - p[2]) +
                   square(r) * (p[2] - q[2])) / determinant,
                  (square(p) * (r[1] - q[1]) + square(q) * (p[1] - r[1]) +
                   square(r) * (q[1] - p[1])) / determinant)
        return (centre, g2_simple_distance(centre, p))
    end
    circle = (first(points), 0.0)
    for i in eachindex(points)
        inside(circle, points[i]) && continue
        circle = (points[i], 0.0)
        for j in 1:(i - 1)
            inside(circle, points[j]) && continue
            circle = diametral(points[i], points[j])
            for k in 1:(j - 1)
                inside(circle, points[k]) || (circle = through(points[i], points[j], points[k]))
            end
        end
    end
    return circle
end

"""
`tex` in the font of the labels of the orbits in the tables of trees, as
`styled_tree_label` sets them: bold, and thickened by an outline, which
needs pdfrender.
"""
g2_square_table_font(tex) = "\\textpdfrender{TextRenderingMode=FillStroke," *
                            "LineWidth=$STYLED_TREE_LABEL_OUTLINE}{\\bfseries\\boldmath $tex}"

"""The label of the orbit of `ray` in the ray catalog, such as `\\gamma_1`."""
function g2_square_ray_type(ray)
    # The eleven orbits of the rays over QQ, and that of the rays 52:54, as in
    # styled_ray_orbit_data.
    orbits = push!(load_index_orbits(RAY_ORBITS_CSV, :ray), collect(52:54))
    return paper_ray_article_label(findfirst(orbit -> ray in orbit, orbits))
end

"""
The number of the orbit of the squares in the tables of maximal trees: the
squares of both fields form one orbit, the only orbit of squares of the
atlas.
"""
function g2_square_orbit_number()
    numbers = unique(
        toric_table_label_number(row["label"])
        for row in toric_table_csv_rows(TORIC_TABLE_ATLAS)
        if row["kind"] == "square" &&
           all(field -> field in ("qq", "gf2"), split(row["fields"], '+'))
    )
    length(numbers) == 1 || error("the atlas has $(length(numbers)) orbits of squares, not one")
    return only(numbers)
end

"""The distances between the taxa of `tree`, as a dictionary by pairs."""
function g2_square_taxon_distances(tree)
    distances = Dict{Tuple{Int,Int},QQFieldElem}()
    for taxon in 1:tree.ntaxa
        reached = Dict(taxon => QQ(0))
        pending = [taxon]
        while !isempty(pending)
            vertex = pop!(pending)
            for (neighbor, length) in tree.adjacency[vertex]
                haskey(reached, neighbor) && continue
                reached[neighbor] = reached[vertex] + length
                push!(pending, neighbor)
            end
        end
        foreach(other -> distances[(taxon, other)] = reached[other], 1:tree.ntaxa)
    end
    return distances
end

"""
The tree of a ray of a cell drawn on `cell_tree`, the tree of the cell: the
same vertices and edges, the pendant lengths of `ray_tree`, and on every
bounded edge the length of its split in the ray, zero for a split of the
cell that the ray lacks, whose ends the ray merges. The splits of a ray of a
cell are splits of the cell; the result is certified to have the metric of
the ray.
"""
function g2_square_ray_on_cell(ray_tree, cell_tree)
    ntaxa = cell_tree.ntaxa
    weights = Dict(edge.cluster => edge.length for edge in rooted_internal_edges(ray_tree))
    pendant(taxon) = last(only(ray_tree.adjacency[taxon]))
    found = Set{Vector{Int}}()
    adjacency = Dict{Int,Vector{Tuple{Int,QQFieldElem}}}()
    for (vertex, neighbors) in cell_tree.adjacency, (neighbor, _) in neighbors
        length = if vertex <= ntaxa || neighbor <= ntaxa
            pendant(min(vertex, neighbor))
        else
            # The side of the edge away from taxon 7, as in rooted_internal_edges.
            side = sort(leaves_on_side(cell_tree, vertex, neighbor))
            7 in side && (side = setdiff(1:ntaxa, side))
            push!(found, side)
            get(weights, side, QQ(0))
        end
        push!(get!(adjacency, vertex, Tuple{Int,QQFieldElem}[]), (neighbor, length))
    end
    issubset(keys(weights), found) || error("a split of the ray is not a split of the cell")
    tree = TreeMetric(ntaxa, adjacency)
    g2_square_taxon_distances(tree) == g2_square_taxon_distances(ray_tree) ||
        error("the ray drawn on the tree of the cell lost its metric")
    return tree
end

"""The four rays of the square `cell` in their cyclic order, from the smallest."""
function g2_square_cycle(source, cell)
    rays = g2_resolve_maximal_cell(source, cell).rays
    length(rays) == 4 || throw(ArgumentError("$rays is a triangle, not a square"))
    neighbours(ray) = sort!([other for other in rays if other != ray &&
                             Tuple(sort!([ray, other])) in source.edges])
    all(ray -> length(neighbours(ray)) == 2, rays) ||
        error("the sides of the square $rays are not edges of the fan")
    cycle = [first(rays), first(neighbours(first(rays)))]
    while length(cycle) < 4
        push!(cycle, only(filter(!=(cycle[end - 1]), neighbours(cycle[end]))))
    end
    return cycle
end

"""
A G2-tree with its taxon shapes, drawn by `g2_tree_commands` with the taxa
in the cyclic `order` of `g2_heptagon_layout`: a function that writes it to an `io`, the centre and
the radius of a circle around it, and points on that circle. The taxa lie on
a heptagon around the origin, so the circle is the same for every tree. The
lengths are set in `weight_font_pt` points, on disks of the TikZ color
`disk`, each on its edge.
"""
function g2_square_shape_tree(tree; order, scale, taxon_scale, bounded_length,
                              weight_font_pt, disk="white")
    radius = G2_TREE_UNIT_CM * G2_TREE_LEAF_RADIUS * scale +
             G2_SQUARE_SYMBOL_REACH_MM / 10 * G2_TREE_SYMBOL_SCALE * taxon_scale
    layout = g2_heptagon_layout(tree, order)
    function write(io, _)
        colors = g2_minimal_tree_styles(io, :color)
        println(io, "\\tikzset{g2tree weight/.append style={minimum size=" *
                    "$(g2_simple_number(weight_font_pt + 1))pt" *
                    (disk == "white" ? "" : ",fill=$disk") * "}}")
        g2_tree_commands(io, tree, colors; scale=scale, taxon_scale=taxon_scale,
                         bounded_length=bounded_length, layout_coordinates=layout,
                         weight_font_pt=weight_font_pt, labels_beside=false)
        println(io, raw"\endgroup")
    end
    points = [(radius * cos(angle), radius * sin(angle)) for angle in range(0, 2pi; length=73)]
    return (write=write, centre=(0.0, 0.0), radius=radius, points=points)
end

"""
A simple tree, as `g2_square_shape_tree` gives a G2-tree: taxon 1 is left of
taxon 4, in a mirror image if need be, and the circle is the smallest around
the drawing for labels of `G2_SQUARE_EM_CM`. The writer takes the prefix of
the names of the nodes.
"""
function g2_square_simple_tree(tree; edge_length, max_lines)
    drawing = g2_simple_tree_drawing(tree; edge_length=edge_length, max_lines=max_lines)
    positions = g2_simple_positions(drawing, G2_SQUARE_EM_CM)
    x(taxon) = positions[findfirst(taxa -> taxon in taxa, drawing.graph.taxa)][1]
    x(1) > x(4) + 1e-9 && (drawing = g2_simple_tree_mirror(drawing))
    points = g2_simple_tree_outline(drawing, G2_SQUARE_EM_CM)
    centre, radius = g2_enclosing_circle(points)
    write(io, prefix) = g2_simple_tree_commands(io, drawing; prefix=prefix)
    return (write=write, centre=centre, radius=radius, points=points)
end

"""
    g2_square_cell_tex(cell=G2_SQUARE_RAYS; simple=false, scale=1.0, taxon_scale=0.75,
                       bounded_length=1.1, ray_scale=0.8, ray_bounded_length=1.1,
                       weight_font_pt=10, edge_length=0.8, max_lines=3, standalone=false)

The figure of the square `cell` of the fan over QQ, given by its index or its
four rays. Its rays go clockwise from the top left corner, in their cyclic
order from the smallest; each is labeled outside its circle by its orbit in
the ray catalog. The trees are G2-trees with the options `scale`,
`taxon_scale`, and `bounded_length` of `g2_tree_tikz`, all with their taxa
where the tree of the cell has them, and with their edge lengths in
`weight_font_pt` points. Every ray is drawn on the vertices and edges of the
cell, by `g2_square_ray_on_cell`. The trees of the rays are smaller by the
factor `ray_scale`, which multiplies their `scale` and `taxon_scale`, and
their bounded edges have the length `ray_bounded_length` instead. With
`simple=true` they are
simple trees with the options `edge_length` and `max_lines` of
`g2_simple_tree_tikz`. Every corner has a circle of the same radius, the
smallest that holds each of the four trees, and the square is the smallest
whose sides between the circles are `G2_SQUARE_SIDE_MIN_CM` long and that
holds the tree of the cell clear of the circles and the sides.
"""
function g2_square_cell_tex(cell=G2_SQUARE_RAYS; simple::Bool=false,
                            scale::Real=1.0, taxon_scale::Real=0.75, bounded_length::Real=1.1,
                            ray_scale::Real=0.8, ray_bounded_length::Real=1.1,
                            weight_font_pt::Real=10, edge_length::Real=0.8, max_lines=3,
                            standalone::Bool=false, source=nothing)
    source === nothing && (source = load_g2_maximal_cell_source(; field=:qq))
    source.field == :qq || throw(ArgumentError("the square cell figure is drawn over QQ"))
    cycle = g2_square_cycle(source, cell)
    selected = g2_maximal_cell_tree(cycle; source=source)
    # The splits of the rays are splits of the cell, so their trees are planar
    # with the taxa of the cell.
    order = topology_compatible_leaf_order(selected.tree)
    # The yellow square and triangle trade places when they form a cherry, so
    # that the yellow taxa mirror the blue ones.
    neighbour(taxon) = first(only(selected.tree.adjacency[taxon]))
    neighbour(4) == neighbour(5) && (order = replace(order, 4 => 5, 5 => 4))
    draw(tree; factor=1.0, bounded=bounded_length, disk="white") = simple ?
        g2_square_simple_tree(tree; edge_length=edge_length, max_lines=max_lines) :
        g2_square_shape_tree(tree; order=order, scale=factor * scale,
                             taxon_scale=factor * taxon_scale, bounded_length=bounded,
                             weight_font_pt=weight_font_pt, disk=disk)
    # Every ray has the vertices and edges of the cell.
    corners = [draw(g2_square_ray_on_cell(g2_ray_tree(ray; source=source).tree, selected.tree);
                    factor=ray_scale, bounded=ray_bounded_length)
               for ray in cycle]
    # The tree of the cell lies on the square, and so do its disks.
    middle = draw(selected.tree; disk="g2SquareFill")

    radius = maximum(corner -> corner.radius, corners) + G2_SQUARE_CIRCLE_PAD_CM
    middle_points = [(x - middle.centre[1], y - middle.centre[2]) for (x, y) in middle.points]
    # Half the side of the square, by steps of a millimetre.
    # The tree is also clear of the label of the cell, at the middle of the
    # bottom side.
    name_width, name_height = G2_SQUARE_CELL_NAME_CM
    clear(half) = all(middle_points) do point
        max(abs(point[1]), abs(point[2])) <= half - G2_SQUARE_CENTRE_GAP_CM &&
            all(g2_simple_distance(point, (x * half, y * half)) >=
                radius + G2_SQUARE_CENTRE_GAP_CM for x in (-1, 1) for y in (-1, 1)) &&
            (abs(point[1]) > name_width / 2 + G2_SQUARE_CENTRE_GAP_CM ||
             point[2] >= -half + G2_SQUARE_NAME_GAP_CM + name_height + G2_SQUARE_CENTRE_GAP_CM)
    end
    half = ceil((radius + G2_SQUARE_SIDE_MIN_CM / 2) * 100) / 100
    while !clear(half)
        half += 0.01
    end
    places = [(-half, half), (half, half), (half, -half), (-half, -half)]
    point((x, y)) = "($(g2_simple_number(x)),$(g2_simple_number(y)))"

    tikz = sprint() do io
        println(io, raw"\begingroup")
        println(io, "\\definecolor{g2RayLabel}{HTML}{$G2_SQUARE_RAY_COLOR}")
        println(io, "\\definecolor{g2CellLabel}{HTML}{$G2_SQUARE_CELL_COLOR}")
        println(io, G2_SQUARE_FILL)
        if simple
            println(io, "\\fontsize{$(G2_SQUARE_FONT_PT)pt}" *
                        "{$(round(1.2G2_SQUARE_FONT_PT; digits=1))pt}\\selectfont")
            println(io, g2_simple_tree_styles())
        end
        println(io, raw"""\tikzset{
  g2square side/.style={draw=g2RayLabel,line width=.8pt},
  g2square circle/.style={draw=g2RayLabel,fill=white,line width=.8pt},
  g2square cell/.style={fill=g2SquareFill},
  g2square name/.style={text=g2RayLabel,inner sep=0pt,font=\boldmath},
  g2square cell name/.style={g2square name,text=g2CellLabel},
}""")
        println(io, raw"\begin{tikzpicture}[x=1cm,y=1cm,baseline=(current bounding box.center)]")
        square = join(point.(places), "--") * "--cycle;"
        println(io, "\\fill[g2square cell] ", square)
        println(io, "\\draw[g2square side] ", square)
        for (ray, place, corner) in zip(cycle, places, corners)
            println(io, "% Ray $ray")
            println(io, "\\draw[g2square circle] $(point(place)) " *
                        "circle[radius=$(g2_simple_number(radius))];")
            # The label touches the circle on the diagonal, away from the square.
            outward = (radius + G2_SQUARE_NAME_GAP_CM) / sqrt(2) .* sign.(place)
            anchor = (place[2] > 0 ? "south" : "north") * " " * (place[1] > 0 ? "west" : "east")
            println(io, "\\node[g2square name,anchor=$anchor] at $(point(place .+ outward)) {" *
                        g2_square_table_font("\\($(g2_square_ray_type(ray))\\)") * "};")
            println(io, "\\begin{scope}[shift={$(point(place .- corner.centre))}]")
            corner.write(io, "r$(ray)n")
            println(io, raw"\end{scope}")
        end
        println(io, "% The cell")
        println(io, "\\node[g2square cell name,anchor=south] at " *
                    "$(point((0.0, -half + G2_SQUARE_NAME_GAP_CM))) {" *
                    g2_square_table_font("\\oldstylenums{$(g2_square_orbit_number())}") * "};")
        println(io, "\\begin{scope}[shift={$(point(.-(middle.centre)))}]")
        middle.write(io, "cn")
        println(io, raw"\end{scope}")
        println(io, raw"\end{tikzpicture}")
        println(io, raw"\endgroup")
    end
    description = selected.description * ", with its rays clockwise from the top left " *
                  "corner; " * (simple ? "simple drawings" : "G2-trees") *
                  ", generated by julia/to_tex/square_cell_figure.jl"
    return g2_wrap_tree_tikz(tikz, description; standalone=standalone)
end

function g2_square_cell_main(args=ARGS)
    script = "julia --project=julia julia/to_tex/square_cell_figure.jl"
    if "--help" in args
        println("Usage: $script [--rays=I,J,K,L] [--scale=FACTOR] [--taxon-scale=FACTOR] " *
                "[--bounded-length=LENGTH] [--ray-scale=FACTOR] " *
                "[--ray-bounded-length=LENGTH] [--weight-font=PT] [--standalone] [--output=FILE]")
        println("       $script --simple [--rays=I,J,K,L] [--edge-length=LENGTH] " *
                "[--max-lines=N] [--standalone] [--output=FILE]")
        return nothing
    end
    simple = "--simple" in args
    valued = simple ? ("rays", "edge-length", "max-lines", "output") :
                      ("rays", "scale", "taxon-scale", "bounded-length", "ray-scale",
                       "ray-bounded-length", "weight-font", "output")
    for arg in args
        arg in ("--simple", "--standalone") ||
            any(name -> startswith(arg, "--$name="), valued) ||
            throw(ArgumentError("unknown option $arg; use --help"))
    end
    g2_tree_option(args, "output") == "" &&
        throw(ArgumentError("output path must not be empty"))
    rays = g2_tree_option(args, "rays")
    cell = rays === nothing ? G2_SQUARE_RAYS : Tuple(parse.(Int, split(rays, ',')))
    number(name, default) = parse(Float64, something(g2_tree_option(args, name), default))
    tex = g2_square_cell_tex(
        cell;
        simple=simple,
        scale=number("scale", "1.0"),
        taxon_scale=number("taxon-scale", "0.75"),
        bounded_length=number("bounded-length", "1.1"),
        ray_scale=number("ray-scale", "0.8"),
        ray_bounded_length=number("ray-bounded-length", "1.1"),
        weight_font_pt=number("weight-font", "10"),
        edge_length=number("edge-length", "0.8"),
        max_lines=parse(Int, something(g2_tree_option(args, "max-lines"), "3")),
        standalone="--standalone" in args,
    )
    return g2_tree_emit(args, tex)
end

if abspath(PROGRAM_FILE) == @__FILE__
    g2_square_cell_main()
end
