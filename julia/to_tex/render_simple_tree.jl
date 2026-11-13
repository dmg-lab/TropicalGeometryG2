###############################################################################
# Simple drawings of the trees of the tropical G2 Grassmannian. Taxa are
# written as numbers. An edge of length n is drawn as n parallel lines, or as
# one line with n beside it when n exceeds max_lines, and an edge of length
# zero is contracted: taxa at distance zero share one label.
#
# The trees are selected as in render_maximal_cell_tree.jl, which this file
# includes. The public entry points are g2_simple_tree_tikz and
# g2_simple_tree_tex; g2_simple_tree_drawing and g2_simple_tree_commands put
# several trees in one picture.
###############################################################################

# Unless a script that includes this file has included the renderer already.
isdefined(@__MODULE__, :g2_tree_tikz) ||
    include(joinpath(@__DIR__, "render_maximal_cell_tree.jl"))

const G2_SIMPLE_LINE_WIDTH_PT = 0.6
# Distance between the centres of two adjacent parallel lines.
const G2_SIMPLE_LINE_GAP_PT = 2.0
# The lines of an edge of negative length, which a lineality star may have.
const G2_SIMPLE_NEGATIVE_DASH = "on 2.4pt off 1.6pt"
# Between a label and the lines that end at it.
const G2_SIMPLE_LABEL_GAP_PT = 0.5
const G2_SIMPLE_PT_PER_CM = 72.27 / 2.54
# The assumed size of a label, in em: the width and the height of a digit,
# and the margin of the label, which the fragment sets.
const G2_SIMPLE_DIGIT_WIDTH_EM = 0.5
const G2_SIMPLE_DIGIT_HEIGHT_EM = 0.65
const G2_SIMPLE_LABEL_SEP_EM = 0.17
# The length written beside an edge drawn as one line: its font, the size of
# that font in em, and the gap between it and the edge.
const G2_SIMPLE_LENGTH_FONT = raw"\scriptsize"
const G2_SIMPLE_LENGTH_SIZE_EM = 0.7
const G2_SIMPLE_LENGTH_GAP_EM = 0.2
# The em for which the order of the edges around the nodes is chosen. The
# fragment itself places the labels by the em of its document.
const G2_SIMPLE_EM_CM = 10 / G2_SIMPLE_PT_PER_CM

"""
The tree with its edges of length zero contracted. `taxa[node]` lists the
taxa on a node, which are at distance zero from each other; it is empty for an
inner node without taxa. `edges` holds the other edges `(node, node, length)`.
"""
function g2_simple_tree_graph(tree::TreeMetric)
    vertices = sort!(collect(keys(tree.adjacency)))
    representative = Dict(vertex => vertex for vertex in vertices)
    function find(vertex)
        while representative[vertex] != vertex
            vertex = representative[vertex]
        end
        return vertex
    end
    for vertex in vertices, (neighbour, edge_length) in tree.adjacency[vertex]
        left, right = find(vertex), find(neighbour)
        iszero(edge_length) && left != right &&
            (representative[max(left, right)] = min(left, right))
    end
    # A node is named by its smallest vertex, so nodes with taxa come first.
    names = sort!(unique!(find.(vertices)))
    node = Dict(name => index for (index, name) in enumerate(names))
    taxa = [Int[] for _ in names]
    for taxon in 1:tree.ntaxa
        push!(taxa[node[find(taxon)]], taxon)
    end
    edges = Tuple{Int,Int,QQFieldElem}[]
    for vertex in vertices, (neighbour, edge_length) in tree.adjacency[vertex]
        vertex < neighbour && !iszero(edge_length) || continue
        push!(edges, (node[find(vertex)], node[find(neighbour)], QQ(edge_length)))
    end
    return (taxa=taxa, edges=edges)
end

g2_simple_distance(p, q) = hypot(p[1] - q[1], p[2] - q[2])

"""Whether the segments `ab` and `cd`, which share no end, meet."""
function g2_simple_segments_meet(a, b, c, d)
    side(p, q, r) = (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])
    opposite(x, y) = (x > 1e-9 && y < -1e-9) || (x < -1e-9 && y > 1e-9)
    return opposite(side(a, b, c), side(a, b, d)) && opposite(side(c, d, a), side(c, d, b))
end

"""
The directions of the edges of a rooted tree: for every node but the root,
the angle of the edge from its parent. `children` lists, for every node, its
neighbours other than `parents[node]`, counterclockwise from the edge to the
parent.

Around a node every edge takes an angular sector proportional to its weight
and points along the bisector of its sector. The weight of an edge to an
outer node is the number of taxa on that node, and it is one for any other
edge. The edges at a node are thus equally spaced when no outer node carries
several taxa.
"""
function g2_simple_angles(root, parents, children, weight)
    angles = Dict{Int,Float64}()
    # A node with the angle of the edge to its parent; the first edge of the
    # root points upwards.
    pending = [(root, pi / 2)]
    while !isempty(pending)
        node, base = pop!(pending)
        before = node == root ? -weight(first(children[node])) / 2 :
                 weight(parents[node]) / 2
        total = sum(weight, children[node]; init=0.0) +
                (node == root ? 0.0 : weight(parents[node]))
        for child in children[node]
            angles[child] = base + 2pi * (before + weight(child) / 2) / total
            before += weight(child)
            isempty(children[child]) || push!(pending, (child, angles[child] + pi))
        end
    end
    return angles
end

"""
How far the label of a node with `taxa` reaches from its centre in the
direction `angle`, in em: the label is a row of digits in a box with the
margin of the fragment. A node without taxa is a point.
"""
function g2_simple_label_reach(taxa, angle)
    isempty(taxa) && return 0.0
    half_width = length(taxa) * G2_SIMPLE_DIGIT_WIDTH_EM / 2 + G2_SIMPLE_LABEL_SEP_EM
    half_height = G2_SIMPLE_DIGIT_HEIGHT_EM / 2 + G2_SIMPLE_LABEL_SEP_EM
    return min(half_width / abs(cos(angle)), half_height / abs(sin(angle)))
end

"""
The angle by which to turn `points` so that their longest direction is
horizontal; zero if they have none, as for a regular star.
"""
function g2_simple_turn(points)
    count = length(points)
    centre = (sum(first, points) / count, sum(last, points) / count)
    xx = sum(point -> (point[1] - centre[1])^2, points)
    yy = sum(point -> (point[2] - centre[2])^2, points)
    xy = sum(point -> (point[1] - centre[1]) * (point[2] - centre[2]), points)
    return hypot(xx - yy, 2xy) > 1e-9 * count ? -atan(2xy, xx - yy) / 2 : 0.0
end

"""
How well a drawing separates its parts: no two edges should cross, then the
smallest distance between a node and another node or an edge it does not end
should be large, then the nodes should be spread out.
"""
function g2_simple_layout_score(positions, edges)
    nodes = sort!(collect(keys(positions)))
    crossings = count(
        !(one[1] in other[1:2]) && !(one[2] in other[1:2]) && g2_simple_segments_meet(
            positions[one[1]], positions[one[2]],
            positions[other[1]], positions[other[2]],
        )
        for (index, one) in enumerate(edges) for other in edges[(index + 1):end]
    )
    clearance = Inf
    spread = 0.0
    for (index, left) in enumerate(nodes), right in nodes[(index + 1):end]
        separation = g2_simple_distance(positions[left], positions[right])
        clearance = min(clearance, separation)
        spread += separation
    end
    for node in nodes, (left, right, _) in edges
        node in (left, right) && continue
        clearance = min(clearance, g2_segment_distance(
            positions[node], positions[left], positions[right],
        ))
    end
    return (-crossings, round(clearance; digits=6), round(spread; digits=6))
end

"""
Places of the nodes of `graph`. The centre of a node lies at
`edge_length * spans[node]`, plus `reaches[node]` in em: every edge is then
`edge_length` long between the labels, or the nodes without taxa, at its ends.

The angles at every node are those of `g2_simple_angles`. Among the orders of
the edges around the nodes, the first one with the best
`g2_simple_layout_score` is taken, for labels of `G2_SIMPLE_EM_CM`. Each
drawing is turned so that its longest direction is horizontal.
"""
function g2_simple_tree_layout(graph, edge_length)
    count = length(graph.taxa)
    count == 1 && return (spans=[(0.0, 0.0)], reaches=[(0.0, 0.0)])
    neighbours = [Int[] for _ in 1:count]
    for (left, right, _) in graph.edges
        push!(neighbours[left], right)
        push!(neighbours[right], left)
    end
    foreach(sort!, neighbours)
    weight(node) = length(neighbours[node]) == 1 ? max(1, length(graph.taxa[node])) : 1

    # Root the tree at a centre: a node from which the farthest node is nearest.
    function depths(source)
        depth = Dict(source => 0)
        queue = [source]
        while !isempty(queue)
            node = popfirst!(queue)
            for neighbour in neighbours[node]
                haskey(depth, neighbour) && continue
                depth[neighbour] = depth[node] + 1
                push!(queue, neighbour)
            end
        end
        return depth
    end
    root = argmin(node -> (maximum(values(depths(node))), node), 1:count)
    depth = depths(root)
    outward = sort!(filter(!=(root), collect(1:count)); by=node -> (depth[node], node))
    parents = Dict(
        node => only(neighbour for neighbour in neighbours[node]
                     if depth[neighbour] < depth[node])
        for node in outward
    )
    inner = [node for node in 1:count if node == root || length(neighbours[node]) > 1]
    # Turning the root does not change the drawing, so its first edge is fixed.
    candidates = [
        node == root ?
        [vcat(first(neighbours[node]), order)
         for order in TropicalGeometryG2.all_permutations(neighbours[node][2:end])] :
        TropicalGeometryG2.all_permutations(filter(!=(parents[node]), neighbours[node]))
        for node in inner
    ]
    step(point, length, angle) = (point[1] + length * cos(angle), point[2] + length * sin(angle))
    best = nothing
    for orders in Iterators.product(candidates...)
        children = Dict(node => Int[] for node in 1:count)
        foreach((node, order) -> children[node] = order, inner, orders)
        angles = g2_simple_angles(root, parents, children, weight)
        # The longest direction is that of the tree with all its nodes as points.
        spans = Dict(root => (0.0, 0.0))
        foreach(node -> spans[node] = step(spans[parents[node]], 1.0, angles[node]), outward)
        turn = g2_simple_turn(collect(values(spans)))
        reaches = Dict(root => (0.0, 0.0))
        for node in outward
            angle = angles[node] + turn
            spans[node] = step(spans[parents[node]], 1.0, angle)
            reaches[node] = step(
                reaches[parents[node]],
                g2_simple_label_reach(graph.taxa[parents[node]], angle) +
                g2_simple_label_reach(graph.taxa[node], angle + pi),
                angle,
            )
        end
        positions = Dict(
            node => (edge_length * spans[node][1] + G2_SIMPLE_EM_CM * reaches[node][1],
                     edge_length * spans[node][2] + G2_SIMPLE_EM_CM * reaches[node][2])
            for node in 1:count
        )
        score = g2_simple_layout_score(positions, graph.edges)
        if best === nothing || score > best.score
            best = (score=score, spans=spans, reaches=reaches)
        end
    end
    best.score[1] == 0 || error("every drawing of this tree has crossing edges")
    return (spans=[best.spans[node] for node in 1:count],
            reaches=[best.reaches[node] for node in 1:count])
end

# Rounded for the fragment; a negative zero would print as -0.
g2_simple_number(value) = string(round(value; digits=4) + 0.0)

"""
    g2_simple_tree_drawing(tree::TreeMetric; edge_length=1.0, max_lines=nothing)

The drawing of `tree` that `g2_simple_tree_tikz` writes: the contracted tree
`graph` of `g2_simple_tree_graph`, the `spans` and `reaches` of its nodes as
in `g2_simple_tree_layout`, `edge_length`, and `max_lines`.
"""
function g2_simple_tree_drawing(tree::TreeMetric; edge_length::Real=1.0, max_lines=nothing)
    scale = Float64(edge_length)
    isfinite(scale) && scale > 0 ||
        throw(ArgumentError("edge_length must be a positive finite number"))
    max_lines === nothing || max_lines isa Integer && !(max_lines isa Bool) && max_lines >= 1 ||
        throw(ArgumentError("max_lines must be a positive integer"))
    graph = g2_simple_tree_graph(tree)
    for (_, _, metric) in graph.edges
        denominator(metric) == 1 || throw(ArgumentError(
            "the simple drawing needs integer edge lengths, " *
            "but an edge has length $(rational_key(metric))",
        ))
    end
    layout = g2_simple_tree_layout(graph, scale)
    return (graph=graph, spans=layout.spans, reaches=layout.reaches,
            edge_length=scale, max_lines=max_lines)
end

"""The mirror image of `drawing` in a vertical line."""
g2_simple_tree_mirror(drawing) = merge(drawing, (
    spans=[(-x, y) for (x, y) in drawing.spans],
    reaches=[(-x, y) for (x, y) in drawing.reaches],
))

"""Whether an edge of length `metric` is one line with its length beside it."""
g2_simple_numbered(drawing, metric) =
    drawing.max_lines !== nothing && abs(metric) > drawing.max_lines

"""The places of the nodes of `drawing`, in centimetres, for labels of `em_cm`."""
function g2_simple_positions(drawing, em_cm::Real=G2_SIMPLE_EM_CM)
    return [(drawing.edge_length * span[1] + em_cm * reach[1],
             drawing.edge_length * span[2] + em_cm * reach[2])
            for (span, reach) in zip(drawing.spans, drawing.reaches)]
end

"""Half the width and half the height, in em, of the length `metric` as written."""
function g2_simple_length_size(metric)
    characters = length(string(numerator(metric)))
    return (characters * G2_SIMPLE_DIGIT_WIDTH_EM * G2_SIMPLE_LENGTH_SIZE_EM / 2,
            G2_SIMPLE_DIGIT_HEIGHT_EM * G2_SIMPLE_LENGTH_SIZE_EM / 2)
end

"""
The lengths written beside the edges of `drawing` that are longer than
`max_lines`: for each, its edge, its length, and the `span` and `reach` of
its centre, as for a node. A length stands at the middle of the visible edge,
on the side where it is farther from the other edges and nodes for labels of
`em_cm`, and above the edge if both sides are as good.
"""
function g2_simple_length_labels(drawing, em_cm::Real=G2_SIMPLE_EM_CM)
    graph = drawing.graph
    positions = g2_simple_positions(drawing, em_cm)
    labels = []
    for (edge, (left, right, metric)) in enumerate(graph.edges)
        g2_simple_numbered(drawing, metric) || continue
        angle = atan(positions[right][2] - positions[left][2],
                     positions[right][1] - positions[left][1])
        # The middle of the visible edge, between the labels at its ends.
        shift = (g2_simple_label_reach(graph.taxa[left], angle) -
                 g2_simple_label_reach(graph.taxa[right], angle + pi)) / 2
        span = ((drawing.spans[left][1] + drawing.spans[right][1]) / 2,
                (drawing.spans[left][2] + drawing.spans[right][2]) / 2)
        middle = ((drawing.reaches[left][1] + drawing.reaches[right][1]) / 2 + shift * cos(angle),
                  (drawing.reaches[left][2] + drawing.reaches[right][2]) / 2 + shift * sin(angle))
        half_width, half_height = g2_simple_length_size(metric)
        best = nothing
        for side in (1, -1)
            normal = angle + side * pi / 2
            away = half_width * abs(cos(normal)) + half_height * abs(sin(normal)) +
                   G2_SIMPLE_LENGTH_GAP_EM
            reach = (middle[1] + away * cos(normal), middle[2] + away * sin(normal))
            centre = (drawing.edge_length * span[1] + em_cm * reach[1],
                      drawing.edge_length * span[2] + em_cm * reach[2])
            clearance = minimum(
                g2_simple_distance(centre, positions[node])
                for node in eachindex(positions) if !(node in (left, right));
                init=Inf,
            )
            for (other, (one, two, _)) in enumerate(graph.edges)
                other == edge && continue
                clearance = min(clearance, g2_segment_distance(
                    centre, positions[one], positions[two],
                ))
            end
            score = (round(clearance; digits=6), round(sin(normal); digits=6))
            (best === nothing || score > best.score) && (best = (score=score, reach=reach))
        end
        push!(labels, (edge=edge, metric=metric, span=span, reach=best.reach))
    end
    return labels
end

"""
Points around which a circle encloses `drawing`, in centimetres, for labels
of `em_cm`: the corners of the labels of the nodes and of the lengths, and
the nodes without taxa.
"""
function g2_simple_tree_outline(drawing, em_cm::Real=G2_SIMPLE_EM_CM)
    positions = g2_simple_positions(drawing, em_cm)
    points = Tuple{Float64,Float64}[]
    corners(centre, half_width, half_height) = [
        (centre[1] + x * em_cm * half_width, centre[2] + y * em_cm * half_height)
        for x in (-1, 1) for y in (-1, 1)
    ]
    for (node, taxa) in enumerate(drawing.graph.taxa)
        if isempty(taxa)
            push!(points, positions[node])
        else
            append!(points, corners(
                positions[node],
                length(taxa) * G2_SIMPLE_DIGIT_WIDTH_EM / 2 + G2_SIMPLE_LABEL_SEP_EM,
                G2_SIMPLE_DIGIT_HEIGHT_EM / 2 + G2_SIMPLE_LABEL_SEP_EM,
            ))
        end
    end
    for label in g2_simple_length_labels(drawing, em_cm)
        centre = (drawing.edge_length * label.span[1] + em_cm * label.reach[1],
                  drawing.edge_length * label.span[2] + em_cm * label.reach[2])
        append!(points, corners(centre, g2_simple_length_size(label.metric)...))
    end
    return points
end

"""
Where the parallel line at `offset` points to the left of an edge ends at
the node `name`, as a TikZ coordinate. `angle` is the direction of the edge
away from the node, and `left` and `right` are the angles to the next edges
around the node on either side.

At a node with taxa the lines end on a front across the edge, just outside
the label. At a node without taxa the edges share the space around the node
along the bisectors of their angles: a line ends on the bisector on its side,
where it meets the line of the neighbouring edge at the same distance from
its own edge, if there is one. The middle line of an edge reaches the node.
Returns the coordinate and how far along the edge the line ends.
"""
function g2_simple_line_end(name, labelled, angle, left, right, offset)
    along = labelled ? G2_SIMPLE_LABEL_GAP_PT :
            offset >= 0 ? offset / tan(left / 2) : -offset / tan(right / 2)
    shift = "(" * g2_simple_number(along * cos(angle) - offset * sin(angle)) * "pt," *
            g2_simple_number(along * sin(angle) + offset * cos(angle)) * "pt)"
    anchor = labelled ? "$name." * g2_simple_number(mod(rad2deg(angle), 360)) : name
    return ("([shift={$shift}]$anchor)", labelled ? 0.0 : along)
end

"""The TikZ styles that the commands of `g2_simple_tree_commands` use."""
function g2_simple_tree_styles()
    return "\\tikzset{\n" *
           "  g2simple taxa/.style={inner sep=$(G2_SIMPLE_LABEL_SEP_EM)em,outer sep=0pt},\n" *
           "  g2simple line/.style={line width=$(G2_SIMPLE_LINE_WIDTH_PT)pt,line cap=round},\n" *
           "  g2simple negative/.style={g2simple line,line cap=butt," *
           "dash pattern=$G2_SIMPLE_NEGATIVE_DASH},\n" *
           "  g2simple length/.style={font=$G2_SIMPLE_LENGTH_FONT,inner sep=0pt},\n}"
end

"""
Write the TikZ commands that draw `drawing` to `io`, for a picture in
centimetres with the styles of `g2_simple_tree_styles`. The nodes are named
`prefix` followed by their number, so that several trees can share a
picture.
"""
function g2_simple_tree_commands(io::IO, drawing; prefix::AbstractString="n")
    graph, scale = drawing.graph, drawing.edge_length
    positions = g2_simple_positions(drawing)
    # A length in centimetres plus a length in the em of the document.
    place(span, reach) = "{" * g2_simple_number(scale * span) * "cm" *
                         (reach < 0 ? "-" : "+") * g2_simple_number(abs(reach)) * "em}"
    at(span, reach) = "($(place(span[1], reach[1])),$(place(span[2], reach[2])))"

    # For every node, the directions of its edges, counterclockwise.
    directions = [Tuple{Float64,Int}[] for _ in positions]
    for (edge, (left, right, _)) in enumerate(graph.edges)
        angle = atan(positions[right][2] - positions[left][2],
                     positions[right][1] - positions[left][1])
        push!(directions[left], (mod2pi(angle), edge))
        push!(directions[right], (mod2pi(angle + pi), edge))
    end
    foreach(sort!, directions)
    # The direction of `edge` at `node` and the angles to its two neighbours.
    function around(node, edge)
        position = findfirst(direction -> direction[2] == edge, directions[node])
        count = length(directions[node])
        angle = directions[node][position][1]
        left = mod2pi(directions[node][mod1(position + 1, count)][1] - angle)
        right = mod2pi(angle - directions[node][mod1(position - 1, count)][1])
        return (angle, count == 1 ? 2pi : left, count == 1 ? 2pi : right)
    end

    for (node, (span, reach)) in enumerate(zip(drawing.spans, drawing.reaches))
        if isempty(graph.taxa[node])
            println(io, "\\coordinate ($prefix$node) at $(at(span, reach));")
        else
            println(io, "\\node[g2simple taxa] ($prefix$node) at $(at(span, reach)) " *
                        "{$(join(graph.taxa[node]))};")
        end
    end
    for (edge, (left, right, metric)) in enumerate(graph.edges)
        lines = g2_simple_numbered(drawing, metric) ? 1 : Int(abs(numerator(metric)))
        style = metric < 0 ? "g2simple negative" : "g2simple line"
        for line in 1:lines
            # Offsets are measured to the left of the edge at each end, so
            # they have opposite signs at the two ends.
            offset = (line - (lines + 1) / 2) * G2_SIMPLE_LINE_GAP_PT
            from, cut_from = g2_simple_line_end(
                "$prefix$left", !isempty(graph.taxa[left]), around(left, edge)..., offset,
            )
            to, cut_to = g2_simple_line_end(
                "$prefix$right", !isempty(graph.taxa[right]), around(right, edge)..., -offset,
            )
            (cut_from + cut_to) / G2_SIMPLE_PT_PER_CM < scale || throw(ArgumentError(
                "edge_length $scale is too short for the parallel lines of this tree",
            ))
            println(io, "\\draw[$style] $from--$to;")
        end
    end
    for label in g2_simple_length_labels(drawing)
        println(io, "\\node[g2simple length] at $(at(label.span, label.reach)) " *
                    "{\$$(numerator(label.metric))\$};")
    end
    return nothing
end

"""
    g2_simple_tree_tikz(tree::TreeMetric; edge_length=1.0, max_lines=nothing)

Return a TikZ fragment that draws `tree`, for a document loading TikZ. Every
edge is drawn with the length `edge_length`, in centimetres, between the
labels or unmarked nodes at its ends, whatever its metric length. An edge of
metric length `n` is `n` parallel lines, dashed when `n` is negative; an edge
of length zero is contracted. If `abs(n)` exceeds `max_lines`, the edge is
one line instead, with `n` written beside it. A node with taxa is the
sequence of their numbers, and a node without taxa is not marked. Edge
lengths must be integers.
"""
function g2_simple_tree_tikz(tree::TreeMetric; edge_length::Real=1.0, max_lines=nothing)
    drawing = g2_simple_tree_drawing(tree; edge_length=edge_length, max_lines=max_lines)
    return sprint() do io
        println(io, raw"\begingroup")
        println(io, g2_simple_tree_styles())
        println(io, raw"\begin{tikzpicture}[x=1cm,y=1cm,baseline=(current bounding box.center)]")
        g2_simple_tree_commands(io, drawing)
        println(io, raw"\end{tikzpicture}")
        println(io, raw"\endgroup")
    end
end

"""
    g2_simple_tree_tex(selected; edge_length=1.0, max_lines=nothing, standalone=false)

Return the input-ready fragment that draws a selected tree as
`g2_simple_tree_tikz` does. `selected` comes from `g2_lineality_tree`,
`g2_ray_tree`, `g2_edge_tree`, or `g2_maximal_cell_tree`. Pass
`standalone=true` for a complete document.
"""
function g2_simple_tree_tex(selected; edge_length::Real=1.0, max_lines=nothing,
                            standalone::Bool=false)
    return g2_wrap_tree_tikz(
        g2_simple_tree_tikz(selected.tree; edge_length=edge_length, max_lines=max_lines),
        selected.description * ", simple drawing";
        standalone=standalone,
    )
end

function g2_simple_tree_main(args=ARGS)
    script = "julia --project=. to_tex/render_simple_tree.jl"
    if "--help" in args
        println("Usage: $script $G2_TREE_SELECTION_USAGE " *
                "[--edge-length=LENGTH] [--max-lines=N] [--standalone] [--output=FILE]")
        println("       $script (--list-rays | --list-edges | --list-cells) " *
                "[--field=QQ|GF2]")
        return nothing
    end
    g2_tree_check_arguments(args, ("--standalone",), ("edge-length", "max-lines"))
    g2_tree_listing(args) && return nothing
    edge_length = parse(Float64, something(g2_tree_option(args, "edge-length"), "1.0"))
    max_lines = g2_tree_option(args, "max-lines")
    max_lines === nothing || (max_lines = parse(Int, max_lines))
    tex = g2_simple_tree_tex(g2_tree_selected(args); edge_length=edge_length,
                             max_lines=max_lines, standalone="--standalone" in args)
    return g2_tree_emit(args, tex)
end

if abspath(PROGRAM_FILE) == @__FILE__
    g2_simple_tree_main()
end
