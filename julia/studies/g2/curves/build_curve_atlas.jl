###############################################################################
# The atlas of a decorated curve graph in the tropical G2-Grassmannian. Its
# first page is a map of the graph: the vertices with their numbers, orange on
# walls of the G2-Grassmannian and green inside a cell, on walls of the
# Bergman fan, and every edge and leg labelled by the maximal cell of the
# coarsest fan holding it.
# The next pages hold the G2-tree of every cell on the map, drawn as in
# main.tex, with its census orbit, its number in the tables of maximal trees,
# and its edges and legs. A graph with at most CURVE_PAGE_CELLS cells has a
# one-page atlas instead. The graph of its cells, adjacent when the curve
# passes from one to the other, is laid out first, and the edges of each cell
# are packed at the place of their cell, in a coloured region that bands join
# where the edges of the cell are in several pieces. Among several such
# layouts, the map is one with the fewest crossings of edges in different
# cells. Its vertices are black on walls of the G2-Grassmannian and white
# inside a cell, and the tree of each cell stands beside the map, pointing to
# its region; --edge-length is then ignored.
#
#   julia --project=julia julia/studies/g2/curves/build_curve_atlas.jl \
#     GRAPH.csv [--drop-pending-trees] [--plain] [--unnumbered] [--output=ATLAS.pdf] [--edge-length=CM] [--no-pdf]
#
# GRAPH.csv is written by decorated_curve_graph.jl. Every edge and leg must
# lie in one maximal cell, as in its reduced graph, whose vertices of degree 2
# are removed only between two edges in one cell. With --drop-pending-trees
# the map is that of the graph without its pending trees and their legs, see
# drop_pending_trees; the core is then drawn with the fewest crossings found,
# from a planar embedding of it or of it without one edge (curve_drawn_core).
# With --unnumbered the vertices are dots without their numbers, and the
# catalog counts the edges and legs of every cell.
# With --plain the atlas is a picture of half a page of main.tex, the map
# alone as on the first page but without vertex numbers and cells, the
# vertices dots of their kind, the core of the graph drawn first and in bold
# and the pending trees small around it; with --drop-pending-trees as well,
# the core alone, thick, its edges labelled by their cells. A segment of the
# core crossing another passes over it with a thin white border. The TeX source is written beside the
# PDF, by default GRAPH_atlas.pdf, or GRAPH_core_atlas.pdf with
# --drop-pending-trees, GRAPH_plain.pdf with --plain, beside GRAPH.csv, and
# the compilation files go to julia/work/tex/curve_atlas/.
###############################################################################

include(joinpath(@__DIR__, "..", "..", "..", "to_tex", "render_maximal_cell_tree.jl"))

import LinearAlgebra
import Random
using TropicalCensus: work_path

const CURVE_GRAPH_HEADER = "edge,tail,head,cells,orbits,tree_numbers,pieces"
# The drawn length of an edge of the map, in centimetres. A smaller map is
# enlarged until it fills an A3 page, its nodes, labels and strokes with it.
const CURVE_EDGE_LENGTH_CM = 1.5
const CURVE_LAYOUT_SWEEPS = 300
const CURVE_SPREAD_SWEEPS = 600
# The least distance of a point from an edge without it, in edge lengths.
const CURVE_CLEARANCE = 0.4
# The width of an entry of the catalog of trees, in centimetres.
const CURVE_ENTRY_WIDTH_CM = 6.2
# An atlas of at most this many cells is one A3 page, its cells coloured
# regions of the map and their trees around it, instead of a catalog.
const CURVE_PAGE_CELLS = 12
# The picture of that page below its header, in centimetres.
const CURVE_PAGE_PICTURE = (40.0, 26.0)
# The colours of its cells, a categorical palette checked for colour
# blindness, in this order; a cell after the eighth takes a colour again, with
# dashed edges.
const CURVE_CELL_COLOURS = ["2A78D6", "EB6834", "1BAF7A", "EDA100", "E87BA4", "008300", "6250D6", "E34948"]
# Two vertices of edges in one cell are laid out at most this far apart.
const CURVE_CELL_DISTANCE = 2.0
# The label of the cell of an edge stands at the first of these places along
# it that is this far from the labels before it, in edge lengths.
const CURVE_LABEL_POSITIONS = (0.5, 0.35, 0.65, 0.25, 0.75)
const CURVE_LABEL_SPACING = 0.3
# Two vertices at most this far apart, in edge lengths, may swap places, and
# a vertex may move on a grid of this step up to 8 steps away.
const CURVE_SWAP_REACH = 3.0
const CURVE_MOVE_STEP = 0.5
# On a one-page atlas, a cell of k edges has radius CURVE_CELL_RADIUS * sqrt(k)
# in the graph of the cells, and every point is drawn with this weight to the
# mean of the places of its cells there.
const CURVE_CELL_RADIUS = 0.9
const CURVE_ANCHOR_WEIGHT = 0.5
# The number of layouts of the graph of the cells tried for a one-page atlas.
const CURVE_LAYOUT_TRIES = 24
# A pending tree drawn from its point on the core first drawn: at most this
# much per edge, the core of largest side 1, and within this angle of its
# direction.
const CURVE_TREE_STEP = 0.036
const CURVE_TREE_DIRECTIONS = 48
# A tree out of the outer cycle of the core counts this many times its radius.
const CURVE_TREE_OUTSIDE = 3
const CURVE_TREE_ANGLE = deg2rad(50)
# The relaxation of the core first drawn: its sweeps, and the least distance
# from a point to a segment of the core more than CURVE_CORE_HOPS edges away
# that it keeps, a chain of one edge of length 1.
const CURVE_CORE_SWEEPS = 300
const CURVE_CORE_CLEARANCE = 0.35
const CURVE_CORE_HOPS = 3
# The plain map fills half a page of main.tex, \\textwidth 6.5in by half of
# \\textheight 8.8in, in centimetres.
const CURVE_PLAIN_SIZE = (16.51, 11.17)
# The core of the plain map, and the white border of a segment of it over another.
const CURVE_CORE_PT = 1.3
const CURVE_BORDER_PT = 0.9
const CURVE_CORE_ONLY_PT = 1.8

"""
The edges and legs of the decorated graph `path`: number, tail, head (0 for a
leg), and the cell, census orbit and tree number of the maximal cell that
holds all its pieces.
"""
function read_curve_graph(path)
    lines = readlines(path)
    first(lines) == CURVE_GRAPH_HEADER || error("$path is not a decorated graph of decorated_curve_graph.jl")
    return map(lines[2:end]) do line
        fields = split(line, ',')
        cells, orbits, numbers = (unique(parse.(Int, split(field, ';'))) for field in fields[4:6])
        length(cells) == 1 ||
            error("edge $(fields[1]) of $path runs through the cells $cells; draw the reduced graph")
        (edge=parse(Int, fields[1]), tail=parse(Int, fields[2]), head=parse(Int, fields[3]),
         cell=only(cells), orbit=only(orbits), number=only(numbers))
    end
end

"The edges of the core of `edges`: the legs are dropped and the leaves removed until none is left."
function curve_core(edges)
    kept = filter(e -> e.head != 0, edges)
    while true
        valence = Dict{Int,Int}()
        for e in kept, v in (e.tail, e.head)
            valence[v] = get(valence, v, 0) + 1
        end
        any(==(1), values(valence)) || return kept
        filter!(e -> valence[e.tail] > 1 && valence[e.head] > 1, kept)
    end
end

"""
The edges of `edges` without the pending trees, the trees hanging from the
rest of the graph, see `curve_core`. A vertex of degree 2 is then replaced by
one edge between its neighbours when its two edges lie in one cell. The
vertices keep their numbers.
"""
function drop_pending_trees(edges)
    kept = curve_core(edges)
    merged = true
    while merged
        merged = false
        incident = Dict{Int,Vector{Int}}()
        for (k, e) in enumerate(kept), v in (e.tail, e.head)
            push!(get!(incident, v, Int[]), k)
        end
        for (v, around) in incident
            length(around) == 2 || continue
            i, j = around
            i != j && kept[i].cell == kept[j].cell || continue
            other(e) = e.tail == v ? e.head : e.tail
            kept[i] = merge(kept[i], (tail=other(kept[i]), head=other(kept[j])))
            deleteat!(kept, j)
            merged = true
            break
        end
    end
    return kept
end

"""
The two ends of every edge and leg: the vertices by their places in `index`,
and the far ends of the legs numbered after them, in the order of `edges`.
"""
function curve_ends(edges, index)
    far = cumsum([e.head == 0 for e in edges])
    return [(index[e.tail], e.head == 0 ? length(index) + far[k] : index[e.head]) for (k, e) in enumerate(edges)]
end

"The point of the segment from `a` to `b` closest to `p`."
function curve_closest(p, a, b)
    d = b .- a
    return a .+ clamp(sum((p .- a) .* d) / max(sum(d .^ 2), 1e-12), 0, 1) .* d
end

"""
The segments joining the pieces of a region, the connected pieces of its
edges `pairs` between the points at `x`, `y`: each joins the closest middles
of edges of two pieces, never a vertex, along a shortest tree on the pieces.
"""
function curve_region_joins(pairs, x, y)
    parent = Dict(p => p for pair in pairs for p in pair)
    root(p) = parent[p] == p ? p : (parent[p] = root(parent[p]))
    foreach(((u, v),) -> parent[root(u)] = root(v), pairs)
    pieces = collect(values(Dict(r => [k for (k, (u, _)) in enumerate(pairs) if root(u) == r]
                                  for r in unique(root(u) for (u, _) in pairs))))
    middle = curve_middles(pairs, x, y)
    length_of((p, q)) = hypot((p .- q)...)
    gap(A, B) = argmin(length_of, [(middle[a], middle[b]) for a in A for b in B])
    tree, joins = [popfirst!(pieces)], Tuple{NTuple{2,Float64},NTuple{2,Float64}}[]
    while !isempty(pieces)
        nearest = argmin(((_, g),) -> length_of(g), [(j, gap(A, pieces[j])) for A in tree for j in eachindex(pieces)])
        push!(joins, nearest[2])
        push!(tree, popat!(pieces, nearest[1]))
    end
    return joins
end

"""
Spread the points `x`, `y`: every edge of `ends` is a spring of length 1, two
points closer than 1 repel each other, and so do a point and an edge without
it closer than `CURVE_CLEARANCE`. The steps shrink at every sweep.
"""
function curve_spread!(x, y, ends)
    N = length(x)
    for sweep in 1:CURVE_SPREAD_SWEEPS
        fx, fy = zeros(N), zeros(N)
        # `i` moves away from `others` by `reach - r` along the offset `d`
        # of length `r`, and each of `others` back by `share` of it.
        function repel!(i, others, d, reach; share=1.0)
            r = max(hypot(d...), 1e-6)
            f = (reach - r) / r
            fx[i] += f * d[1]
            fy[i] += f * d[2]
            foreach(k -> (fx[k] -= share * f * d[1]; fy[k] -= share * f * d[2]), others)
        end
        for (u, v) in ends
            u == v || repel!(u, (v,), (x[u] - x[v], y[u] - y[v]), 1.0)
        end
        for i in 1:N, j in i+1:N
            hypot(x[i] - x[j], y[i] - y[j]) < 1 && repel!(i, (j,), (x[i] - x[j], y[i] - y[j]), 1.0)
        end
        for p in 1:N, (a, b) in ends
            (p == a || p == b || a == b) && continue
            d = (x[p], y[p]) .- curve_closest((x[p], y[p]), (x[a], y[a]), (x[b], y[b]))
            hypot(d...) < CURVE_CLEARANCE && repel!(p, (a, b), d, CURVE_CLEARANCE; share=0.5)
        end
        step = 0.15 * (1 - sweep / CURVE_SPREAD_SWEEPS) + 0.01
        for i in 1:N
            m = hypot(fx[i], fy[i])
            m > 0 && (x[i] += min(m, step) / m * fx[i]; y[i] += min(m, step) / m * fy[i])
        end
    end
    return x, y
end

"""
Whether the edges `a` and `b` of `ends` between the points at `x`, `y` cross,
without sharing an end and, given the `colours` of the edges, of different
colours.
"""
function curve_crossing(x, y, ends, a, b; colours=nothing)
    isnothing(colours) || colours[a] != colours[b] || return false
    (p, q), (r, s) = ends[a], ends[b]
    allunique((p, q, r, s)) || return false
    turn(i, j, k) = sign((x[j] - x[i]) * (y[k] - y[i]) - (y[j] - y[i]) * (x[k] - x[i]))
    return turn(p, q, r) * turn(p, q, s) < 0 && turn(r, s, p) * turn(r, s, q) < 0
end

"The number of crossings of the edges `ends` between the points at `x`, `y`, as `curve_crossing` counts them."
curve_crossings(x, y, ends; colours=nothing) =
    count(curve_crossing(x, y, ends, a, b; colours) for a in eachindex(ends) for b in a+1:length(ends))

"""
Swap the places of two vertices at most `CURVE_SWAP_REACH` apart, or move a
vertex to a point around the middle of its neighbours or of the grid of step
`CURVE_MOVE_STEP` around it, at least half an edge from the other points, the
far ends of their legs, the points of degree 1, moving with them, whenever
that lowers the number of crossings of the edges `ends`, of different
`colours` if given, without bringing more
points closer than `CURVE_CLEARANCE` to edges without them at the moved
vertices, until nothing does. A move goes to the nearest of the points with
the fewest crossings.
"""
function curve_uncross!(x, y, ends; colours=nothing)
    N = length(x)
    incident = [Int[] for _ in 1:N]
    for (k, (u, v)) in enumerate(ends)
        push!(incident[u], k)
        u == v || push!(incident[v], k)
    end
    far = [length(incident[p]) == 1 for p in 1:N]
    legs = [[q for k in incident[p] for q in ends[k] if q != p && far[q]] for p in 1:N]
    point(p) = (x[p], y[p])
    crossings(S) = count(curve_crossing(x, y, ends, a, b; colours) for a in S for b in eachindex(ends) if !(b in S) || b > a)
    function shift!(p, d)
        for q in [p; legs[p]]
            x[q] += d[1]
            y[q] += d[2]
        end
    end
    function swap!(u, v)
        d = point(v) .- point(u)
        shift!(u, d)
        shift!(v, .-d)
    end
    # The pairs of a point and an edge without it closer than `CURVE_CLEARANCE`
    # at `v`: `v` and the other edges, and the edges of `v` and the other points.
    gap(q, (a, b)) = hypot((point(q) .- curve_closest(point(q), point(a), point(b)))...)
    tight(v) = sum(enumerate(ends)) do (k, (a, b))
        k in incident[v] ? count(q -> q != a && q != b && gap(q, (a, b)) < CURVE_CLEARANCE, 1:N) :
                           v != a && v != b && gap(v, (a, b)) < CURVE_CLEARANCE
    end
    improved = true
    while improved
        improved = false
        for u in 1:N, v in u+1:N
            (far[u] || far[v] || hypot((point(u) .- point(v))...) > CURVE_SWAP_REACH) && continue
            S = union(incident[u], incident[v])
            before, tight_before = crossings(S), tight(u) + tight(v)
            swap!(u, v)
            crossings(S) < before && tight(u) + tight(v) <= tight_before ? (improved = true) : swap!(u, v)
        end
        for v in 1:N
            far[v] && continue
            around = [q for k in incident[v] for q in ends[k] if q != v && !far[q]]
            isempty(around) && continue
            middle = reduce((p, q) -> p .+ q, point.(around)) ./ length(around)
            targets = [[middle .+ r .* (cos(a), sin(a)) for r in (0.0, 0.6, 1.2) for a in 0:π/4:7π/4];
                       [point(v) .+ CURVE_MOVE_STEP .* (i, j) for i in -8:8 for j in -8:8 if 0 < hypot(i, j) <= 8]]
            before, tight_before = crossings(incident[v]), tight(v)
            best = (before, 0.0, nothing)
            for target in targets
                all(q -> q == v || q in legs[v] || hypot((point(q) .- target)...) >= 0.5, 1:N) || continue
                d = target .- point(v)
                shift!(v, d)
                c, m = crossings(incident[v]), hypot(d...)
                (c, m) < best[1:2] && tight(v) <= tight_before && (best = (c, m, target))
                shift!(v, .-d)
            end
            best[1] < before && (shift!(v, best[3] .- point(v)); improved = true)
        end
    end
    return x, y
end

"""
Points in the plane whose distances approach `distance`, by stress
majorization with weights `1/d^2`, from classical scaling or from the places
`start`. A point is also drawn to its place in `anchors`, if given, with
weight `CURVE_ANCHOR_WEIGHT`.
"""
function curve_stress(distance; start=nothing, anchors=nothing)
    N = size(distance, 1)
    N == 1 && return [0.0], [0.0]
    if isnothing(start)
        centring = LinearAlgebra.I(N) - fill(1 / N, N, N)
        scaling = LinearAlgebra.eigen(LinearAlgebra.Symmetric(-centring * distance .^ 2 * centring / 2), N-1:N)
        x, y = (scaling.vectors[:, k] .* sqrt(max(scaling.values[k], 0)) for k in (2, 1))
    else
        x, y = first.(start), last.(start)
    end
    for _ in 1:CURVE_LAYOUT_SWEEPS, i in 1:N
        sx = sy = total = 0.0
        for j in 1:N
            j == i && continue
            d = distance[i, j]
            r = max(hypot(x[i] - x[j], y[i] - y[j]), 1e-9)
            sx += (x[j] + d * (x[i] - x[j]) / r) / d^2
            sy += (y[j] + d * (y[i] - y[j]) / r) / d^2
            total += 1 / d^2
        end
        if !isnothing(anchors)
            sx += CURVE_ANCHOR_WEIGHT * anchors[i][1]
            sy += CURVE_ANCHOR_WEIGHT * anchors[i][2]
            total += CURVE_ANCHOR_WEIGHT
        end
        x[i], y[i] = sx / total, sy / total
    end
    return x, y
end

"""
The places of the cells `cells` of `edges` on the map, from their graph: two
cells are adjacent when an edge of each ends at one point, at the sum of
their radii, `CURVE_CELL_RADIUS` times the square root of their numbers of
edges, and laid out by `curve_stress`, from classical scaling for `seed` 0
and from random places drawn with `seed` otherwise. Returns, for each of the
`N` points of `ends`, the mean of the places of the cells of its edges.
"""
function curve_cell_anchors(edges, ends, cells, N; seed=0)
    index = Dict(c => i for (i, c) in enumerate(cells))
    m = length(cells)
    around = [Set{Int}() for _ in 1:N]
    for (e, (u, v)) in zip(edges, ends)
        push!(around[u], index[e.cell])
        push!(around[v], index[e.cell])
    end
    radius = [CURVE_CELL_RADIUS * sqrt(count(e -> e.cell == c, edges)) for c in cells]
    distance = [i == j ? 0.0 : Inf for i in 1:m, j in 1:m]
    for s in around, i in s, j in s
        i == j || (distance[i, j] = min(distance[i, j], radius[i] + radius[j]))
    end
    for k in 1:m, i in 1:m, j in 1:m
        distance[i, j] = min(distance[i, j], distance[i, k] + distance[k, j])
    end
    distance[isinf.(distance)] .= maximum(filter(isfinite, distance)) + 2maximum(radius)
    rng = Random.Xoshiro(seed)
    start = seed == 0 ? nothing : [Tuple(maximum(distance) .* randn(rng, 2)) for _ in 1:m]
    X, Y = curve_stress(distance; start)
    return [(sum(X[i] for i in s) / length(s), sum(Y[i] for i in s) / length(s)) for s in around]
end

"""
The layout of a one-page atlas of the cells `cells` of `edges`: of the
layouts of `curve_layout` drawn to the places of `curve_cell_anchors` for the
seeds `0:CURVE_LAYOUT_TRIES-1`, one with the fewest crossings of edges in
different cells; edges of one cell may cross.
"""
function curve_cell_layout(edges, ends, cells, N)
    groups = [unique(p for (e, pair) in zip(edges, ends) if e.cell == c for p in pair) for c in cells]
    colours = [e.cell for e in edges]
    layouts = map(0:CURVE_LAYOUT_TRIES-1) do seed
        x, y = curve_layout(ends, N; groups, anchors=curve_cell_anchors(edges, ends, cells, N; seed), colours)
        (curve_crossings(x, y, ends; colours), x, y)
    end
    _, x, y = argmin(first, layouts)
    return x, y
end

"""
Positions of the `N` points joined by the pairs `ends`, by `curve_stress`:
the distances in the plane approach those in the graph, every edge of
length 1. Points in different components are one more than the diameter
apart, and two points of one of `groups` at most `CURVE_CELL_DISTANCE`.
Given `anchors`, the points start there, a little apart, and are drawn to
them. They are then spread by `curve_spread!`, swapped by `curve_uncross!`,
which counts only the crossings of edges of different `colours` if given,
spread again and swapped again, and the longer side is horizontal.
"""
function curve_layout(ends, N; groups=(), anchors=nothing, colours=nothing)
    neighbours = [Int[] for _ in 1:N]
    for (u, v) in ends
        u == v && continue
        push!(neighbours[u], v)
        push!(neighbours[v], u)
    end
    distance = fill(-1.0, N, N)
    for s in 1:N
        distance[s, s] = 0
        queue = [s]
        for u in queue, v in neighbours[u]
            distance[s, v] < 0 && (distance[s, v] = distance[s, u] + 1; push!(queue, v))
        end
    end
    distance[distance .< 0] .= maximum(distance) + 1
    for group in groups, i in group, j in group
        i == j || (distance[i, j] = min(distance[i, j], CURVE_CELL_DISTANCE))
    end
    # The points of one anchor start around it at golden angles, a little apart.
    start = isnothing(anchors) ? nothing : [a .+ 0.3 .* (cos(2.39996p), sin(2.39996p)) for (p, a) in enumerate(anchors)]
    x, y = curve_stress(distance; start, anchors)
    curve_spread!(x, y, ends)
    curve_uncross!(x, y, ends; colours)
    curve_spread!(x, y, ends)
    curve_uncross!(x, y, ends; colours)
    return maximum(x) - minimum(x) >= maximum(y) - minimum(y) ? (x, y) : (y, x)
end

"""
The chains of the core, the edges `bold` of `ends`: the paths between its
points of degree at least 3 through points of degree 2. Returns those
points and the chains, each as its points, from one of them to another, and
its edges.
"""
function curve_core_chains(ends, bold)
    incident = Dict{Int,Vector{Int}}()
    foreach(k -> foreach(p -> push!(get!(incident, p, Int[]), k), ends[k]), bold)
    branch = sort!([p for (p, around) in incident if length(around) >= 3])
    chains = Tuple{Vector{Int},Vector{Int}}[]
    seen = Set{Int}()
    for b in branch, k in incident[b]
        k in seen && continue
        points, path, p, e = [b], Int[], b, k
        while true
            push!(seen, e)
            push!(path, e)
            u, v = ends[e]
            p = u == p ? v : u
            push!(points, p)
            length(incident[p]) == 2 || break
            e = only(filter(!=(e), incident[p]))
        end
        push!(chains, (points, path))
    end
    return branch, chains
end

"""
The planar embeddings of the connected graph of `m` vertices and the edges
`pairs`, of degrees at most 3, by brute force over its rotation systems up to
mirror images: those with `E - V + 2` faces. Each is a list of faces, lists
of darts, the dart `2c - 1` running along `pairs[c]` and `2c` back. A vertex
of `crossings` has degree 4, the crossing of two edges: its darts, in their
order, alternate between them. With `stop`, only the first embedding found.
"""
function curve_embeddings(m, pairs; crossings=(), stop=false)
    E = length(pairs)
    tail(d) = isodd(d) ? pairs[(d+1)÷2][1] : pairs[d÷2][2]
    twin(d) = isodd(d) ? d + 1 : d - 1
    out = [Int[] for _ in 1:m]
    foreach(d -> push!(out[tail(d)], d), 1:2E)
    all(v -> length(out[v]) <= (v in crossings ? 4 : 3), 1:m) || error("the core has a vertex of degree above 3")
    free = findall(o -> length(o) >= 3, out)
    isempty(free) || pop!(free)
    target = E - count(!isempty, out) + 2
    embeddings = Vector{Vector{Int}}[]
    following = zeros(Int, 2E)
    for flip in 0:(1<<length(free))-1
        for v in 1:m
            k = findfirst(==(v), free)
            o = isnothing(k) || iseven(flip >> (k - 1)) ? out[v] : reverse(out[v])
            foreach(i -> following[o[i]] = o[mod1(i + 1, length(o))], eachindex(o))
        end
        seen, faces = falses(2E), Vector{Int}[]
        for d in 1:2E
            seen[d] && continue
            face, e = Int[], d
            while !seen[e]
                seen[e] = true
                push!(face, e)
                e = following[twin(e)]
            end
            push!(faces, face)
        end
        length(faces) == target && (push!(embeddings, faces); stop && break)
    end
    return embeddings
end

"""
A shortest way from the vertex `s` to the vertex `t` through the `faces` of
an embedding of the edges `pairs`, see `curve_embeddings`: the faces it
passes and the darts it crosses, one between two consecutive faces.
"""
function curve_route(pairs, faces, s, t)
    tail(d) = isodd(d) ? pairs[(d+1)÷2][1] : pairs[d÷2][2]
    twin(d) = isodd(d) ? d + 1 : d - 1
    face_of = zeros(Int, 2length(pairs))
    foreach(((f, face),) -> face_of[face] .= f, enumerate(faces))
    # The faces by their distance from those at `s`, each with the face and the dart it is reached from.
    distance, back = Dict{Int,Int}(), Dict{Int,Tuple{Int,Int}}()
    queue = [f for f in eachindex(faces) if any(d -> tail(d) == s, faces[f])]
    foreach(f -> (distance[f] = 0; back[f] = (0, 0)), queue)
    for f in queue, d in faces[f]
        g = face_of[twin(d)]
        haskey(distance, g) || (distance[g] = distance[f] + 1; back[g] = (f, d); push!(queue, g))
    end
    goal = argmin(f -> distance[f], [f for f in eachindex(faces) if any(d -> tail(d) == t, faces[f])])
    path, crossed, f = [goal], Int[], goal
    while back[f][1] != 0
        f, d = back[f]
        pushfirst!(path, f)
        pushfirst!(crossed, d)
    end
    return path, crossed
end

"""
Tutte's barycentric drawing of the graph of `m` vertices and the edges
`pairs` with the weights `weights`: the cycle `order` on the unit circle,
its edges taking arcs proportional to `lengths`, every other vertex at the
weighted mean of its neighbours. Crossing-free, with convex faces, when the
graph is planar and 3-connected and the cycle bounds a face.
"""
function curve_tutte(m, pairs, weights, order, lengths)
    L = zeros(m, m)
    for ((u, v), w) in zip(pairs, weights)
        L[u, u] += w
        L[v, v] += w
        L[u, v] -= w
        L[v, u] -= w
    end
    angles = 2π .* [0; cumsum(lengths)[1:end-1]] ./ sum(lengths)
    X = zeros(m, 2)
    X[order, :] = [cos.(angles) sin.(angles)]
    inner = setdiff(1:m, order)
    isempty(inner) || (X[inner, :] = -L[inner, inner] \ (L[inner, order] * X[order, :]))
    return X
end

"""
`X`, the places of the branch points of the chains `pairs` of the core,
relaxed by stress majorization towards their distances in the graph of
chains, a chain of `k` edges of length `sqrt(k)`, after the best scaling:
sweep after sweep, a branch point moves to its place of least stress, or
a half, a quarter, a tenth or a twentieth of the way, only when the core `placed(X)` keeps its
number of crossings and its clearance, the least distance from a point of
the core to a segment more than `CURVE_CORE_HOPS` edges away along it,
stays at least the smaller of its value and `CURVE_CORE_CLEARANCE`.
"""
function curve_relax_core(X, pairs, sizes, placed, crossings_of, clearance_of)
    m = size(X, 1)
    d = fill(Inf, m, m)
    foreach(i -> d[i, i] = 0, 1:m)
    for ((u, v), k) in zip(pairs, sizes)
        d[u, v] = d[v, u] = min(d[u, v], sqrt(k))
    end
    for k in 1:m, i in 1:m, j in 1:m
        d[i, j] = min(d[i, j], d[i, k] + d[k, j])
    end
    r(X, i, j) = max(hypot((X[i, :] .- X[j, :])...), 1e-9)
    pairs_off = [(i, j) for i in 1:m for j in i+1:m]
    X = X .* (sum(r(X, i, j) / d[i, j] for (i, j) in pairs_off) / sum(r(X, i, j)^2 / d[i, j]^2 for (i, j) in pairs_off))
    crossings, clearance = crossings_of(placed(X)), clearance_of(placed(X))
    for _ in 1:CURVE_CORE_SWEEPS, i in 1:m
        weights = [j == i ? 0.0 : 1 / d[i, j]^2 for j in 1:m]
        target = sum(weights[j] .* (X[j, :] .+ d[i, j] .* (X[i, :] .- X[j, :]) ./ r(X, i, j)) for j in 1:m if j != i) ./ sum(weights)
        for share in (1.0, 0.5, 0.25, 0.1, 0.05)
            Y = copy(X)
            Y[i, :] = X[i, :] .+ share .* (target .- X[i, :])
            place = placed(Y)
            crossings_of(place) == crossings || continue
            c = clearance_of(place)
            c >= min(clearance, CURVE_CORE_CLEARANCE) || continue
            X, clearance = Y, c
            break
        end
    end
    return X
end

"""
A drawing of the core, the edges `bold` of `ends`, with the fewest crossings
this finds, or `nothing`. The graph of its chains, see `curve_core_chains`,
is drawn by `curve_tutte` from a planar embedding, see `curve_embeddings`,
with every face as the outer cycle in turn, a chain of `k` edges of weight
`1/k` and arc `k`. A nonplanar graph is drawn without the chain whose
`curve_route` back crosses the fewest chains, then drawn along that route,
square across the middles of the chains it crosses. For each chain drawn
last, the drawing whose outer cycle has the most `weight` on its points,
then whose clearance, see `curve_relax_core`, over its largest side, is
largest, is relaxed by `curve_relax_core`, and the
relaxed drawing best in the same order is kept. The points inside a chain
divide it evenly. Returns the places of the points of the core, the drawing
of largest side 1, the number of crossings, and the points of the outer
cycle in their order.
"""
function curve_drawn_core(ends, bold; weight=Dict{Int,Int}())
    branch, chains = curve_core_chains(ends, bold)
    m = length(branch)
    index = Dict(p => i for (i, p) in enumerate(branch))
    pairs = [(index[first(points)], index[last(points)]) for (points, _) in chains]
    # Tutte's drawing needs a simple graph.
    (m >= 4 && allunique(minmax(p...) for p in pairs) && all(((u, v),) -> u != v, pairs)) || return nothing
    sizes = [length(path) for (_, path) in chains]
    # (removed chain or 0, the other chains, an embedding of them, the faces and crossed darts back)
    options = [(0, collect(eachindex(pairs)), faces, ([0], Int[])) for faces in curve_embeddings(m, pairs)]
    if isempty(options)
        for k in eachindex(pairs)
            rest = setdiff(eachindex(pairs), k)
            append!(options, [(k, rest, faces, curve_route(pairs[rest], faces, pairs[k]...))
                              for faces in curve_embeddings(m, pairs[rest])])
        end
        isempty(options) && return nothing
        fewest = minimum(length(o[4][2]) for o in options)
        filter!(o -> length(o[4][2]) == fewest, options)
    end
    core_ends = [ends[k] for k in bold]
    function placed(X, removed, rest, crossed)
        place = Dict{Int,Tuple{Float64,Float64}}()
        along(c, t) = Tuple(X[pairs[c][1], :] .+ t .* (X[pairs[c][2], :] .- X[pairs[c][1], :]))
        for c in rest, (i, p) in enumerate(chains[c][1])
            place[p] = along(c, (i - 1) / sizes[c])
        end
        removed == 0 && return place
        # Square across the middle of the middle segment of every chain crossed, between corners
        # half that segment away on either side.
        corners, ends_at = [Tuple(X[pairs[removed][1], :])], Tuple(X[pairs[removed][2], :])
        for (i, d) in enumerate(crossed)
            c = rest[(d+1)÷2]
            middle = along(c, (fld(sizes[c], 2) + 0.5) / sizes[c])
            e = X[pairs[c][2], :] .- X[pairs[c][1], :]
            normal = (-e[2], e[1]) ./ hypot(e...)
            following = i == length(crossed) ? ends_at : along(rest[(crossed[i+1]+1)÷2], 0.5)
            reach = 0.5hypot(e...) / sizes[c]
            sum(normal .* (following .- middle)) < 0 && (normal = .-normal)
            append!(corners, [middle .- reach .* normal, middle, middle .+ reach .* normal])
        end
        push!(corners, ends_at)
        run = [0; cumsum([hypot((corners[i+1] .- corners[i])...) for i in 1:length(corners)-1])]
        points = chains[removed][1]
        share = [run[end] * (i - 1) / sizes[removed] for i in eachindex(points)]
        # The inner points nearest to the corners beside a crossing move onto them, when there are enough.
        if length(points) - 2 >= 2length(crossed)
            for j in 2:length(corners)-1
                mod(j - 2, 3) == 1 && continue
                i = argmin(i -> abs(share[i] - run[j]), 2:length(points)-1)
                share[i] = run[j]
            end
        end
        for (i, p) in enumerate(points)
            j = clamp(searchsortedlast(run, share[i]), 1, length(corners) - 1)
            place[p] = corners[j] .+ ((share[i] - run[j]) / max(run[j+1] - run[j], 1e-12)) .* (corners[j+1] .- corners[j])
        end
        return place
    end
    crossings_of(place) = curve_crossings(Dict(p => q[1] for (p, q) in place), Dict(p => q[2] for (p, q) in place), core_ends)
    # The least distance from a point of the core to a segment more than `CURVE_CORE_HOPS` edges away along it.
    hops = Dict{Int,Dict{Int,Int}}()
    neighbours = Dict{Int,Vector{Int}}()
    foreach(((a, b),) -> (push!(get!(neighbours, a, Int[]), b); push!(get!(neighbours, b, Int[]), a)), core_ends)
    for p in keys(neighbours)
        hops[p] = Dict(p => 0)
        queue = [p]
        for u in queue, v in neighbours[u]
            haskey(hops[p], v) || (hops[p][v] = hops[p][u] + 1; push!(queue, v))
        end
    end
    far = [(p, a, b) for p in keys(neighbours) for (a, b) in core_ends if min(hops[p][a], hops[p][b]) > CURVE_CORE_HOPS]
    clearance_of(place) = minimum(hypot((place[p] .- curve_closest(place[p], place[a], place[b]))...) for (p, a, b) in far)
    extent(place) = maximum(maximum(q[k] for q in values(place)) - minimum(q[k] for q in values(place)) for k in 1:2)
    # For every chain drawn last, or none, the best Tutte drawing, then relaxed.
    best = Dict{Int,Any}()
    for (removed, rest, faces, (route, crossed)) in options
        tail(d) = isodd(d) ? pairs[rest[(d+1)÷2]][1] : pairs[rest[d÷2]][2]
        for (f, face) in enumerate(faces)
            order = tail.(face)
            # The outer face is a cycle off the route back.
            allunique(order) && !(f in route) || continue
            X = curve_tutte(m, pairs[rest], 1 ./ sizes[rest], order, [sizes[rest[(d+1)÷2]] for d in face])
            place = placed(X, removed, rest, crossed)
            crossings_of(place) == length(crossed) || continue
            outer = reduce(vcat, (isodd(d) ? chains[rest[(d+1)÷2]][1] : reverse(chains[rest[d÷2]][1]))[1:end-1] for d in face)
            score = (sum(p -> get(weight, p, 0), outer), clearance_of(place) / extent(place))
            (!haskey(best, removed) || score > best[removed][1]) && (best[removed] = (score, X, rest, crossed, outer))
        end
    end
    isempty(best) && return nothing
    relaxed = map(collect(best)) do (removed, ((outside, _), X, rest, crossed, outer))
        at(X) = placed(X, removed, rest, crossed)
        place = at(curve_relax_core(X, pairs, sizes, at, crossings_of, clearance_of))
        ((outside, clearance_of(place) / extent(place)), place, length(crossed), outer)
    end
    _, place, crossings, outer = argmax(first, relaxed)
    return Dict(p => q ./ extent(place) for (p, q) in place), crossings, outer
end

"""
The distance from `p` along the unit vector `d` to the first of the
`segments` without the point `skip`, or `Inf`.
"""
function curve_ray(p, d, segments, skip)
    t = Inf
    for (a, b, pa, pb) in segments
        (a == skip || b == skip) && continue
        e = pb .- pa
        denominator = d[1] * e[2] - d[2] * e[1]
        abs(denominator) < 1e-12 && continue
        w = pa .- p
        s, r = (w[1] * e[2] - w[2] * e[1]) / denominator, (w[1] * d[2] - w[2] * d[1]) / denominator
        s > 1e-9 && -1e-9 <= r <= 1 + 1e-9 && (t = min(t, s))
    end
    return t
end

"""
Places for the `N` points of `ends`: the core, the edges `bold`, drawn by
`curve_drawn_core`, and every pending tree drawn radially from the point of
the core where it hangs, in sectors of angle in proportion to the numbers of
leaves. The outer cycle of the core is the one that holds the most leaves
of trees. The trees are placed from the largest on: each in the circular
sector about its point of the largest radius, up to `CURVE_TREE_STEP` of the
core per edge, that meets neither the core nor the sectors of the trees
before it, among `CURVE_TREE_DIRECTIONS` directions and its angle, within
`CURVE_TREE_ANGLE`, or two thirds or a third of it; a sector out of the outer cycle
counts `CURVE_TREE_OUTSIDE` times its radius. Returns `nothing` when
`curve_drawn_core` does.
"""
function curve_core_first_layout(ends, N, bold)
    around = [Int[] for _ in 1:N]
    for (k, (u, v)) in enumerate(ends)
        k in bold && continue
        push!(around[u], v)
        push!(around[v], u)
    end
    leaves(v, parent) = (children = filter(!=(parent), around[v]); isempty(children) ? 1 : sum(c -> leaves(c, v), children))
    # The leaves of the trees at every point of the core, for an outer cycle that holds many.
    weight = Dict(p => sum(q -> leaves(q, p), around[p]; init=0) for k in bold for p in ends[k])
    drawn = curve_drawn_core(ends, bold; weight)
    isnothing(drawn) && return nothing
    place, crossings, outer = drawn
    println("The core is drawn with $crossings crossing", crossings == 1 ? "" : "s", "; ",
            sum(p -> weight[p], outer), " of the ", sum(values(weight)), " leaves of the trees hang from its outer cycle")
    x, y = zeros(N), zeros(N)
    foreach(((p, q),) -> (x[p], y[p]) = q, place)
    core = [(place[a], place[b]) for (a, b) in (ends[k] for k in bold)]
    depth(v, parent) = 1 + maximum((depth(c, v) for c in around[v] if c != parent); init=0)
    function hang!(v, parent, root, step, low, high, level)
        middle = (low + high) / 2
        x[v], y[v] = place[root] .+ (level * step) .* (cos(middle), sin(middle))
        children = filter(!=(parent), around[v])
        total, start = sum(c -> leaves(c, v), children; init=0), low
        for c in children
            width = (high - low) * leaves(c, v) / total
            hang!(c, v, root, step, start, start + width, level + 1)
            start += width
        end
    end
    turn(a, b, c) = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
    crossing((a, b), (c, d)) = turn(a, b, c) * turn(a, b, d) < 0 && turn(c, d, a) * turn(c, d, b) < 0
    function inside(q, polygon)
        within = false
        for (a, b) in zip(polygon, circshift(polygon, -1))
            (a[2] > q[2]) != (b[2] > q[2]) && q[1] < a[1] + (q[2] - a[2]) / (b[2] - a[2]) * (b[1] - a[1]) && (within = !within)
        end
        return within
    end
    sides(polygon) = collect(zip(polygon, circshift(polygon, -1)))
    boundary = [place[p] for p in outer]
    # A sector as a polygon, its apex a little off the point of the core, which it may then share.
    sector(p, direction, spread, radius) =
        [place[p] .+ 0.04radius .* (cos(direction), sin(direction));
         [place[p] .+ radius .* (cos(a), sin(a)) for a in range(direction - spread, direction + spread; length=7)]]
    placed = Tuple{Vector{NTuple{2,Float64}},NTuple{2,Float64},Float64}[]
    function free(polygon, centre, radius)
        for segment in core
            any(side -> crossing(side, segment), sides(polygon)) && return false
            any(q -> inside(q, polygon), segment) && return false
        end
        for (other, c, r) in placed
            hypot((c .- centre)...) > r + radius && continue
            any(side -> any(s -> crossing(side, s), sides(other)), sides(polygon)) && return false
            (any(q -> inside(q, other), polygon) || any(q -> inside(q, polygon), other)) && return false
        end
        return true
    end
    branches = sort!([(p, q) for p in keys(place) for q in around[p]]; by=((p, q),) -> -leaves(q, p))
    for (p, q) in branches
        levels, fan = depth(q, p), leaves(q, p)
        reach = (levels + 0.3) * CURVE_TREE_STEP
        best = (0.0, 0.0, 0.0, 0.0)
        widest = clamp(fan * deg2rad(7), deg2rad(14), CURVE_TREE_ANGLE)
        for direction in range(0, 2π; length=CURVE_TREE_DIRECTIONS + 1)[1:end-1], spread in (widest, 2widest / 3, widest / 3)
            low, high = 0.0, reach
            free(sector(p, direction, spread, 0.05reach), place[p], 0.05reach) || continue
            for _ in 1:9
                middle = (low + high) / 2
                free(sector(p, direction, spread, middle), place[p], middle) ? (low = middle) : (high = middle)
            end
            # Out of the outer cycle of the core rather than in one of its faces.
            outward = !inside(place[p] .+ 0.5low .* (cos(direction), sin(direction)), boundary)
            score = low * sqrt(spread / widest) * (outward ? CURVE_TREE_OUTSIDE : 1)
            score > best[1] && (best = (score, direction, spread, low))
        end
        _, direction, spread, radius = best
        # A tree with no room at all is drawn tiny, towards the first direction.
        radius = max(radius, 0.05reach)
        push!(placed, (sector(p, direction, spread, radius), place[p], radius))
        hang!(q, p, p, radius / (levels + 0.3), direction - spread, direction + spread, 1)
    end
    return maximum(x) - minimum(x) >= maximum(y) - minimum(y) ? (x, y) : (y, x)
end

"`E - V + C` of the `n` vertices and the edges among `ends`; the legs, whose far ends follow the vertices, are not counted."
function curve_genus(ends, n)
    parent = collect(1:n)
    root(v) = parent[v] == v ? v : (parent[v] = root(parent[v]))
    bounded = filter(((_, v),) -> v <= n, ends)
    foreach(((u, v),) -> parent[root(u)] = root(v), bounded)
    return length(bounded) - n + count(v -> root(v) == v, 1:n)
end

"""
The TikZ commands of the map, for a picture or a scope: the points at `x`,
`y`, the edges with their cells, multiple edges bent apart, the legs as
arrows, and on top the vertices with their numbers `vertices`, as walls of
the G2-Grassmannian where `on_walls`. With `regions`, every edge has the
colour `cellC` of its cell `C`, dashed for the cells `dashed`, over a
translucent band of that colour, the region of the cell. With `plain`, there
are no numbers and no cells: the vertices are dots, and the edges `bold`
are drawn thick, over the others, one after the other, each over those
before it that it crosses with a `curve core border` of `CURVE_BORDER_PT`
on each side, `unit` centimetres to a unit of the picture, the core
`core_pt` thick; with `cell_labels`, the edges still have their cells,
`label_spacing` apart. Unless `numbered`, the vertices are dots without
their numbers.
"""
function curve_map_commands(io::IO, edges, ends, vertices, on_walls, x, y; regions=false, dashed=Set{Int}(),
                            plain=false, bold=Set{Int}(), unit=1.0, core_pt=CURVE_CORE_PT, cell_labels=false,
                            label_spacing=CURVE_LABEL_SPACING, numbered=true)
    for k in eachindex(x)
        @printf(io, "\\coordinate (c%d) at (%.4f,%.4f);\n", k, x[k], y[k])
    end
    multiplicity = Dict{Tuple{Int,Int},Int}()
    foreach(p -> multiplicity[minmax(p...)] = get(multiplicity, minmax(p...), 0) + 1, ends)
    seen = Dict{Tuple{Int,Int},Int}()
    connectors = map(ends) do (u, v)
        u == v && return "to[out=60,in=120,looseness=12]"
        k = seen[minmax(u, v)] = get(seen, minmax(u, v), 0) + 1
        bend = 40 * (k - (multiplicity[minmax(u, v)] + 1) / 2) * (u < v ? 1 : -1)
        iszero(bend) ? "--" : "to[bend left=$bend]"
    end
    paths = ["(c$u) $connector (c$v)" for ((u, v), connector) in zip(ends, connectors)]
    if regions
        for cell in unique(e.cell for e in edges)
            mine = [e.cell == cell for e in edges]
            joins = [@sprintf("(%.4f,%.4f) -- (%.4f,%.4f)", p..., q...) for (p, q) in curve_region_joins(ends[mine], x, y)]
            println(io, "\\draw[curve band,draw=cell$cell] ", join([paths[mine]; joins], " "), ";")
        end
    end
    for (k, (e, path)) in enumerate(zip(edges, paths))
        k in bold && continue
        style = e.head == 0 ? "curve leg" : "curve edge"
        regions && (style *= ",draw=cell$(e.cell)" * (e.cell in dashed ? ",dashed" : ""))
        println(io, "\\draw[$style] $path;")
    end
    drawn = Int[]
    for k in sort!(collect(bold))
        (u, v), d = ends[k], (x[ends[k][2]] - x[ends[k][1]], y[ends[k][2]] - y[ends[k][1]])
        for j in drawn
            curve_crossing(x, y, ends, j, k) || continue
            (a, b) = ends[j]
            e = (x[b] - x[a], y[b] - y[a])
            # The crossing at `t` along `k`, and a stretch of `k` wide enough to clear `j` and its border.
            t = ((x[a] - x[u]) * e[2] - (y[a] - y[u]) * e[1]) / (d[1] * e[2] - d[2] * e[1])
            sine = abs(d[1] * e[2] - d[2] * e[1]) / (hypot(d...) * hypot(e...))
            half = (core_pt / 2 + 2CURVE_BORDER_PT) / 28.4528 / unit / max(sine, 0.3) / hypot(d...)
            @printf(io, "\\draw[curve core border] (%.4f,%.4f) -- (%.4f,%.4f);\n",
                    x[u] + (t - half) * d[1], y[u] + (t - half) * d[2], x[u] + (t + half) * d[1], y[u] + (t + half) * d[2])
        end
        println(io, "\\draw[curve core] $(paths[k]);")
        push!(drawn, k)
    end
    placed = Tuple{Float64,Float64}[]
    for (e, (u, v), connector) in zip(edges, ends, connectors)
        plain && !cell_labels && break
        at(t) = (x[u] + t * (x[v] - x[u]), y[u] + t * (y[v] - y[u]))
        free = [t for t in CURVE_LABEL_POSITIONS if all(q -> hypot((at(t) .- q)...) >= label_spacing, placed)]
        t = isempty(free) ? 0.5 : first(free)
        push!(placed, at(t))
        println(io, "\\path (c$u) $connector node[curve cell,pos=$t] {$(e.cell)} (c$v);")
    end
    numbers = numbered && !plain
    for (k, (v, wall)) in enumerate(zip(vertices, on_walls))
        println(io, "\\node[", numbers ? "curve vertex" : "curve dot", ",", wall ? "curve g2 wall" : "curve bergman",
                "] at (c$k) {", numbers ? v : "", "};")
    end
end

"The middles of the edges and legs `ends` between the points at `x`, `y`."
curve_middles(ends, x, y) = [((x[u] + x[v]) / 2, (y[u] + y[v]) / 2) for (u, v) in ends]

"""
The colour of each of `cells`, by its index in `CURVE_CELL_COLOURS`: the
first cells take them in order, and every later one the colour whose cells
are farthest from it on the map, between the middles `middle` of `edges`.
"""
function curve_cell_colours(cells, edges, middle)
    members(c) = findall(e -> e.cell == c, edges)
    apart(a, b) = minimum(hypot((middle[k] .- middle[l])...) for k in members(a) for l in members(b))
    colour = Dict{Int,Int}()
    for (i, c) in enumerate(cells)
        colour[c] = i <= length(CURVE_CELL_COLOURS) ? i :
            argmax(j -> minimum(apart(c, d) for d in keys(colour) if colour[d] == j), eachindex(CURVE_CELL_COLOURS))
    end
    return colour
end

"""
The places on the one-page atlas of `k` cells, in the picture
`CURVE_PAGE_PICTURE`: two columns of frames of width `w` and height `h` for
the trees of the cells, drawn at `scale`, and between them the map at `x`,
`y`, `unit` centimetres to an edge, shifted by `offset`.
"""
function curve_page_geometry(x, y, k)
    W, H = CURVE_PAGE_PICTURE
    h = min(6.4, H / cld(k, 2))
    scale = min(1.0, (h - 2.2) / 2G2_TREE_LEAF_RADIUS)
    w = 2G2_TREE_LEAF_RADIUS * scale + 1.4
    left, right = w + 1.5, W - w - 1.5
    unit = min((right - left) / maximum(x), H / maximum(y))
    offset = (left + (right - left - unit * maximum(x)) / 2, (H - unit * maximum(y)) / 2)
    return (; unit, offset, w, h, scale)
end

"""
The picture of the one-page atlas: the map with its cells as coloured
regions, and the G2-tree of every cell of `cells` in a frame of its colour,
dashed for the cells `dashed`. The cells whose middle is left of that of the
map have their trees in the column on the left, the others on the right, in
the order of the heights of their middles, each with a line from beside its
title, clear of the tree, to the nearest middle of an edge of its region.
"""
function curve_page_tikz(io::IO, edges, ends, vertices, on_walls, x, y, cells, dashed, geometry, source)
    W, H = CURVE_PAGE_PICTURE
    (; unit, offset, w, h, scale) = geometry
    middle = curve_middles(ends, x, y)
    members(c) = findall(e -> e.cell == c, edges)
    centre(c) = reduce((p, q) -> p .+ q, middle[members(c)]) ./ length(members(c))
    order = sort(cells; by=c -> centre(c)[1])
    frames = []
    for (side, column) in ((-1, order[1:cld(end, 2)]), (1, order[cld(end, 2)+1:end]))
        for (i, c) in enumerate(sort(column; by=c -> -centre(c)[2]))
            fx, fy = side < 0 ? w / 2 : W - w / 2, H - (i - 0.5) * H / length(column)
            foot = (fx - side * w / 2, fy + h / 2 - .45)
            anchor = argmin(p -> hypot((p .- foot)...), [offset .+ unit .* middle[k] for k in members(c)])
            push!(frames, (cell=c, x=fx, y=fy, foot, anchor))
        end
    end
    colours = g2_minimal_tree_styles(io, :color)
    println(io, "\\begin{tikzpicture}[x=1cm,y=1cm]")
    println(io, "\\useasboundingbox (0,0) rectangle ($W,$H);")
    for f in frames
        @printf(io, "\\draw[curve leader,draw=cell%d] (%.3f,%.3f) -- (%.3f,%.3f);\n", f.cell, f.foot..., f.anchor...)
        @printf(io, "\\fill[cell%d] (%.3f,%.3f) circle (2.4pt);\n", f.cell, f.anchor...)
    end
    @printf(io, "\\begin{scope}[shift={(%.3fcm,%.3fcm)},x=%.4fcm,y=%.4fcm]\n", offset..., unit, unit)
    curve_map_commands(io, edges, ends, vertices, on_walls, x, y; regions=true, dashed)
    println(io, "\\end{scope}")
    for f in frames
        e = edges[first(members(f.cell))]
        @printf(io, "\\draw[curve frame,draw=cell%d%s] (%.3f,%.3f) rectangle (%.3f,%.3f);\n", f.cell,
                f.cell in dashed ? ",dashed" : "", f.x - w / 2, f.y - h / 2 + .12, f.x + w / 2, f.y + h / 2 - .12)
        @printf(io, "\\node[curve frame title] at (%.3f,%.3f) {\\textbf{cell %d}\\\\orbit %d, Type %d};\n",
                f.x, f.y + h / 2 - .22, f.cell, e.orbit, e.number)
        @printf(io, "\\begin{scope}[shift={(%.3f,%.3f)}]\n", f.x, f.y - .5)
        g2_tree_commands(io, g2_maximal_cell_tree(f.cell; source).tree, colours; label=true, scale, taxon_scale=0.7)
        println(io, "\\end{scope}")
    end
    println(io, "\\end{tikzpicture}")
    println(io, "\\endgroup")
end

"""
The entry of the catalog for `cell`: its numbers, its G2-tree, and the edges
and legs `members` that lie in it, by their vertices if `numbered`, or else
their numbers.
"""
function curve_catalog_entry(io::IO, cell, members, source; numbered=true)
    edges = ["$(e.tail)--$(e.head)" for e in members if e.head != 0]
    legs = [string(e.tail) for e in members if e.head == 0]
    if !numbered
        e = first(members)
        println(io, "\\begin{minipage}[t]{$(CURVE_ENTRY_WIDTH_CM)cm}\\centering")
        println(io, "{\\large\\bfseries cell $cell}\\\\[1pt]{\\small orbit $(e.orbit), Type $(e.number)}\\\\[3pt]")
        print(io, g2_tree_tikz(g2_maximal_cell_tree(cell; source).tree; palette=:color, label=true, taxon_scale=0.8))
        println(io, "\\\\[3pt]\\scriptsize $(length(edges)) edge", length(edges) == 1 ? "" : "s",
                isempty(legs) ? "" : ", $(length(legs)) leg" * (length(legs) == 1 ? "" : "s"), ".")
        println(io, "\\end{minipage}\\hspace{6mm plus 4mm}%")
        return
    end
    e = first(members)
    println(io, "\\begin{minipage}[t]{$(CURVE_ENTRY_WIDTH_CM)cm}\\centering")
    println(io, "{\\large\\bfseries cell $cell}\\\\[1pt]{\\small orbit $(e.orbit), Type $(e.number)}\\\\[3pt]")
    print(io, g2_tree_tikz(g2_maximal_cell_tree(cell; source).tree; palette=:color, label=true,
                           taxon_scale=0.8))
    println(io, "\\\\[3pt]\\raggedright\\scriptsize")
    isempty(edges) || print(io, "$(length(edges)) edge", length(edges) == 1 ? "" : "s", ": ", join(edges, ", "), ". ")
    isempty(legs) || print(io, "$(length(legs)) leg", length(legs) == 1 ? "" : "s", " at ", join(legs, ", "), ".")
    println(io, "\n\\end{minipage}\\hspace{6mm plus 4mm}%")
end

"""
The LaTeX source of the atlas of the decorated graph `path`, without its
pending trees unless `pending_trees`. With at most `CURVE_PAGE_CELLS` cells
it is one page, see `curve_page_tikz`; otherwise a map, with edges of `unit`
centimetres or, by default, of the length that `CURVE_EDGE_LENGTH_CM`
describes, and the catalog of the trees of its cells. Without the pending
trees, the map is drawn by `curve_core_first_layout` when it can. Unless
`numbered`, its vertices are dots without their numbers.
"""
function curve_atlas_tex(path; unit=nothing, pending_trees=true, numbered=true)
    edges = read_curve_graph(path)
    # A vertex whose edges and legs all lie in one cell is inside it, on a wall
    # of the Bergman fan of the linear space; the others are on walls of the
    # G2-Grassmannian. Dropping the pending trees keeps their kind.
    cells_at = Dict{Int,Set{Int}}()
    for e in edges, v in (e.tail, e.head)
        v == 0 || push!(get!(cells_at, v, Set{Int}()), e.cell)
    end
    pending_trees || (edges = drop_pending_trees(edges))
    vertices = sort!(unique(v for e in edges for v in (e.tail, e.head) if v != 0))
    on_walls = [length(cells_at[v]) > 1 for v in vertices]
    n = length(vertices)
    ends = curve_ends(edges, Dict(v => k for (k, v) in enumerate(vertices)))
    cells = sort!(unique(e.cell for e in edges); by=c -> (only(unique(e.number for e in edges if e.cell == c)), c))
    one_page = length(cells) <= CURVE_PAGE_CELLS
    N = n + count(e -> e.head == 0, edges)
    core_first = one_page || pending_trees ? nothing : curve_core_first_layout(ends, N, Set(eachindex(edges)))
    # The map takes edges of length about 1, as `curve_layout` draws them.
    isnothing(core_first) ||
        (core_first = core_first ./ sort([hypot(core_first[1][u] - core_first[1][v], core_first[2][u] - core_first[2][v])
                                          for (u, v) in ends])[cld(length(ends), 2)])
    x, y = one_page ? curve_cell_layout(edges, ends, cells, N) :
           isnothing(core_first) ? curve_layout(ends, N) : core_first
    x .-= minimum(x)
    y .-= minimum(y)
    geometry = one_page ? curve_page_geometry(x, y, length(cells)) : nothing
    colour = one_page ? curve_cell_colours(cells, edges, curve_middles(ends, x, y)) : Dict{Int,Int}()
    unit = one_page ? geometry.unit :
        something(unit, max(CURVE_EDGE_LENGTH_CM, min((42 - 4) / maximum(x), (29.7 - 7) / maximum(y))))
    size(points) = round(points * max(1, unit / CURVE_EDGE_LENGTH_CM); digits=2)
    dashed = Set(cells[length(CURVE_CELL_COLOURS)+1:end])
    width = one_page ? 42.0 : max(42.0, unit * maximum(x) + 4)
    height = one_page ? 29.7 : max(29.7, unit * maximum(y) + 7)
    source = load_g2_maximal_cell_source(; field=:qq)
    return sprint() do io
        print(io, """
        \\documentclass{article}
        \\usepackage[paperwidth=$(width)cm,paperheight=$(height)cm,margin=1cm]{geometry}
        \\usepackage[T1]{fontenc}
        \\usepackage{lmodern}
        \\usepackage{xcolor}
        \\usepackage{tikz}
        \\usetikzlibrary{arrows.meta}
        \\renewcommand{\\familydefault}{\\sfdefault}
        \\pagestyle{empty}
        \\setlength{\\parindent}{0pt}
        \\tikzset{
          curve edge/.style={draw=black!75,line width=$(size(.6))pt},
          curve leg/.style={draw=black!40,line width=$(size(.6))pt,-{Stealth[length=$(size(4))pt]}},
          curve cell/.style={fill=white,inner sep=$(size(.7))pt,text=$(one_page ? "black!80" : "blue!55!black"),font=\\sffamily\\fontsize{$(size(5.5))}{$(size(6))}\\selectfont},
          curve vertex/.style={circle,draw=black,line width=$(size(.4))pt,inner sep=$(size(.4))pt,minimum size=$(size(8))pt,font=\\sffamily\\fontsize{$(size(4.5))}{$(size(5))}\\selectfont},
          curve key/.style={circle,draw=black,line width=.4pt,minimum size=$(numbered ? 8 : 5)pt},
          curve dot/.style={circle,draw=black,line width=$(size(.3))pt,inner sep=0pt,minimum size=$(size(4))pt},
          curve g2 wall/.style={$(one_page ? "fill=black!80,text=white" : "fill=orange!55")},
          curve bergman/.style={$(one_page ? "fill=white" : "fill=green!45!cyan!40")},
          curve band/.style={line width=$(round(0.36unit; digits=3))cm,line cap=round,line join=round,draw opacity=.25},
          curve leader/.style={line width=.8pt},
          curve frame/.style={line width=1.4pt,rounded corners=4pt,fill=white},
          curve frame title/.style={anchor=north,align=center,font=\\sffamily\\small},
        }
        $(join(("\\definecolor{cell$c}{HTML}{$(CURVE_CELL_COLOURS[j])}" for (c, j) in colour), "\n"))
        \\begin{document}
        {\\LARGE\\bfseries $(replace(basename(path), "_" => "\\_"))}\\par\\medskip
        """)
        # Only the atlas with a catalog explains its map; the one-page atlas has its title alone.
        one_page || print(io, """
        $(pending_trees ? "" : "Without the pending trees and their legs; a vertex of degree 2 joins edges in different cells. ")
        $n vertices, $(count(e -> e.head != 0, edges)) edges, $(count(e -> e.head == 0, edges)) legs (arrows),
        genus $(curve_genus(ends, n)). $(numbered ? "Vertices are numbered as in the graph; each" : "Each") edge and leg is labelled
        in blue by the maximal cell of the coarsest tropical G\\(_2\\)-Grassmannian that holds it,
        \\texttt{data/g2/grassmannian/02\\_coarsest\\_fan.mrdi}. The trees of its $(length(cells)) cells follow the map.
        \\tikz[baseline=-.6ex]\\node[curve key,curve g2 wall] {};~$(count(on_walls)) vertices on walls of the
        G\\(_2\\)-Grassmannian, whose edges and legs lie in several cells;
        \\tikz[baseline=-.6ex]\\node[curve key,curve bergman] {};~$(count(!, on_walls)) inside a cell, on walls of the
        Bergman fan, whose edges and legs all lie in that cell$(pending_trees ? "" : ", counting those of the pending trees").\\par
        """)
        println(io, "\\vfill\\centering")
        if one_page
            curve_page_tikz(io, edges, ends, vertices, on_walls, x, y, cells, dashed, geometry, source)
            println(io, "\\vfill\n\\end{document}")
            return
        end
        println(io, "\\begin{tikzpicture}[x=$(unit)cm,y=$(unit)cm]")
        curve_map_commands(io, edges, ends, vertices, on_walls, x, y; numbered)
        println(io, raw"\end{tikzpicture}")
        print(io, """
        \\vfill\\clearpage
        {\\LARGE\\bfseries The $(length(cells)) maximal cells}\\par\\medskip
        The G\\(_2\\)-tree of each cell, of the sum of the balanced generators of its rays, as in the tables of maximal
        trees of \\texttt{main.tex}; cells ordered by their number there. Below each tree, its edges and legs$(numbered ? "" : ", counted").\\par\\bigskip
        {\\raggedright\\setlength{\\lineskip}{7mm}\\lineskiplimit=0pt
        """)
        foreach(c -> curve_catalog_entry(io, c, filter(e -> e.cell == c, edges), source; numbered), cells)
        println(io, "\\par}\n\\end{document}")
    end
end

"""
The LaTeX source of the plain map of the decorated graph `path`, a picture
of `CURVE_PLAIN_SIZE` for half a page of main.tex: no numbers and no cells,
the vertices dots, orange on walls of the G2-Grassmannian and green inside
a cell, the legs short arrows, and the core in bold, drawn first by
`curve_core_first_layout`, the pending trees hanging from it. Without
`pending_trees`, it is the map of the core alone, see `drop_pending_trees`,
`CURVE_CORE_ONLY_PT` thick, each edge labelled by its cell.
"""
function curve_plain_tex(path; pending_trees=true)
    edges = read_curve_graph(path)
    cells_at = Dict{Int,Set{Int}}()
    for e in edges, v in (e.tail, e.head)
        v == 0 || push!(get!(cells_at, v, Set{Int}()), e.cell)
    end
    pending_trees || (edges = drop_pending_trees(edges))
    vertices = sort!(unique(v for e in edges for v in (e.tail, e.head) if v != 0))
    on_walls = [length(cells_at[v]) > 1 for v in vertices]
    n = length(vertices)
    ends = curve_ends(edges, Dict(v => k for (k, v) in enumerate(vertices)))
    N = n + count(e -> e.head == 0, edges)
    core_numbers = Set(e.edge for e in curve_core(edges))
    bold = Set(findall(e -> e.edge in core_numbers, edges))
    layout = curve_core_first_layout(ends, N, bold)
    isnothing(layout) && @warn "The core could not be drawn first; the whole graph is laid out at once"
    x, y = isnothing(layout) ? curve_layout(ends, N) : layout
    x .-= minimum(x)
    y .-= minimum(y)
    W, H = CURVE_PLAIN_SIZE
    unit = round(min(W / maximum(x), H / maximum(y)); digits=4)
    core_pt = pending_trees ? CURVE_CORE_PT : CURVE_CORE_ONLY_PT
    return sprint() do io
        print(io, """
        % $(basename(path)): $n vertices, $(count(e -> e.head != 0, edges)) edges, $(count(e -> e.head == 0, edges)) legs,
        % genus $(curve_genus(ends, n)), $(pending_trees ? "" : "without the pending trees, ")core of $(length(bold)) edges in bold; $(count(on_walls)) vertices on walls of the
        % G2-Grassmannian (orange), $(count(!, on_walls)) inside a cell (green). Written by build_curve_atlas.jl --plain.
        \\documentclass[tikz,border=1pt]{standalone}
        \\usetikzlibrary{arrows.meta}
        \\tikzset{
          curve edge/.style={draw=black!70,line width=.4pt},
          curve leg/.style={draw=black!45,line width=.35pt,-{Stealth[length=2.4pt,width=1.9pt]}},
          curve core/.style={draw=black,line width=$(core_pt)pt,line cap=round},
          curve core border/.style={draw=white,line width=$(core_pt + 2CURVE_BORDER_PT)pt,line cap=butt},
          curve dot/.style={circle,draw=black,line width=$(pending_trees ? .2 : .35)pt,inner sep=0pt,minimum size=$(pending_trees ? 2.3 : 4.2)pt},
          curve cell/.style={fill=white,inner sep=.5pt,text=blue!55!black,font=\\sffamily\\fontsize{6.5}{7}\\selectfont},
          curve g2 wall/.style={fill=orange!55},
          curve bergman/.style={fill=green!45!cyan!40},
        }
        \\begin{document}
        \\begin{tikzpicture}[x=$(unit)cm,y=$(unit)cm]
        """)
        curve_map_commands(io, edges, ends, vertices, on_walls, x, y; plain=true, bold, unit, core_pt,
                           cell_labels=!pending_trees, label_spacing=0.35 / unit)
        println(io, "\\end{tikzpicture}\n\\end{document}")
    end
end

"Compile `tex` with pdflatex in the work directory, to `pdf`."
function compile_curve_atlas(tex, pdf)
    pdflatex = Sys.which("pdflatex")
    pdflatex === nothing && error("pdflatex is required")
    build = work_path("tex", "curve_atlas")
    mkpath(build)
    tex = abspath(tex)
    # From the folder of the file, so that it may include the files beside it by their names.
    run(Cmd(Cmd([pdflatex, "-interaction=nonstopmode", "-halt-on-error", "-file-line-error",
                 "-output-directory=$build", tex]); dir=dirname(tex)))
    cp(joinpath(build, splitext(basename(tex))[1] * ".pdf"), pdf; force=true)
    return pdf
end

if abspath(PROGRAM_FILE) == @__FILE__
    graphs = filter(arg -> !startswith(arg, "--"), ARGS)
    length(graphs) == 1 || error("usage: build_curve_atlas.jl GRAPH.csv [--drop-pending-trees] [--plain] " *
                                 "[--unnumbered] [--output=ATLAS.pdf] [--edge-length=CM] [--no-pdf]")
    for arg in filter(arg -> startswith(arg, "--"), ARGS)
        arg in ("--no-pdf", "--drop-pending-trees", "--plain", "--unnumbered") || startswith(arg, "--output=") ||
            startswith(arg, "--edge-length=") || error("unknown option $arg")
    end
    graph = abspath(only(graphs))
    pending_trees = !("--drop-pending-trees" in ARGS)
    plain = "--plain" in ARGS
    suffix = (pending_trees ? "" : "_core") * (plain ? "_plain.pdf" : "_atlas.pdf")
    pdf = something(g2_tree_option(ARGS, "output"), splitext(graph)[1] * suffix)
    edge_length = g2_tree_option(ARGS, "edge-length")
    unit = isnothing(edge_length) ? nothing : parse(Float64, edge_length)
    tex = splitext(abspath(pdf))[1] * ".tex"
    numbered = !("--unnumbered" in ARGS)
    write(tex, plain ? curve_plain_tex(graph; pending_trees) : curve_atlas_tex(graph; unit, pending_trees, numbered))
    println("Wrote $tex")
    "--no-pdf" in ARGS || println("Wrote ", compile_curve_atlas(tex, abspath(pdf)))
end
