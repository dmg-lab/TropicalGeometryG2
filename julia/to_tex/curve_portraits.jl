###############################################################################
# Two pictures of a decorated curve graph, drawn from one drawing of its core,
# the graph without its legs and pending trees (see curve_core in
# build_curve_atlas.jl), every vertex kept, planar or with one crossing:
#
# - PREFIX_berkovich.pdf, the whole curve without labels: the core in bold,
#   its vertices small dots, those of degree 3 in
#   red (for a trivalent core of genus g, the 2g - 2 points of its canonical
#   divisor), and every pending tree a bushy fan from its point of the core,
#   each level shorter and thinner than the one before, its legs tapering to a
#   point; the fans push each other apart. With POINTS, the points of
#   point_places.csv of curve_details.jl are red instead, at their places on
#   the core or on the trees (a point on a leg at its middle).
# - PREFIX_core_regions.pdf, an atlas of the core: the convex hull of its
#   drawing cut into one region for every maximal cell, which holds the edges
#   and vertices of that cell, drawn in its colour; a vertex between two cells
#   is split between their colours, the core as large as in the first picture.
#   The G2-tree of every cell stands beside the map in a frame of its colour.
#   With more cells than colours, two cells whose regions meet differ; with
#   --connected-regions every cell has a colour of its own, up to 13.
#
#   julia --project=. -t 8 to_tex/curve_portraits.jl \
#     GRAPH.csv [--output-prefix=PREFIX] [--outer-face=K] [--seed=S] [--no-pdf]
#     [--fragments=FIGURE [--standalone]] [--points=POINTS] [--crossing=C] [--side-by-side]
#     [--chain-power=P] [--tree-scale=F] [--frames] [--trivalent-nodes]
#     [--connected-regions] [--region-starts=K] [--region-moves=M] [--region-drawing=FILE] [--region-smoothing=S]
#     [--atlas-scale=A] [--overlapping-trees]
#
# GRAPH.csv is a reduced graph of decorated_curve_graph.jl; PREFIX is GRAPH by
# default. With FIGURE, the two pictures are also written as FIGURE_berkovich.tex
# and FIGURE_regions.tex (with --standalone, as the documents of
# g2_standalone_figure, whose PDFs main.tex includes), for main.tex, which they fit on one page
# together; with --side-by-side, in one row of the width of main.tex, the
# types of the cells, coloured, in their regions instead of the frames, one
# without room in its region beside it with a line to it (with --frames the
# frames, smaller, their columns as high as the Berkovich picture), and
# PREFIX_side_by_side.pdf shows the row. A core that is not planar is drawn
# with one crossing, the C-th of `portrait_crossings`, by default the best of
# all, where one edge passes over the other with a thin border in the colour below; the
# crossing stays where the two segments meet. The drawings are made on all
# threads. A chain of k edges between points of degree 3 is drawn as long as
# k^P edges, by default k. With --trivalent-nodes, only the vertices of degree
# 3 of the core are drawn. With F, the fans are at most F times as long and
# push each other apart within F times the distance, by default 1; with
# --overlapping-trees, they may cross each other, never the core. The outer face of the core is one that meets the most cells, the
# K-th of them (default 1) in the order of `portrait_outer_faces`. The core is
# drawn by Tutte's method from it, then relaxed by stress majorization towards
# its graph distances, shortened within a cell and lengthened between cells,
# and relaxed again with a point in each bounded face drawn to its corners, so
# that every face leaves room for the trees at its corners; no move makes two
# edges cross. The regions cut every face as `portrait_partition_face` says.
#
# With --connected-regions, both pictures draw the core as the atlas wants, see
# `portrait_region_drawing`, from K starts annealed by M moves each (default 12
# and 50000), in which the region of every cell is connected and the edges of
# two cells cross as little as possible, those of one cell more freely: a point
# of the hull, or of a margin around it, belongs to the cell of the nearest
# edge; where edges cross, the first of `portrait_all_crossings` passes over
# the other in both. The core of the atlas is then A times as large as in the
# Berkovich picture (default 1, as large), and side by side the two pictures
# stand in boxes of one height, their cores at the same height.
# The drawing is written to
# PREFIX_region_drawing.csv, the place of every vertex of the graph; with FILE
# it is read from there instead, lines starting with # aside. With S, it is
# then smoothed by S steps of the flow of `portrait_region_smoothing`, which
# straightens the edges and makes the regions rounder, every vertex and edge
# staying in its region, and the smoothed drawing is written to
# PREFIX_region_drawing.csv.
###############################################################################

include(joinpath(@__DIR__, "build_curve_atlas.jl"))

# The distances in the drawing of the core between vertices of one cell, and
# of no common cell, over their distances in the graph.
const PORTRAIT_SAME_CELL = 0.8
const PORTRAIT_OTHER_CELL = 1.15
# A chain of k edges between points of degree 3 is drawn as long as k^this
# many edges, by default; see --chain-power.
const PORTRAIT_CHAIN_POWER = 1.0
# The point of a bounded face is drawn to its corners from this many edges
# times the square root of its trees, with weight this many times the stress.
const PORTRAIT_FACE_REACH = 0.5
const PORTRAIT_FACE_PULL = 10.0
const PORTRAIT_SWEEPS = (400, 600)
const PORTRAIT_CLEARANCE = 0.35
# The second relaxation sweeps the points in the order of these seeds, the
# first in their order and the others shuffled, and the drawing whose
# smallest fan is largest is kept.
const PORTRAIT_SEEDS = 0:7
# A fan: its first edge at most this long, core edges of median length 1, a
# level this much shorter than the one before, a leg this much shorter than an
# edge at its depth; the children of a vertex spread over at least the first
# angle each, at most the second in all; a fan at most the half angle, in
# degrees, and the room of a gap between core edges shared by its fans as
# (fans + 1)^PORTRAIT_GAP_SHARE.
const PORTRAIT_FIRST_EDGE = 0.42
const PORTRAIT_LEVEL = 0.88
const PORTRAIT_LEG = 0.9
const PORTRAIT_SPREAD = (45, 140)
const PORTRAIT_FAN_HALF = 70
const PORTRAIT_GAP_SHARE = 1.5
# The fans of the Berkovich picture push each other apart within this
# distance, core edges of median length 1, see `portrait_repel!`, rewarded
# this much for the log of their first edges, in at most this many sweeps of
# moves turning or opening by at most this many degrees.
const PORTRAIT_REPEL = (reach=0.25, reward=1.0, sweeps=120, turn=4)
# Strokes of the Berkovich picture, in pt: the core, the first level of a fan
# (thinner by 0.92 at every level, at least .3), the dots of the core, and
# those of the divisor, its vertices of degree 3, in red.
const PORTRAIT_CORE_PT, PORTRAIT_EDGE_PT, PORTRAIT_DOT_PT, PORTRAIT_DIVISOR_PT = 1.5, 0.62, 1.55, 2.1
const PORTRAIT_DIVISOR_COLOUR = "D62728"
# The largest width and height of the Berkovich picture, in cm. With the
# height of the region atlas, the two pictures and their captions fill at
# most a page of main.tex (6.5in x 8.8in).
const PORTRAIT_BERKOVICH_SIZE = (16.5, 9.6)
# The region atlas, in cm, its core at the scale of the Berkovich picture: its
# height, that of the columns of frames unless the map is higher, its frames
# and their gap to the map, its width that of the map and the frames; in pt
# its margin around the hull and its vertices.
const PORTRAIT_ATLAS_HEIGHT = 8.4
const PORTRAIT_FRAME = (2.45, 2.45)
const PORTRAIT_FRAME_GAP = 0.3
const PORTRAIT_TREE_SCALE, PORTRAIT_TAXON_SCALE = 0.45, 0.42
# The regions in colours of the figures of main.tex, g2TaxonBlue, muted,
# g2HeartRed, g2TaxonYellow and g2AppleBody, in the order of the cells of
# `portrait_curve`.
const PORTRAIT_REGION_COLOURS = ["1676A3", "647584", "D7443E", "E39A24", "7CB342"]
# With --connected-regions, every cell its own colour, of those and eight more,
# no two closer than 15 in CIELAB as lines or as the tints of the regions.
const PORTRAIT_MORE_COLOURS = ["7B1FA2", "1A237E", "D81B60", "4DB6AC", "8B5E3C", "004D40", "ADA200", "D96A9A"]
const PORTRAIT_REGION_TINT = 0.22
const PORTRAIT_MARGIN_PT = 5.0
const PORTRAIT_NODE_PT = 2.5
const PORTRAIT_PT = 1 / 28.4528
# The row of --side-by-side, in cm: its width, at most its height, the gap
# between the two pictures, with --frames the least gap between two frames of
# a column, and with --connected-regions at most the height of the boxes in
# which the two cores stand at one height; the least room, in cm, of a region for its label, the
# spacing of the points searched for it, and how far beyond the hull, in cm,
# a label without that room stands.
const PORTRAIT_SIDE_BY_SIDE = (width=16.5, height=9.0, gap=0.5, frame_gap=0.15, aligned=9.66)
const PORTRAIT_LABEL_ROOM, PORTRAIT_LABEL_GRID = 0.17, 0.05
const PORTRAIT_LABEL_BESIDE = (0.6, 1.0)
# The region atlas of --connected-regions: its points along the edges this far
# apart while the core is drawn and then, as shares of the largest side of the
# core, the first and last this far from the ends of an edge; the regions this
# far beyond the hull. Between two components of the edges of a cell, the
# distance is at most this many times the square root of its edges. The
# weights of the terms of `portrait_region_energy`, the starts and moves of
# the annealing, its first temperature and step, and the core of the atlas
# this many times as large as in the Berkovich picture.
const PORTRAIT_REGION_SAMPLES = (spacing=0.025, final=0.008, ends=0.006)
const PORTRAIT_REGION_MARGIN = 0.03
const PORTRAIT_REGION_PULL = 0.35
const PORTRAIT_REGION_WEIGHTS = (split=1000.0, mixed=100.0, within=12.0, stress=1.0, clearance=30.0, bend=3.0, angle=10.0,
                                 border=20.0)
const PORTRAIT_REGION_STARTS, PORTRAIT_REGION_MOVES = 12, 50_000
const PORTRAIT_REGION_ANNEAL = (temperature=5.0, step=0.5)
const PORTRAIT_REGION_SCALE = 1.0
# The flow of --region-smoothing, lengths as shares of the largest side of the
# drawing it starts from: the weights of the terms of
# `portrait_smoothing_terms`, the step of a point, falling from the first to
# the second, the reach of the clearance, the difference of the gradient, and
# the steps between two looks at the drawing.
const PORTRAIT_SMOOTHING_WEIGHTS = (bend=10.0, even=2.0, angle=1.0, clearance=1.0, anchor=20.0, convex=15.0, round=25.0)
const PORTRAIT_SMOOTHING_STEP = (0.002, 0.0002)
const PORTRAIT_SMOOTHING_REACH, PORTRAIT_SMOOTHING_DIFFERENCE, PORTRAIT_SMOOTHING_LOOK = 0.02, 5e-4, 20
# A crossing lies this far inside both its segments, as a share of them, and
# a red point of the Berkovich picture this many pt from it.
const PORTRAIT_CROSSING_INSIDE = 0.15
const PORTRAIT_CROSSING_CLEARANCE_PT = 12.0
# The border of an edge over another in the atlas, in pt on either side.
const PORTRAIT_CROSSING_GAP_PT = 0.45

portrait_turn(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])
portrait_area(poly) = sum(portrait_turn((0.0, 0.0), poly[i], poly[mod1(i + 1, length(poly))]) for i in eachindex(poly)) / 2
"The CIELAB coordinates of the colour `hex` mixed with white, `tint` of it."
function portrait_lab(hex; tint=1.0)
    rgb = [tint * parse(Int, hex[i:i+1]; base=16) / 255 + 1 - tint for i in (1, 3, 5)]
    linear = [c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055)^2.4 for c in rgb]
    xyz = [0.4124 0.3576 0.1805; 0.2126 0.7152 0.0722; 0.0193 0.1192 0.9505] * linear ./ [0.9505, 1.0, 1.089]
    f = [t > 0.008856 ? cbrt(t) : 7.787t + 16 / 116 for t in xyz]
    return (116f[2] - 16, 500(f[1] - f[2]), 200(f[2] - f[3]))
end
function portrait_centroid(poly)
    A, cx, cy = 0.0, 0.0, 0.0
    for (p, q) in zip(poly, circshift(poly, -1))
        w = p[1] * q[2] - q[1] * p[2]
        A, cx, cy = A + w, cx + (p[1] + q[1]) * w, cy + (p[2] + q[2]) * w
    end
    return (cx / 3A, cy / 3A)
end
portrait_crossing((a, b), (c, d)) =
    portrait_turn(a, b, c) * portrait_turn(a, b, d) < 0 && portrait_turn(c, d, a) * portrait_turn(c, d, b) < 0

"""
The distances in the graph of the edges `ends` on `1:n`, an edge of a chain
of k edges between points of degree at least 3 of length k^(power - 1), see
`curve_core_chains`.
"""
function portrait_graph_distances(n, ends; power=1.0)
    D = fill(Inf, n, n)
    foreach(s -> D[s, s] = 0, 1:n)
    _, chains = curve_core_chains(ends, Set(eachindex(ends)))
    for (_, path) in chains, k in path
        u, v = ends[k]
        D[u, v] = D[v, u] = min(D[u, v], length(path)^(power - 1))
    end
    for k in 1:n, i in 1:n, j in 1:n
        D[i, j] = min(D[i, j], D[i, k] + D[k, j])
    end
    return D
end

"""
The faces of the planar embedding of the core with the chains `chains`
between its branch points, see `curve_core_chains`, as cycles of its points,
with the number of cells of `colour` on each, the edges through a `crossing`
passing through it. Returns those meeting the most
cells, then with the most points of degree 3, then by their sorted points.
"""
function portrait_outer_faces(ends, colour, n; crossing=nothing)
    branch, chains = curve_core_chains(ends, Set(eachindex(ends)))
    if !isnothing(crossing)
        # The chains at the crossing first, by their next points around it, so that their darts alternate there.
        around = [crossing.over[1], crossing.under[1], crossing.over[2], crossing.under[2]]
        next(points) = first(points) == crossing.x ? points[2] : last(points) == crossing.x ? points[end-1] : 0
        sort!(chains; by=((points, _),) -> something(findfirst(==(next(points)), around), 5))
    end
    index = Dict(p => i for (i, p) in enumerate(branch))
    pairs = [(index[first(points)], index[last(points)]) for (points, _) in chains]
    embeddings = curve_embeddings(length(branch), pairs; crossings=isnothing(crossing) ? () : (index[crossing.x],))
    isempty(embeddings) && error("the core is not planar")
    degree = zeros(Int, n)
    foreach(((u, v),) -> (degree[u] += 1; degree[v] += 1), ends)
    cycles = [reduce(vcat, (isodd(d) ? chains[(d+1)÷2][1] : reverse(chains[d÷2][1]))[1:end-1] for d in face)
              for face in first(embeddings)]
    cells(cycle) = length(unique(colour[minmax(cycle[k], cycle[mod1(k + 1, length(cycle))])] for k in eachindex(cycle)))
    filter!(allunique, cycles)
    most = maximum(cells, cycles)
    return sort(filter(c -> cells(c) == most, cycles); by=c -> (-count(v -> degree[v] == 3, c), sort(c)))
end

"""
The weighted Tutte drawing of the graph `neighbours` on `1:n`: the cycle
`outer` evenly on the unit circle, every other vertex at the mean of its
neighbours.
"""
function portrait_tutte(n, neighbours, outer)
    L = zeros(n, n)
    for u in 1:n, v in neighbours[u]
        L[u, u] += 1
        L[u, v] -= 1
    end
    X = zeros(n, 2)
    angles = 2π .* (0:length(outer)-1) ./ length(outer)
    X[outer, :] = [cos.(angles) sin.(angles)]
    inner = setdiff(1:n, outer)
    X[inner, :] = -L[inner, inner] \ (L[inner, outer] * X[outer, :])
    return [(X[k, 1], X[k, 2]) for k in 1:n]
end

"""
The point where the segments from `a` to `b` and from `p` to `q` meet,
`inside` inside both, or nothing.
"""
function portrait_meeting(a, b, p, q; inside=PORTRAIT_CROSSING_INSIDE)
    d, e, r = b .- a, q .- p, p .- a
    det = d[1] * e[2] - d[2] * e[1]
    abs(det) < 1e-12 && return nothing
    t, s = (r[1] * e[2] - r[2] * e[1]) / det, (r[1] * d[2] - r[2] * d[1]) / det
    lo, hi = inside, 1 - inside
    return lo <= t <= hi && lo <= s <= hi ? a .+ t .* d : nothing
end

"The point `crossing.x` moved to the meeting of the segments through it, or false when they do not meet."
function portrait_follow!(pos, crossing)
    (a, b), (p, q) = crossing.over, crossing.under
    y = portrait_meeting(pos[a], pos[b], pos[p], pos[q])
    isnothing(y) && return false
    pos[crossing.x] = y
    return true
end

"""
`pos` relaxed by stress majorization towards the distances `D`, with the
weights `W`, after the best scaling: point after point, in their order or
shuffled by `rng` at every sweep, a move to its place
of least stress, or part of the way, is kept only when no edge of `segments`
at the point then crosses another, and the least distance from a point to a
segment without it stays at least the smaller of its value and
`PORTRAIT_CLEARANCE`. The point of a `crossing` does not move itself but
follows the ends of its edges, see `portrait_follow!`.
"""
function portrait_relax(pos, segments, D, W; sweeps, rng=nothing, crossing=nothing)
    n = length(pos)
    incident = [Int[] for _ in 1:n]
    foreach(((k, (u, v)),) -> (push!(incident[u], k); push!(incident[v], k)), enumerate(segments))
    r(p, q) = max(hypot((p .- q)...), 1e-9)
    off = [(i, j) for i in 1:n for j in i+1:n]
    s = sum(r(pos[i], pos[j]) / D[i, j] for (i, j) in off) / sum(r(pos[i], pos[j])^2 / D[i, j]^2 for (i, j) in off)
    pos = [p .* s for p in pos]
    function crosses(k, l)
        (a, b), (c, d) = segments[k], segments[l]
        allunique((a, b, c, d)) && portrait_crossing((pos[a], pos[b]), (pos[c], pos[d]))
    end
    gap(p, (a, b)) = hypot((pos[p] .- curve_closest(pos[p], pos[a], pos[b]))...)
    function clearance(i)
        c = Inf
        for (k, (a, b)) in enumerate(segments)
            if k in incident[i]
                foreach(p -> p == a || p == b || (c = min(c, gap(p, (a, b)))), 1:n)
            elseif i != a && i != b
                c = min(c, gap(i, (a, b)))
            end
        end
        return c
    end
    x = isnothing(crossing) ? 0 : crossing.x
    through = isnothing(crossing) ? () : (crossing.over..., crossing.under...)
    for _ in 1:sweeps, i in (isnothing(rng) ? (1:n) : Random.randperm(rng, n))
        i == x && continue
        total = sum(W[i, j] for j in 1:n if j != i)
        target = reduce(.+, (W[i, j] .* (pos[j] .+ D[i, j] .* (pos[i] .- pos[j]) ./ r(pos[i], pos[j])) for j in 1:n if j != i)) ./ total
        moving = i in through ? [i, x] : [i]
        before, old = clearance.(moving), pos[moving]
        for share in (1.0, 0.5, 0.25, 0.1)
            pos[i] = old[1] .+ share .* (target .- old[1])
            (i in through ? portrait_follow!(pos, crossing) : true) &&
                !any(crosses(k, l) for m in moving for k in incident[m] for l in eachindex(segments) if l != k) &&
                all(clearance(m) >= min(b, PORTRAIT_CLEARANCE) for (m, b) in zip(moving, before)) && break
            pos[moving] = old
        end
    end
    return pos
end

"""
The faces of the straight-line planar drawing `pos` of the graph
`neighbours` on `1:n`, as cycles of vertices, counterclockwise but the outer
one.
"""
function portrait_faces(pos, neighbours)
    dir(v, w) = atan(pos[w][2] - pos[v][2], pos[w][1] - pos[v][1])
    ring = [sort(neighbours[v]; by=w -> dir(v, w)) for v in eachindex(pos)]
    seen = Set{Tuple{Int,Int}}()
    faces = Vector{Int}[]
    for u in eachindex(pos), v in neighbours[u]
        (u, v) in seen && continue
        face, a, b = Int[], u, v
        while !((a, b) in seen)
            push!(seen, (a, b))
            push!(face, a)
            a, b = b, ring[b][mod1(findfirst(==(a), ring[b]) - 1, length(ring[b]))]
        end
        push!(faces, face)
    end
    return faces
end

"""
The planar drawing of the core, the edges `ends` on `1:n` with the cells
`colour`, from the outer face `outer`, see the header, the second relaxation
in the order of `seed`, the chains as long as `chain_power` says, see
`portrait_graph_distances`; with largest side 1. With a `crossing`, its point
moves to the meeting of its segments in the drawing of Tutte, or else after
a first relaxation with the point free; nothing when they do not meet there
either, or the move makes two edges cross.
"""
function portrait_core_drawing(ends, colour, n, outer; seed=0, crossing=nothing, chain_power=PORTRAIT_CHAIN_POWER)
    neighbours = [Int[] for _ in 1:n]
    foreach(((u, v),) -> (push!(neighbours[u], v); push!(neighbours[v], u)), ends)
    cells = [Set(colour[minmax(k, w)] for w in neighbours[k]) for k in 1:n]
    G = portrait_graph_distances(n, ends; power=chain_power)
    D = [i == j ? 0.0 : G[i, j] * (isempty(intersect(cells[i], cells[j])) ? PORTRAIT_OTHER_CELL : PORTRAIT_SAME_CELL)
         for i in 1:n, j in 1:n]
    pos = portrait_tutte(n, neighbours, outer)
    weights = [i == j ? 0.0 : D[i, j]^-2 for i in 1:n, j in 1:n]
    if !isnothing(crossing)
        crossed(pos) = any(portrait_crossing((pos[a], pos[b]), (pos[c], pos[d])) for (a, b) in ends for (c, d) in ends
                           if allunique((a, b, c, d)))
        if !portrait_follow!(pos, crossing) || crossed(pos)
            pos = portrait_relax(portrait_tutte(n, neighbours, outer), ends, D, weights; sweeps=PORTRAIT_SWEEPS[1])
            portrait_follow!(pos, crossing) && !crossed(pos) || return nothing
        end
    end
    pos = portrait_relax(pos, ends, D, weights; sweeps=PORTRAIT_SWEEPS[1], crossing)
    # A point in every bounded face, drawn to its corners from a distance that grows with the
    # trees there: one at a vertex of degree 3, half of one at a vertex of degree 2, none at a crossing.
    faces = filter(f -> portrait_area(pos[f]) > 0, portrait_faces(pos, neighbours))
    trees(v) = !isnothing(crossing) && v == crossing.x ? 0.0 : length(neighbours[v]) == 3 ? 1.0 : 0.5
    reach = [PORTRAIT_FACE_REACH * sqrt(sum(trees, f)) for f in faces]
    F = length(faces)
    E = zeros(n + F, n + F)
    E[1:n, 1:n] = D
    for (a, f) in enumerate(faces)
        foreach(j -> E[n+a, j] = E[j, n+a] = reach[a] + minimum(D[v, j] for v in f), 1:n)
        for (b, g) in enumerate(faces)
            a == b || (E[n+a, n+b] = reach[a] + reach[b] + minimum(D[v, w] for v in f for w in g))
        end
    end
    W = [i == j ? 0.0 : E[i, j]^-2 for i in 1:n+F, j in 1:n+F]
    for (a, f) in enumerate(faces), v in f
        W[n+a, v] = W[v, n+a] = PORTRAIT_FACE_PULL * W[v, n+a]
    end
    pos = portrait_relax([pos; [portrait_centroid(pos[f]) for f in faces]], ends, E, W; sweeps=PORTRAIT_SWEEPS[2],
                         rng=seed == 0 ? nothing : Random.Xoshiro(seed), crossing)[1:n]
    lo = (minimum(first.(pos)), minimum(last.(pos)))
    side = maximum(maximum(p[k] for p in pos) - lo[k] for k in 1:2)
    return [(p .- lo) ./ side for p in pos]
end

function portrait_on_segment(p, a, b; tol=1e-9)
    abs(portrait_turn(a, b, p)) <= tol * max(1.0, hypot((b .- a)...)) || return false
    return min(a[1], b[1]) - tol <= p[1] <= max(a[1], b[1]) + tol && min(a[2], b[2]) - tol <= p[2] <= max(a[2], b[2]) + tol
end

"Whether `p` lies in the polygon `poly` or on its boundary."
function portrait_inside(p, poly)
    sides = zip(poly, circshift(poly, -1))
    any(portrait_on_segment(p, a, b) for (a, b) in sides) && return true
    within = false
    for (a, b) in sides
        (a[2] > p[2]) != (b[2] > p[2]) && p[1] < a[1] + (p[2] - a[2]) / (b[2] - a[2]) * (b[1] - a[1]) && (within = !within)
    end
    return within
end

"The corners of a shortest path inside the simple polygon `poly` from its corner `s` to its corner `t`."
function portrait_geodesic(poly, s, t)
    n = length(poly)
    function visible(i, j)
        (j == mod1(i + 1, n) || i == mod1(j + 1, n)) && return true
        a, b = poly[i], poly[j]
        for k in 1:n
            c, d = poly[k], poly[mod1(k + 1, n)]
            portrait_turn(a, b, c) * portrait_turn(a, b, d) < -1e-14 && portrait_turn(c, d, a) * portrait_turn(c, d, b) < -1e-14 &&
                return false
        end
        return all(portrait_inside(a .+ t .* (b .- a), poly) for t in 0.05:0.05:0.95)
    end
    distance, back, done = fill(Inf, n), zeros(Int, n), falses(n)
    distance[s] = 0
    while true
        u = argmin(k -> done[k] ? Inf : distance[k], 1:n)
        (done[u] || isinf(distance[u]) || u == t) && break
        done[u] = true
        for v in 1:n
            d = distance[u] + hypot((poly[u] .- poly[v])...)
            !done[v] && d < distance[v] - 1e-15 && visible(u, v) && (distance[v] = d; back[v] = u)
        end
    end
    path = [t]
    while path[1] != s
        pushfirst!(path, back[path[1]])
    end
    return path
end

"""
The pieces of the polygon `poly`, counterclockwise, whose side from
`poly[k]` has the cell `cells[k]`, or 0 for a side of the hull: each run of
sides of one cell takes the part between it and the shortest path inside
the polygon joining its ends. The rest, inside all these paths, goes to a
cell of several runs, the longest, so that their pieces join; otherwise it is
cut by spokes from its centroid, a side of the hull halved between the runs
beside it.
"""
function portrait_partition_face(poly, cells; tiny=1e-10)
    n = length(poly)
    starts = [k for k in 1:n if cells[k] != cells[mod1(k - 1, n)]]
    isempty(starts) && return [(cells[1], poly)]
    runs = [(s, starts[mod1(r + 1, length(starts))]) for (r, s) in enumerate(starts)]
    pieces = Tuple{Int,Vector{NTuple{2,Float64}}}[]
    paths = Vector{Int}[]
    run_length = Dict{Int,Float64}()
    for (s, e) in runs
        path = portrait_geodesic(poly, s, e)
        push!(paths, path)
        along = s < e ? collect(s:e) : [collect(s:n); collect(1:e)]
        run_length[cells[s]] = get(run_length, cells[s], 0.0) + sum(hypot((poly[i] .- poly[j])...) for (i, j) in zip(along, along[2:end]))
        piece = [poly[along]; poly[reverse(path[2:end-1])]]
        cells[s] != 0 && abs(portrait_area(piece)) > tiny && push!(pieces, (cells[s], piece))
    end
    middle = poly[reduce(vcat, (path[1:end-1] for path in paths))]
    length(middle) >= 3 && abs(portrait_area(middle)) > tiny || return pieces
    owner = reduce(vcat, (fill(r, length(path) - 1) for (r, path) in enumerate(paths)))
    run_cells = [cells[s] for (s, _) in runs]
    repeated = [c for c in unique(run_cells) if c != 0 && count(==(c), run_cells) >= 2]
    isempty(repeated) || return push!(pieces, (argmax(c -> run_length[c], repeated), middle))
    centre = portrait_centroid(middle)
    spokes = Tuple{Int,Vector{NTuple{2,Float64}}}[]
    for (j, r) in enumerate(owner)
        p, q = middle[j], middle[mod1(j + 1, length(middle))]
        if run_cells[r] != 0
            push!(spokes, (run_cells[r], [centre, p, q]))
        else
            half = (p .+ q) ./ 2
            push!(spokes, (run_cells[mod1(r - 1, length(runs))], [centre, p, half]),
                  (run_cells[mod1(r + 1, length(runs))], [centre, half, q]))
        end
    end
    # A middle that is not star-shaped from its centroid goes whole to its longest run.
    all(portrait_area(t) > -tiny for (_, t) in spokes) && return append!(pieces, spokes)
    return push!(pieces, (argmax(c -> c == 0 ? -1.0 : run_length[c], unique(run_cells)), middle))
end

"The corners of the convex hull of the points `pos[keys]`, counterclockwise."
function portrait_hull(pos, keys)
    order = sort(collect(keys); by=k -> pos[k])
    chain(points) = foldl(points; init=Int[]) do hull, k
        while length(hull) >= 2 && portrait_turn(pos[hull[end-1]], pos[hull[end]], pos[k]) <= 1e-12
            pop!(hull)
        end
        push!(hull, k)
    end
    return [chain(order)[1:end-1]; chain(reverse(order))[1:end-1]]
end

"""
The pieces of the regions of the convex hull of the planar drawing `pos` of
the core, its edges `ends` with the cells `colour`: those of
`portrait_partition_face` of its bounded faces, and of the pockets between
its outer face and its hull, each a side of the hull and the path of the
outer face below it.
"""
function portrait_region_pieces(pos, ends, colour)
    neighbours = [Int[] for _ in eachindex(pos)]
    foreach(((u, v),) -> (push!(neighbours[u], v); push!(neighbours[v], u)), ends)
    side(a, b) = colour[minmax(a, b)]
    faces = portrait_faces(pos, neighbours)
    outer = argmin(f -> portrait_area(pos[f]), faces)
    pieces = Tuple{Int,Vector{NTuple{2,Float64}}}[]
    for f in faces
        f === outer || append!(pieces, portrait_partition_face(pos[f], [side(a, b) for (a, b) in zip(f, circshift(f, -1))]))
    end
    O = reverse(outer)
    at = Dict(v => k for (k, v) in enumerate(O))
    hull = portrait_hull(pos, O)
    for (a, b) in zip(hull, circshift(hull, -1))
        path = at[a] < at[b] ? O[at[a]:at[b]] : [O[at[a]:end]; O[1:at[b]]]
        pocket = [path[1]; reverse(path[2:end])]
        length(path) > 2 && abs(portrait_area(pos[pocket])) > 1e-10 || continue
        append!(pieces, portrait_partition_face(pos[pocket], [0; [side(p, q) for (p, q) in zip(pocket[2:end], circshift(pocket, -1)[2:end])]]))
    end
    return pieces
end

"""
The points along the segments `ends` of the drawing `pos`, `spacing` apart,
the first and last `gap` from its ends, so that the border of the regions of
two segments at a vertex passes through it, and the segment of every point;
with `counts`, `counts[k]` points on the k-th segment, however long.
"""
function portrait_region_samples(pos, ends; spacing, gap, counts=nothing)
    δ = gap
    points, owner = NTuple{2,Float64}[], Int[]
    for (k, (u, v)) in enumerate(ends)
        a, b = pos[u], pos[v]
        L = hypot((b .- a)...)
        if isnothing(counts) && L < 2δ + 1e-9
            push!(points, (a .+ b) ./ 2)
            push!(owner, k)
            continue
        end
        m = isnothing(counts) ? max(1, round(Int, (L - 2δ) / spacing)) + 1 : counts[k]
        isnothing(counts) || (L = max(L, 3δ))
        for t in range(δ / L, 1 - δ / L; length=m)
            push!(points, a .+ t .* (b .- a))
            push!(owner, k)
        end
    end
    return points, owner
end

"The convex polygon `poly`, its corners with the tags of their sides, cut into `out` to the side of `p` of the bisector of `p` and `q`, a new side tagged `j`."
function portrait_clip!(out, poly, p, q, j)
    mid, d = (p .+ q) ./ 2, q .- p
    side(x) = (x[1] - mid[1]) * d[1] + (x[2] - mid[2]) * d[2]
    empty!(out)
    isempty(poly) && return out
    sb = side(poly[1][1])
    for i in eachindex(poly)
        (a, t), (b, _) = poly[i], poly[mod1(i + 1, length(poly))]
        sa, sb = sb, side(b)
        sa <= 0 && push!(out, (a, t))
        (sa <= 0) != (sb <= 0) && push!(out, (a .+ (sa / (sa - sb)) .* (b .- a), sa <= 0 ? j : t))
    end
    return out
end

"""
The Voronoi cells of `points` in the convex polygon `domain`,
counterclockwise: each cell its corners, each with the point beyond its side
from there, or 0 on the side of the domain. The points nearby are found ring
after ring of a grid of `bucket`.
"""
function portrait_voronoi(points, domain; bucket)
    lo = (minimum(first.(domain)), minimum(last.(domain)))
    nx, ny = (floor(Int, (maximum(p[k] for p in domain) - lo[k]) / bucket) + 1 for k in 1:2)
    at(p) = (clamp(floor(Int, (p[1] - lo[1]) / bucket) + 1, 1, nx), clamp(floor(Int, (p[2] - lo[2]) / bucket) + 1, 1, ny))
    grid = [Int[] for _ in 1:nx, _ in 1:ny]
    foreach(i -> push!(grid[at(points[i])...], i), eachindex(points))
    cells = Vector{Vector{Tuple{NTuple{2,Float64},Int}}}(undef, length(points))
    for (i, p) in enumerate(points)
        poly, spare = [(q, 0) for q in domain], Tuple{NTuple{2,Float64},Int}[]
        (gx, gy), r = at(p), 0
        while true
            for x in max(gx - r, 1):min(gx + r, nx), y in max(gy - r, 1):min(gy + r, ny), j in grid[x, y]
                max(abs(x - gx), abs(y - gy)) == r && j != i || continue
                # A point beyond twice the radius of the cell does not cut it, nor one at the same place.
                d2 = sum(abs2, points[j] .- p)
                d2 > (1e-9 * bucket)^2 && any(((c, _),) -> sum(abs2, c .- p) > d2 / 4, poly) || continue
                portrait_clip!(spare, poly, p, points[j], j)
                poly, spare = spare, poly
            end
            # The points beyond ring r are at least r buckets away.
            r * bucket >= 2maximum(((c, _),) -> hypot((c .- p)...), poly; init=0.0) || r > max(nx, ny) ? break : (r += 1)
        end
        cells[i] = poly
    end
    return cells
end

"""
The map of the regions of the drawing `pos` of the edges `ends` of the cells
`cell`: the Voronoi cells of the points of `portrait_region_samples`, each
with the cell of its edge, within the hull of `pos` `PORTRAIT_REGION_MARGIN`
times its largest side beyond it, with round corners. Returns the cells, the
cell of every point, the number of components of the region of every cell,
cells of one colour joined when they share a side, the length of the
borders between regions, as shares of the largest side of `pos`, and the
convexity and roundness of the region of every cell, its area over that of
its hull and 4π its area over its perimeter squared. Lengths are shares of
`side`, by default the largest side of `pos`; with `counts`, so many points
on every edge.
"""
function portrait_region_map(pos, ends, cell; spacing=PORTRAIT_REGION_SAMPLES.spacing, counts=nothing,
                             side=maximum(maximum(p[k] for p in pos) - minimum(p[k] for p in pos) for k in 1:2))
    points, owner = portrait_region_samples(pos, ends; spacing=spacing * side, gap=PORTRAIT_REGION_SAMPLES.ends * side, counts)
    hull = pos[portrait_hull(pos, eachindex(pos))]
    m = PORTRAIT_REGION_MARGIN * side
    normal(a, b) = (d = b .- a; (d[2], -d[1]) ./ hypot(d...))
    domain = NTuple{2,Float64}[]
    for (i, a) in enumerate(hull)
        θ1 = atan(reverse(normal(hull[mod1(i - 1, end)], a))...)
        θ2 = atan(reverse(normal(a, hull[mod1(i + 1, end)]))...)
        θ2 < θ1 && (θ2 += 2π)
        append!(domain, [a .+ m .* (cos(θ), sin(θ)) for θ in range(θ1, θ2; length=6)])
    end
    cells = portrait_voronoi(points, domain; bucket=2spacing * side)
    colour = cell[owner]
    parent = collect(eachindex(points))
    root(x) = parent[x] == x ? x : (parent[x] = root(parent[x]))
    colours = sort!(unique(cell))
    border, area, perimeter = 0.0, Dict(c => 0.0 for c in colours), Dict(c => 0.0 for c in colours)
    corners = Dict(c => NTuple{2,Float64}[] for c in colours)
    for (i, poly) in enumerate(cells)
        length(poly) >= 3 && (area[colour[i]] += portrait_area(first.(poly)))
        append!(corners[colour[i]], first.(poly))
        for (k, (a, j)) in enumerate(poly)
            ℓ = hypot((poly[mod1(k + 1, end)][1] .- a)...)
            (j == 0 || colour[j] != colour[i]) && (perimeter[colour[i]] += ℓ)
            j == 0 || ℓ < 1e-9 * side || (colour[j] == colour[i] ? (parent[root(i)] = root(j)) : (border += ℓ / 2))
        end
    end
    roots = Dict{Int,Set{Int}}()
    foreach(i -> push!(get!(roots, colour[i], Set{Int}()), root(i)), eachindex(points))
    convexity = Dict(c => area[c] / portrait_area(corners[c][portrait_hull(corners[c], eachindex(corners[c]))]) for c in colours)
    roundness = Dict(c => 4π * area[c] / perimeter[c]^2 for c in colours)
    return (; cells, colour, components=Dict(c => length(r) for (c, r) in roots), border=border / side, convexity, roundness)
end

"""
The paths through the segments `segments`, joined where their ends meet
within `tol`: with `directed`, each from the end of one to the start of the
next, every path closed; otherwise each to the end of the next either way,
a path ending where not two segments meet.
"""
function portrait_join(segments; directed, tol)
    key(p) = round.(Int, p ./ tol)
    at = Dict{NTuple{2,Int},Vector{Tuple{Int,Bool}}}()
    for (s, (a, b)) in enumerate(segments)
        push!(get!(at, key(a), []), (s, false))
        directed || push!(get!(at, key(b), []), (s, true))
    end
    meeting(p) = [h for x in -1:1, y in -1:1 for h in get(at, key(p) .+ (x, y), ()) if hypot((p .- (h[2] ? segments[h[1]][2] : segments[h[1]][1]))...) <= tol]
    used = falses(length(segments))
    paths = Vector{NTuple{2,Float64}}[]
    # Open paths start at an end where not two segments meet, closed ones anywhere.
    starts = directed ? collect(eachindex(segments)) :
             [[s for (s, (a, b)) in enumerate(segments) if length(meeting(a)) != 2 || length(meeting(b)) != 2]; eachindex(segments)]
    for s in starts
        used[s] && continue
        a, b = segments[s]
        !directed && length(meeting(b)) != 2 && length(meeting(a)) == 2 && ((a, b) = (b, a))
        path = [a, b]
        used[s] = true
        while true
            next = filter(h -> !used[h[1]], meeting(path[end]))
            (isempty(next) || (!directed && length(meeting(path[end])) != 2)) && break
            t, flipped = first(next)
            used[t] = true
            push!(path, flipped ? segments[t][1] : segments[t][2])
        end
        directed && hypot((path[end] .- path[1])...) <= tol && pop!(path)
        # Corners on a straight line are dropped.
        keep = [k for k in eachindex(path) if (!directed && k in (1, length(path))) ||
                abs(portrait_turn(path[mod1(k - 1, end)], path[k], path[mod1(k + 1, end)])) > 1e-12]
        push!(paths, path[keep])
    end
    return paths
end

"""
The components of the edges of every cell, of the edges `ends` on `1:n` of
the cells `cell`, as (vertex, cell) => a vertex of its component.
"""
function portrait_cell_components(n, ends, cell)
    component = Dict{Tuple{Int,Int},Int}()
    for c in unique(cell)
        parent = collect(1:n)
        root(x) = parent[x] == x ? x : (parent[x] = root(parent[x]))
        foreach(k -> cell[k] == c && (parent[root(ends[k][1])] = root(ends[k][2])), eachindex(ends))
        foreach(k -> cell[k] == c && foreach(v -> component[(v, c)] = root(v), ends[k]), eachindex(ends))
    end
    return component
end

"The crossings of the segments `ends` of the drawing `pos` without a common end, as (k, l, point)."
function portrait_all_crossings(pos, ends)
    crossings = Tuple{Int,Int,NTuple{2,Float64}}[]
    for k in eachindex(ends), l in k+1:length(ends)
        (a, b), (p, q) = ends[k], ends[l]
        allunique((a, b, p, q)) && portrait_crossing((pos[a], pos[b]), (pos[p], pos[q])) || continue
        x = portrait_meeting(pos[a], pos[b], pos[p], pos[q]; inside=0.0)
        push!(crossings, (k, l, something(x, (pos[a] .+ pos[b]) ./ 2)))
    end
    return crossings
end

"""
The energy of the drawing `pos` of the core with the edges `ends` of the
cells `cell`, its vertices with the `neighbours`, for --connected-regions, its
terms weighted by `PORTRAIT_REGION_WEIGHTS`: the components of regions of
`portrait_region_map` beyond one a cell, the crossings of edges of two cells
and of one, the stress towards the distances `D` with the weights `W`, the
vertices closer than a quarter of the shortest edge to an edge without them,
the bends of the vertices of degree 2, (1 + cos)/2 of their angle, the angles
at the others below 50 degrees, and the length of the borders. Returns it
and its terms.
"""
function portrait_region_energy(pos, ends, cell, D, W, neighbours)
    regions = portrait_region_map(pos, ends, cell)
    split = sum(values(regions.components)) - length(regions.components)
    crossings = portrait_all_crossings(pos, ends)
    mixed = count(((k, l, _),) -> cell[k] != cell[l], crossings)
    n = length(pos)
    stress = sum(W[i, j] * (hypot((pos[i] .- pos[j])...) - D[i, j])^2 for i in 1:n for j in i+1:n)
    reach = minimum(D[u, v] for (u, v) in ends) / 4
    clearance = sum(max(0.0, 1 - hypot((pos[p] .- curve_closest(pos[p], pos[u], pos[v]))...) / reach)^2
                    for (u, v) in ends for p in 1:n if p != u && p != v)
    bend, angle = 0.0, 0.0
    for v in 1:n
        if length(neighbours[v]) == 2
            a, b = (pos[w] .- pos[v] for w in neighbours[v])
            bend += (1 + sum(a .* b) / (hypot(a...) * hypot(b...))) / 2
        else
            θ = sort([atan(reverse(pos[w] .- pos[v])...) for w in neighbours[v]])
            angle += sum(max(0.0, 1 - g / deg2rad(50))^2 for g in [diff(θ); θ[1] + 2π - θ[end]])
        end
    end
    terms = (; split, mixed, within=length(crossings) - mixed, stress, clearance, bend, angle, border=regions.border)
    return sum(PORTRAIT_REGION_WEIGHTS[k] * terms[k] for k in keys(terms)), terms
end

"`X` relaxed by `sweeps` steps of stress majorization towards the distances `D` with the weights `W`."
function portrait_smacof(X, D, W; sweeps)
    n = size(X, 1)
    V = -copy(W)
    foreach(i -> V[i, i] = sum(W[i, :]) - W[i, i], 1:n)
    inverse = LinearAlgebra.pinv(V)
    for _ in 1:sweeps
        B = [i == j ? 0.0 : -W[i, j] * D[i, j] / max(hypot((X[i, :] .- X[j, :])...), 1e-12) for i in 1:n, j in 1:n]
        foreach(i -> B[i, i] = -sum(B[i, :]), 1:n)
        X = inverse * (B * X)
    end
    return X
end

"""
The drawing of the core of `curve` for its region atlas with
--connected-regions, in which every cell has one region, largest side 1:
edges may cross, best those of one cell. The distances are those of
`portrait_graph_distances` with the chain power `chain_power`, shortened
within a cell and lengthened between cells as in `portrait_core_drawing`, and
between vertices of two components of the edges of one cell at most
`PORTRAIT_REGION_PULL` times the square root of its edges. From `starts`
starts, the first the classical scaling of the distances, the others random
by their numbers, each relaxed by stress majorization and then annealed by
`moves` moves of single points, see `portrait_region_energy`: the drawing of
least energy whose regions stay connected with points `PORTRAIT_REGION_SAMPLES.final` apart.
"""
function portrait_region_drawing(curve; starts=PORTRAIT_REGION_STARTS, moves=PORTRAIT_REGION_MOVES, chain_power=PORTRAIT_CHAIN_POWER)
    (; n, ends, colour) = curve
    cell = [colour[minmax(e...)] for e in ends]
    neighbours = [Int[] for _ in 1:n]
    foreach(((u, v),) -> (push!(neighbours[u], v); push!(neighbours[v], u)), ends)
    cells = [Set(cell[k] for k in eachindex(ends) if v in ends[k]) for v in 1:n]
    component = portrait_cell_components(n, ends, cell)
    G = portrait_graph_distances(n, ends; power=chain_power)
    D = [i == j ? 0.0 : G[i, j] * (isempty(intersect(cells[i], cells[j])) ? PORTRAIT_OTHER_CELL : PORTRAIT_SAME_CELL)
         for i in 1:n, j in 1:n]
    for i in 1:n, j in 1:n, c in intersect(cells[i], cells[j])
        component[(i, c)] == component[(j, c)] || (D[i, j] = min(D[i, j], PORTRAIT_REGION_PULL * sqrt(count(==(c), cell))))
    end
    for k in 1:n, i in 1:n, j in 1:n
        D[i, j] = min(D[i, j], D[i, k] + D[k, j])
    end
    W = [i == j ? 0.0 : D[i, j]^-2 for i in 1:n, j in 1:n]
    energy(pos) = portrait_region_energy(pos, ends, cell, D, W, neighbours)
    drawings = Vector{Any}(undef, starts)
    Threads.@threads for s in 1:starts
        rng = Random.Xoshiro(s)
        if s == 1
            # The classical scaling: the two largest eigenvectors of the doubly centred squared distances.
            J = LinearAlgebra.I - fill(1 / n, n, n)
            F = LinearAlgebra.eigen(LinearAlgebra.Symmetric(-J * (D .^ 2) * J / 2))
            top = sortperm(F.values; rev=true)[1:2]
            X = F.vectors[:, top] .* sqrt.(max.(F.values[top], 0))'
        else
            X = randn(rng, n, 2)
        end
        X = portrait_smacof(X, D, W; sweeps=PORTRAIT_SWEEPS[1])
        pos = [(X[k, 1], X[k, 2]) for k in 1:n]
        current, _ = energy(pos)
        best, kept = current, copy(pos)
        for move in 1:moves
            # The temperature and the step fall linearly.
            t = move / moves
            k = rand(rng, 1:n)
            old = pos[k]
            pos[k] = old .+ (PORTRAIT_REGION_ANNEAL.step * (1 - t) + 0.02) .* (randn(rng), randn(rng))
            e, _ = energy(pos)
            if e <= current || rand(rng) < exp((current - e) / (PORTRAIT_REGION_ANNEAL.temperature * (1 - t) + 1e-3))
                current = e
                e < best && ((best, kept) = (e, copy(pos)))
            else
                pos[k] = old
            end
        end
        split = portrait_region_map(kept, ends, cell; spacing=PORTRAIT_REGION_SAMPLES.final).components
        drawings[s] = (split=sum(values(split)) - length(split), energy=best, terms=energy(kept)[2], pos=kept, start=s)
    end
    for d in drawings
        @printf("region start %d: energy %.1f, %d split regions at the final spacing, %s\n", d.start, d.energy, d.split,
                join(("$k $(round(v; digits=2))" for (k, v) in pairs(d.terms)), ", "))
    end
    # Nothing outside the loop above is named like its variables, which would make them shared by the threads.
    chosen = argmin(d -> (d.split, d.energy), drawings)
    println("Kept region start $(chosen.start)")
    lo = (minimum(first.(chosen.pos)), minimum(last.(chosen.pos)))
    side = maximum(maximum(p[k] for p in chosen.pos) - lo[k] for k in 1:2)
    return [(p .- lo) ./ side for p in chosen.pos]
end

"""
The terms of the energy of `portrait_region_smoothing` of the drawing `pos`
of the core with the edges `ends` of the cells `cell`, its vertices in the
cells `cells` with the `neighbours`, from the drawing `start`, its regions
with `counts` points on every edge: the bends of the vertices of degree 2,
(1 + cos)/2 of their angle, and the unevenness of their two edges, the angles
at the others below 90 degrees, the vertices closer than
`PORTRAIT_SMOOTHING_REACH` to an edge without them, a barrier when an edge
of theirs is of another cell, the squared distances from `start`, and how
far the regions of `portrait_region_map` are from convex and from round.
"""
function portrait_smoothing_terms(pos, start, ends, cell, cells, neighbours, counts)
    bend = even = angle = 0.0
    for v in eachindex(pos)
        if length(neighbours[v]) == 2
            a, b = (pos[w] .- pos[v] for w in neighbours[v])
            la, lb = hypot(a...), hypot(b...)
            bend += (1 + sum(a .* b) / (la * lb)) / 2
            even += ((la - lb) / (la + lb))^2
        else
            θ = sort([atan(reverse(pos[w] .- pos[v])...) for w in neighbours[v]])
            angle += sum(max(0.0, 1 - g / deg2rad(90))^2 for g in [diff(θ); θ[1] + 2π - θ[end]])
        end
    end
    clearance, reach = 0.0, PORTRAIT_SMOOTHING_REACH
    for (k, (u, v)) in enumerate(ends), p in eachindex(pos)
        (p == u || p == v) && continue
        d = hypot((pos[p] .- curve_closest(pos[p], pos[u], pos[v]))...)
        d < reach || continue
        clearance += cells[p] == Set([cell[k]]) ? (1 - d / reach)^2 : (reach / max(d, 1e-9) - 1)^2
    end
    anchor = sum(sum(abs2, p .- q) for (p, q) in zip(pos, start))
    regions = portrait_region_map(pos, ends, cell; counts, side=1.0)
    return (; bend, even, angle, clearance, anchor, convex=sum(1 - x for x in values(regions.convexity)),
            round=sum(1 - x for x in values(regions.roundness)))
end

"The gradient of `energy` at the points `pos` by central differences `ε` apart, on all threads."
function portrait_gradient(energy, pos, ε)
    gradient = Vector{NTuple{2,Float64}}(undef, length(pos))
    Threads.@threads for k in eachindex(pos)
        p = copy(pos)
        f(d) = (p[k] = pos[k] .+ d; e = energy(p); p[k] = pos[k]; e)
        gradient[k] = ((f((ε, 0.0)) - f((-ε, 0.0))) / 2ε, (f((0.0, ε)) - f((0.0, -ε))) / 2ε)
    end
    return gradient
end

"""
The drawing `pos` of the core of `curve` for --connected-regions smoothed by
`steps` steps of a flow, largest side 1: at every step every point moves
against the gradient of the terms of `portrait_smoothing_terms`, weighted by
`PORTRAIT_SMOOTHING_WEIGHTS`, by at most a step that falls linearly over
`PORTRAIT_SMOOTHING_STEP`, a point with less than the median force less. The
regions have as many points on every edge as at the start, so that they move
continuously. A step is halved, at most four times and then not taken, while
it makes edges of two cells cross, more edges cross than at the start, or a
region split, so that every vertex and edge stays in its region. Returns the
drawing of least energy, looked at every `PORTRAIT_SMOOTHING_LOOK` steps,
whose regions stay connected with points `PORTRAIT_REGION_SAMPLES.final` apart.
"""
function portrait_region_smoothing(curve, pos; steps)
    (; n, ends, colour) = curve
    cell = [colour[minmax(e...)] for e in ends]
    neighbours = [Int[] for _ in 1:n]
    foreach(((u, v),) -> (push!(neighbours[u], v); push!(neighbours[v], u)), ends)
    cells = [Set(cell[k] for k in eachindex(ends) if v in ends[k]) for v in 1:n]
    unit(drawing) = (lo = (minimum(first.(drawing)), minimum(last.(drawing)));
                     side = maximum(maximum(q[k] for q in drawing) - lo[k] for k in 1:2); [(q .- lo) ./ side for q in drawing])
    start = unit(pos)
    counts = [max(1, round(Int, (hypot((start[v] .- start[u])...) - 2PORTRAIT_REGION_SAMPLES.ends) / PORTRAIT_REGION_SAMPLES.spacing)) + 1
              for (u, v) in ends]
    within = length(portrait_all_crossings(start, ends))
    terms(drawing) = portrait_smoothing_terms(drawing, start, ends, cell, cells, neighbours, counts)
    energy(drawing) = (t = terms(drawing); sum(PORTRAIT_SMOOTHING_WEIGHTS[k] * t[k] for k in keys(PORTRAIT_SMOOTHING_WEIGHTS)))
    connected(drawing; spacing=PORTRAIT_REGION_SAMPLES.spacing) =
        all(==(1), values(portrait_region_map(drawing, ends, cell; spacing, counts=spacing == PORTRAIT_REGION_SAMPLES.spacing ? counts : nothing,
                                              side=1.0).components))
    function admissible(drawing)
        crossings = portrait_all_crossings(drawing, ends)
        return !any(((k, l, _),) -> cell[k] != cell[l], crossings) && length(crossings) <= within && connected(drawing)
    end
    report(drawing) = join(("$k $(round(v; digits=3))" for (k, v) in pairs(terms(drawing))), ", ")
    println("smoothing from: ", report(start))
    now, best = start, (energy(start), start)
    for step in 1:steps
        gradient = portrait_gradient(energy, now, PORTRAIT_SMOOTHING_DIFFERENCE)
        median = sort([hypot(x...) for x in gradient])[end ÷ 2]
        σ = PORTRAIT_SMOOTHING_STEP[1] + (PORTRAIT_SMOOTHING_STEP[2] - PORTRAIT_SMOOTHING_STEP[1]) * step / steps
        for _ in 1:5
            trial = [q .- (σ / max(hypot(x...), median)) .* x for (q, x) in zip(now, gradient)]
            admissible(trial) && ((now = trial); break)
            σ /= 2
        end
        if step % PORTRAIT_SMOOTHING_LOOK == 0 || step == steps
            current = energy(now)
            current < best[1] && connected(now; spacing=PORTRAIT_REGION_SAMPLES.final) && (best = (current, copy(now)))
        end
    end
    println("smoothed to: ", report(best[2]))
    return unit(best[2])
end

"""
The data of a curve: its reduced graph `path`, the core by `curve_core`
and its edges `ends` on `1:n` with their cells `colour`, and the branches
hanging from every vertex; `core_ends` keeps the ends of the edges of the
core when `portrait_crossings` adds a crossing.
"""
function portrait_curve(path)
    edges = read_curve_graph(path)
    core = curve_core(edges)
    vertices = sort!(unique(v for e in core for v in (e.tail, e.head)))
    ends = curve_ends(core, Dict(v => k for (k, v) in enumerate(vertices)))
    colour = Dict(minmax(u, v) => e.cell for (e, (u, v)) in zip(core, ends))
    on_core = Set(vertices)
    around, legs = Dict{Int,Vector{Int}}(), Dict{Int,Int}()
    for e in edges
        e.head == 0 && (legs[e.tail] = get(legs, e.tail, 0) + 1; continue)
        e.tail in on_core && e.head in on_core && continue
        push!(get!(around, e.tail, Int[]), e.head)
        push!(get!(around, e.head, Int[]), e.tail)
    end
    cells = sort!(unique(e.cell for e in core); by=c -> (only(unique(e.number for e in core if e.cell == c)), c))
    return (; edges, vertices, ends, colour, around, legs, cells, n=length(vertices),
            number=Dict(e.cell => e.number for e in core), core, core_ends=ends, crossing=nothing)
end

"Whether the core of `curve` is planar."
function portrait_planar(curve)
    branch, chains = curve_core_chains(curve.ends, Set(eachindex(curve.ends)))
    index = Dict(p => i for (i, p) in enumerate(branch))
    return !isempty(curve_embeddings(length(branch), [(index[first(points)], index[last(points)]) for (points, _) in chains]; stop=true))
end

"""
The ways to draw the core of `curve`, not planar, with one crossing: for
any two of its chains between points of degree 3, see `curve_core_chains`,
without a common end, whose crossing makes it planar, `curve` with a new
point `n + 1`, of the vertex 0, where an edge of the first chain crosses one
of the second, the middle ones of their edges in the cell they share most
if they share one, else their middle ones, among those not in `avoid` if
there are any. Its `crossing` holds the point and the ends of the two edges,
the first drawn over the second.
"""
function portrait_crossings(curve; avoid=())
    (; ends, colour, n) = curve
    branch, chains = curve_core_chains(ends, Set(eachindex(ends)))
    index = Dict(p => i for (i, p) in enumerate(branch))
    pairs = [(index[first(points)], index[last(points)]) for (points, _) in chains]
    x = length(branch) + 1
    cell(k) = colour[minmax(ends[k]...)]
    crossed = []
    for i in eachindex(chains), j in i+1:length(chains)
        (a, b), (p, q) = pairs[i], pairs[j]
        isempty(intersect((a, b), (p, q))) || continue
        # The darts from the crossing alternate, as `curve_embeddings` wants.
        isempty(curve_embeddings(x, [[(x, a), (x, p), (x, b), (x, q)]; pairs[setdiff(eachindex(pairs), (i, j))]]; crossings=(x,), stop=true)) &&
            continue
        first_path, second_path = chains[i][2], chains[j][2]
        shared = intersect(cell.(first_path), cell.(second_path))
        best = isempty(shared) ? 0 : argmax(c -> min(count(k -> cell(k) == c, first_path), count(k -> cell(k) == c, second_path)), shared)
        function middle(path)
            ks = best == 0 ? path : filter(k -> cell(k) == best, path)
            free = filter(k -> !(curve.core[k].edge in avoid), ks)
            ks = isempty(free) ? ks : free
            return ks[cld(length(ks), 2)]
        end
        k, l = middle(first_path), middle(second_path)
        y = n + 1
        split = [(ends[k][1], y), (y, ends[k][2]), (ends[l][1], y), (y, ends[l][2])]
        new_ends = [ends[setdiff(eachindex(ends), (k, l))]; split]
        new_colour = merge(colour, Dict(minmax(e...) => cell(m) for (e, m) in zip(split, (k, k, l, l))))
        push!(crossed, merge(curve, (; ends=new_ends, colour=new_colour, n=y, vertices=[curve.vertices; 0],
                                     crossing=(x=y, over=ends[k], under=ends[l]))))
    end
    return crossed
end

# A branch of a fan is a vertex off the core or a leg, -1, hanging from `v`.
portrait_children(curve, v, parent) = [[w for w in get(curve.around, v, Int[]) if w != parent]; fill(-1, get(curve.legs, v, 0))]
portrait_leaves(curve, b, parent) = b == -1 ? 1 : sum(c -> portrait_leaves(curve, c, b), portrait_children(curve, b, parent); init=0)
portrait_height(curve, b, parent) =
    b == -1 ? 0 : 1 + maximum((portrait_height(curve, c, b) for c in portrait_children(curve, b, parent)); init=0)

"""
The segments of the fan of the branch `b` of `v` from `p`, as (from, to,
level, is a leg, vertex, branch): its first edge of length `ℓ` along the middle of the
sector `lo..hi`, every later one `PORTRAIT_LEVEL` times its parent, the
children sharing the sector of their parent equally, widened or narrowed as
`PORTRAIT_SPREAD` says, the biggest in the middle and, of two, the bigger on
the left at odd levels and on the right at even ones, so that a stem with a
leg at every vertex zigzags.
"""
function portrait_fan(curve, b, v, p, lo, hi, ℓ, level=1)
    φ = (lo + hi) / 2
    q = p .+ ℓ .* (cos(φ), sin(φ))
    out = [(p, q, level, b == -1, v, b)]
    b == -1 && return out
    children = portrait_children(curve, b, v)
    leaves(c) = portrait_leaves(curve, c, b)
    width = clamp(hi - lo, deg2rad(PORTRAIT_SPREAD[1]) * length(children), deg2rad(PORTRAIT_SPREAD[2]))
    start = φ - width / 2
    arranged = Int[]
    foreach(((i, c),) -> isodd(i) ? push!(arranged, c) : pushfirst!(arranged, c), enumerate(sort(children; by=c -> -leaves(c))))
    iseven(level) && reverse!(arranged)
    for c in arranged
        w = width / length(children)
        append!(out, portrait_fan(curve, c, b, q, start, start + w, ℓ * PORTRAIT_LEVEL * (c == -1 ? PORTRAIT_LEG : 1.0), level + 1))
        start += w
    end
    return out
end

"""
The fans of the pending trees and legs of `curve` around its core drawn at
`pos`: at every vertex of the core its branches, the tallest first, each
into the gap between core edges of largest room, its angle times the free
ray along its middle, shared as `PORTRAIT_GAP_SHARE` says; the branches of
one gap share it by their numbers of legs. Every fan starts with a first
edge of `scale` times `PORTRAIT_FIRST_EDGE`, core edges of median length 1, and the fan
of the larger first edge of two that cross, or a fan crossing the core,
shrinks by 5% until none do; with `repel`, they then push each other apart
as `portrait_repel!` says. With `overlap`, fans may cross each other, not the
core. Returns the core so scaled, the segments, and the smallest first edge.
"""
function portrait_fans(curve, pos; repel=false, scale=1.0, overlap=false)
    (; n, ends, vertices) = curve
    median = sort([hypot((pos[u] .- pos[v])...) for (u, v) in ends])[cld(length(ends), 2)]
    pos = [p ./ median for p in pos]
    neighbours = [Int[] for _ in 1:n]
    foreach(((u, v),) -> (push!(neighbours[u], v); push!(neighbours[v], u)), ends)
    segments = [(u, v, pos[u], pos[v]) for (u, v) in ends]
    fans = []
    for k in 1:n
        p, v = pos[k], vertices[k]
        dir(w) = atan(pos[w][2] - p[2], pos[w][1] - p[1])
        ws = sort(neighbours[k]; by=dir)
        gaps = [(dir(ws[i]), i == length(ws) ? dir(ws[1]) + 2π : dir(ws[i+1])) for i in eachindex(ws)]
        room = [(hi - lo) * min(curve_ray(p, (cos((lo + hi) / 2), sin((lo + hi) / 2)), segments, k), 3.0) for (lo, hi) in gaps]
        held = [Int[] for _ in gaps]
        for b in sort(portrait_children(curve, v, 0); by=b -> (-portrait_height(curve, b, v), -portrait_leaves(curve, b, v)))
            push!(held[argmax(g -> room[g] / (length(held[g]) + 1)^PORTRAIT_GAP_SHARE, eachindex(gaps))], b)
        end
        for ((lo, hi), bs) in zip(gaps, held)
            isempty(bs) && continue
            margin = 0.2 * (hi - lo) / length(bs)
            lo, hi = lo + margin, hi - margin
            widths = [portrait_leaves(curve, b, v) + 1.0 for b in bs]
            start = lo
            for (b, w) in zip(bs, widths)
                sector = (hi - lo) * w / sum(widths)
                half = min(sector / 2 * (length(bs) == 1 ? 0.8 : 0.9), deg2rad(PORTRAIT_FAN_HALF))
                push!(fans, (k=k, v=v, b=b, φ=start + sector / 2, half=half))
                start += sector
            end
        end
    end
    first_edge = fill(scale * PORTRAIT_FIRST_EDGE, length(fans))
    drawn(a) = (f = fans[a]; portrait_fan(curve, f.b, f.v, pos[f.k], f.φ - f.half, f.φ + f.half,
                                          first_edge[a] * (f.b == -1 ? PORTRAIT_LEG : 1.0)))
    core = [(pos[u], pos[v]) for (u, v) in ends]
    for round in 1:401
        round == 401 && error("the fans still cross after 400 rounds")
        lines = [[(p, q) for (p, q, _, _) in drawn(a)] for a in eachindex(fans)]
        hits = Set(a for a in eachindex(fans) if any(portrait_crossing(s, t) for s in lines[a] for t in core))
        for a in eachindex(fans), c in a+1:length(fans)
            overlap && break
            any(portrait_crossing(s, t) for s in lines[a] for t in lines[c]) && push!(hits, first_edge[a] >= first_edge[c] ? a : c)
        end
        isempty(hits) && break
        first_edge[collect(hits)] .*= 0.95
    end
    repel && portrait_repel!(curve, pos, fans, first_edge; scale, overlap)
    return pos, reduce(vcat, drawn.(eachindex(fans))), minimum(first_edge)
end

"""
The fans `fans` with first edges `first_edge` around the core at `pos`,
none crossing another or the core, pushing each other apart: one fan after
another turns, opens, closes, grows or shrinks by a step, whichever lowers
its energy most without a crossing, the steps getting smaller sweep after
sweep. The energy of a fan is, for each of its points and segment of another
fan or of the core, and each point of another fan or of the core and segment
of it, (r / d - 1)^2 where their distance d is below r = `scale` times
`PORTRAIT_REPEL.reach`, less `PORTRAIT_REPEL.reward` times the log of its
first edge, which stays at most `scale` times `PORTRAIT_FIRST_EDGE`. With
`overlap`, a fan may cross another, not the core.
"""
function portrait_repel!(curve, pos, fans, first_edge; scale=1.0, overlap=false)
    (; reach, reward, sweeps, turn) = PORTRAIT_REPEL
    reach *= scale
    drawn(f, ℓ) = [(p, q) for (p, q, _, _) in portrait_fan(curve, f.b, f.v, pos[f.k], f.φ - f.half, f.φ + f.half,
                                                              ℓ * (f.b == -1 ? PORTRAIT_LEG : 1.0))]
    lines = [drawn(f, ℓ) for (f, ℓ) in zip(fans, first_edge)]
    core = [(pos[u], pos[v]) for (u, v) in curve.ends]
    near(p, (a, b)) = (d = hypot((p .- curve_closest(p, a, b))...); d < reach ? (reach / max(d, 1e-9) - 1)^2 : 0.0)
    function energy(a, L)
        e = -reward * log(hypot((L[1][2] .- L[1][1])...))
        for (c, M) in enumerate(lines)
            c == a && continue
            e += sum(near(q, s) for (_, q) in L for s in M) + sum(near(q, s) for (_, q) in M for s in L)
        end
        e += sum(near(q, t) for (_, q) in L for t in core)
        return e + sum(near(p, s) for (k, p) in enumerate(pos) if k != fans[a].k for s in L)
    end
    clear(a, L) = (overlap || !any(portrait_crossing(s, t) for (c, M) in enumerate(lines) if c != a for t in M for s in L)) &&
                  !any(portrait_crossing(s, t) for s in L for t in core)
    for sweep in 1:sweeps
        step = 1 - 0.85 * (sweep - 1) / (sweeps - 1)
        moved = false
        for (a, f) in enumerate(fans)
            best, choice = energy(a, lines[a]) - 1e-9, nothing
            for (dφ, dh, grow) in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1), (1, 0, 1), (-1, 0, 1), (0, -1, 1))
                g = (; f..., φ=f.φ + dφ * step * deg2rad(turn),
                     half=clamp(f.half + dh * step * deg2rad(turn), deg2rad(12), deg2rad(PORTRAIT_FAN_HALF)))
                ℓ = min(first_edge[a] * 1.06^(grow * step), scale * PORTRAIT_FIRST_EDGE)
                L = drawn(g, ℓ)
                e = energy(a, L)
                e < best && clear(a, L) && ((best, choice) = (e, (g, ℓ, L)))
            end
            isnothing(choice) && continue
            fans[a], first_edge[a], lines[a] = choice
            moved = true
        end
        moved || break
    end
    return fans
end

portrait_point(p) = @sprintf("(%.4f,%.4f)", p...)

"The middle of the box of the points `pos`."
portrait_middle(pos) = Tuple((minimum(p[k] for p in pos) + maximum(p[k] for p in pos)) / 2 for k in 1:2)

"""
The TeX line that makes the box of a picture, whatever it draws, run from
`left` to `right` cm across from `middle` and from `down` cm below it to `up`
cm above.
"""
portrait_box(middle, (left, right, down, up)) = "\\pgfresetboundingbox\\path[use as bounding box] " *
                                                portrait_point(middle .+ (left, -down)) * " rectangle " * portrait_point(middle .+ (right, up)) * ";"

"How far the points `ink` reach from `middle`, `margin` beyond them: leftwards, as a negative, rightwards, downwards and upwards."
portrait_reach(ink, middle, margin) = (minimum(p[1] - middle[1] for p in ink) - margin, maximum(p[1] - middle[1] for p in ink) + margin,
                                       maximum(middle[2] - p[2] for p in ink) + margin, maximum(p[2] - middle[2] for p in ink) + margin)

"""
How far the Berkovich picture of `curve` of the layout `layout`, with the red
points `points` if given, reaches from the middle of the box of its core, in
cm, as `portrait_reach` says.
"""
function portrait_berkovich_reach(curve, layout, points)
    (; pos, segments, unit) = layout
    ink = [first.(segments); [s[2] for s in segments]; pos; something(points, NTuple{2,Float64}[])]
    return portrait_reach([p .* unit for p in ink], portrait_middle(pos) .* unit, PORTRAIT_DIVISOR_PT * PORTRAIT_PT)
end

"""
The layout of the Berkovich picture of `curve` with its core drawn at `pos`:
the core scaled by `portrait_fans` and its fans pushed apart, the lowest
corner of their box, the cm of a unit of the scaled core in a picture of
`PORTRAIT_BERKOVICH_SIZE`, the cm of a unit of `pos`, and the scale of the
frames of the region atlas and its height. With `side_by_side`, the unit is
the largest for which the picture and the region atlas beside it, its core at
the same scale, fit in the row of `PORTRAIT_SIDE_BY_SIDE`, or with `atlas`,
its core that many times as large, fit in it side by side in boxes of one
height in which the cores stand at one height; with `frames`, the atlas keeps its frames, at most of
`PORTRAIT_FRAME`, its two columns as high as the picture. With `overlap`, the
fans may cross each other.
"""
function portrait_berkovich_layout(curve, pos; side_by_side=false, tree_scale=1.0, frames=false, overlap=false, atlas=nothing)
    scaled, segments, _ = portrait_fans(curve, pos; repel=true, scale=tree_scale, overlap)
    points = [first.(segments); [s[2] for s in segments]; scaled]
    lo = (minimum(first.(points)), minimum(last.(points)))
    extent(ps, k) = maximum(p[k] for p in ps) - minimum(p[k] for p in ps)
    frame, atlas_height = PORTRAIT_FRAME[1], PORTRAIT_ATLAS_HEIGHT
    if side_by_side && frames
        (; width, height, gap, frame_gap) = PORTRAIT_SIDE_BY_SIDE
        # Both cores and the two columns of frames fill the width, and the k frames of a column the height of the picture.
        k, across, high = cld(length(curve.cells), 2), extent(points, 1) + extent(scaled, 1), extent(points, 2)
        room = width - gap - 0.2 - 2PORTRAIT_FRAME_GAP
        unit = (room + 2(k - 1) * frame_gap / k) / (across + 2high / k)
        frame = (unit * high - (k - 1) * frame_gap) / k
        frame > PORTRAIT_FRAME[1] && ((frame, unit) = (PORTRAIT_FRAME[1], (room - 2PORTRAIT_FRAME[1]) / across))
        unit = min(unit, height / high, (height - 0.6) / extent(scaled, 2))
        frame = min(frame, (unit * high - (k - 1) * frame_gap) / k)
        atlas_height = unit * high
    elseif side_by_side && !isnothing(atlas)
        (; width, gap, aligned) = PORTRAIT_SIDE_BY_SIDE
        # Side by side, the atlas 0.6 cm wider than its core and its labels beside it; up and down, boxes of one height
        # from the lowest to the highest point of the two from the middles of the cores, its dots 0.2 cm more.
        middle = portrait_middle(scaled)
        down, up = (max(maximum(s * (p[2] - middle[2]) for p in points), atlas * extent(scaled, 2) / 2) for s in (-1, 1))
        unit = min((width - gap - 0.6 - PORTRAIT_LABEL_BESIDE[2]) / (extent(points, 1) + atlas * extent(scaled, 1)), (aligned - 0.2) / (down + up))
    elseif side_by_side
        (; width, height, gap) = PORTRAIT_SIDE_BY_SIDE
        # The atlas takes 0.6 cm beyond its core both ways, see `portrait_regions_tex`, and the borders 0.2 cm.
        unit = min((width - gap - 0.8) / (extent(points, 1) + extent(scaled, 1)), height / extent(points, 2),
                   (height - 0.6) / extent(scaled, 2))
    else
        W, H = PORTRAIT_BERKOVICH_SIZE
        unit = min(W / extent(points, 1), H / extent(points, 2))
    end
    u, v = first(curve.ends)
    return (pos=scaled, segments, lo, unit, core_unit=unit * hypot((scaled[u] .- scaled[v])...) / hypot((pos[u] .- pos[v])...),
            frame_scale=frame / PORTRAIT_FRAME[1], atlas_height)
end

"""
The places in the layout `layout` of `portrait_berkovich_layout` of the
points `places`, each an edge or leg of the graph of `curve` and its lattice
length from the tail, and the length of the edge (nothing for a leg): along
its segment of the core or of a fan, a point on a leg at its middle. The
points of a segment keep their order but are moved apart, and away from its
ends and a crossing on it, until their dots and those of the ends do not
touch, see `portrait_spread`.
"""
function portrait_points(curve, layout, places)
    (; pos, segments, unit) = layout
    on_core = Dict(e.edge => k for (k, e) in enumerate(curve.core))
    # The segments of the legs at every vertex, taken by its legs in their order.
    leg_segments = Dict{Int,Vector{Int}}()
    foreach(s -> segments[s][4] && push!(get!(leg_segments, segments[s][5], Int[]), s), eachindex(segments))
    leg_number = Dict{Int,Int}()
    for e in curve.edges
        e.head == 0 && (leg_number[e.edge] = count(f -> f.head == 0 && f.tail == e.tail && f.edge <= e.edge, curve.edges))
    end
    edge_of = Dict(e.edge => e for e in curve.edges)
    # Every point as its segment, from the tail of its edge, and its share of it.
    on = map(places) do (edge, from_tail, len)
        e = edge_of[edge]
        if haskey(on_core, edge)
            u, v = curve.core_ends[on_core[edge]]
            return (pos[u], pos[v], from_tail / len)
        elseif e.head == 0
            p, q = segments[leg_segments[e.tail][leg_number[edge]]]
            return (p, q, 0.5)
        end
        s = only(s for s in segments if !s[4] && minmax(s[5], s[6]) == minmax(e.tail, e.head))
        # The segment runs from its vertex, the tail of the edge or its head.
        return s[5] == e.tail ? (s[1], s[2], from_tail / len) : (s[2], s[1], from_tail / len)
    end
    shares = last.(on)
    # A crossing on a segment parts it like an end.
    crossings = [isnothing(curve.crossing) ? [] : [pos[curve.crossing.x]]; last.(portrait_all_crossings(pos, curve.ends))]
    for (p, q) in unique((p, q) for (p, q, _) in on)
        long = hypot((q .- p)...) * unit / PORTRAIT_PT
        cuts = [0.0; sort([hypot((x .- p)...) / hypot((q .- p)...) for x in crossings if portrait_on_segment(x, p, q; tol=1e-6)]); 1.0]
        for (a, b) in zip(cuts, cuts[2:end])
            mine = findall(k -> on[k][1] == p && on[k][2] == q && a <= shares[k] < b + (b == 1), eachindex(on))
            isempty(mine) && continue
            part = long * (b - a)
            margin(cut) = (cut in (0, 1) ? PORTRAIT_DIVISOR_PT + PORTRAIT_DOT_PT + 0.4 : PORTRAIT_CROSSING_CLEARANCE_PT) / part
            shares[mine] = a .+ (b - a) .* portrait_spread((shares[mine] .- a) ./ (b - a), (margin(a), margin(b)),
                                                           (2PORTRAIT_DIVISOR_PT + 0.4) / part)
        end
    end
    return [p .+ f .* (q .- p) for ((p, q, _), f) in zip(on, shares)]
end

"""
The shares `shares` of a segment moved apart in their order, as little as a
pass each way does it, until they are `gap` apart and `margins` from its
ends; evenly spread when the segment is too short.
"""
function portrait_spread(shares, margins, gap)
    order = sortperm(shares)
    f = shares[order]
    k = length(f)
    if sum(margins) + (k - 1) * gap > 1
        f = collect(range(1 / (k + 1), k / (k + 1); length=k))
    else
        f[1] = max(f[1], margins[1])
        foreach(i -> f[i] = max(f[i], f[i-1] + gap), 2:k)
        f[k] = min(f[k], 1 - margins[2])
        foreach(i -> f[i] = min(f[i], f[i+1] - gap), k-1:-1:1)
    end
    spread = similar(f)
    spread[order] = f
    return spread
end

"""
The TeX lines that draw the edge `over` of a crossing at `x` over the other,
`under`, between the points `pos`, `at` giving their TeX points and `unit`
the cm of a unit: a border along it at the crossing in the colour `border`,
that of the background, `gap` pt on either side, then the edge in `style`.
"""
function portrait_over(pos, x, over, under, at, unit, width, style; border="white", gap=CURVE_BORDER_PT)
    d, e = pos[over[2]] .- pos[over[1]], pos[under[2]] .- pos[under[1]]
    d = d ./ hypot(d...)
    sine = abs(d[1] * e[2] - d[2] * e[1]) / hypot(e...)
    half = (width / 2 + 2gap) * PORTRAIT_PT / unit / max(sine, 0.3)
    return string(@sprintf("\\draw[%s,line width=%.2fpt,line cap=butt] %s -- %s;\n", border, width + 2gap,
                           at(pos[x] .- half .* d), at(pos[x] .+ half .* d)),
                  "\\draw[$style,line width=$(width)pt,line cap=round] ", at(pos[over[1]]), " -- ", at(pos[x]), " -- ",
                  at(pos[over[2]]), ";\n")
end

"""
The TeX source of the Berkovich picture of `curve` of the layout `layout` of
`portrait_berkovich_layout`, see the header, with the red points `points` of
`portrait_points` if given; with `trivalent`, no dots at the vertices of
degree 2 of the core; with `box`, its box reaches so far from the middle
of the box of the core, see `portrait_box`.
"""
function portrait_berkovich_tex(curve, layout, name; points=nothing, trivalent=false, box=nothing)
    (; pos, segments, lo, unit) = layout
    at(p) = portrait_point((p .- lo) .* unit)
    return sprint() do io
        println(io, "% $name: the curve, its core in bold, ", isnothing(points) ? "the vertices of degree 3 of the core" :
                "the $(length(points)) points of its divisor", " in red, every pending tree a fan from it. Written by curve_portraits.jl.")
        println(io, "\\documentclass[tikz,border=2pt]{standalone}\n\\begin{document}")
        println(io, "\\definecolor{g2DivisorRed}{HTML}{$PORTRAIT_DIVISOR_COLOUR}\n\\begin{tikzpicture}")
        for (p, q, level, leg) in segments
            width = max(PORTRAIT_EDGE_PT * 0.92^(level - 1), 0.3)
            if leg
                # A leg tapers from the width of its edge to a point.
                d = q .- p
                ν = (-d[2], d[1]) .* (width * PORTRAIT_PT / unit / 2 / hypot(d...))
                println(io, "\\fill[black!60] ", at(p .+ ν), " -- ", at(q), " -- ", at(p .- ν), " -- cycle;")
            else
                @printf(io, "\\draw[black!85,line width=%.2fpt,line cap=round] %s -- %s;\n", width, at(p), at(q))
            end
        end
        for (u, v) in curve.ends
            println(io, "\\draw[black,line width=$(PORTRAIT_CORE_PT)pt,line cap=round] ", at(pos[u]), " -- ", at(pos[v]), ";")
        end
        crossing = curve.crossing
        isnothing(crossing) || print(io, portrait_over(pos, crossing.x, crossing.over, crossing.under, at, unit, PORTRAIT_CORE_PT, "black"))
        for (k, l, x) in portrait_all_crossings(pos, curve.ends)
            # The other crossings, the first edge over the second as in the region atlas, drawn again only near them.
            (a, b) = curve.ends[k]
            d = pos[b] .- pos[a]
            near = clamp.((hypot((x .- pos[a])...) .+ (-1, 1) .* (0.4 / unit)) ./ hypot(d...), 0, 1)
            print(io, portrait_over([pos; [pos[a] .+ near[1] .* d, x, pos[a] .+ near[2] .* d]], curve.n + 2, (curve.n + 1, curve.n + 3),
                                    curve.ends[l], at, unit, PORTRAIT_CORE_PT, "black"))
        end
        degree = zeros(Int, curve.n)
        foreach(((u, v),) -> (degree[u] += 1; degree[v] += 1), curve.ends)
        vertex(k) = (isnothing(crossing) || k != crossing.x) && !(trivalent && degree[k] == 2)
        red = isnothing(points) ? findall(==(3), degree) : Int[]
        foreach(k -> println(io, "\\fill[black] ", at(pos[k]), " circle ($(PORTRAIT_DOT_PT)pt);"), filter(k -> vertex(k) && !(k in red), 1:curve.n))
        foreach(p -> println(io, "\\fill[g2DivisorRed] ", at(p), " circle ($(PORTRAIT_DIVISOR_PT)pt);"), isnothing(points) ? pos[red] : points)
        isnothing(box) || println(io, portrait_box((portrait_middle(pos) .- lo) .* unit, box))
        println(io, "\\end{tikzpicture}\n\\end{document}")
    end
end

"""
The TeX source of the region atlas of `curve` with its core drawn at `pos`,
`unit` cm for a unit of `pos`, as in the Berkovich picture:
the regions of `portrait_region_pieces`, `PORTRAIT_MARGIN_PT` beyond the hull
with round corners, white lines between regions, the edges and vertices in
the colours of their cells, and the trees of the cells, without labels, in
frames in two columns with their types as `\\celltype` of main.tex writes
them, as near the heights of their regions as they fit, each with a line to
its region. With `labels`, the types instead, in the colours of their cells,
in every region with room for one, see `portrait_label_places`, and the
picture 0.6 cm larger than the core both ways. The frames and their trees are
`frame_scale` times as large, and the columns at least `atlas_height` high.
With `trivalent`, only the vertices of degree 3 are drawn.
With `connected`, for a drawing of `portrait_region_drawing`, the regions are
those of `portrait_region_map` instead, and the edges that cross, of one cell
or of two, pass one over the other. An edge over another has a thin border
in the colour of the region below the crossing. With `box`, the box of the
picture reaches so far from the middle of the box of the core, see
`portrait_box`; with `measure`, returns instead how far the picture reaches
from there, see `portrait_reach`, a label taken 0.5 cm across and 0.3 cm high.
"""
function portrait_regions_tex(curve, pos, name; unit, labels=false, frame_scale=1.0, atlas_height=PORTRAIT_ATLAS_HEIGHT,
                              trivalent=false, connected=false, box=nothing, measure=false)
    (; n, ends, colour, cells, crossing) = curve
    fw, fh = frame_scale .* PORTRAIT_FRAME
    width, height = maximum(first.(pos)), maximum(last.(pos))
    H = labels ? unit * height + 0.6 : max(atlas_height, unit * height + 0.6)
    W = unit * width + (labels ? 0.6 : 2 * (fw + PORTRAIT_FRAME_GAP))
    offset = ((W - unit * width) / 2, (H - unit * height) / 2)
    at(p) = offset .+ unit .* p
    place = at.(pos)
    if connected
        regions = portrait_region_map(pos, ends, [colour[minmax(e...)] for e in ends]; spacing=PORTRAIT_REGION_SAMPLES.final)
        split = filter(((_, k),) -> k > 1, regions.components)
        isempty(split) || println("Split regions of the atlas: ", split)
        pieces = [(c, at.(first.(poly))) for (c, poly) in zip(regions.colour, regions.cells)]
        # The sides of the cells between regions and on the margin, counterclockwise around their regions.
        outline = Dict{Int,Vector{NTuple{2,NTuple{2,Float64}}}}()
        walls = NTuple{2,NTuple{2,Float64}}[]
        for (i, poly) in enumerate(regions.cells), (k, (a, j)) in enumerate(poly)
            j == 0 || regions.colour[j] != regions.colour[i] || continue
            side = (at(a), at(poly[mod1(k + 1, end)][1]))
            hypot((side[2] .- side[1])...) > 1e-9 || continue
            push!(get!(outline, regions.colour[i], []), side)
            j > i && push!(walls, side)
        end
    else
        pieces = [(c, at.(poly)) for (c, poly) in portrait_region_pieces(pos, ends, colour)]
    end
    key(p) = round.(p; digits=5)
    neighbours = [Int[] for _ in 1:n]
    foreach(((u, v),) -> (push!(neighbours[u], v); push!(neighbours[v], u)), ends)
    palette = Dict(c => PORTRAIT_REGION_COLOURS[mod1(j, end)] for (j, c) in enumerate(cells))
    colours = connected ? [PORTRAIT_REGION_COLOURS; PORTRAIT_MORE_COLOURS] : PORTRAIT_REGION_COLOURS
    distinct = connected && length(cells) <= length(colours)
    if distinct || length(cells) > length(colours)
        # The cells every cell meets, at a corner of their regions or at a vertex of the core, by its regions or its edges.
        meeting = Dict{NTuple{2,Float64},Set{Int}}()
        foreach(((c, poly),) -> foreach(p -> push!(get!(meeting, key(p), Set{Int}()), c), poly), pieces)
        foreach(v -> union!(get!(meeting, key(place[v]), Set{Int}()), (colour[minmax(v, w)] for w in neighbours[v])), 1:n)
        meets = Dict(c => union(Set{Int}(), (s for s in values(meeting) if c in s)...) for c in cells)
        area = Dict(c => sum((abs(portrait_area(poly)) for (d, poly) in pieces if d == c); init=0.0) for c in cells)
        tint = [portrait_lab(h; tint=PORTRAIT_REGION_TINT) for h in colours]
        taken = Dict{Int,Int}()
        for c in sort(cells; by=c -> -area[c])
            if distinct
                # The largest regions first, every cell takes the colour left whose tint differs most from those of the cells it meets.
                free = filter(j -> !(j in values(taken)), eachindex(colours))
                taken[c] = argmax(j -> minimum((hypot((tint[j] .- tint[taken[d]])...) for d in meets[c] if haskey(taken, d)); init=Inf), free)
            else
                # The largest regions first, every cell takes the first colour of none of the cells it meets.
                used = Set(taken[d] for d in meets[c] if haskey(taken, d))
                taken[c] = something(findfirst(j -> !(j in used), eachindex(colours)), mod1(length(taken) + 1, length(colours)))
            end
        end
        palette = Dict(c => colours[taken[c]] for c in cells)
    end
    shade = round(Int, 100PORTRAIT_REGION_TINT)
    # The margin: a strip beyond every side of a piece on the hull, and wedges at its corners.
    hull = portrait_hull(place, 1:n)
    m = PORTRAIT_MARGIN_PT * PORTRAIT_PT
    normal(a, b) = (d = b .- a; (d[2], -d[1]) ./ hypot(d...))
    on_hull = [Tuple{Int,NTuple{2,Float64},NTuple{2,Float64}}[] for _ in hull]
    margin = Tuple{Int,Vector{NTuple{2,Float64}}}[]
    # The regions of a connected atlas reach beyond the hull themselves.
    for (i, (a, b)) in enumerate(zip(hull, circshift(hull, -1))), (c, poly) in (connected ? () : pieces), (p, q) in zip(poly, circshift(poly, -1))
        A, B = place[a], place[b]
        hypot((p .- q)...) > 1e-9 && portrait_on_segment(p, A, B; tol=1e-7) && portrait_on_segment(q, A, B; tol=1e-7) || continue
        push!(margin, (c, [p, q, q .+ m .* normal(A, B), p .+ m .* normal(A, B)]))
        push!(on_hull[i], (c, p, q))
    end
    for (i, a) in enumerate(connected ? () : hull)
        A, before = place[a], mod1(i - 1, length(hull))
        touching(j) = first(c for (c, p, q) in on_hull[j] if min(hypot((p .- A)...), hypot((q .- A)...)) < 1e-7)
        θ1 = atan(reverse(normal(place[hull[before]], A))...)
        θ2 = atan(reverse(normal(A, place[hull[mod1(i + 1, end)]]))...)
        θ2 < θ1 && (θ2 += 2π)
        arc(s, t) = [A .+ m .* (cos(θ), sin(θ)) for θ in range(s, t; length=8)]
        push!(margin, (touching(before), [A; arc(θ1, (θ1 + θ2) / 2)]), (touching(i), [A; arc((θ1 + θ2) / 2, θ2)]))
    end
    # The sectors of the regions and margin around the point `p`, as (cell, from, to) in degrees:
    # the cells on a small circle around it, every change of cell found to a hundredth of a degree.
    region_at(p) = (i = findfirst(((_, poly),) -> portrait_inside(p, poly), [pieces; margin]); isnothing(i) ? 0 : [pieces; margin][i][1])
    function region_sectors(p; radius=0.02)
        towards(θ) = region_at(p .+ radius .* (cosd(θ), sind(θ)))
        θs = 0.5:1.0:359.5
        seen = towards.(θs)
        0 in seen && error("a point of the core at $p outside every region")
        changes = [k for k in eachindex(θs) if seen[k] != seen[mod1(k - 1, end)]]
        isempty(changes) && return [(seen[1], 0.0, 360.0)]
        borders_at = map(changes) do k
            below, above = θs[mod1(k - 1, end)] - (k == 1 ? 360 : 0), θs[k]
            for _ in 1:20
                middle = (below + above) / 2
                towards(middle) == seen[k] ? (above = middle) : (below = middle)
            end
            (below + above) / 2
        end
        return [(seen[k], borders_at[j], borders_at[mod1(j + 1, end)] + (j == length(changes) ? 360 : 0))
                for (j, k) in enumerate(changes)]
    end
    side_cell = Dict((key(p), key(q)) => c for (c, poly) in (connected ? () : pieces) for (p, q) in zip(poly, circshift(poly, -1)))
    borders = connected ? walls : [(p, q) for ((p, q), c) in side_cell if p < q && get(side_cell, (q, p), c) != c]
    # Frames in two columns, of the split of the cells into halves whose lines are shortest in all.
    function column_frames(column, side)
        column = sort(column; by=c -> -centre[c][2])
        k = length(column)
        gap = (H - k * fh) / max(k - 1, 1)
        y = [clamp(centre[c][2], fh / 2, H - fh / 2) for c in column]
        for _ in 1:50, i in 1:k-1
            short = fh + gap - (y[i] - y[i+1])
            short > 0 && (y[i] += short / 2; y[i+1] -= short / 2)
            y .= clamp.(y, fh / 2, H - fh / 2)
        end
        return map(zip(column, y)) do (c, fy)
            fx = side < 0 ? fw / 2 : W - fw / 2
            foot = (fx - side * fw / 2, fy)
            # The line ends a little inside the point of the region nearest to the frame.
            near = argmin(q -> hypot((q .- foot)...), [curve_closest(foot, p, q) for (d, poly) in pieces if d == c
                                                        for (p, q) in zip(poly, circshift(poly, -1))])
            (cell=c, x=fx, y=fy, foot, tip=near .+ 0.25 .* (centre[c] .- near) ./ hypot((centre[c] .- near)...))
        end
    end
    frames, source, written = [], nothing, []
    if labels
        written = portrait_label_places(pieces, borders, place, ends, colour)
    else
        centre = Dict(c => portrait_centroid(argmax(poly -> abs(portrait_area(poly)), [poly for (d, poly) in pieces if d == c])) for c in cells)
        splits = [left for mask in 0:(1<<length(cells))-1 for left in (cells[[isodd(mask >> (i - 1)) for i in eachindex(cells)]],)
                  if abs(2length(left) - length(cells)) <= 1]
        frames = argmin(fs -> sum(hypot((f.tip .- f.foot)...) for f in fs),
                        [[column_frames(left, -1); column_frames(setdiff(cells, left), 1)] for left in splits])
        source = load_g2_maximal_cell_source(; field=:qq)
    end
    core_middle = offset .+ unit .* portrait_middle(pos)
    if measure
        ink = [[p for (_, poly) in [pieces; margin] for p in poly]; place;
               [p .+ (x, y) for (_, p, _) in written for x in (-0.25, 0.25) for y in (-0.15, 0.15)];
               [(f.x + x, f.y + y) for f in frames for x in (-fw / 2, fw / 2) for y in (-fh / 2, fh / 2)]]
        return portrait_reach(ink, core_middle, 0.0)
    end
    polygon(poly) = join(portrait_point.(poly), " -- ") * " -- cycle"
    return sprint() do io
        println(io, "% $name: the core and the regions of its cells, with ", labels ? "their types" : "their G2-trees and types",
                ". Written by curve_portraits.jl.")
        print(io, "\\documentclass[12pt,tikz,border=1pt]{standalone}\n\\usepackage{amsmath,pdfrender}\n",
              "% The types as main.tex writes them.\n",
              "\\providecommand{\\celltype}[1]{\\text{\\textpdfrender{TextRenderingMode=Fill,LineWidth=.05pt}",
              "{\\normalfont\\bfseries\\oldstylenums{#1}}}}\n\\begin{document}\n")
        # The styles open a group, which keeps the colours of the cells too.
        styles = g2_minimal_tree_styles(io, :color)
        foreach(c -> println(io, "\\definecolor{cell$c}{HTML}{$(palette[c])}"), cells)
        # Labels beside the hull may pass the box of the picture, which then holds what it draws.
        println(io, "\\begin{tikzpicture}[x=1cm,y=1cm]", labels ? "" : "\n\\useasboundingbox (0,0) rectangle ($W,$H);")
        if connected
            # A region as one path of its outlines, so that its cells leave no seams; its holes run clockwise.
            for c in cells
                loops = portrait_join(outline[c]; directed=true, tol=1e-6)
                println(io, "\\fill[cell$(c)!$shade] ", join(polygon.(loops), " "), ";")
            end
            for line in portrait_join(borders; directed=false, tol=1e-6)
                println(io, "\\draw[white,line width=.8pt,line cap=round,line join=round] ", join(portrait_point.(line), " -- "), ";")
            end
        else
            for (c, poly) in [pieces; margin]
                println(io, "\\fill[cell$(c)!$shade,draw=cell$(c)!$shade,line width=.25pt,line join=round] ", polygon(poly), ";")
            end
            foreach(((p, q),) -> println(io, "\\draw[white,line width=.8pt,line cap=round] ", portrait_point(p), " -- ", portrait_point(q), ";"), borders)
        end
        foreach(f -> println(io, "\\draw[cell$(f.cell),line width=.7pt] ", portrait_point(f.foot), " -- ", portrait_point(f.tip), ";"), frames)
        for (u, v) in ends
            println(io, "\\draw[cell$(colour[minmax(u, v)]),line width=1.5pt,line cap=round] ", portrait_point(place[u]), " -- ", portrait_point(place[v]), ";")
        end
        # The border of an edge over another has the colour of the region below the crossing.
        tint(p) = (c = region_at(p); c == 0 ? "white" : "cell$(c)!$shade")
        if !isnothing(crossing)
            over = colour[minmax(crossing.over[1], crossing.x)]
            print(io, portrait_over(place, crossing.x, crossing.over, crossing.under, portrait_point, 1.0, 1.5, "cell$over";
                                    border=tint(place[crossing.x]), gap=PORTRAIT_CROSSING_GAP_PT))
        end
        for (k, l, x) in (connected ? portrait_all_crossings(place, ends) : ())
            # The first edge over the second, drawn again only near the crossing.
            (a, b) = ends[k]
            d = place[b] .- place[a]
            ℓ = hypot(d...)
            near = clamp.((hypot((x .- place[a])...) .+ (-1, 1) .* 0.4) ./ ℓ, 0, 1)
            print(io, portrait_over([place; [place[a] .+ near[1] .* d, x, place[a] .+ near[2] .* d]], n + 2, (n + 1, n + 3), ends[l],
                                    portrait_point, 1.0, 1.5, "cell$(colour[minmax(a, b)])"; border=tint(x), gap=PORTRAIT_CROSSING_GAP_PT))
        end
        for v in 1:n
            isnothing(crossing) || v != crossing.x || continue
            around = [(w, colour[minmax(v, w)]) for w in neighbours[v]]
            here = unique(last.(around))
            trivalent && length(around) == 2 && continue
            if length(here) == 1
                println(io, "\\fill[cell$(only(here))] ", portrait_point(place[v]), " circle ($(PORTRAIT_NODE_PT)pt);")
            end
            # A vertex between cells: split as the regions split the plane around it.
            for (c, from, to) in (length(here) > 1 ? region_sectors(place[v]) : ())
                @printf(io, "\\fill[cell%d] %s -- ++(%.2f:%gpt) arc[start angle=%.2f,end angle=%.2f,radius=%gpt] -- cycle;\n",
                        c, portrait_point(place[v]), from, PORTRAIT_NODE_PT, from, to, PORTRAIT_NODE_PT)
            end
            println(io, "\\draw[white,line width=.5pt] ", portrait_point(place[v]), " circle ($(PORTRAIT_NODE_PT)pt);")
        end
        for f in frames
            @printf(io, "\\draw[cell%d,line width=1.2pt,rounded corners=4pt,fill=white] (%.3f,%.3f) rectangle (%.3f,%.3f);\n",
                    f.cell, f.x - fw / 2, f.y - fh / 2, f.x + fw / 2, f.y + fh / 2)
            @printf(io, "\\node[anchor=north west,inner sep=2pt] at (%.3f,%.3f) {\\celltype{%d}};\n",
                    f.x - fw / 2 + 0.05, f.y + fh / 2 - 0.05, curve.number[f.cell])
            @printf(io, "\\begin{scope}[shift={(%.3f,%.3f)}]\n", f.x, f.y - 0.06)
            g2_tree_commands(io, g2_maximal_cell_tree(f.cell; source).tree, styles; scale=frame_scale * PORTRAIT_TREE_SCALE,
                             taxon_scale=frame_scale * PORTRAIT_TAXON_SCALE, edge_labels=:none)
            println(io, "\\end{scope}")
        end
        for (c, p, tip) in written
            # A label outside its regions points to them with a line of the colour of its edges.
            println(io, "\\node[text=cell$c,inner sep=0pt,outer sep=1.5pt,font=\\small] (label$c) at ", portrait_point(p),
                    " {\\celltype{$(curve.number[c])}};")
            isnothing(tip) || println(io, "\\draw[cell$c,line width=.6pt,line cap=round] (label$c) -- ", portrait_point(tip), ";")
        end
        isnothing(box) || println(io, portrait_box(core_middle, box))
        println(io, "\\end{tikzpicture}\n\\endgroup\n\\end{document}")
    end
end

"""
The labels of the regions `pieces` of the cells of `colour`, cut by the
segments `borders`, around the core with the vertices at `place` and the
edges `ends`, as (cell, point, tip): in every region, a connected part of
the pieces of a cell on a grid of `PORTRAIT_LABEL_GRID`, its point of the
grid farthest from the borders, the hull and the core, when that room is at
least `PORTRAIT_LABEL_ROOM`, and the tip nothing; a cell without such a
region outside the hull, `PORTRAIT_LABEL_BESIDE` from it, where it is nearest
to its regions and apart from the other labels if it can be, and the tip of
its line to them halfway from their nearest point to the nearest edge of the
cell, which lies in them too.
"""
function portrait_label_places(pieces, borders, place, ends, colour)
    hull = place[portrait_hull(place, eachindex(place))]
    walls = [borders; collect(zip(hull, circshift(hull, -1)))]
    core = [(place[u], place[v]) for (u, v) in ends]
    gap(p, (a, b)) = hypot((p .- curve_closest(p, a, b))...)
    room(p) = min(minimum(s -> gap(p, s), walls), minimum(s -> gap(p, s), core) - 0.75PORTRAIT_PT,
                  minimum(q -> hypot((p .- q)...), place) - PORTRAIT_NODE_PT * PORTRAIT_PT)
    xs = minimum(first.(hull)):PORTRAIT_LABEL_GRID:maximum(first.(hull))
    ys = minimum(last.(hull)):PORTRAIT_LABEL_GRID:maximum(last.(hull))
    grid = [(x, y) for x in xs, y in ys]
    owner = map(grid) do p
        i = findfirst(((_, poly),) -> portrait_inside(p, poly), pieces)
        isnothing(i) ? 0 : pieces[i][1]
    end
    # The regions, flooded over the grid, each with its point of most room.
    seen = falses(size(grid))
    best = Dict{Int,Vector{Tuple{Float64,NTuple{2,Float64},Int}}}()
    for start in CartesianIndices(grid)
        (seen[start] || owner[start] == 0) && continue
        c, queue = owner[start], [start]
        seen[start] = true
        for i in queue, d in (CartesianIndex(1, 0), CartesianIndex(-1, 0), CartesianIndex(0, 1), CartesianIndex(0, -1))
            j = i + d
            checkbounds(Bool, grid, j) && !seen[j] && owner[j] == c && (seen[j] = true; push!(queue, j))
        end
        push!(get!(best, c, []), (argmax(first, [(room(grid[i]), grid[i]) for i in queue])..., length(queue)))
    end
    cells = sort!(unique(colour[minmax(u, v)] for (u, v) in ends))
    labels = Tuple{Int,NTuple{2,Float64},Union{Nothing,NTuple{2,Float64}}}[(c, p, nothing) for c in cells
                                                                           for (r, p, _) in get(best, c, []) if r >= PORTRAIT_LABEL_ROOM]
    free(p) = min(minimum(s -> gap(p, s), core) - 0.75PORTRAIT_PT, minimum(q -> hypot((p .- q)...), place) - PORTRAIT_NODE_PT * PORTRAIT_PT,
                  minimum((hypot((p .- q)...) - 2PORTRAIT_LABEL_ROOM for (_, q, _) in labels); init=Inf))
    # The places beside the hull for the labels of the cells without room, nearest to their regions.
    outside(p) = !portrait_inside(p, hull) && PORTRAIT_LABEL_BESIDE[1] <= minimum(s -> gap(p, s), zip(hull, circshift(hull, -1))) <= PORTRAIT_LABEL_BESIDE[2]
    beside = filter(outside, [(x, y) for x in xs[1]-PORTRAIT_LABEL_BESIDE[2]:PORTRAIT_LABEL_GRID:xs[end]+PORTRAIT_LABEL_BESIDE[2]
                              for y in ys[1]-PORTRAIT_LABEL_BESIDE[2]:PORTRAIT_LABEL_GRID:ys[end]+PORTRAIT_LABEL_BESIDE[2]])
    for c in cells
        any(l -> l[1] == c, labels) && continue
        sides = [(a, b) for (d, poly) in pieces if d == c for (a, b) in zip(poly, circshift(poly, -1))]
        apart = filter(p -> free(p) >= 0, beside)
        p = argmin(p -> minimum(s -> gap(p, s), sides), isempty(apart) ? beside : apart)
        near = argmin(q -> hypot((q .- p)...), [curve_closest(p, a, b) for (a, b) in sides])
        tip = (near .+ argmin(q -> hypot((q .- near)...), [curve_closest(near, a, b) for ((a, b), (u, v)) in zip(core, ends)
                                                           if colour[minmax(u, v)] == c])) ./ 2
        push!(labels, (c, p, tip))
    end
    return labels
end

"A number of point_places.csv, an integer or a fraction `a//b`."
portrait_number(s) = (parts = parse.(Int, split(s, "//")); length(parts) == 1 ? Float64(only(parts)) : parts[1] / parts[2])

"The places of the points of `path`, a point_places.csv of curve_details.jl, on the reduced graph, see `portrait_points`."
function portrait_read_places(path)
    lines = readlines(path)
    header = split(first(lines), ',')
    columns = [findfirst(==(name), header) for name in ("reduced_edge", "from_tail", "length")]
    any(isnothing, columns) && error("$path has no places on the reduced graph; rerun curve_details.jl")
    return map(lines[2:end]) do line
        edge, from_tail, len = split(line, ',')[columns]
        (parse(Int, edge), portrait_number(from_tail), isempty(len) ? nothing : portrait_number(len))
    end
end

function curve_portraits_main(args=ARGS)
    graphs = filter(arg -> !startswith(arg, "--"), args)
    usage = "usage: curve_portraits.jl GRAPH.csv [--output-prefix=PREFIX] [--outer-face=K] [--seed=S] [--no-pdf] " *
            "[--fragments=FIGURE [--standalone]] [--points=POINTS] [--crossing=C] [--side-by-side] [--chain-power=P] [--tree-scale=F] [--frames] " *
            "[--trivalent-nodes] [--connected-regions] [--region-starts=K] [--region-moves=M] [--region-drawing=FILE] " *
            "[--region-smoothing=S] [--atlas-scale=F] [--overlapping-trees]"
    length(graphs) == 1 || error(usage)
    for arg in filter(arg -> startswith(arg, "--"), args)
        arg in ("--no-pdf", "--side-by-side", "--frames", "--trivalent-nodes", "--connected-regions", "--overlapping-trees",
                "--standalone") ||
            any(startswith(arg, "--$name=") for name in ("output-prefix", "outer-face", "seed", "fragments", "points", "crossing",
                                                          "chain-power", "tree-scale", "region-starts", "region-moves",
                                                          "region-drawing", "region-smoothing", "atlas-scale")) ||
            error("unknown option $arg; $usage")
    end
    fragments = g2_tree_option(args, "fragments")
    "--standalone" in args && isnothing(fragments) && error("--standalone writes the fragments as documents; pass --fragments=FIGURE")
    side_by_side, frames, trivalent = "--side-by-side" in args, "--frames" in args, "--trivalent-nodes" in args
    connected = "--connected-regions" in args
    graph = abspath(only(graphs))
    prefix = abspath(something(g2_tree_option(args, "output-prefix"), splitext(graph)[1]))
    mkpath(dirname(prefix))
    curve = portrait_curve(graph)
    core_curve = curve
    places = g2_tree_option(args, "points")
    places = isnothing(places) ? nothing : portrait_read_places(places)
    chosen, crossing = g2_tree_option(args, "outer-face"), g2_tree_option(args, "crossing")
    seeds = isnothing(g2_tree_option(args, "seed")) ? PORTRAIT_SEEDS : [parse(Int, g2_tree_option(args, "seed"))]
    chain_power = parse(Float64, something(g2_tree_option(args, "chain-power"), string(PORTRAIT_CHAIN_POWER)))
    tree_scale = parse(Float64, something(g2_tree_option(args, "tree-scale"), "1"))
    option(name, default) = parse(typeof(default), something(g2_tree_option(args, name), string(default)))
    overlap, share = "--overlapping-trees" in args, option("atlas-scale", PORTRAIT_REGION_SCALE)
    if connected
        # The drawing, by vertex of the graph, read from FILE, its lines after # aside, or annealed, then smoothed;
        # one made here is written beside the pictures.
        saved, smoothing = g2_tree_option(args, "region-drawing"), option("region-smoothing", 0)
        if isnothing(saved)
            drawn = portrait_region_drawing(core_curve; starts=option("region-starts", PORTRAIT_REGION_STARTS),
                                            moves=option("region-moves", PORTRAIT_REGION_MOVES), chain_power)
        else
            rows = split.(filter(line -> !startswith(line, '#'), readlines(saved))[2:end], ',')
            at = Dict(parse(Int, r[1]) => (parse(Float64, r[2]), parse(Float64, r[3])) for r in rows)
            drawn = [at[v] for v in core_curve.vertices]
        end
        smoothing > 0 && (drawn = portrait_region_smoothing(core_curve, drawn; steps=smoothing))
        if isnothing(saved) || smoothing > 0
            write(prefix * "_region_drawing.csv",
                  join(["vertex,x,y"; [@sprintf("%d,%.6f,%.6f", v, p...) for (v, p) in zip(core_curve.vertices, drawn)]], "\n") * "\n")
            println("Wrote ", prefix * "_region_drawing.csv")
        end
        # Both pictures draw the core so.
        pos = drawn
    else
        variants = portrait_planar(curve) ? [curve] : portrait_crossings(curve; avoid=isnothing(places) ? () : first.(places))
        isempty(variants) && error("the core of $graph needs more than one crossing")
        # A crossing by the edges of the graph that cross.
        edge_at(variant, ends) = curve.core[findfirst(==(ends), curve.core_ends)].edge
        describe(variant) = isnothing(variant.crossing) ? "" :
            "crossing $(findfirst(==(variant), variants)), edge $(edge_at(variant, variant.crossing.over)) over " *
            "edge $(edge_at(variant, variant.crossing.under)), "
        isnothing(first(variants).crossing) || println(length(variants), " crossings make the core planar")
        jobs = []
        for variant in (isnothing(crossing) ? variants : [variants[parse(Int, crossing)]])
            candidates = portrait_outer_faces(variant.ends, variant.colour, variant.n; variant.crossing)
            isnothing(variant.crossing) && println(length(candidates), " faces of the core meet the most cells")
            foreach(((f, seed),) -> push!(jobs, (variant, candidates, f, seed)),
                    Iterators.product(isnothing(chosen) ? eachindex(candidates) : [parse(Int, chosen)], seeds))
        end
        # The drawings, on all threads; one whose fans cannot be drawn is dropped.
        drawings, failures = Vector{Any}(nothing, length(jobs)), fill("the crossing leaves its segments", length(jobs))
        Threads.@threads for j in eachindex(jobs)
            # Local to each job, not shared by the threads.
            local variant, candidates, f, seed, drawn
            variant, candidates, f, seed = jobs[j]
            drawn = portrait_core_drawing(variant.ends, variant.colour, variant.n, candidates[f]; seed, variant.crossing, chain_power)
            isnothing(drawn) && continue
            try
                drawings[j] = (portrait_fans(variant, drawn; scale=tree_scale)[3], variant, f, seed, drawn)
            catch failure
                failure isa ErrorException || rethrow()
                failures[j] = failure.msg
            end
        end
        for ((variant, _, f, seed), drawing, failure) in zip(jobs, drawings, failures)
            isnothing(drawing) ? println(describe(variant), "outer face $f, seed $seed: ", failure) :
                @printf("%souter face %d, seed %d: smallest first edge of a fan %.3f\n", describe(variant), f, seed, drawing[1])
        end
        _, curve, f, seed, pos = argmax(first, filter(!isnothing, drawings))
        println("Kept ", describe(curve), "outer face $f, seed $seed")
    end
    name = basename(graph)
    layout = portrait_berkovich_layout(curve, pos; side_by_side, tree_scale, frames, overlap, atlas=connected ? share : nothing)
    points = isnothing(places) ? nothing : portrait_points(curve, layout, places)
    atlas = (curve=curve, pos=pos, unit=layout.core_unit)
    # With --connected-regions, the core of the atlas is A times as large as in the Berkovich picture.
    connected && (atlas = (; atlas..., unit=share * layout.core_unit))
    regions_tex(; box=nothing, measure=false) = portrait_regions_tex(atlas.curve, atlas.pos, name; atlas.unit, labels=side_by_side && !frames,
                                                                     layout.frame_scale, layout.atlas_height, trivalent, connected, box, measure)
    # Side by side, the two pictures in boxes of one height, their cores at the same height.
    reach = connected && side_by_side ? [portrait_berkovich_reach(curve, layout, points), regions_tex(; measure=true)] : nothing
    box(k) = isnothing(reach) ? nothing : (reach[k][1], reach[k][2], maximum(r[3] for r in reach), maximum(r[4] for r in reach))
    for (suffix, figure, tex) in (("_berkovich", "_berkovich", portrait_berkovich_tex(curve, layout, name; points, trivalent, box=box(1))),
                                  ("_core_regions", "_regions", regions_tex(; box=box(2))))
        write(prefix * suffix * ".tex", tex)
        println("Wrote ", prefix * suffix * ".tex")
        "--no-pdf" in args || println("Wrote ", compile_curve_atlas(prefix * suffix * ".tex", prefix * suffix * ".pdf"))
        isnothing(fragments) && continue
        # The fragment: the body of the document in a group, after its first comment and the command.
        lines = split(tex, '\n')
        body = lines[findfirst(==("\\begin{document}"), lines)+1:findfirst(==("\\end{document}"), lines)-1]
        first(body) == "\\begingroup" || (body = ["\\begingroup"; body; "\\endgroup"])
        # Outside the picture every line ends with %, so that the fragment adds no space in a row of pictures.
        inside = cumsum(startswith.(body, "\\begin{tikzpicture}")) .- cumsum(startswith.(body, "\\end{tikzpicture}"))
        body = [inside[i] > 0 && !startswith(line, "\\begin{tikzpicture}") ? line : line * "%" for (i, line) in enumerate(body)]
        command = join(["% Generated by julia --project=julia julia/to_tex/curve_portraits.jl"; args], " ")
        fragment = join([lines[1]; command; body], "\n") * "\n"
        write(fragments * figure * ".tex", "--standalone" in args ? g2_standalone_figure(fragment) : fragment)
        println("Wrote ", fragments * figure * ".tex")
    end
    if side_by_side && !("--no-pdf" in args)
        # The row of the two pictures, centred on one line.
        row = prefix * "_side_by_side.tex"
        write(row, "% The Berkovich picture and the region atlas of $name side by side. Written by curve_portraits.jl.\n" *
                   "\\documentclass[border=2pt]{standalone}\n\\usepackage{graphicx}\n\\begin{document}\n" *
                   "\$\\vcenter{\\hbox{\\includegraphics{$(basename(prefix))_berkovich.pdf}}}\\hspace{$(PORTRAIT_SIDE_BY_SIDE.gap)cm}" *
                   "\\vcenter{\\hbox{\\includegraphics{$(basename(prefix))_core_regions.pdf}}}\$\n\\end{document}\n")
        println("Wrote ", compile_curve_atlas(row, prefix * "_side_by_side.pdf"))
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    curve_portraits_main()
end
