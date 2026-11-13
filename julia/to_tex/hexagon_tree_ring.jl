# Writes ../figures/hexagon_tree_ring.tex of the repository: a tropical hexagon, in the notation of Section 9 of main.tex, drawn as its six G2-trees
# p^(12), p^(23), ..., p^(61) in the plane, such that consecutive trees p^(k-1,k) and p^(k,k+1) cross in their common point
# x^(k). The points x^(1), ..., x^(6) are apples at the corners of a regular hexagon. In the tree p^(k,k+1) the path from
# x^(k) to x^(k+1) is a chain from corner to corner, kinked only at its inner vertices, its pieces in proportion to their
# lengths, and the rest of the tree hangs outwards. An edge through an apple is drawn straight. Each tree has its own
# color; the taxa are the shapes of the paper, kept clear of each other. No lengths are written (`tree_tikz` can with
# `lengths=true`); pendant edges are not to scale.
#
# The hexagon is completed by tropical_complete_hexagon from the flags INPUT_1 = (x^(1), p^(12)) and INPUT_4 = (x^(4), p^(45)).
#
# Run from julia/:   julia --project=. to_tex/hexagon_tree_ring.jl [--complete]
#   [--output=FILE]   (write FILE instead of ../figures/hexagon_tree_ring.tex)
#   [--standalone]    (with --output: the document of g2_standalone_figure, whose PDF main.tex includes)

include(joinpath(@__DIR__, "..", "studies", "hexagon", "tropical_complete_hexagon.jl"))
include(joinpath(@__DIR__, "decorated_trees_figure.jl"))   # taxa, apples, trees.jl

const OUTPUT = normpath(joinpath(@__DIR__, "..", "..", "figures", "hexagon_tree_ring.tex"))   # the figure of main.tex

# The example of the search for presentable hexagons: (p_1, L_1) and (p_4, L_4), min convention. All six trees are
# maximal and of different topologies.
const INPUT_1 = ([4, 3, 11, 0, 15, 1, 2], [4, 4, 1, 6, 2, 3, 3, 8, 5, 13, 6, 0, 13, 1, 2, 2, 6, 3, 3, 4, 4])
const INPUT_4 = ([6, 3, 10, 0, 3, 4, 15], [7, 6, 4, 7, 8, 9, 3, 13, 18, 9, 6, 0, 3, 4, 13, 13, 6, 3, 9, 6, 7])
# The hexagon completed from them by tropical_complete_hexagon (each vector shifted to minimum 0): the points x^(1), ...,
# x^(6) and the lines p^(12), ..., p^(61). Stored, as the completion takes about 80 s; `--complete` computes it again and
# checks these values.
const HEXAGON = (
    points=[[4, 3, 11, 0, 15, 1, 2], [6, 3, 3, 0, 5, 1, 2], [3, 0, 2, 1, 4, 5, 5],
            [6, 3, 10, 0, 3, 4, 15], [6, 3, 3, 0, 3, 4, 6], [3, 0, 2, 1, 4, 5, 3]],
    lines=[[4, 4, 1, 6, 2, 3, 3, 8, 5, 13, 6, 0, 13, 1, 2, 2, 6, 3, 3, 4, 4],
           [6, 6, 3, 8, 4, 5, 3, 0, 5, 1, 2, 2, 9, 3, 4, 4, 2, 3, 5, 6, 6],
           [9, 8, 3, 6, 7, 11, 5, 0, 3, 4, 8, 2, 5, 6, 15, 15, 8, 5, 11, 8, 9],
           [7, 6, 4, 7, 8, 9, 3, 13, 18, 9, 6, 0, 3, 4, 13, 13, 6, 3, 9, 6, 7],
           [7, 6, 3, 6, 7, 9, 3, 0, 3, 4, 6, 2, 5, 6, 6, 6, 6, 3, 9, 6, 7],
           [4, 6, 3, 8, 4, 5, 5, 0, 7, 1, 2, 2, 15, 3, 4, 4, 2, 5, 5, 6, 4]])

# Drawing parameters, in cm and degrees.
const HEX_FLATTEN = 0.72               # the corners lie on an ellipse, flattened vertically by this factor
const HEX_RADIUS = 3.7              # the corners p_k of the hexagon; its sides have the same length
const BOUNDED_CM = 0.38                # drawn length of a hanging bounded edge per unit of length
const BOUNDED_MIN_CM = 0.55
const PENDANT_CM = 1.15                # drawn length of a pendant edge (or of its part beyond a point)
const BOW = -15                      # angle of the path to the side at p_k (positive: outwards)
const SPREAD_PER_LEAF = 22             # angle per taxon of the hanging subtrees
const MAX_FAN = 120                    # largest angle of the subtrees hanging at an inner vertex of the path
const FAN_MARGIN = 25                  # their angle to the path
const LEAN = 30                        # how much they lean away from the nearer corner, where the neighbour continues
const CORNER_FAN = 70                  # largest angle of the subtrees hanging at a corner
const TAXON_SCALE = 0.55
const APPLE_SCALE = 0.75
# Clearances kept by the taxa, in cm, from each other, from the apples and from the edges; the pendant edges turn by
# steps of RELAX_STEP degrees, at most MAX_TURN, to keep them.
const TAXON_GAP, APPLE_GAP, EDGE_GAP = 0.55, 0.5, 0.28
const RELAX_STEP, MAX_TURN = 3, 50
const CROSSING_PENALTY = 10.0          # for each two edges crossing
const FAN_TURN, FAN_STEP = 60, 5       # the fans of hanging subtrees turn by at most FAN_TURN, in steps of FAN_STEP
# p^(12), ..., p^(61): purple, teal, black, magenta, slate, brown; neighbours differ in hue and lightness
const TREE_COLORS = ("7B1FA2", "00897B", "212121", "D81B60", "78909C", "8D5524")

# ---------------------------------------------------------------------------------------------
# Positions of the points on the trees (as in hexagon_example_figure.jl)
# ---------------------------------------------------------------------------------------------

"Distance between two vertices of the tree."
function vertex_distance(T, a, b)
    a == b && return QQ(0)
    stack = [(a, 0, QQ(0))]
    while !isempty(stack)
        v, parent, d = pop!(stack)
        for (w, l) in T.adjacency[v]
            w == parent && continue
            w == b && return d + l
            push!(stack, (w, v, d + l))
        end
    end
    error("disconnected tree")
end

"""
The position of the point with valuations `x` on the tree `T`, as `(u, v, t)`: the point lies on the edge `(u, v)` at
distance `t` from `u`; for a leaf edge `u` is its inner vertex, and `t` may exceed the pendant length, as leaves are rays.
"""
function point_position(T, x)
    found = []
    for u in sort!(collect(keys(T.adjacency))), (v, l) in T.adjacency[u]
        u < v || continue
        A, B = leaves_on_side(T, u, v), leaves_on_side(T, v, u)
        alpha = unique([QQ(x[i]) + vertex_distance(T, i, u) for i in A])
        beta = unique([QQ(x[i]) + vertex_distance(T, i, v) + l for i in B])
        length(alpha) == 1 && length(beta) == 1 || continue
        t = (beta[1] - alpha[1]) / 2
        if u <= T.ntaxa                      # leaf edge: measure from the inner vertex v
            l - t >= 0 && push!(found, (v, u, l - t))
        elseif 0 <= t <= l
            push!(found, (u, v, t))
        end
    end
    isempty(found) && error("the point does not lie on the tree")
    return first(found)
end

# ---------------------------------------------------------------------------------------------
# A tree with its two points as vertices, and the path between them
# ---------------------------------------------------------------------------------------------

"""
The tree `T` as a mutable graph: `adj[v]` the neighbours, `len[(v, w)]` the lengths, `bounded[(v, w)]` whether the edge
lies on a bounded edge of `T`. Pendant edges get length `big`, as they are rays.
"""
function tree_graph(T, big)
    adj = Dict(v => Int[w for (w, _) in T.adjacency[v]] for v in keys(T.adjacency))
    len, bounded = Dict{Tuple{Int,Int},QQFieldElem}(), Dict{Tuple{Int,Int},Bool}()
    for v in keys(T.adjacency), (w, l) in T.adjacency[v]
        b = v > T.ntaxa && w > T.ntaxa
        len[(v, w)] = b ? l : QQ(big)
        bounded[(v, w)] = b
    end
    return (adj=adj, len=len, bounded=bounded)
end

"The path from `a` to `b` in the graph `G`."
function graph_path(G, a, b)
    stack = [[a]]
    while !isempty(stack)
        path = pop!(stack)
        path[end] == b && return path
        for w in G.adj[path[end]]
            length(path) > 1 && w == path[end - 1] && continue
            push!(stack, vcat(path, w))
        end
    end
    error("no path")
end

"Insert the vertex `p` at distance `t` from `u` on the path from `u` to `v` (an edge of the original tree); returns the vertex."
function insert_point!(G, u, v, t, p)
    path = graph_path(G, u, v)
    d = QQ(0)
    for (a, b) in zip(path[1:end-1], path[2:end])
        t == d && return a
        l = G.len[(a, b)]
        if t < d + l
            s = t - d
            filter!(!=(b), G.adj[a]); filter!(!=(a), G.adj[b])
            push!(G.adj[a], p); push!(G.adj[b], p); G.adj[p] = [a, b]
            for (x, y, m) in ((a, p, s), (p, b, l - s))
                G.len[(x, y)] = G.len[(y, x)] = m
                G.bounded[(x, y)] = G.bounded[(y, x)] = G.bounded[(a, b)]
            end
            return p
        end
        d += l
    end
    t == d && return path[end]
    error("the point lies beyond the edge")
end

"Taxa in the component of `v` after removing the edge to `parent`."
function graph_taxa(G, v, parent, ntaxa)
    v <= ntaxa && return [v]
    return reduce(vcat, [graph_taxa(G, w, v, ntaxa) for w in G.adj[v] if w != parent]; init=Int[])
end

# ---------------------------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------------------------

corner(k) = (a = 90 - 60 * (mod1(k, 6) - 1); (HEX_RADIUS * cosd(a), HEX_FLATTEN * HEX_RADIUS * sind(a)))
direction(θ) = (cosd(θ), sind(θ))
plus(a, b, s=1.0) = (a[1] + s * b[1], a[2] + s * b[2])

"Drawn length of an edge hanging off the path."
drawn_length(G, a, b) = G.bounded[(a, b)] ? max(BOUNDED_MIN_CM, BOUNDED_CM * Float64(G.len[(a, b)])) : PENDANT_CM

"""
Lay out the subtree at `v`, reached from `parent`, in the direction `θ` within a sector of width `width`, widened to
SPREAD_PER_LEAF per taxon.
"""
function layout_subtree!(pos, G, v, parent, θ, width, ntaxa)
    pos[v] = plus(pos[parent], direction(θ), drawn_length(G, parent, v))
    children = [w for w in G.adj[v] if w != parent]
    isempty(children) && return
    sizes = [length(graph_taxa(G, w, v, ntaxa)) for w in children]
    total = sum(sizes)
    width = max(width, SPREAD_PER_LEAF * (total - 1))
    start = θ + width / 2
    for (w, s) in zip(children, sizes)
        sub = width * s / total
        layout_subtree!(pos, G, w, v, start - sub / 2, sub, ntaxa)
        start -= sub
    end
end

"""
The drawing of the line `L_k` with its points `x_a = p_k`, `x_b = p_{k+1}`: coordinates of all vertices, the edges with
their lengths, and the vertices of the two points.
"""
function tree_layout(L, xa, xb, k)
    T = tree_from_min_pluecker(QQ.(L))
    ua, va, ta = point_position(T, xa)
    ub, vb, tb = point_position(T, xb)
    big = maximum(vcat([ta, tb], [l for v in keys(T.adjacency) for (_, l) in T.adjacency[v]])) + 1
    G = tree_graph(T, big)
    A = insert_point!(G, ua, va, ta, 101)
    B = insert_point!(G, ub, vb, tb, 102)
    A == B && error("p_$k and p_$(mod1(k + 1, 6)) lie at the same place of L_$k")
    path = graph_path(G, A, B)
    # The path is a chain from the corner p_k to the corner p_{k+1}, with kinks only at its inner vertices: it leaves p_k
    # at the angle BOW to the side and turns by equal angles, so that it arrives at the angle -BOW. A path of one segment
    # is straight. The chain is scaled and rotated to end at p_{k+1}; its pieces are in proportion to their lengths.
    P, Q = corner(k), corner(k + 1)
    side = atand(Q[2] - P[2], Q[1] - P[1])
    lengths = [Float64(G.len[(a, b)]) for (a, b) in zip(path[1:end-1], path[2:end])]
    m = length(lengths)
    turns = m == 1 ? [0.0] : [BOW * (1 - 2 * (i - 1) / (m - 1)) for i in 1:m]
    V = reduce(plus, [ℓ .* direction(side + d) for (ℓ, d) in zip(lengths, turns)])
    scale = hypot(Q[1] - P[1], Q[2] - P[2]) / hypot(V...)
    rotation = side - atand(V[2], V[1])
    angles = side .+ turns .+ rotation             # the directions of the pieces of the path
    pos = Dict{Int,Tuple{Float64,Float64}}(A => P)
    for (i, (a, b)) in enumerate(zip(path[1:end-1], path[2:end]))
        pos[b] = plus(pos[a], direction(angles[i]), scale * lengths[i])
    end
    pos[B] = Q                                     # exactly, so that the trees meet in the apple
    # The subtrees hanging at a vertex of the path share a fan, of SPREAD_PER_LEAF per taxon, each a sector in proportion
    # to its taxa. At a corner it points straight on along the path, so that an edge through an apple is straight, and is
    # at most CORNER_FAN wide. At an inner vertex it points into the outer angle between the two pieces of the path,
    # leaning away from the nearer corner. `rotate_fans!` turns them afterwards.
    drawn = cumsum(vcat(0.0, lengths))
    fans = Tuple{Int,Vector{Int}}[]                # vertex of the path => the vertices of its hanging subtrees
    for (i, v) in enumerate(path)
        children = [w for w in G.adj[v] if !(w in path)]
        isempty(children) && continue
        sizes = [length(graph_taxa(G, w, v, T.ntaxa)) for w in children]
        if v == A || v == B
            θ = v == A ? angles[1] + 180 : angles[end]
            width = min(CORNER_FAN, SPREAD_PER_LEAF * (sum(sizes) - 1))
        else
            forward, back = angles[i], angles[i - 1] + 180
            outer = mod(back - forward, 360)       # the outer side is to the left of the path
            width = min(MAX_FAN, outer - 2 * FAN_MARGIN, SPREAD_PER_LEAF * (sum(sizes) - 1))
            f = drawn[i] / drawn[end]              # near p_k lean forwards, near p_{k+1} backwards
            θ = clamp(forward + outer / 2 - LEAN * (1 - 2f), forward + FAN_MARGIN + width / 2, back - FAN_MARGIN - width / 2)
        end
        start = θ + width / 2
        for (w, s) in zip(children, sizes)
            sub = width * s / sum(sizes)
            layout_subtree!(pos, G, w, v, start - sub / 2, sub, T.ntaxa)
            start -= sub
        end
        # A fan at an apple inside an edge continues that edge, so it does not turn.
        v in (101, 102) || push!(fans, (v, reduce(vcat, [subtree_vertices(G, w, v) for w in children])))
    end
    edges = [(a, b, G.len[(a, b)], G.bounded[(a, b)], a in path && b in path) for a in keys(G.adj) for b in G.adj[a] if a < b]
    # The middle of the path, for the label of the line.
    j = findfirst(>=(drawn[end] / 2), drawn)
    middle = plus(pos[path[j - 1]], direction(angles[j - 1]), scale * (drawn[end] / 2 - drawn[j - 1]))
    return (tree=T, pos=pos, edges=edges, A=A, B=B, path=path, ntaxa=T.ntaxa, parent=Dict(t => only(G.adj[t]) for t in 1:T.ntaxa),
            fans=fans, start_angle=angles[1], end_angle=angles[end], middle=middle)
end

"The vertices of the subtree at `v`, reached from `parent`."
subtree_vertices(G, v, parent) = vcat([v], [u for w in G.adj[v] if w != parent for u in subtree_vertices(G, w, v)])

"The direction in which the apple at p_k hangs: into the inner angle between the paths of p^(k-1,k) and p^(k,k+1)."
function apple_angle(layouts, k)
    v = plus(direction(layouts[k].start_angle), direction(layouts[mod1(k - 1, 6)].end_angle + 180))
    return atand(v[2], v[1])
end

"The centre of the body of the apple at p_k."
apple_centre(layouts, k) = plus(corner(k), direction(apple_angle(layouts, k)), 0.35 * APPLE_SCALE)

"Do the segments pq and ab cross in a point inside both? Touching, as at a shared apple, does not count."
function segments_cross(p, q, a, b)
    orientation(u, v, w) = sign((v[1] - u[1]) * (w[2] - u[2]) - (v[2] - u[2]) * (w[1] - u[1]))
    return orientation(p, q, a) * orientation(p, q, b) < 0 && orientation(a, b, p) * orientation(a, b, q) < 0
end

dist(p, q) = hypot(p[1] - q[1], p[2] - q[2])
"`p` turned by `δ` degrees about `c`."
turn(p, c, δ) = (s = sind(δ); co = cosd(δ); (c[1] + co * (p[1] - c[1]) - s * (p[2] - c[2]), c[2] + s * (p[1] - c[1]) + co * (p[2] - c[2])))

"""
How badly the drawing is cluttered: CROSSING_PENALTY for each two edges crossing, and the squares of the shortfalls of the
taxa from the clearances TAXON_GAP to each other, APPLE_GAP to the apples and EDGE_GAP to the edges.
"""
function clutter(layouts, apples)
    taxa = [(k, t) for k in 1:6 for t in 1:layouts[k].ntaxa]
    segments = [(k, a, b, layouts[k].pos[a], layouts[k].pos[b]) for k in 1:6 for (a, b) in layouts[k].edges]
    P = [layouts[k].pos[t] for (k, t) in taxa]
    s = 0.0
    for i in eachindex(P)
        for j in (i + 1):length(P)
            s += max(0, TAXON_GAP - dist(P[i], P[j]))^2
        end
        for a in apples
            s += max(0, APPLE_GAP - dist(P[i], a))^2
        end
        k, t = taxa[i]
        for (k2, a, b, pa, pb) in segments
            k2 == k && t in (a, b) && continue
            s += max(0, EDGE_GAP - g2_segment_distance(P[i], pa, pb))^2
        end
    end
    for i in eachindex(segments), j in (i + 1):length(segments)
        segments_cross(segments[i][4], segments[i][5], segments[j][4], segments[j][5]) && (s += CROSSING_PENALTY)
    end
    return s
end

"""
Turn the fans of subtrees hanging at the vertices of the paths about these vertices, each in turn to its best angle
within ±FAN_TURN, as long as this reduces the `clutter`. Fans at an apple inside an edge stay, so the edge stays straight.
"""
function rotate_fans!(layouts)
    apples = [apple_centre(layouts, k) for k in 1:6]
    current = clutter(layouts, apples)
    turned = Dict((k, v) => 0 for k in 1:6 for (v, _) in layouts[k].fans)
    for _ in 1:10
        improved = false
        for k in 1:6, (v, vertices) in layouts[k].fans
            pos, c = layouts[k].pos, layouts[k].pos[v]
            original = Dict(u => pos[u] for u in vertices)
            best, bestδ = current, 0
            for δ in -FAN_TURN:FAN_STEP:FAN_TURN
                abs(turned[(k, v)] + δ) <= FAN_TURN || continue
                for u in vertices
                    pos[u] = turn(original[u], c, δ)
                end
                new = clutter(layouts, apples)
                if new < best - 1e-9
                    best, bestδ = new, δ
                end
            end
            for u in vertices
                pos[u] = turn(original[u], c, bestδ)
            end
            if bestδ != 0
                current, improved = best, true
                turned[(k, v)] += bestδ
            end
        end
        improved || break
    end
    return layouts
end

"""
Turn the pendant edges of all trees about their inner vertices, by steps of RELAX_STEP and at most MAX_TURN, as long as
this reduces the `clutter`. Pendant edges through an apple stay straight.
"""
function relax_taxa!(layouts)
    taxa = [(k, t) for k in 1:6 for t in 1:layouts[k].ntaxa]
    movable = [(k, t) for (k, t) in taxa if !(layouts[k].parent[t] in (101, 102))]   # 101, 102: apples inside edges
    apples = [apple_centre(layouts, k) for k in 1:6]
    penalty() = clutter(layouts, apples)
    turned = Dict(kt => 0 for kt in movable)
    current = penalty()
    for _ in 1:200
        improved = false
        for (k, t) in movable
            D = layouts[k]
            p, c = D.pos[t], D.pos[D.parent[t]]
            for δ in (RELAX_STEP, -RELAX_STEP)
                abs(turned[(k, t)] + δ) <= MAX_TURN || continue
                φ = atand(p[2] - c[2], p[1] - c[1]) + δ
                D.pos[t] = plus(c, direction(φ), dist(p, c))
                new = penalty()
                if new < current - 1e-12
                    current, improved = new, true
                    turned[(k, t)] += δ
                    break
                end
                D.pos[t] = p
            end
        end
        improved || break
    end
    if current >= CROSSING_PENALTY
        segments = [(k, a, b) for k in 1:6 for (a, b) in layouts[k].edges]
        at(k, v) = round.(layouts[k].pos[v]; digits=2)
        crossing = ["p^($(LINE_NAMES[k])) $(at(k, a))--$(at(k, b)) and p^($(LINE_NAMES[l])) $(at(l, c))--$(at(l, d))"
                    for (i, (k, a, b)) in enumerate(segments) for (l, c, d) in segments[(i + 1):end]
                    if segments_cross(layouts[k].pos[a], layouts[k].pos[b], layouts[l].pos[c], layouts[l].pos[d])]
        @warn "edges still cross: " * join(crossing, "; ")
    end
    current > 0 && @warn "some taxa are still too close to each other, an apple or an edge (clutter $current)"
    return layouts
end

# ---------------------------------------------------------------------------------------------
# TikZ
# ---------------------------------------------------------------------------------------------

function number_tex(q)
    q = QQ(q)
    denominator(q) == 1 && return string(numerator(q))
    return "\\tfrac{$(numerator(q))}{$(denominator(q))}"
end

xy(p) = @sprintf("(%.4f,%.4f)", p[1], p[2])
const LINE_NAMES = ("12", "23", "34", "45", "56", "61")

function tree_tikz(io, D, k, colors; lengths::Bool=false)
    c = "hexTree$k"
    for (a, b, l, bounded, onpath) in D.edges
        style = bounded ? "hex internal" : "hex pendant"
        @printf(io, "\\draw[%s,draw=%s] %s--%s;\n", style, c, xy(D.pos[a]), xy(D.pos[b]))
    end
    # Lengths: the bounded edges off the path, and every piece of the path, on the side away from the center (the apples
    # hang towards it).
    for (a, b, l, bounded, onpath) in D.edges
        lengths && (bounded || onpath) || continue
        p, q = D.pos[a], D.pos[b]
        m = ((p[1] + q[1]) / 2, (p[2] + q[2]) / 2)
        d = hypot(q[1] - p[1], q[2] - p[2])
        n = (-(q[2] - p[2]) / d, (q[1] - p[1]) / d)
        s = sign(m[1] * n[1] + m[2] * n[2] + 1e-9)
        @printf(io, "\\node[hex length,text=%s] at %s {\\(%s\\)};\n", c, xy(plus(m, n, 0.2s)), number_tex(l))
    end
    for v in keys(D.pos)
        v > D.ntaxa && v < 100 && @printf(io, "\\node[hex vertex,fill=%s] at %s {};\n", c, xy(D.pos[v]))
    end
    for t in 1:D.ntaxa
        println(io, "\\coordinate (v$t) at $(xy(D.pos[t]));")
        g2_minimal_tree_taxon(io, t, colors; label=false, taxon_scale=TAXON_SCALE)
    end
end

function document_body(H)
    io = IOBuffer()
    println(io, "% Generated by julia/to_tex/hexagon_tree_ring.jl")
    colors = g2_minimal_tree_styles(io, :color)
    print(io, g2_apple_colors_tex())
    for k in 1:6
        println(io, "\\definecolor{hexTree$k}{HTML}{$(TREE_COLORS[k])}")
    end
    println(io, "\\tikzset{")
    println(io, g2_apple_pic())
    println(io, raw"""  hex pendant/.style={line width=.8pt,line cap=round},
  hex internal/.style={line width=2.2pt,line cap=round},
  hex vertex/.style={circle,draw=white,line width=.4pt,inner sep=0pt,minimum size=4.2pt},
  hex length/.style={font=\sffamily\fontsize{7}{8}\selectfont\bfseries,inner sep=.5pt},
  hex line label/.style={font=\bfseries},
  hex point label/.style={font=\small,text=g2AppleDark},
}""")
    println(io, raw"\begin{tikzpicture}[x=1cm,y=1cm]")
    layouts = relax_taxa!(rotate_fans!([tree_layout(H.lines[k], H.points[k], H.points[mod1(k + 1, 6)], k) for k in 1:6]))
    for k in 1:6
        println(io, "% p^($(LINE_NAMES[k]))")
        tree_tikz(io, layouts[k], k, colors)
    end
    for k in 1:6
        # The apple at x^(k), its label beyond it, and the label of p^(k,k+1) inside the hexagon next to its path.
        a = apple_angle(layouts, k)
        @printf(io, "\\pic[rotate=%.1f,scale=%.3f] at %s {g2apple};\n", a + 90, APPLE_SCALE, xy(corner(k)))
        @printf(io, "\\node[hex point label] at %s {\\(a^{(%d)}\\)};\n", xy(plus(corner(k), direction(a), 0.85)), k)
        m = layouts[k].middle
        @printf(io, "\\node[hex line label,text=hexTree%d] at %s {\\(d^{(%s)}\\)};\n", k,
                xy(plus(m, m, -0.65 / hypot(m...))), LINE_NAMES[k])
    end
    println(io, raw"""% Center the hexagon of apples, whose center is the origin, on the page: extend the bounding box
% symmetrically about x=0, so that \centering in main.tex centers the origin and not the pendant edges.
\newdimen\hexbbwest \newdimen\hexbbeast
\pgfextractx{\hexbbwest}{\pgfpointanchor{current bounding box}{west}}
\pgfextractx{\hexbbeast}{\pgfpointanchor{current bounding box}{east}}
\path (-\hexbbeast,0cm) -- (-\hexbbwest,0cm);""")
    println(io, raw"\end{tikzpicture}")
    println(io, raw"\endgroup")
    return String(take!(io))
end

function hexagon_tree_ring_main(args=ARGS)
    if "--complete" in args
        t = @elapsed C = tropical_complete_hexagon(INPUT_1, INPUT_4)
        println("completed the hexagon in $(round(t; digits=1)) s"); flush(stdout)
        C.points == HEXAGON.points && C.lines == HEXAGON.lines ||
            error("the completion differs from HEXAGON: points $(C.points), lines $(C.lines)")
    end
    outputs = [String(split(arg, '='; limit=2)[2]) for arg in args if startswith(arg, "--output=")]
    standalone = "--standalone" in args
    standalone && isempty(outputs) && error("--standalone writes a document, which needs --output=FILE")
    output = isempty(outputs) ? OUTPUT : abspath(only(outputs))
    open(io -> print(io, standalone ? g2_standalone_figure(document_body(HEXAGON)) : document_body(HEXAGON)), output, "w")
    println("wrote ", output)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    hexagon_tree_ring_main()
end
