###############################################################################
# The fins of the bounded part of the tropical K3 surface X18 of
# k3_surface_polygons.jl, and the 2-sphere that is left without them.
#
# The bounded part of X18 is a homology 2-sphere, but not a surface: some of
# its edges lie in 1 or in 3 to 5 bounded polygons. Its second homology over QQ
# is spanned by one cycle z. The polygons off the support of z are the fins.
# The script checks that the support of z is a 2-sphere: every edge lies in
# exactly two of its polygons, the link of every vertex is one circle, and the
# Euler characteristic is 2. It then checks that the core of the curve X17 runs
# on the sphere, and that the core cuts the sphere into as many discs as the
# planar core has faces, E - V + 2, so that the core is the edge graph of the
# sphere, cut into these faces. Last, the core without its vertices of valence
# 2 is simple, trivalent and 3-connected, so by Steinitz's theorem it is the
# edge graph of a simple 3-polytope.
#
#   julia --project=. studies/g2/curves/k3_fins.jl DIR [K] [--input=FILE]
#
# DIR and K as for k3_surface_polygons.jl, whose output FILE is read, by
# default results/curves/genus_10/k3/DIR_x18.mrdi. The polygons are
# exact; only the cyclic order of the vertices of a polygon is found in
# floating point, and the edges it gives are checked against the exact edge
# keys of the polygon.
###############################################################################

include(joinpath(@__DIR__, "k3_surface_polygons.jl"))
import Combinatorics
using Printf

"The edges of a convex polygon with the exact vertices `verts`, as pairs of consecutive vertices."
function cyclic_edges(verts)
    m = sum(verts) ./ length(verts)
    u = Float64.(verts[1] .- m)
    # The second axis: the vertex farthest from the line through m and the first vertex.
    across = [(d = Float64.(v .- m); d .- (sum(d .* u) / sum(u .* u)) .* u) for v in verts]
    w = across[argmax(sum.(abs2, across))]
    angle(v) = (d = Float64.(v .- m); atan(sum(d .* w), sum(d .* u)))
    order = sortperm(verts; by=angle)
    return [(verts[order[i]], verts[order[mod1(i + 1, length(order))]]) for i in eachindex(order)]
end

"The edge key of k3_surface_polygons.jl of the edge from `a` to `b` of the polygon `(c, S)`, `nothing` on a wall."
function curve_key(c, S, a, b, t)
    x = (a .+ b) ./ 2
    s = first(S)
    tight = [k for k in eachindex(x) if !(k in S) && x[k] - t[k] == x[s] - t[s]]
    length(tight) <= 1 || error("an edge of ($c, $S) attains the minimum $(length(tight) + 4) times")
    return isempty(tight) ? nothing : (:curve, c, sort([S; only(tight)]))
end

"The indices of the edges of `G` in its core: the legs are dropped, then the leaves removed until none is left."
function core_edges(G)
    kept = [k for (k, e) in enumerate(G.edges) if !is_leg(e)]
    while true
        valence = zeros(Int, length(G.points))
        foreach(k -> (valence[G.edges[k].tail] += 1; valence[G.edges[k].head] += 1), kept)
        any(k -> valence[G.edges[k].tail] == 1 || valence[G.edges[k].head] == 1, kept) || return kept
        filter!(k -> valence[G.edges[k].tail] > 1 && valence[G.edges[k].head] > 1, kept)
    end
end

function find_root!(parent, x)
    while parent[x] != x
        parent[x] = parent[parent[x]]
        x = parent[x]
    end
    return x
end

"The number of classes of `1:n` joined by the pairs `pairs`."
function n_classes(n, pairs)
    parent = collect(1:n)
    foreach(((a, b),) -> parent[find_root!(parent, a)] = find_root!(parent, b), pairs)
    return count(x -> find_root!(parent, x) == x, 1:n)
end

function k3_fins(x18, G)
    t = x18.t
    bounded = [(cs, p) for (cs, p) in x18.polygons if p.bounded]
    # Vertices, edges and the oriented boundary of each bounded polygon.
    vertex_index = Dict{Vector{QQFieldElem},Int}()
    vertex(v) = get!(vertex_index, v, length(vertex_index) + 1)
    edge_index = Dict{Tuple{Int,Int},Int}()
    edge_curve_key = Dict{Int,Any}()
    boundaries = Vector{Vector{Tuple{Int,Int}}}()   # (edge, sign) for each polygon
    for ((c, S), p) in bounded
        pairs = cyclic_edges(p.verts)
        length(pairs) == length(p.keys) ||
            error("polygon ($c, $S): $(length(pairs)) sides, $(length(p.keys)) edge keys")
        found = [curve_key(c, S, a, b, t) for (a, b) in pairs]
        Set(filter(!isnothing, found)) == Set(k for k in p.keys if k[1] == :curve) ||
            error("polygon ($c, $S): the sides on the curve differ from its edge keys")
        push!(boundaries, map(zip(pairs, found)) do ((a, b), key)
            i, j = vertex(a), vertex(b)
            e = get!(edge_index, minmax(i, j), length(edge_index) + 1)
            isnothing(key) || (get!(edge_curve_key, e, key) == key || error("two keys for one edge"))
            (e, i < j ? 1 : -1)
        end)
    end
    nv, ne, nf = length(vertex_index), length(edge_index), length(bounded)
    @printf("bounded part: %d vertices, %d edges, %d polygons, Euler characteristic %d\n", nv, ne, nf, nv - ne + nf)
    polygons_of_edge = [Int[] for _ in 1:ne]
    foreach(f -> foreach(((e, _),) -> push!(polygons_of_edge[e], f), boundaries[f]), 1:nf)
    println("edges by the number of bounded polygons containing them: ", histogram(length.(polygons_of_edge)))

    # The second homology over QQ: the kernel of the boundary map from polygons to edges.
    D = zero_matrix(QQ, ne, nf)
    foreach(f -> foreach(((e, s),) -> D[e, f] += s, boundaries[f]), 1:nf)
    b2, Z = nullspace(D)
    println("b2 = $b2")
    b2 == 1 || error("the second homology is not one-dimensional")
    z = [Z[f, 1] for f in 1:nf]
    support = findall(!iszero, z)
    println("the cycle z: coefficients ", histogram(abs.(z[support] ./ z[first(support)])), " (relative to the first),",
            " on $(length(support)) polygons; $(nf - length(support)) polygons are fins")

    # The support of z as a surface.
    on = falses(nf); on[support] .= true
    sphere_polygons_of_edge = [filter(f -> on[f], polygons_of_edge[e]) for e in 1:ne]
    sphere_edges = findall(!isempty, sphere_polygons_of_edge)
    println("edges of the support by the number of its polygons containing them: ",
            histogram(length.(sphere_polygons_of_edge[sphere_edges])))
    ends = Dict(e => ij for (ij, e) in edge_index)
    sphere_vertices = sort!(unique!([v for e in sphere_edges for v in ends[e]]))
    χ = length(sphere_vertices) - length(sphere_edges) + length(support)
    # The link of a vertex: its edges, joined by the polygons around it.
    edges_at = Dict(v => Int[] for v in sphere_vertices)
    foreach(e -> foreach(v -> push!(edges_at[v], e), ends[e]), sphere_edges)
    links = map(sphere_vertices) do v
        local_number = Dict(e => k for (k, e) in enumerate(edges_at[v]))
        pairs = Tuple{Int,Int}[]
        for f in support
            here = [local_number[e] for (e, _) in boundaries[f] if haskey(local_number, e)]
            isempty(here) && continue
            length(here) == 2 || error("a polygon meets the vertex $v in $(length(here)) edges")
            push!(pairs, (here[1], here[2]))
        end
        degrees = zeros(Int, length(local_number))
        foreach(((a, b),) -> (degrees[a] += 1; degrees[b] += 1), pairs)
        all(==(2), degrees) && n_classes(length(local_number), pairs) == 1
    end
    surface = all(e -> length(sphere_polygons_of_edge[e]) == 2, sphere_edges) && all(links)
    @printf("support: %d vertices, %d edges, %d polygons, Euler characteristic %d; %s\n", length(sphere_vertices),
            length(sphere_edges), length(support), χ,
            surface ? "every edge in 2 polygons and every vertex link one circle" :
                      "NOT a surface ($(count(!, links)) vertex links are not one circle)")
    position = Dict(f => k for (k, f) in enumerate(support))
    adjacent = Tuple{Int,Int}[]
    for e in sphere_edges
        fs = sphere_polygons_of_edge[e]
        foreach(f -> push!(adjacent, (position[fs[1]], position[f])), fs[2:end])
    end
    connected = n_classes(length(support), adjacent) == 1
    sphere = surface && connected && χ == 2
    sphere && println("the support of z is a connected closed surface with Euler characteristic 2: a 2-sphere")

    # The core of the curve on the sphere, and the pieces of the sphere cut along it.
    core_indices = core_edges(G)
    core_keys = Set(x18.curve_keys[k] for k in core_indices)
    curve_keys = Set(values(x18.curve_keys))
    sphere_keys = Set(edge_curve_key[e] for e in sphere_edges if haskey(edge_curve_key, e))
    fin_keys = setdiff(Set(edge_curve_key[e] for e in 1:ne if haskey(edge_curve_key, e)), sphere_keys)
    println("core: $(length(core_indices)) edges of the subdivided curve, $(length(core_keys)) edges of X18; ",
            "$(length(core_keys ∩ sphere_keys)) of them on the sphere")
    println("bounded edges of the curve: $(length(curve_keys)) edges of X18, $(length(curve_keys ∩ sphere_keys)) on ",
            "the sphere, $(length(curve_keys ∩ fin_keys)) on fins only")
    core_vertices = unique!([v for k in core_indices for v in (G.edges[k].tail, G.edges[k].head)])
    faces = length(core_indices) - length(core_vertices) + 2
    joins = Tuple{Int,Int}[]
    for e in sphere_edges
        get(edge_curve_key, e, nothing) in core_keys && continue
        fs = sphere_polygons_of_edge[e]
        foreach(f -> push!(joins, (position[fs[1]], position[f])), fs[2:end])
    end
    pieces = n_classes(length(support), joins)
    println("the sphere cut along the core: $pieces pieces; the planar core has E - V + 2 = $faces faces")
    # Steinitz: the core, without its vertices of valence 2, is the edge graph of a simple 3-polytope when it is
    # simple, trivalent and 3-connected.
    C = desubdivide(core(G))
    n = length(C.points)
    ends_of_core = [(e.tail, e.head) for e in C.edges]
    simple = all(((a, b),) -> a != b, ends_of_core) && allunique(minmax.(first.(ends_of_core), last.(ends_of_core)))
    valences = [count(e -> v in e, ends_of_core) for v in 1:n]
    three_connected = all(Combinatorics.combinations(1:n, 2)) do removed
        rest = setdiff(1:n, removed)
        number = Dict(v => k for (k, v) in enumerate(rest))
        kept = [(number[a], number[b]) for (a, b) in ends_of_core if haskey(number, a) && haskey(number, b)]
        n_classes(length(rest), kept) == 1
    end
    println("the core: $n vertices of valences ", histogram(valences), ", $(length(ends_of_core)) edges, ",
            simple ? "simple" : "NOT simple", ", ", three_connected ? "3-connected" : "NOT 3-connected")
    return (; sphere, pieces, faces, core_on_sphere=core_keys ⊆ sphere_keys,
            polytope_graph=simple && all(==(3), valences) && three_connected)
end

if abspath(PROGRAM_FILE) == @__FILE__
    positional = filter(arg -> !startswith(arg, "--"), ARGS)
    length(positional) in (1, 2) || error("usage: k3_fins.jl DIR [K] [--input=FILE]")
    dir = abspath(positional[1])
    k = length(positional) == 2 ? parse(Int, positional[2]) : 1
    given = [String(split(arg, '='; limit=2)[2]) for arg in ARGS if startswith(arg, "--input=")]
    input = isempty(given) ? project_path("results", "curves", "genus_10", "k3", basename(dir) * "_x18.mrdi") :
                             only(given)
    G = load_decorated_graphs(joinpath(dir, "decorated_graph_$k.mrdi")).subdivided
    result = k3_fins(load_k3_polygons(input), G)
    result.sphere && result.core_on_sphere && result.pieces == result.faces && result.polytope_graph ||
        error("the core is not the edge graph of a 2-sphere of X18: $result")
    println("the core is the edge graph of the 2-sphere left by removing the fins, and of a simple 3-polytope")
end
