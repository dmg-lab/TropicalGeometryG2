###############################################################################
# How the region atlas of curve_portraits.jl of a planar core matches the
# tropical K3 surface X18 of k3_surface_polygons.jl, in two checks.
#
# 1. The bounded part of X18 and the faces of the drawing: the cells it meets
#    and their areas; its pieces once cut along the core, the polygons glued
#    along all their edges but those of the core, each with the faces of the
#    drawing whose borders hold all of its core edges; and the polygons ring
#    after ring from the core (ring 0 holds a core edge, ring r + 1 shares an
#    edge with ring r), with the share of each in the cells of the core.
# 2. The parts of the bounded part in the cells of the core against the
#    regions of the atlas: whether each part is connected and holds the core
#    edges of its cell, its share of their area, which parts meet along walls
#    and their shares of those walls, against the regions, their areas and the
#    lengths of their borders.
#
# The drawing is that of curve_portraits.jl with --outer-face=F --seed=S
# (default 1 and 0) and the default chain power on the reduced graph of
# experiment K of DIR, the regions those of `portrait_region_pieces`. Areas
# and lengths on X18 are Euclidean, in the slice sum(x) = 0.
#
#   julia --project=. studies/g2/curves/k3_core_regions.jl DIR [K] [--polygons=FILE]
#     [--outer-face=F] [--seed=S]
#
# FILE is the output of k3_surface_polygons.jl, by default where it writes it.
# Prints the two checks in Markdown.
###############################################################################

include(joinpath(@__DIR__, "k3_surface_polygons.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "to_tex", "curve_portraits.jl"))
using LinearAlgebra

"The area of the convex polygon with the vertices `verts`, in floating point."
function k3_area(verts)
    V = [Float64.(v) for v in verts]
    m = sum(V) / length(V)
    a = normalize(V[1] - m)
    b = normalize(argmax(norm, [(v - m) - dot(v - m, a) * a for v in V]))
    corners = sort([(dot(v - m, a), dot(v - m, b)) for v in V]; by=p -> atan(p[2], p[1]))
    return abs(sum(p[1] * q[2] - q[1] * p[2] for (p, q) in zip(corners, circshift(corners, -1)))) / 2
end

"The classes of the union-find `parent` after joining the elements of every group of `groups`."
function k3_classes(elements, groups)
    parent = Dict(i => i for i in elements)
    root(x) = parent[x] == x ? x : (parent[x] = root(parent[x]))
    foreach(g -> foreach(i -> parent[root(i)] = root(first(g)), g), groups)
    classes = Dict{Int,Vector{Int}}()
    foreach(i -> push!(get!(classes, root(i), Int[]), i), elements)
    return collect(values(classes))
end

percent(x) = @sprintf("%.1f%%", 100x)

"Print the two checks of the header for experiment `k` of `dir`, its polygons `data` of k3_surface_polygons.jl."
function k3_core_regions(dir, k, data; outer_face=1, seed=0)
    G = load_decorated_graphs(joinpath(dir, "decorated_graph_$k.mrdi")).subdivided
    _, number = cell_labels()
    cell_type(c) = "$c | $(number[c])"

    # The drawing of the atlas, its faces as sets of edges of the reduced graph, and its regions.
    curve = portrait_curve(joinpath(dir, "decorated_graph_$(k)_reduced.csv"))
    @req portrait_planar(curve) "the core of $dir is not planar"
    pos = portrait_core_drawing(curve.ends, curve.colour, curve.n, portrait_outer_faces(curve.ends, curve.colour, curve.n)[outer_face]; seed)
    neighbours = [Int[] for _ in 1:curve.n]
    foreach(((u, v),) -> (push!(neighbours[u], v); push!(neighbours[v], u)), curve.ends)
    faces = portrait_faces(pos, neighbours)
    outer = argmin(f -> portrait_area(pos[f]), faces)
    edge_of = Dict(minmax(u, v) => curve.core[i].edge for (i, (u, v)) in enumerate(curve.ends))
    face_edges = [Set(edge_of[minmax(a, b)] for (a, b) in zip(f, circshift(f, -1))) for f in faces]
    regions = portrait_region_pieces(pos, curve.ends, curve.colour)

    # The edges of X18 on the core: those of the pieces of the edges of the reduced graph.
    piece_edge = Dict(only(e.pieces) => i for (i, e) in enumerate(G.edges) if !is_leg(e))
    core_key_edge = Dict{Any,Int}()
    for line in readlines(joinpath(dir, "decorated_graph_$(k)_reduced.csv"))[2:end]
        fields = split(line, ',')
        edge = parse(Int, fields[1])
        any(es -> edge in es, face_edges) || continue
        foreach(p -> core_key_edge[data.curve_keys[piece_edge[p]]] = edge, parse.(Int, split(fields[7], ';')))
    end
    core_cells = sort!(unique(key[2] for key in keys(core_key_edge)))

    bounded = [(c, S, p) for ((c, S), p) in data.polygons if p.bounded]
    areas = [k3_area(p.verts) for (_, _, p) in bounded]
    total = sum(areas)
    holders = Dict{Any,Vector{Int}}()
    foreach(((i, (_, _, p)),) -> foreach(key -> push!(get!(holders, key, Int[]), i), p.keys), enumerate(bounded))
    vertices_count = length(unique(v for (_, _, p) in bounded for v in p.verts))
    println("X18 has $(length(data.polygons)) polygons. Its bounded part has $vertices_count vertices, $(length(holders)) edges ",
            "and $(length(bounded)) polygons.\n")

    println("## Check 1: the bounded part and the faces of the drawing\n")
    by_cell = Dict{Int,Vector{Int}}()
    foreach(i -> push!(get!(by_cell, bounded[i][1], Int[]), i), eachindex(bounded))
    share(is) = sum(areas[is]) / total
    in_core(is) = sum((areas[i] for i in is if bounded[i][1] in core_cells); init=0.0) / sum(areas[is])
    println("The bounded part meets $(length(by_cell)) maximal cells; the $(length(core_cells)) cells of the core hold ",
            percent(in_core(eachindex(bounded))), " of its area. The cells holding at least 1% of it:\n")
    println("| cell | type | polygons | area | |\n|---:|---:|---:|---:|---|")
    for (c, is) in sort(collect(by_cell); by=x -> -share(x[2]))
        share(is) >= 0.01 && println("| $(cell_type(c)) | $(length(is)) | $(percent(share(is))) | $(c in core_cells ? "core" : "") |")
    end
    pieces = k3_classes(eachindex(bounded), [is for (key, is) in holders if !haskey(core_key_edge, key)])
    println("\nCut along the $(length(core_key_edge)) edges of the core, it falls into $(length(pieces)) pieces. The drawing has ",
            "$(length(faces)) faces, the outer one $(findfirst(f -> f === outer, faces)).\n")
    println("| polygons | area | core edges on its border | faces holding them all | area in the cells of the core |\n|---:|---:|---:|---|---:|")
    for is in sort(pieces; by=is -> -share(is))
        border = Set(core_key_edge[key] for i in is for key in bounded[i][3].keys if haskey(core_key_edge, key))
        holding = [f for f in eachindex(faces) if issubset(border, face_edges[f])]
        println("| $(length(is)) | $(percent(share(is))) | $(length(border)) | $(isempty(holding) ? "none" : join(holding, ", ")) | ",
                percent(in_core(is)), " |")
    end
    sheets = [count(i -> bounded[i][1] == key[2], get(holders, key, Int[])) for key in keys(core_key_edge)]
    println("\nAlong every edge of the core, 5 polygons of X18 meet; the bounded ones among them: ",
            join(("$s^$(count(==(s), sheets))" for s in sort!(unique(sheets))), " "), ".\n")
    println("| ring | polygons | area | in the cells of the core |\n|---:|---:|---:|---:|")
    ring, frontier, depth = fill(-1, length(bounded)), [i for i in eachindex(bounded) if any(key -> haskey(core_key_edge, key), bounded[i][3].keys)], 0
    ring[frontier] .= 0
    while !isempty(frontier)
        println("| $depth | $(length(frontier)) | $(percent(share(frontier))) | $(percent(in_core(frontier))) |")
        next = unique([j for i in frontier for key in bounded[i][3].keys for j in holders[key] if ring[j] == -1])
        ring[next] .= depth + 1
        frontier, depth = next, depth + 1
    end

    println("\n## Check 2: the parts in the cells of the core and the regions of the atlas\n")
    region_area = Dict(c => sum(abs(portrait_area(poly)) for (d, poly) in regions if d == c) for c in core_cells)
    core_area = sum(sum(areas[by_cell[c]]) for c in core_cells)
    println("| cell | type | bounded polygons | pieces | core edges, in one piece | area on X18 | area of the region |\n|---:|---:|---:|---:|---|---:|---:|")
    for c in core_cells
        is = by_cell[c]
        parts = k3_classes(is, [filter(i -> bounded[i][1] == c, js) for js in values(holders)])
        edges = count(key -> key[2] == c, keys(core_key_edge))
        together = count(part -> any(i -> any(key -> haskey(core_key_edge, key), bounded[i][3].keys), part), parts) == 1
        println("| $(cell_type(c)) | $(length(is)) | $(length(parts)) | $edges, $(together ? "yes" : "no") | ",
                percent(sum(areas[is]) / core_area), " | ", percent(region_area[c] / sum(values(region_area))), " |")
    end
    # The length of a wall edge, from the two vertices it shares with the polygons on it.
    function edge_length(key)
        is = holders[key]
        common = collect(intersect((Set(bounded[i][3].verts) for i in is)...))
        length(common) == 2 ? norm(Float64.(common[1]) - Float64.(common[2])) : 0.0
    end
    walls = Dict{Tuple{Int,Int},Float64}()
    for (key, is) in holders
        key[1] == :wall || continue
        cs = filter(in(core_cells), unique(bounded[i][1] for i in is))
        length(cs) >= 2 || continue
        ℓ = edge_length(key)
        foreach(pair -> walls[pair] = get(walls, pair, 0.0) + ℓ, ((a, b) for a in cs for b in cs if a < b))
    end
    rounded(p) = round.(p; digits=7)
    sides = Dict{Any,Vector{Int}}()
    foreach(((c, poly),) -> foreach(((a, b),) -> push!(get!(sides, Set([rounded(a), rounded(b)]), Int[]), c), zip(poly, circshift(poly, -1))),
            regions)
    borders = Dict{Tuple{Int,Int},Float64}()
    for (side, cs) in sides
        length(unique(cs)) == 2 || continue
        a, b = collect(side)
        pair = Tuple(sort(unique(cs)))
        borders[pair] = get(borders, pair, 0.0) + hypot((a .- b)...)
    end
    println("\n| cells | wall on X18 | border in the atlas |\n|---|---:|---:|")
    for pair in sort(unique([collect(keys(walls)); collect(keys(borders))]); by=p -> -get(walls, p, 0.0))
        println("| $(pair[1]), $(pair[2]) | ", haskey(walls, pair) ? percent(walls[pair] / sum(values(walls))) : "none", " | ",
                haskey(borders, pair) ? percent(borders[pair] / sum(values(borders))) : "none", " |")
    end
    apart = [(a, b) for a in core_cells for b in core_cells if a < b && !haskey(walls, (a, b)) && !haskey(borders, (a, b))]
    println("\nMeeting neither on X18 nor in the atlas: ", isempty(apart) ? "none" : join(("$a and $b" for (a, b) in apart), ", "), ".")
end

if abspath(PROGRAM_FILE) == @__FILE__
    positional = filter(arg -> !startswith(arg, "--"), ARGS)
    length(positional) in (1, 2) || error("usage: k3_core_regions.jl DIR [K] [--polygons=FILE] [--outer-face=F] [--seed=S]")
    dir = abspath(positional[1])
    k = length(positional) == 2 ? parse(Int, positional[2]) : 1
    option(name) = (given = [String(split(arg, '='; limit=2)[2]) for arg in ARGS if startswith(arg, "--$name=")]; isempty(given) ? nothing : only(given))
    polygons = something(option("polygons"), project_path("results", "curves", "genus_10", "k3", basename(dir) * "_x18.mrdi"))
    k3_core_regions(dir, k, load_k3_polygons(polygons); outer_face=parse(Int, something(option("outer-face"), "1")),
                    seed=parse(Int, something(option("seed"), "0")))
end
