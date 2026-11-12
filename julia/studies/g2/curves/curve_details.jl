###############################################################################
# curve_details: the core and the points of a curve of
# hyperplane_intersection_census.jl, in the directory of its experiment.
#
# The curve is the stored curve.mrdi, in the slice x_1 = 0 of the lineality
# (1,...,1), whose lattice is ZZ^21 there; lengths are lattice lengths. Every
# edge of the core (curve_core.csv) gets its length, the number of its
# segments and their multiplicities, and the maximal cells of the coarsest
# tropical G2-Grassmannian along it; every vertex of the core its point and
# the cells containing it. Every point of the fifth cut (points.mrdi) gets its
# coordinates, its multiplicity, the two coordinates where the minimum of
# x - t over the fifth hyperplane, of shift t, is attained, and its place on
# the curve: on the core, at a distance from the tail of an edge of the core,
# or on a pending tree or its legs, at a distance from the point of the core
# where the tree hangs, and its edge or leg of curve_reduced.csv, with the
# distance from its tail and its length (none for a leg). The graph of the
# core gets its girth, its diameter, its connectivity and the order of its
# automorphism group.
#
#   julia --project=julia julia/studies/g2/curves/curve_details.jl EXPERIMENT_DIR
#
# writes core_edges.csv, core_vertices.csv and point_places.csv in
# EXPERIMENT_DIR and prints a summary.
###############################################################################

include(joinpath(@__DIR__, "hyperplane_intersection_census.jl"))

"The lattice length of the segment from `a` to `b` in the slice x_1 = 0, whose lattice is ZZ^n there."
function lattice_length(a, b)
    d = b .- a
    all(iszero, d) && return QQ(0)
    p = primitive_direction(d)
    i = findfirst(!iszero, p)
    return d[i] // p[i]
end

"The pieces of the edge `e` of the graph record `G`, in their order from its tail, with the point where each starts."
function ordered_pieces(G, e, list)
    at, remaining, ordered = G.points[G.tails[e]], copy(G.pieces[e]), Tuple{Int,Vector{QQFieldElem}}[]
    while !isempty(remaining)
        i = findfirst(k -> any(==(at), collect.(vertices(list[k]))), remaining)
        k = popat!(remaining, i)
        push!(ordered, (k, at))
        others = filter(!=(at), collect.(vertices(list[k])))
        isempty(others) || (at = only(others))
    end
    return ordered
end

"`p_ij` for the coordinate `k` of QQ^21, in the order 12, 13, ..., 67."
coordinate_name(k) = "p" * join(PLUCKER_PAIRS[k])

"Girth, diameter, vertex connectivity (at most 3) and automorphism group order of the graph of `n` vertices and the edges `pairs`."
function graph_invariants(n, pairs)
    around = [Int[] for _ in 1:n]
    foreach(((u, v),) -> (push!(around[u], v); push!(around[v], u)), pairs)
    function distances(s; without=())
        d = fill(-1, n)
        s in without && return d
        d[s] = 0
        queue = [s]
        for u in queue, v in around[u]
            v in without || d[v] >= 0 || (d[v] = d[u] + 1; push!(queue, v))
        end
        return d
    end
    # The shortest cycle through an edge (u, v) has length 1 + the distance from u to v without that edge.
    girth = minimum(pairs) do (u, v)
        d = fill(-1, n)
        d[u] = 0
        queue = [u]
        for a in queue, b in around[a]
            (a, b) == (u, v) && continue
            d[b] >= 0 || (d[b] = d[a] + 1; push!(queue, b))
        end
        d[v] < 0 ? typemax(Int) : d[v] + 1
    end
    diameter = maximum(maximum(distances(s)) for s in 1:n)
    connected_without(removed) = (rest = setdiff(1:n, removed); all(>=(0), distances(first(rest); without=removed)[rest]))
    connectivity = any(v -> !connected_without([v]), 1:n) ? 1 :
                   any(!connected_without([u, v]) for u in 1:n for v in u+1:n) ? 2 : 3
    group = automorphism_group(graph_from_edges(Undirected, [[u, v] for (u, v) in pairs]))
    return (; girth, diameter, connectivity, automorphisms=order(group))
end

function curve_details(directory)
    file(name) = joinpath(directory, name)
    curve, points, hyperplanes = load(file("curve.mrdi")), load(file("points.mrdi")), load(file("hyperplanes.mrdi"))
    list, S, C = curve.polyhedra, curve.subdivided, curve.core
    cells = collect(maximal_cones(load(data_path("g2", "grassmannian", "02_coarsest_fan.mrdi"))))
    orbit, number = cell_labels()
    cell_names(cs) = join(("$c ($(orbit[c])/$(number[c]))" for c in cs), ' ')

    # The edges of the core.
    core_piece = Dict(k => e for e in eachindex(C.tails) for k in C.pieces[e])
    rows = ["edge,tail,head,length,segments,multiplicities,cells"]
    lengths = map(eachindex(C.tails)) do e
        ordered = ordered_pieces(C, e, list)
        len = sum(lattice_length(collect.(vertices(list[k]))...) for (k, _) in ordered)
        push!(rows, "$e,$(C.tails[e]),$(C.heads[e]),$len,$(length(ordered))," *
                    "$(join(unique(curve.multiplicities[k] for (k, _) in ordered), ';')),$(join(unique(C.cells[e]), ';'))")
        len
    end
    write(file("core_edges.csv"), join(rows, '\n') * '\n')

    # The vertices of the core.
    rows = ["vertex,cells,orbits/types," * join(coordinate_name.(1:21), ',')]
    for (v, p) in enumerate(C.points)
        around = containing_cells(convex_hull(permutedims(p)), cells)
        push!(rows, "$v,$(join(around, ';')),$(join(("$(orbit[c])/$(number[c])" for c in around), ';')),$(join(p, ','))")
    end
    write(file("core_vertices.csv"), join(rows, '\n') * '\n')

    # The places of the points on the curve: the subdivided graph without the core is a forest of pending trees.
    vertex_of = Dict(p => i for (i, p) in enumerate(S.points))
    on_core = Set(vertex_of[collect(v)] for k in keys(core_piece) for v in vertices(list[k]))
    edge_of = Dict(only(S.pieces[e]) => e for e in eachindex(S.tails))
    hanging = Dict{Int,Vector{Tuple{Int,Int}}}()
    for e in eachindex(S.tails)
        only(S.pieces[e]) in keys(core_piece) && continue
        for (a, b) in ((S.tails[e], S.heads[e]), (S.heads[e], S.tails[e]))
            a == 0 || push!(get!(hanging, a, Tuple{Int,Int}[]), (b, e))
        end
    end
    "The point of the core where the tree through the vertex `w` of `S` hangs, from `w` away from the edge `from`, and the length to it."
    function to_core(w, from)
        w in on_core && return w, QQ(0)
        for (u, e) in hanging[w]
            (e == from || u == 0) && continue
            found = to_core(u, e)
            isnothing(found) || return found[1], found[2] + lattice_length(S.points[w], S.points[u])
        end
        return nothing
    end
    "The place of the point `p` of the core: a vertex of the core, or a distance from the tail of an edge of it."
    function core_place(p)
        v = findfirst(==(p), C.points)
        isnothing(v) || return "core vertex $v"
        for e in eachindex(C.tails), (k, start) in ordered_pieces(C, e, list)
            ends = collect.(vertices(list[k]))
            p in ends || continue
            ordered = ordered_pieces(C, e, list)
            before = sum((lattice_length(collect.(vertices(list[j]))...) for (j, _) in ordered[1:findfirst(==((k, start)), ordered)-1]); init=QQ(0))
            return "core edge $e at $(before + lattice_length(start, p)) of $(lengths[e]) from vertex $(C.tails[e])"
        end
    end
    R = curve.reduced
    reduced_edge = Dict(k => e for e in eachindex(R.tails) for k in R.pieces[e])
    "The edge or leg of `R` through the point `x` of the piece `k`, the length to `x` from its tail, and its length."
    function reduced_place(x, k)
        e = reduced_edge[k]
        ordered = ordered_pieces(R, e, list)
        j = findfirst(o -> o[1] == k, ordered)
        before = sum((lattice_length(collect.(vertices(list[o[1]]))...) for o in ordered[1:j-1]); init=QQ(0))
        total = R.heads[e] == 0 ? "" : sum(lattice_length(collect.(vertices(list[o[1]]))...) for o in ordered)
        return "$e,$(before + lattice_length(ordered[j][2], x)),$total"
    end
    t = hyperplanes.shifts[N_HYPERPLANES]
    rows = ["point,multiplicity,minimum_at,place,reduced_edge,from_tail,length," * join(coordinate_name.(1:21), ',')]
    summary = String[]
    for (i, x) in enumerate(points.points)
        # The pieces of the curve the point lies on, by the cut that made it.
        k = something(findfirst(k -> x in list[k], points.curve_cells[i]), Some(nothing))
        @req !isnothing(k) "point $i is on none of the pieces $(points.curve_cells[i]) of the curve"
        k = points.curve_cells[i][k]
        values = x .- t
        least = findall(==(minimum(values)), values)
        @req length(least) >= 2 "point $i is not on the fifth hyperplane"
        place = if haskey(core_piece, k)
            e = core_piece[k]
            ordered = ordered_pieces(C, e, list)
            j = findfirst(o -> o[1] == k, ordered)
            before = sum((lattice_length(collect.(vertices(list[o[1]]))...) for o in ordered[1:j-1]); init=QQ(0))
            "on core edge $e at $(before + lattice_length(ordered[j][2], x)) of $(lengths[e]) from vertex $(C.tails[e])"
        else
            e = edge_of[k]
            tail, head = S.tails[e], S.heads[e]
            # Towards the core through the tail, or else through the head of the piece.
            found = to_core(tail, e)
            through = tail
            if isnothing(found) || head != 0 && !(tail in on_core) && isnothing(to_core(tail, e))
                found, through = to_core(head, e), head
            end
            w, len = found
            kind = head == 0 ? "leg" : "pending tree"
            "on a $kind at $(len + lattice_length(S.points[through], x)) from $(core_place(S.points[w]))"
        end
        push!(rows, "$i,$(points.multiplicities[i]),$(join(coordinate_name.(least), ';')),$place,$(reduced_place(x, k)),$(join(x, ','))")
        push!(summary, "point $i (multiplicity $(points.multiplicities[i]), min at $(join(coordinate_name.(least), '='))): $place")
    end
    write(file("point_places.csv"), join(rows, '\n') * '\n')

    invariants = graph_invariants(length(C.points), collect(zip(C.tails, C.heads)))
    println("Core: $(length(C.points)) vertices, $(length(C.tails)) edges, total length $(sum(lengths)); ",
            "edge lengths $(join(sort(lengths), ' '))")
    println("Core graph: girth $(invariants.girth), diameter $(invariants.diameter), ",
            "vertex connectivity $(invariants.connectivity), $(invariants.automorphisms) automorphisms")
    foreach(println, summary)
    println("Wrote core_edges.csv, core_vertices.csv and point_places.csv in $directory")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: curve_details.jl EXPERIMENT_DIR")
    curve_details(abspath(only(ARGS)))
end
