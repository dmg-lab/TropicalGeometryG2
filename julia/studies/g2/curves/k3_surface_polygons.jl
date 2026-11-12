###############################################################################
# The tropical K3 surface X18 = trop(G2/P2) ∩ (U(18,21) + t) that contains the
# genus-10 curve X17 = trop(G2/P2) ∩ (U(17,21) + t) of experiment K of an
# archive DIR of the Method 1 runs (decorated_graph_K.mrdi of
# decorated_curve_graph.jl, translation_K.mrdi), polygon by polygon, exactly.
#
# A polygon is C ∩ (B_S + t), C a maximal cell of the coarsest fan and
# B_S = {x : x_S constant, x_j ≥ x_S off S}, |S| = 4, taken in the slice
# sum(x) = 0. The search starts from the five sheets C ∩ (B_S + t), S ⊂ T,
# along every bounded edge C ∩ (B_T + t), |T| = 5, of the curve, and crosses
# every edge of a polygon: along a wall of C to the other cells around the
# wall, with the same S, and along an edge of the curve, where another
# coordinate k attains the minimum of x - t, to the other 4-subsets of
# S ∪ {k} in C. X18 is connected, so the search finds all of it.
#
#   julia --project=julia julia/studies/g2/curves/k3_surface_polygons.jl DIR [K] [--output=FILE]
#
# Writes FILE, by default julia/results/curves/genus_10/k3/DIR_x18.mrdi: every
# polygon (cell, S) with its edges as keys, (:wall, wall, S) or
# (:curve, cell, T), whether it is bounded, and its vertices; the key of every
# bounded edge of the subdivided graph of the curve; and the walls of the
# fan, numbered, with the cells around them.
###############################################################################

include(joinpath(@__DIR__, "decorated_curve_graph.jl"))
import Combinatorics
using Printf

"""
The polygons of X18 for the curve `G` (a subdivided graph of
decorated_curve_graph.jl) and the translation `t`, see the header.
"""
function k3_polygons(G, t)
    fan = load(data_path("g2", "grassmannian", "02_coarsest_fan.mrdi"))
    cells = collect(maximal_cones(fan))
    n = ambient_dim(fan)
    inequalities = [Matrix{QQFieldElem}(linear_inequality_matrix(facets(C))) for C in cells]
    equations = [Matrix{QQFieldElem}(linear_equation_matrix(linear_span(C))) for C in cells]
    # The cells are simplicial modulo a lineality of dimension 3, without rays to name their walls: a wall is a
    # facet of a cell, found again among those of the same linear span.
    facet_cones = [faces(C, dim(C) - 1) for C in cells]
    walls, by_span = Cone{QQFieldElem}[], Dict{Matrix{QQFieldElem},Vector{Int}}()
    wall_of = map(facet_cones) do Fs
        map(Fs) do F
            span = Matrix{QQFieldElem}(rref(matrix(QQ, Matrix{QQFieldElem}(linear_equation_matrix(linear_span(F)))))[2])
            bucket = get!(by_span, span, Int[])
            found = findfirst(w -> walls[w] == F, bucket)
            isnothing(found) || return bucket[found]
            push!(walls, F)
            push!(bucket, length(walls))
            return length(walls)
        end
    end
    cells_of_wall = Dict{Int,Vector{Int}}()
    foreach(c -> foreach(w -> push!(get!(cells_of_wall, w, Int[]), c), wall_of[c]), eachindex(cells))
    println(length(walls), " walls, cells around a wall: ",
            join(("$k^$(count(v -> length(v) == k, values(cells_of_wall)))" for k in sort!(unique(length.(values(cells_of_wall))))), " "))

    # The polygon C ∩ (B_S + t) in the slice, or nothing when it is not 2-dimensional.
    function polygon(c, S)
        s = first(S)
        others = [k for k in 1:n if !(k in S)]
        rows = [inequalities[c]; vcat(([j == s ? 1 : j == k ? -1 : 0 for j in 1:n]' for k in others)...)]
        rhs = [zeros(QQ, size(inequalities[c], 1)); [t[s] - t[k] for k in others]]
        eq = [equations[c]; vcat(([j == i ? 1 : j == s ? -1 : 0 for j in 1:n]' for i in S if i != s)...); ones(Int, n)']
        eqrhs = [zeros(QQ, size(equations[c], 1)); [t[i] - t[s] for i in S if i != s]; 0]
        P = polyhedron((Matrix{QQFieldElem}(rows), rhs), (Matrix{QQFieldElem}(eq), eqrhs))
        return dim(P) == 2 ? P : nothing
    end
    # The edges of a polygon by what holds them: another coordinate at the minimum, or a wall of its cell.
    function edge_keys(c, S, P)
        s = first(S)
        return map(filter(F -> n_vertices(F) > 0, faces(P, 1))) do F
            # A point inside the edge: the middle of a segment, or a ray one step from its vertex. Faces
            # without vertices are those at infinity between two rays.
            ends = [collect(v) for v in vertices(F)]
            x = length(ends) == 2 ? (ends[1] .+ ends[2]) ./ 2 : ends[1] .+ collect(first(rays(F)))
            tight = [k for k in 1:n if !(k in S) && x[k] - t[k] == x[s] - t[s]]
            on = [j for (j, W) in enumerate(facet_cones[c]) if x in W]
            @req length(tight) + length(on) == 1 "an edge of ($c, $S) lies on the coordinates $tight and the facets $on"
            isempty(tight) ? (:wall, wall_of[c][only(on)], S) : (:curve, c, sort([S; only(tight)]))
        end
    end

    curve_keys = Dict{Int,Any}()
    queue = Tuple{Int,Vector{Int}}[]
    for (k, e) in enumerate(G.edges)
        is_leg(e) && continue
        d = (G.points[e.tail] .+ G.points[e.head]) ./ 2 .- t
        T = findall(==(minimum(d)), d)
        @req length(T) == 5 "edge $k of the curve attains the minimum of x - t $(length(T)) times"
        curve_keys[k] = (:curve, only(e.cells), T)
        foreach(S -> push!(queue, (only(e.cells), sort(S))), Combinatorics.combinations(T, 4))
    end
    polygons, flat = Dict{Tuple{Int,Vector{Int}},Any}(), Set{Tuple{Int,Vector{Int}}}()
    started = time()
    while !isempty(queue)
        c, S = pop!(queue)
        (haskey(polygons, (c, S)) || (c, S) in flat) && continue
        P = polygon(c, S)
        isnothing(P) && (push!(flat, (c, S)); continue)
        keys = edge_keys(c, S, P)
        polygons[(c, S)] = (; keys, bounded=is_bounded(P), verts=[collect(v) for v in vertices(P)])
        for key in keys
            key[1] == :curve ? foreach(S2 -> push!(queue, (c, sort(S2))), Combinatorics.combinations(key[3], 4)) :
                               foreach(c2 -> push!(queue, (c2, S)), cells_of_wall[key[2]])
        end
        length(polygons) % 500 == 0 && @printf("%d polygons, %.0f s\n", length(polygons), time() - started)
    end
    @printf("%d polygons, %d bounded; %d pairs met across an edge were not 2-dimensional; %.0f s\n", length(polygons),
            count(p -> p.bounded, values(polygons)), length(flat), time() - started)
    return (; polygons=collect(polygons), curve_keys, cells_of_wall)
end

if abspath(PROGRAM_FILE) == @__FILE__
    positional = filter(arg -> !startswith(arg, "--"), ARGS)
    length(positional) in (1, 2) || error("usage: k3_surface_polygons.jl DIR [K] [--output=FILE]")
    dir = abspath(positional[1])
    k = length(positional) == 2 ? parse(Int, positional[2]) : 1
    given = [String(split(arg, '='; limit=2)[2]) for arg in ARGS if startswith(arg, "--output=")]
    output = isempty(given) ? project_path("results", "curves", "genus_10", "k3", basename(dir) * "_x18.mrdi") : only(given)
    G = load_decorated_graphs(joinpath(dir, "decorated_graph_$k.mrdi")).subdivided
    t = QQ.(load(joinpath(dir, "translation_$k.mrdi")))
    mkpath(dirname(abspath(output)))
    save_k3_polygons(output, (; k3_polygons(G, t)..., t))
    println("Wrote ", output)
end
