###############################################################################
# Is the core of a genus-10 curve trivalent? The curve is the intersection of
# the tropical G2-Grassmannian with a translate B + t of the uniform Bergman
# fan, computed exactly cell by cell with `exact_piece` of
# bounded_cell_valences.jl, after a pass in floating point that keeps, with a
# margin, the subsets T whose piece may be more than a point. Every vertex
# inside a cell has 6 branches, so a vertex of the core is trivalent when 3 of
# its branches are pending trees.
#
# Kinds of translations, from a given translation t0, by default that of
# julia/curves_g10/uniform_shifts/translation_3.mrdi (experiment 3 of the
# uniform shifts), and the span W of the fan:
#
#   normal    the part of t0 in W plus a random part orthogonal to W;
#   wpart     the part of t0 orthogonal to W plus a random part in W;
#   fixed     a fixed random part in W plus a random part orthogonal to W;
#   uniform   uniformly random, as in mukai_curves/run_experiments.jl.
#
#   julia --project=julia -t 24 julia/studies/g2/curves/core_trivalence.jl KIND TRIALS
#         [--translation=FILE] [--batch=N] [--output=DIR]
#   julia --project=julia -t 3 julia/studies/g2/curves/core_trivalence.jl --check=FILE
#
# Trial k draws its random part from Xoshiro(k), so a run can be resumed or
# extended: the trials already in DIR/core_trivalence_KIND.csv are skipped.
# The curves with a trivalent core or of genus 10 are saved as
# DIR/KIND_K/decorated_graph_K.mrdi with the CSV files of
# decorated_curve_graph.jl, for build_curve_atlas.jl. DIR defaults to
# results/curves/genus_10/core_trivalence. The second form only describes the
# curve of the translation in FILE.
###############################################################################

include("bounded_cell_valences.jl")
include("decorated_curve_graph.jl")
import LinearAlgebra

const FLOAT_MARGIN = 1e-6
const CORE_KINDS = ("normal", "wpart", "fixed", "uniform")

# polymake, behind `exact_piece` for the subsets of deficient rank, is not thread-safe.
const POLYMAKE_LOCK = ReentrantLock()

"""
For every translation of `shifts`, the subsets of `subsets` whose piece of the
cell of span `Mf` and facets `AMf`, in floating point, may be more than a
point, with a margin relative to the size of the translation.
"""
function candidate_subsets(Mf, AMf, subsets, shifts)
    found = [Vector{Int}[] for _ in shifts]
    margins = [FLOAT_MARGIN * (1 + maximum(abs, t)) for t in shifts]
    for T in subsets
        k = T[1]
        E = vcat(Mf[T[2:end], :] .- Mf[k:k, :], Mf[k:k, :])
        F = LinearAlgebra.svd(E; full=true)
        r = count(>(1e-9 * F.S[1]), F.S)
        P = F.V[:, 1:r] * LinearAlgebra.Diagonal(1 ./ F.S[1:r]) * F.U[:, 1:r]'
        xd, ad = Mf * F.V[:, 6], AMf * F.V[:, 6]
        for (n, t) in enumerate(shifts)
            f = [t[T[2:end]] .- t[k]; t[k]]
            y0 = P * f
            if r < 5
                LinearAlgebra.norm(E * y0 - f) <= margins[n] && push!(found[n], T)
                continue
            end
            x0, a0 = Mf * y0, AMf * y0
            lo, hi, possible = -Inf, Inf, true
            # The constraint c + s g ≥ 0 on the parameter s of the line, relaxed by the margin.
            function constrain!(c, g)
                c += margins[n]
                if abs(g) < 1e-12
                    possible &= c >= 0
                elseif g > 0
                    lo = max(lo, -c / g)
                else
                    hi = min(hi, -c / g)
                end
            end
            for j in 1:21
                j in T || constrain!((x0[j] - t[j]) - (x0[k] - t[k]), xd[j] - xd[k])
            end
            foreach(i -> constrain!(-a0[i], -ad[i]), eachindex(a0))
            possible && lo <= hi && push!(found[n], T)
        end
    end
    return found
end

"The curve of the translation `t` over `cells` from the candidate subsets `found[c]` of each cell `c`, and its number of pieces of dimension at least 2."
function exact_curve(cells, found, t::Vector{QQFieldElem})
    points, index = Vector{QQFieldElem}[], Dict{Vector{QQFieldElem},Int}()
    vertex!(v) = get!(() -> (push!(points, v); length(points)), index, v)
    edges, seen, higher = DecoratedEdge[], Set{Any}(), 0
    for (c, (M, A)) in enumerate(cells)
        AM = A * M
        for T in found[c]
            # Only the subsets of deficient rank reach polymake; they are few.
            piece = rank(vcat([M[i:i, :] - M[T[1]:T[1], :] for i in T[2:end]]..., M[T[1]:T[1], :])) == 5 ?
                exact_piece(M, A, AM, t, T) : lock(() -> exact_piece(M, A, AM, t, T), POLYMAKE_LOCK)
            isnothing(piece) && continue
            first(piece) == :higher && (higher += 1; continue)
            a = vertex!(piece[2])
            b = first(piece) == :segment ? vertex!(piece[3]) : 0
            key = b == 0 ? (a, piece[3]) : minmax(a, b)
            key in seen && continue
            push!(seen, key)
            push!(edges, DecoratedEdge((a, b, [c], [length(edges) + 1])))
        end
    end
    return (graph=DecoratedGraph((points, edges)), higher=higher)
end

"""
The curves `C ∩ (B + t)` over the cells `cells`, pairs of a span and facets,
for the translations `shifts`, as decorated graphs, each with the number of its
pieces of dimension at least 2.
"""
function curves(cells, shifts::Vector{Vector{QQFieldElem}}; approximations=nothing)
    subsets = collect(Combinatorics.combinations(1:21, 5))
    floats = [Float64.(t) for t in shifts]
    approximations = something(approximations, [(Float64.(Matrix(M)), Float64.(Matrix(A * M))) for (M, A) in cells])
    found = Vector{Vector{Vector{Vector{Int}}}}(undef, length(cells))
    Threads.@threads for c in eachindex(cells)
        found[c] = candidate_subsets(approximations[c]..., subsets, floats)
    end
    results = Vector{Any}(undef, length(shifts))
    Threads.@threads for n in eachindex(shifts)
        results[n] = exact_curve(cells, [found[c][n] for c in eachindex(cells)], shifts[n])
    end
    return results
end

"The valences of the vertices of `G`, legs included."
function graph_valences(G::DecoratedGraph)
    valence = zeros(Int, length(G.points))
    foreach(e -> (valence[e.tail] += 1; is_leg(e) || (valence[e.head] += 1)), G.edges)
    return valence
end

"""
A row of the table for the curve `G`: its size and genus, the valences of its
core, and those of the vertices of the core inside a cell of a bounded orbit.
"""
function describe_curve(G::DecoratedGraph, higher, orbit)
    C = desubdivide(core(G))
    where = Dict(p => i for (i, p) in enumerate(G.points))
    cells_at = [Int[] for _ in G.points]
    foreach(e -> (push!(cells_at[e.tail], only(e.cells)); is_leg(e) || push!(cells_at[e.head], only(e.cells))), G.edges)
    valence = graph_valences(C)
    bounded = [valence[v] for (v, p) in enumerate(C.points)
               if allequal(cells_at[where[p]]) && orbit[first(cells_at[where[p]])] in BOUNDED_ORBITS]
    legs = count(is_leg, G.edges)
    return (genus=genus(G), vertices=length(G.points), edges=length(G.edges) - legs, legs=legs,
            core="$(length(C.points))/$(length(C.edges))", core_valences=histogram(valence),
            trivalent=!isempty(valence) && all(==(3), valence), bounded=histogram(bounded),
            curve_valences=histogram(graph_valences(G)), higher=higher)
end

const CORE_COLUMNS = (:genus, :vertices, :edges, :legs, :core, :core_valences, :trivalent, :bounded, :curve_valences, :higher)

"Saves the graphs of the curve `G` of trial `k` as decorated_curve_graph.jl does, in `dir`."
function save_curve(dir, k, G::DecoratedGraph, orbit, number)
    mkpath(dir)
    graphs = (subdivided=G, desubdivided=desubdivide(G), reduced=desubdivide(G; within_cells=true),
              core=desubdivide(core(G)))
    save_decorated_graphs(joinpath(dir, "decorated_graph_$(k).mrdi"), graphs)
    foreach(name -> write_decorated_graph(joinpath(dir, "decorated_graph_$(k)_$(name).csv"), graphs[name], orbit, number),
            (:desubdivided, :reduced, :core))
end

"The translation of trial `k` of `kind` from the translation `t0`, as a function of `k`, with the span of `fan`."
function trial_translations(fan, kind, t0)
    span = vcat(matrix(QQ, reduce(vcat, permutedims.(collect.(rays_modulo_lineality(fan).rays_modulo_lineality)))),
                matrix(QQ, lineality_space(fan)))
    normals = kernel(span; side=:right)
    project(v) = v - normals * solve(transpose(normals) * normals, transpose(normals) * v; side=:right)
    t0 = matrix(QQ, 21, 1, QQ.(t0))
    fixed = project(matrix(QQ, 21, 1, rand(Random.Xoshiro(0), SHIFT_RANGE, 21)))
    function translation(k)
        rng = Random.Xoshiro(k)
        normal() = normals * matrix(QQ, 6, 1, rand(rng, SHIFT_RANGE, 6))
        t = kind == "normal" ? project(t0) + normal() :
            kind == "wpart" ? (t0 - project(t0)) + project(matrix(QQ, 21, 1, rand(rng, SHIFT_RANGE, 21))) :
            kind == "fixed" ? fixed + normal() :
            matrix(QQ, 21, 1, rand(rng, SHIFT_RANGE, 21))
        return Vector{QQFieldElem}(t[:, 1])
    end
    return translation
end

function main(args)
    option(name, default) = last([[default]; [split(a, '='; limit=2)[2] for a in args if startswith(a, "--$name=")]])
    positional = filter(a -> !startswith(a, "--"), args)
    fan = load(data_path("g2", "grassmannian", "02_coarsest_fan.mrdi"))
    cells = [cell_span(C) for C in maximal_cones(fan)]
    orbit, number = cell_labels()

    check = option("check", "")
    if !isempty(check)
        t = Vector{QQFieldElem}(QQ.(load(String(check))))
        d = describe_curve(only(curves(cells, [t])).graph, 0, orbit)
        println(join(("$c=$(d[c])" for c in CORE_COLUMNS), ", "))
        return
    end

    length(positional) == 2 && positional[1] in CORE_KINDS ||
        error("usage: core_trivalence.jl KIND TRIALS [--translation=FILE] [--batch=N] [--output=DIR], KIND in $CORE_KINDS")
    kind, trials = positional[1], parse(Int, positional[2])
    batch = parse(Int, option("batch", string(max(Threads.nthreads(), 1))))
    file = option("translation", project_path("curves_g10", "uniform_shifts", "translation_3.mrdi"))
    output = abspath(option("output", project_path("results", "curves", "genus_10", "core_trivalence")))
    mkpath(output)

    translation = trial_translations(fan, kind, load(String(file)))

    table = joinpath(output, "core_trivalence_$(kind).csv")
    done = isfile(table) ? Set(parse(Int, first(split(l, ','))) for l in readlines(table)[2:end]) : Set{Int}()
    isfile(table) || open(io -> println(io, join(["trial"; collect(string.(CORE_COLUMNS)); "translation"], ",")), table, "w")
    todo = [k for k in 1:trials if !(k in done)]
    println("$kind: $(length(todo)) of $trials trials to run, translation $file, $(Threads.nthreads()) threads, output $table")
    approximations = [(Float64.(Matrix(M)), Float64.(Matrix(A * M))) for (M, A) in cells]
    for range in Iterators.partition(todo, batch)
        shifts = translation.(range)
        elapsed = @elapsed results = curves(cells, shifts; approximations)
        open(table, "a") do io
            for (k, t, result) in zip(range, shifts, results)
                d = describe_curve(result.graph, result.higher, orbit)
                println(io, join([k; [d[c] for c in CORE_COLUMNS]; join(string.(t), ";")], ","))
                (d.trivalent || d.genus == 10) && save_curve(joinpath(output, "$(kind)_$(k)"), k, result.graph, orbit, number)
                d.trivalent && println("trial $k: trivalent core, genus $(d.genus), core $(d.core)")
                d.genus == 10 && println("trial $k: genus 10, core $(d.core), valences $(d.core_valences)")
            end
        end
        println("trials $(first(range))-$(last(range)): $(round(elapsed; digits=1)) s")
        flush(stdout)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
