###############################################################################
# The unbounded rays of a genus-10 curve, with their multiplicities, and its
# degree in every coordinate: for every coordinate i, the sum over the rays
# of m (u_i - min u), with u the primitive direction and m the multiplicity,
# is the number of points of the curve on {p_i = 0}, 18 for a curve section
# of G2/P2. Directions are modulo the diagonal, in the min convention.
#
# trop(G2/P2) lies in x12 = x67, x13 = x57, x17 = x56, x23 = x47, x27 = x46
# and x37 = x45 (checked on its rays and lineality), so a ray points along
# e_ab + e_cd for these pairs, never along one of them alone: a curve has at
# most 15 x 18 = 270 rays counted with multiplicity.
#
# Two kinds of DIR:
#   - an experiment of hyperplane_intersection_census.jl: the rays and the
#     multiplicities stored in DIR/curve.mrdi;
#   - a curve trop(G2/P2) ∩ (U(17,21) + t) of a certified run of
#     run_experiments.jl: the curve is recomputed exactly from
#     DIR/translation_K.mrdi as in core_trivalence.jl, every piece of
#     C ∩ (B_S + t) of multiplicity m_C [ZZ^21 : L_C + L_S], with m_C over
#     QQbar, L_C = span(C) ∩ ZZ^21 and L_S = {x_S constant}. The curve is
#     checked balanced and its rays against DIR/decorated_graph_K.mrdi.
#
#   julia --project=julia -t 8 julia/studies/g2/curves/curve_rays.jl DIR [K]
###############################################################################

include(joinpath(@__DIR__, "core_trivalence.jl"))
using TropicalCensus: study_path, csv_rows

"The primitive vector of the direction `d` in `ZZ^n`."
function primitive(d)
    d = d .* lcm(denominator.(d))
    return ZZ.(d .// gcd(numerator.(d)))
end

"The pairs of coordinates equal on every ray and on the lineality of `fan`."
function equal_coordinates(fan)
    vectors = [[collect(r) for r in rays_modulo_lineality(fan).rays_modulo_lineality];
               [collect(r) for r in lineality_space(fan)]]
    return [(i, j) for i in 1:21 for j in i+1:21 if all(v -> v[i] == v[j], vectors)]
end

"The multiplicities over QQbar of the maximal cells of the coarsest fan."
function cell_weights()
    raw = Dict(row["cone"] => parse(Int, row["mult_QQbar"]) for row in csv_rows(study_path(:g2, "weights", "weights.csv")))
    return [only(unique(raw[c] for c in split(row["original_cells"], ';')))
            for row in csv_rows(data_path("g2", "grassmannian", "02_coarsest_cells.csv"))]
end

"A basis of the lattice `span(M) ∩ ZZ^21`, as columns."
function saturated_lattice(M)
    N = Matrix(kernel(transpose(M); side=:right))
    return kernel(transpose(matrix(ZZ, ZZ.(N .* lcm(denominator.(N))))); side=:right)
end

"`[ZZ^21 : L + L_S]` with `L_S = {x_S constant}`: the index of the image of `L` under `x -> (x_s - x_S[1])`."
function lattice_index(L, S)
    D = matrix(ZZ, [j == S[i+1] ? 1 : j == S[1] ? -1 : 0 for i in 1:length(S)-1, j in 1:21])
    B = D * L
    @req rank(B) == length(S) - 1 "the cell is not transverse to B_S"
    s = snf(B)
    return Int(prod(s[i, i] for i in 1:length(S)-1))
end

"The rays of the census curve of `dir`, as pairs (primitive direction, multiplicity)."
function census_rays(dir)
    curve = load(joinpath(dir, "curve.mrdi"))
    println("stored curve: genus $(curve.genus), balanced $(curve.balanced)")
    return [(primitive(collect(only(rays(P)))), Int(m))
            for (P, m) in zip(curve.polyhedra, curve.multiplicities) if is_feasible(P) && !is_bounded(P)]
end

"The rays of the curve of translation `k` of `dir`, recomputed exactly with their multiplicities."
function bergman_rays(dir, k, fan)
    cells = [cell_span(C) for C in maximal_cones(fan)]
    weights = cell_weights()
    @req length(weights) == length(cells) "one weight per cell"
    t = Vector{QQFieldElem}(QQ.(load(joinpath(dir, "translation_$k.mrdi"))))
    subsets = collect(Combinatorics.combinations(1:21, 5))
    found = Vector{Vector{Vector{Int}}}(undef, length(cells))
    pieces = [Any[] for _ in cells]
    Threads.@threads for c in eachindex(cells)
        M, A = cells[c]
        found[c] = only(candidate_subsets(Float64.(Matrix(M)), Float64.(Matrix(A * M)), subsets, [Float64.(t)]))
        L, AM = saturated_lattice(M), A * M
        for T in found[c]
            piece = lock(() -> exact_piece(M, A, AM, t, T), POLYMAKE_LOCK)
            isnothing(piece) || push!(pieces[c], (c, piece, first(piece) == :higher ? 0 : weights[c] * lattice_index(L, T)))
        end
    end
    pieces = reduce(vcat, pieces)
    @req !any(p -> first(p[2]) == :higher, pieces) "a piece of dimension at least 2"
    labels = [first(p) == :leg ? (p[2], primitive(p[3])) : minmax(p[2], p[3]) for (_, p, _) in pieces]
    @req allunique(labels) "a piece lies in two cells"
    total = Dict{Vector{QQFieldElem},Vector{ZZRingElem}}()
    add!(v, d, m) = (total[v] = get(total, v, zeros(ZZ, 21)) + m .* d)
    for (_, p, m) in pieces
        first(p) == :leg ? add!(p[2], primitive(p[3]), m) :
            (add!(p[2], primitive(p[3] - p[2]), m); add!(p[3], primitive(p[2] - p[3]), m))
    end
    println("recomputed curve: $(count(p -> first(p[2]) == :segment, pieces)) segments, $(length(total)) vertices, ",
            "balanced $(all(iszero, values(total))), piece multiplicities $(histogram(last.(pieces)))")
    G = load_decorated_graphs(joinpath(dir, "decorated_graph_$k.mrdi")).subdivided
    stored = sort!([(normalized(G.points[e.tail]), only(e.cells)) for e in G.edges if is_leg(e)])
    legs = sort!([(p[2], c) for (c, p, _) in pieces if first(p) == :leg])
    println("same rays as decorated_graph_$k.mrdi (vertex and cell): $(stored == legs)")
    return [(primitive(p[3]), m) for (_, p, m) in pieces if first(p) == :leg]
end

function describe_rays(rays)
    println("rays: $(length(rays)), with multiplicity $(sum(last, rays)), multiplicities $(histogram(last.(rays)))")
    by_direction = Dict{Vector{ZZRingElem},Int}()
    foreach(((u, m),) -> by_direction[u] = get(by_direction, u, 0) + m, rays)
    println("directions: $(length(by_direction)), with their total multiplicity:")
    for (u, m) in sort!(collect(by_direction); by=((u, _),) -> findall(>(minimum(u)), u))
        println("  ", join(("e" * join(PLUCKER_PAIRS[i]) for i in findall(>(minimum(u)), u)), " + "),
                all(u[i] - minimum(u) in (0, 1) for i in 1:21) ? "" : " (not 0/1)", ": $m")
    end
    println("degree in every coordinate: ", join((sum(m * (u[i] - minimum(u)) for (u, m) in rays) for i in 1:21), ' '))
    println("balanced at infinity: ", iszero(normalized(sum(m .* u for (u, m) in rays))))
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (1, 2) || error("usage: curve_rays.jl DIR [K]")
    dir, k = abspath(ARGS[1]), length(ARGS) == 2 ? parse(Int, ARGS[2]) : 1
    fan = load(data_path("g2", "grassmannian", "02_coarsest_fan.mrdi"))
    println("equal coordinates on trop(G2/P2): ",
            join(("x" * join(PLUCKER_PAIRS[i]) * " = x" * join(PLUCKER_PAIRS[j]) for (i, j) in equal_coordinates(fan)), ", "))
    describe_rays(isfile(joinpath(dir, "curve.mrdi")) ? census_rays(dir) : bergman_rays(dir, k, fan))
end
