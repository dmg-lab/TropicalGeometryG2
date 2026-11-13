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
#     DIR/translation_K.mrdi by bergman_curve.jl, every piece of
#     C ∩ (B_S + t) of multiplicity m_C [ZZ^21 : L_C + L_S], with m_C over
#     QQbar, L_C = span(C) ∩ ZZ^21 and L_S = {x_S constant}. The curve is
#     checked balanced and its rays against DIR/decorated_graph_K.mrdi.
#
#   julia --project=. -t 8 studies/g2/curves/curve_rays.jl DIR [K]
###############################################################################

include(joinpath(@__DIR__, "bergman_curve.jl"))

"The pairs of coordinates equal on every ray and on the lineality of `fan`."
function equal_coordinates(fan)
    vectors = [[collect(r) for r in rays_modulo_lineality(fan).rays_modulo_lineality];
               [collect(r) for r in lineality_space(fan)]]
    return [(i, j) for i in 1:21 for j in i+1:21 if all(v -> v[i] == v[j], vectors)]
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
    weights = coarsest_cell_weights()
    @req length(weights) == length(cells) "one weight per cell"
    t = Vector{QQFieldElem}(QQ.(load(joinpath(dir, "translation_$k.mrdi"))))
    found = [only(f) for f in cell_candidates(cells, [t])]
    pieces = bergman_curve_pieces(cells, found, t; weights)
    @req !any(p -> first(p.piece) == :higher, pieces) "a piece of dimension at least 2"
    labels = [first(p.piece) == :leg ? (p.piece[2], primitive(p.piece[3])) : minmax(p.piece[2], p.piece[3]) for p in pieces]
    @req allunique(labels) "a piece lies in two cells"
    curve = balancing_and_degree([p.piece for p in pieces], [p.multiplicity for p in pieces])
    println("recomputed curve: $(count(p -> first(p.piece) == :segment, pieces)) segments, $(curve.vertices) vertices, ",
            "balanced $(curve.balanced), piece multiplicities $(histogram([p.multiplicity for p in pieces]))")
    G = load_decorated_graphs(joinpath(dir, "decorated_graph_$k.mrdi")).subdivided
    stored = sort!([(normalized(G.points[e.tail]), only(e.cells)) for e in G.edges if is_leg(e)])
    legs = sort!([(p.piece[2], p.cell) for p in pieces if first(p.piece) == :leg])
    println("same rays as decorated_graph_$k.mrdi (vertex and cell): $(stored == legs)")
    return [(primitive(p.piece[3]), p.multiplicity) for p in pieces if first(p.piece) == :leg]
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
    println("degree in every coordinate: ", join(curve_degree(rays), ' '))
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
