# Section 10 of the paper: run by paper_results.jl, which defines the helpers used here.

"""
The Mukai model of genus 10: G2/P2 ⊂ P^13 by the 4x4 Pfaffians of a skew 7x7 matrix. The matrix is a deliberate copy of
that of studies/g2/curves/mukai_models.jl, a script with no function to call.
"""
function mukai_model(K=QQ)
    R, (a, b, c, d, e, f, g, h, i, j, k, l, m, n) = polynomial_ring(K, [string(Char('a' + s)) for s in 0:13])
    M = matrix(R, [0 -f e g h i a; f 0 -d j k l b; -e d 0 m n -g-k c; -g -j -m 0 c -b d;
                   -h -k -n -c 0 a e; -i -l g+k b -a 0 f; -a -b -c -d -e -f 0])
    return R, ideal(R, pfaffians(M, 4))
end

"The core of a graph (leaves removed until none is left), then without its vertices of valence 2."
function core_graph(n, edges)
    kept = collect(edges)
    while true
        valence = zeros(Int, n)
        foreach(((u, v),) -> (valence[u] += 1; valence[v] += 1), kept)
        any(((u, v),) -> valence[u] == 1 || valence[v] == 1, kept) || break
        filter!(((u, v),) -> valence[u] > 1 && valence[v] > 1, kept)
    end
    subdivided = copy(kept)
    # Merge the edges at vertices of valence 2.
    while true
        valence = zeros(Int, n)
        foreach(((u, v),) -> (valence[u] += 1; valence[v] += 1), kept)
        w = findfirst(k -> valence[k] == 2 && count(e -> e[1] == k && e[2] == k, kept) == 0, 1:n)
        w === nothing && break
        a, b = [e for e in kept if w in e]
        other(e) = e[1] == w ? e[2] : e[1]
        filter!(e -> !(w in e), kept)
        push!(kept, (other(a), other(b)))
    end
    vertices = unique(vcat(first.(kept), last.(kept)))
    valence = Dict(v => count(e -> e[1] == v, kept) + count(e -> e[2] == v, kept) for v in vertices)
    return (; vertices, edges=kept, valence, subdivided)
end

"The primitive integer vector of a rational direction."
s10_primitive(d) = ZZ.(primitive_integer(d))

"A direction in R^21/R1 as its primitive representative with minimum 0."
s10_direction(d) = s10_primitive(d .- minimum(d))

"The lattice length of the segment from `a` to `b` in R^21/R1, whose lattice is Z^21/Z1."
function s10_lattice_length(a, b)
    d = b .- a
    d = d .- minimum(d)
    u = s10_primitive(d)
    k = findfirst(!iszero, u)
    return d[k] // u[k]
end

"A direction with minimum 0 named by its support, as e12+e67."
s10_name(u) = join(("e" * join(PAIRS[i]) for i in findall(>(minimum(u)), u)), "+")

"Whether `x` lies in the open segment from `a` to `b`."
function s10_inside(x, a, b)
    d = b .- a
    k = findfirst(!iszero, d)
    λ = (x[k] - a[k]) // d[k]
    return 0 < λ < 1 && x == a .+ λ .* d
end

"How often min(x - q) is attained: x lies in q + B_c (min convention of the files) when it is at least c + 1."
s10_ties(x, q) = (y = x .- q; count(==(minimum(y)), y))

"Whether `x` lies in the relative interior of a codimension-1 cone: two closed maximal cones, one tight facet in each."
function s10_on_wall(L, x)
    cones = locate_closed(L, x)
    return length(cones) >= 2 && all(c -> count(iszero, L[c].facets * x) == 1, cones)
end

"The connected pieces of a list of segments (pairs of vertices; a ray from u is (u, u)), as lists of segment indices."
function s10_pieces(segments)
    root = Dict{Int,Int}()
    function top(v)
        r = get!(root, v, v)
        r == v && return v
        t = top(r)
        root[v] = t
        return t
    end
    for (u, v) in segments
        a, b = top(u), top(v)
        root[a] = b
    end
    groups = Dict{Int,Vector{Int}}()
    for (k, (u, _)) in enumerate(segments)
        push!(get!(groups, top(u), Int[]), k)
    end
    return sort!(collect(values(groups)); by=length, rev=true)
end

"Whether a connected graph with vertices of valence at most 3 is planar: some rotation system has E - V + 2 faces."
function s10_planar_subcubic(edges)
    vs = unique(vcat(first.(edges), last.(edges)))
    index = Dict(v => k for (k, v) in enumerate(vs))
    E = length(edges)
    # Dart 2k - 1 runs along edge k from its first vertex, dart 2k back.
    tail(d) = index[isodd(d) ? edges[(d + 1) ÷ 2][1] : edges[d ÷ 2][2]]
    twin(d) = isodd(d) ? d + 1 : d - 1
    out = [Int[] for _ in vs]
    foreach(d -> push!(out[tail(d)], d), 1:2E)
    all(o -> length(o) <= 3, out) || error("a vertex of valence above 3")
    flips = findall(o -> length(o) == 3, out)
    flip = Dict(v => k for (k, v) in enumerate(flips))
    following = zeros(Int, 2E)
    seen = falses(2E)
    # The rotation at the last vertex of valence 3 stays fixed: the other choices give mirror images.
    for mask in 0:(1 << max(length(flips) - 1, 0)) - 1
        for (v, o) in enumerate(out)
            k = get(flip, v, 0)
            r = (k == 0 || iseven(mask >> (k - 1))) ? o : reverse(o)
            for i in eachindex(r)
                following[r[i]] = r[mod1(i + 1, length(r))]
            end
        end
        fill!(seen, false)
        faces = 0
        for d in 1:2E
            seen[d] && continue
            faces += 1
            e = d
            while !seen[e]
                seen[e] = true
                e = following[twin(e)]
            end
        end
        faces == E - length(vs) + 2 && return true
    end
    return false
end

"The edges of the skeleton as chains of segments of `core.subdivided`, from one vertex of valence 3 to another."
function s10_chains(core)
    at = Dict{Int,Vector{Int}}()
    for (k, (u, v)) in enumerate(core.subdivided)
        push!(get!(at, u, Int[]), k)
        push!(get!(at, v, Int[]), k)
    end
    used = falses(length(core.subdivided))
    chains = NamedTuple{(:tail, :head, :segments),Tuple{Int,Int,Vector{Int}}}[]
    for v in core.vertices, start in at[v]
        used[start] && continue
        k, w, segments = start, v, Int[]
        while true
            used[k] = true
            push!(segments, k)
            a, b = core.subdivided[k]
            w = a == w ? b : a
            w in core.vertices && break
            k = only(j for j in at[w] if !used[j])
        end
        push!(chains, (tail=v, head=w, segments=segments))
    end
    return chains
end

"Print a list as lines of at most `width` characters."
function s10_list(items; width=104)
    line = ""
    for item in string.(items)
        if !isempty(line) && length(line) + length(item) + 2 > width
            info("  " * line * ",")
            line = ""
        end
        line = isempty(line) ? item : line * ", " * item
    end
    isempty(line) || info("  " * line)
end

"""
A stored curve of genus 10 (Method 1 or 2): its graph, valences, skeleton, the maximal cone of trop(G) of each bounded
edge, ray and skeleton segment (located at a midpoint, or at u + r on a ray), and for Method 2 its rays with their
primitive directions and multiplicities.
"""
function s10_curve(method, L)
    stored_curve = load(data_path("mukai_curves", "genus_10_method_$method.mrdi"))
    graph = curve_graph(stored_curve.complex)
    n = length(graph.points)
    valence = zeros(Int, n)
    foreach(((u, v),) -> (valence[u] += 1; valence[v] += 1), graph.edges)
    foreach(((u, _),) -> valence[u] += 1, graph.rays)
    core = core_graph(n, graph.edges)
    edge_cones = [only(locate(L, (graph.points[u] + graph.points[v]) / 2)) for (u, v) in graph.edges]
    ray_cones = [only(locate(L, graph.points[u] + r)) for (u, r) in graph.rays]
    cone_of = Dict(zip(graph.edges, edge_cones))
    skeleton_cones = [cone_of[e] for e in core.subdivided]
    legs = hasproperty(stored_curve, :multiplicities) ?
           [(s10_primitive(QQ.(collect(only(Oscar.rays(P))))), Int(m))
            for (P, m) in zip(maximal_polyhedra(stored_curve.complex), stored_curve.multiplicities) if !is_bounded(P)] :
           nothing
    return (; method, stored_curve, graph, n, valence, core, edge_cones, ray_cones, skeleton_cones, legs)
end

"The vertices of valence 2 on the skeleton whose two segments lie in different maximal cones: the broken lines."
function s10_broken_lines(C)
    at = Dict{Int,Vector{Int}}()
    for (k, (u, v)) in enumerate(C.core.subdivided)
        push!(get!(at, u, Int[]), k)
        push!(get!(at, v, Int[]), k)
    end
    two = sort([v for (v, ks) in at if length(ks) == 2])
    return two, [v for v in two if length(unique(C.skeleton_cones[at[v]])) == 2]
end

"Print the regions of a right panel: for each maximal cone met by the skeleton, its segments and connected pieces."
function s10_regions(C, types)
    cones = sort(unique(C.skeleton_cones); by=c -> (parse(Int, types[("cell", c)]), c))
    info("Method $(C.method): the regions, one for each maximal cone met by the skeleton (cone, type: skeleton segments, " *
         "and the segments of each connected piece):")
    for c in cones
        segments = findall(==(c), C.skeleton_cones)
        pieces = length.(s10_pieces(C.core.subdivided[segments]))
        info("  cone $c, type $(types[("cell", c)]): $(length(segments)) segment$(length(segments) == 1 ? "" : "s"), " *
             "piece$(length(pieces) == 1 ? "" : "s") of $(join(pieces, ", "))")
    end
end

"Print the points of (d) by the maximal cone of trop(G) containing them."
function s10_points_by_cone(method, cones, types)
    info("Method $method: the 18 points of (d) by maximal cone (cone, type: points):")
    s10_list(["$c, type $(types[("cell", c)]): $(count(==(c), cones))"
              for c in sort(unique(cones); by=c -> (-count(==(c), cones), c))])
end

function section_10()
    subheading("The linear sections of G (here in the Mukai model of genus 10, G2/P2 ⊂ P^13)")
    R, I = mukai_model()
    monomials = unique(vcat([collect(AbstractAlgebra.monomials(f)) for f in gens(I)]...))
    quadrics = rank_of([[coeff(f, m) for m in monomials] for f in gens(I)])
    check("the ideal of the model: independent quadrics among the 35 Pfaffians; dimension (projective) and degree",
          [28, 5, 18], [quadrics, dim(I) - 1, Int(degree(I))])
    # The same sections in a graded ring, for their Hilbert polynomials.
    Rg, xg = graded_polynomial_ring(QQ, [string(Char('a' + s)) for s in 0:13])
    φ = hom(R, Rg, xg)
    Ig = ideal(Rg, [φ(f) for f in gens(I)])
    hilbert = Dict(2 => [1, 7//2, 9//2, 3], 3 => [2, 0, 9, 0], 4 => [-9, 18, 0, 0], 5 => [18, 0, 0, 0])
    Random.seed!(10)
    names = ["(a) Fano threefold in P^11", "(b) K3 surface in P^10", "(c) canonical curve in P^9", "(d) 18 points in P^8"]
    for (c, name) in zip(2:5, names)
        forms = [sum(rand(-9:9) * v for v in gens(R)) for _ in 1:c]
        section = I + ideal(R, forms)
        check("$name: dimension and degree", [5 - c, 18], [dim(section) - 1, Int(degree(section))])
        A, _ = quo(Rg, Ig + ideal(Rg, [φ(f) for f in forms]))
        P = hilbert_polynomial(A)
        coefficients = [coeff(P, k) for k in 0:3]
        if c < 5
            genus = Int(factorial(5 - c) * coefficients[6-c]) ÷ 2 + 1
            check("$name: Hilbert polynomial (coefficients of 1, t, t^2, t^3) and genus g, by H^$(5 - c) = 2g - 2",
                  (QQ.(hilbert[c]), 10), (coefficients, genus))
        else
            check("$name: Hilbert polynomial (coefficients of 1, t, t^2, t^3)", QQ.(hilbert[c]), coefficients)
        end
    end
    cited("analogous constructions for g = 7, 8, 9 yield birational models of M_g (Kap1, Kap2)")
    if QUICK
        skipped("the Betti table eq:bettitable and the Gorenstein property of (d)")
    else
        # The 18 points: a random linear section of the model in P^8, over GF(32003) (over QQ the coefficients grow
        # beyond several GB); the Betti numbers of a generic section are those over QQ for all but finitely many p.
        K = GF(32003)
        Rp, Ip = mukai_model(K)
        S, y = graded_polynomial_ring(K, ["y$s" for s in 0:8])
        Random.seed!(10)
        restriction = hom(Rp, S, [sum(K(rand(1:32002)) * v for v in y) for _ in 1:14])
        points = ideal(S, [restriction(f) for f in gens(Ip)])
        check("the section over GF(32003): dimension and degree", [0, 18], [dim(points) - 1, Int(degree(points))])
        B = minimal_betti_table(points)
        info("the minimal Betti table of the ideal of the 18 points (d):")
        foreach(line -> info("  " * line), split(sprint(show, B), '\n'))
        # Keys (homological step, degree) of the resolution of the ideal; its last module S(-11) is the socle.
        D = Oscar.as_dictionary(B)
        last_step = maximum(first, keys(D))
        check("(d) is arithmetically Gorenstein: the resolution of its ideal ends with one S(-11) (step, degree, rank)",
              [(7, 11, 1)], [(i, Int(d[1]), b) for ((i, d), b) in D if i == last_step])
        table = split(sprint(show, B), '\n')
        betti_row(r) = [x == "-" ? 0 : parse(Int, x)
                        for x in split(strip(split(only(l for l in table if startswith(strip(l), "$r:")), ':')[2]))][1:7]
        check("the Betti table eq:bettitable, its two rows", [[28, 105, 162, 84, 0, 0, 0], [0, 0, 0, 84, 162, 105, 28]],
              [betti_row(2), betti_row(3)])
    end

    subheading("The tropical moduli space trop(M_10) and the tropical linear spaces L_4")
    cited("trop(M_10) is a stacky fan of dimension 27 with 2317172 maximal cones, one for each trivalent graph of genus " *
          "10, with 18 vertices and 27 edges; 41301 of these graphs are simple")
    # Speyer's bound on the bounded i-faces of a tropical linear space of a (d, n) matroid. A face of a generic L_4 is
    # a bounded face of a stratum x_T = ∞ (a generic space of U(17 - t, 21 - t), t = |T|) plus the cone of the e_i, i in T.
    bounded_faces(i, d, n) = (i < 0 || d - i - 1 < 0 || n - 2i - 2 < 0) ? 0 :
                             binomial(n - i - 2, i) * binomial(n - 2i - 2, d - i - 1)
    info("Speyer's formula: bounded faces of a generic L_4 of dimensions 0, 1, 2, 3: " *
         "$(join([bounded_faces(i, 17, 21) for i in 0:3], ", "))")
    f_L4 = [sum(binomial(21, t) * bounded_faces(k - t, 17 - t, 21 - t) for t in 0:k) for k in 0:16]
    check("Speyer's formula: L_4 can have up to this many maximal cells", 735471, last(f_L4))
    info("the expected f-vector of L_4 (no longer printed in the paper): $(show_value(f_L4))")

    subheading("Proposition prop:mukai: the two curves of genus 10 (data/mukai_curves/genus_10_method_{1,2}.mrdi)")
    L = locator(fan("grassmannian"))
    types = orbit_types("grassmannian_orbits.csv", "qq")
    curves = [s10_curve(method, L) for method in 1:2]
    M1, M2 = curves
    for C in curves
        (; method, graph, n, valence, core) = C
        info("Method $method: $(length(graph.edges) + length(graph.rays)) maximal cells, $(length(graph.rays)) rays, " *
             "$(length(graph.edges)) bounded edges on $n vertices, valences $(sort([(k, count(==(k), valence)) for k in unique(valence)]))")
        check("Method $method: genus of the bounded part", 10, length(graph.edges) - n + 1)
        check("Method $method: the core without vertices of valence 2 is trivalent with (vertices, edges)", (18, 27, true),
              (length(core.vertices), length(core.edges), all(==(3), values(core.valence))))
    end
    q_paper = [204958, 1849487, 1813412, -1389486, 611356, 155086, 128988, 1698926, -1863394, 602010, 2239681,
               1904628, -124320, 1025736, 87111, 1360439, -705715, 619752, -2004970, -1344973, -232336]
    q = Int.(load(project_path("curves_g10", "overnight_20261028", "wpart_2722", "translation_1.mrdi")))
    # main.tex prints q and q_1, ..., q_4 in its max convention: minus the shifts stored in the min convention of the files.
    check("Method 1: the shift vector q, minus the stored one (translation_1.mrdi)", q_paper, -q)
    experiment = project_path("curves_g10", "overnight_20261028", "experiment_100128")
    planes = load(joinpath(experiment, "hyperplanes.mrdi"))
    support = [findall(!iszero, planes.normals[h]) for h in 1:5]
    check("Method 2: supports of the normal vectors of the five hyperplanes (q_5 of full support)", [3, 8, 9, 15, 21],
          length.(support))
    check("Method 2: the shift vectors q_1, ..., q_4 on the supports of their normals, minus the stored ones (hyperplanes.mrdi)",
          [["p13" => 217936, "p37" => 84861, "p46" => -34008],
           ["p14" => -1013667, "p15" => -132815, "p25" => -744071, "p26" => 794071, "p37" => -721948, "p46" => 60165,
            "p47" => -360635, "p67" => -705771],
           ["p15" => 806183, "p26" => 732664, "p27" => 107361, "p34" => 721205, "p35" => 846796, "p36" => 183135,
            "p46" => -542064, "p47" => -278271, "p67" => -78144],
           ["p12" => -722953, "p14" => -1004831, "p15" => 757760, "p16" => 217944, "p17" => -733597, "p23" => -589026,
            "p24" => -213921, "p26" => -650364, "p34" => 890415, "p35" => 891447, "p36" => 636797, "p37" => -914965,
            "p46" => -315068, "p56" => 862815, "p57" => -692897]],
          [[pair_name(k) => -Int(planes.shifts[h][k]) for k in support[h]] for h in 1:4])
    stage = only(rows(joinpath(experiment, "summary.csv")))
    stored("Method 2: pieces of trop(G) ∩ L_c for c = 1, ..., 5, the last the 18 points of (d) (sample_genus_10.jl, summary.csv)",
           replace(stage["pieces"], ';' => ", "))

    subheading("The tropical curves (left panels): degree, rays, Berkovich skeleton, the 18 points of (d)")
    grays, glin = fan_rays(fan("grassmannian"))
    equal = [(pair_name(i), pair_name(j)) for i in 1:21 for j in i+1:21 if all(v -> v[i] == v[j], vcat(grays, glin))]
    check("trop(G) lies in x_i = x_j for 6 pairs (on its rays and lineality), which span a space of dimension " *
          "21 - 6 = 15: the curves live in dimension 15",
          ([("p12", "p67"), ("p13", "p57"), ("p17", "p56"), ("p23", "p47"), ("p27", "p46"), ("p37", "p45")], 15),
          (equal, rank_of(vcat(grays, glin))))
    claim("Method 1: at most 15 × 18 = 270 rays", length(M1.graph.rays) <= 270, "$(length(M1.graph.rays)) rays")
    record = TOML.parsefile(project_path("curves_g10", "overnight_20261028", "wpart_2722", "wpart_2722_record.toml"))
    row = record["rays"]
    stored("Method 1: rays, rays with multiplicity (the bound 270 is reached), multiplicities (curve_rays.jl, 94 s on " *
           "8 threads; wpart_2722_record.toml)", "$(row["rays"]), $(row["rays_with_multiplicity"]), $(row["multiplicities"])")
    directions = Dict{Vector{ZZRingElem},Int}()
    foreach(((_, r),) -> directions[s10_direction(r)] = get(directions, s10_direction(r), 0) + 1, M1.graph.rays)
    order(u) = (count(>(0), u), findall(>(0), u))
    info("Method 1: the rays by direction (direction: rays; genus_10_method_1.mrdi stores no multiplicities):")
    s10_list(["$(s10_name(u)): $k" for (u, k) in sort!(collect(directions); by=order ∘ first)])
    legs = M2.legs
    check("Method 2: multiplicities of the maximal cells", [1, 2, 4, 5], sort(unique(Int.(M2.stored_curve.multiplicities))))
    check("Method 2: degree 18 in every coordinate (the sum over the rays of m (u_i - min u), u primitive)", fill(18, 21),
          [Int(sum(m * (u[i] - minimum(u)) for (u, m) in legs)) for i in 1:21])
    claim("Method 2: at most 15 × 18 = 270 rays, counted with multiplicity", sum(last, legs) <= 270,
          "$(length(legs)) rays, $(sum(last, legs)) with multiplicity")
    by_direction = Dict{Vector{ZZRingElem},Tuple{Int,Int}}()
    for (u, m) in legs
        k, total = get(by_direction, u .- minimum(u), (0, 0))
        by_direction[u .- minimum(u)] = (k + 1, total + m)
    end
    info("Method 2: the rays by direction (direction: rays/total multiplicity):")
    s10_list(["$(s10_name(u)): $k/$m" for (u, (k, m)) in sort!(collect(by_direction); by=order ∘ first)])
    for C in curves
        chains = s10_chains(C.core)
        label = Dict(v => k for (k, v) in enumerate(sort(C.core.vertices)))
        lengths = [sum(s10_lattice_length(C.graph.points[C.core.subdivided[k][1]], C.graph.points[C.core.subdivided[k][2]])
                       for k in chain.segments) for chain in chains]
        crossed(chain) = (cs = C.skeleton_cones[chain.segments];
                          join((types[("cell", c)] for (i, c) in enumerate(cs) if i == 1 || c != cs[i-1]), " "))
        info("Method $(C.method): the metric graph, a point of trop(M_10): its $(length(chains)) edges (vertices numbered " *
             "1 to 18) with their lattice lengths and the types of the cones they cross:")
        s10_list(["$(label[c.tail])-$(label[c.head]): $ℓ ($(crossed(c)))" for (c, ℓ) in zip(chains, lengths)])
        if C.method == 2
            csv = [QQ(parse(Int, r["length"])) for r in rows(joinpath(experiment, "core_edges.csv"))]
            claim("Method 2: the lattice lengths of the 27 edges equal those of experiment_100128/core_edges.csv " *
                  "(curve_details.jl)", sort(lengths) == sort(csv))
        end
    end
    m1_ties = [s10_ties(x, q) for x in M1.graph.points]
    check("Method 1: vertices of the curve with min(x - q) attained 5 times (in L_4) and 6 times (in L_5)", (183, 18),
          (count(==(5), m1_ties), count(==(6), m1_ties)))
    claim("Method 1: the 18 points of (d) = trop(G) ∩ L_5 are the 18 vertices of the skeleton, those of valence 6",
          Set(findall(>=(6), m1_ties)) == Set(M1.core.vertices) == Set(findall(==(6), M1.valence)))
    divisor = load(joinpath(experiment, "points.mrdi")).points
    # In the min convention of the files, x lies on the hyperplane h when min(x_i - s_i) over the support of its normal is
    # attained twice, s the stored shift; the paper's q_h is -s and its points -x (max convention).
    on_plane(x, h) = (gap = x[support[h]] .- planes.shifts[h][support[h]]; count(==(minimum(gap)), gap) >= 2)
    curve_points = vcat(M2.graph.points, [(M2.graph.points[u] + M2.graph.points[v]) / 2 for (u, v) in M2.graph.edges],
                        [M2.graph.points[u] + r for (u, r) in M2.graph.rays])
    claim("Method 2: the curve (vertices, midpoints of the edges, points u + r of the rays) lies on the hyperplanes " *
          "q_h + B_{S_h} for h = 1, ..., 4 and the 18 points of (d) on all five",
          all(x -> all(h -> on_plane(x, h), 1:4), curve_points) && all(x -> all(h -> on_plane(x, h), 1:5), divisor))
    inside = [count(((u, v),) -> s10_inside(x, M2.graph.points[u], M2.graph.points[v]), M2.graph.edges) for x in divisor]
    check("Method 2: the points of (d) (points.mrdi): their number, the bounded edges containing each in their interior, " *
          "and whether one is a vertex", (18, fill(1, 18), false), (length(divisor), inside, any(in(M2.graph.points), divisor)))
    places = rows(joinpath(experiment, "point_places.csv"))
    info("Method 2: $(count(r -> startswith(r["place"], "on core edge"), places)) of them on the skeleton, " *
         "$(count(r -> startswith(r["place"], "on a pending tree"), places)) on pending trees (curve_details.jl, point_places.csv)")

    subheading("The traces on trop(G) (right panels): the maximal cones met by L_4")
    check("trop(G) has 248 maximal cones", 248, length(L))
    claim("Method 1: the trivalent skeleton is planar (Figure fig:planar_curve)", s10_planar_subcubic(M1.core.edges))
    by_number(ts) = sort(ts; by=t -> parse(Int, t))
    met = sort(unique(M1.skeleton_cones))
    check("Method 1: the Berkovich skeleton meets 5 maximal cones of trop(G), of types", ["11", "11", "12", "12", "7"],
          sort([types[("cell", c)] for c in met]))
    # main.tex: "Only five of the $248$ maximal cones of $\trop(\linespace)$ meet the skeleton".
    info("Method 1: the whole curve trop(G) ∩ L_4, trees and rays included, meets " *
         "$(length(unique(vcat(M1.edge_cones, M1.ray_cones)))) maximal cones of trop(G)")
    s10_regions(M1, types)
    seven = only(c for c in met if types[("cell", c)] == "7")
    # main.tex: "the cone of type 7 intersects $\mathfrak{L}_4$ in two pieces on the skeleton: one edge and one tree with
    # five edges".
    check("Method 1: the cone of type 7 meets the skeleton in two pieces (segments): one edge and a tree with five edges",
          [5, 1], length.(s10_pieces(M1.core.subdivided[findall(==(seven), M1.skeleton_cones)])))
    in_seven = vcat(M1.graph.edges[findall(==(seven), M1.edge_cones)],
                    [(u, u) for (u, _) in M1.graph.rays[findall(==(seven), M1.ray_cones)]])
    info("Method 1: the cone of type 7 meets the whole curve trop(G) ∩ L_4, trees and rays included, in pieces of " *
         "$(join(length.(s10_pieces(in_seven)), ", ")) segments and rays")
    divisor_cones = [only(locate(L, M1.graph.points[v])) for v in M1.core.vertices]
    check("Method 1: points of (d) in the cone of type 7, the two blue vertices", 2, count(==(seven), divisor_cones))
    s10_points_by_cone(1, divisor_cones, types)
    two, crossings = s10_broken_lines(M1)
    info("Method 1: $(length(two)) vertices of valence 2 on the skeleton, $(length(crossings)) of them between two regions")
    # main.tex: "Broken lines between two colored regions are $3$-valent vertices of the curve that become $2$-valent in the
    # Berkovich skeleton".
    check("Method 1: valences in the curve of the broken lines, the vertices of valence 2 on the skeleton between two " *
          "regions (the paper: 3-valent vertices of the curve)", [3], sort(unique(M1.valence[crossings])))
    claim("Method 1: the broken lines lie in the relative interior of codimension-1 cones of trop(G)",
          all(v -> s10_on_wall(L, M1.graph.points[v]), crossings))
    claim("Method 2: the trivalent skeleton is not planar (Figure fig:nonplanar_curve)", !s10_planar_subcubic(M2.core.edges))
    met = unique(M2.skeleton_cones)
    check("Method 2: the skeleton meets more than 5 maximal cones of trop(G), among them types 17, 12 and 16", (true, true),
          (length(met) > 5, issubset(["17", "12", "16"], [types[("cell", c)] for c in met])))
    m2_types = by_number([types[("cell", c)] for c in met])
    info("Method 2: the skeleton meets $(length(met)) maximal cones of trop(G), of types $(join(m2_types, ", ")); the " *
         "whole curve, trees and rays included, meets $(length(unique(vcat(M2.edge_cones, M2.ray_cones))))")
    s10_regions(M2, types)
    s10_points_by_cone(2, [only(locate(L, x)) for x in divisor], types)
    two, crossings = s10_broken_lines(M2)
    claim("Method 2: the broken lines lie in the relative interior of codimension-1 cones of trop(G)",
          all(v -> s10_on_wall(L, M2.graph.points[v]), crossings),
          "$(length(crossings)) of the $(length(two)) vertices of valence 2 on the skeleton, of valences " *
          "$(join(sort(unique(M2.valence[crossings])), ", ")) in the curve")

    subheading("The K3 surface (b) of Method 1: trop(G) ∩ L_3")
    statement = "the K3 surface of Method 1: its bounded part collapses onto a 2-sphere, cut by the core into its " *
                "E - V + 2 faces; the core is the edge graph of a simple 3-polytope"
    # The polygons of k3_surface_polygons.jl, written by the stage zenodo, or by the script into results/.
    candidates = [data_path("mukai_curves", "genus_10_method_1_k3.mrdi"),
                  project_path("results", "curves", "genus_10", "k3", "wpart_2722_x18.mrdi")]
    found = findfirst(isfile, candidates)
    if !isnothing(found)
        x18 = candidates[found]
        m = run_in_module(study_path("g2", "curves", "k3_fins.jl"); name=:K3Fins)
        G = call(m, :load_decorated_graphs,
                 project_path("curves_g10", "overnight_20261028", "wpart_2722", "decorated_graph_1.mrdi")).subdivided
        polygons = call(m, :load_k3_polygons, x18)
        # k3_fins prints its own report; only the claim goes into this one.
        result = redirect_stdout(devnull) do
            call(m, :k3_fins, polygons, G)
        end
        claim(statement, result.sphere && result.collapses && result.core_on_sphere &&
              result.pieces == result.faces && result.polytope_graph, "$(result.pieces) pieces, $(result.faces) faces")
    else
        # Without the polygons, the result of k3_fins.jl as recorded in wpart_2722_record.toml.
        k3 = record["k3_surface"]
        k3["core_simple_and_3_connected"] || error("wpart_2722_record.toml: the core is not recorded simple and 3-connected")
        stored(statement * " (k3_fins.jl, wpart_2722_record.toml; data/mukai_curves/genus_10_method_1_k3.mrdi of " *
               "the stage zenodo is absent)",
               "a sphere of $(k3["sphere_polygons"]) polygons ($(k3["sphere_vertices"]) vertices, $(k3["sphere_edges"]) " *
               "edges), $(k3["fins"]) fins removed by $(k3["polygon_collapses"]) + $(k3["edge_collapses"]) " *
               "collapses, $(k3["core_discs"]) discs; the core simple and 3-connected")
    end
end
