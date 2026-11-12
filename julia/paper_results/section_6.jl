# Section 6 of the paper: run by julia/paper_results.jl, which defines the helpers used here.

"""
The tree of a row of 04_trees.csv, as `s6_drawn_trees` gives a drawn tree: its internal splits (sides without 7, sorted),
the pendant lengths of 1-7 and the splits => length.
"""
function s6_tree(r)
    lengths = Dict{Vector{Int},Int}()
    for s in split(r["internal_splits"], ';'; keepempty=false)
        side, len = split(s, ':')
        A = Set(c - '0' for c in first(split(side, '|')))
        lengths[sort(collect(7 in A ? setdiff(Set(1:7), A) : A))] = parse(Int, len)
    end
    return (; splits=sort(collect(keys(lengths))), pendant=parse.(Int, split(r["pendant_lengths"], ';')), lengths)
end

"A tree of `s6_tree` or `s6_drawn_trees` as text: the pendant lengths of 1-7, then the splits (side without 7) with their lengths."
s6_tree_text(t) = (t.pendant === nothing ? "" : "pendant " * join(t.pendant, ";") * ", ") * "splits " *
                  join([join(S) * (t.lengths === nothing ? "" : ":$(get(t.lengths, S, 0))") for S in t.splits], " ")

"""
The G2-trees drawn in a figure, in their order, each named by the first capture of the last match of `label` before it
(a tree starts at its leaf v1; v1-v7 are the taxa): its internal splits (sides without 7, sorted) and, when weights are
drawn, the pendant lengths of 1-7 and the splits => length, without length 0. Each weight goes to its nearest edge, since
some are drawn beside it.
"""
function s6_drawn_trees(tex, label)
    marks = [(m.offset, String(m[1])) for m in eachmatch(label, tex)]
    starts = [m.offset for m in eachmatch(r"\\coordinate \(v1\) at", tex)]
    trees = Pair{String,Any}[]
    for (n, start) in enumerate(starts)
        mark = findlast(m -> first(m) < start, marks)
        mark === nothing && continue
        name = marks[mark][2]
        block = SubString(tex, start, n < length(starts) ? prevind(tex, starts[n+1]) : lastindex(tex))
        point = Dict(parse(Int, m[1]) => parse.(Float64, (m[2], m[3]))
                     for m in eachmatch(r"\\coordinate \(v(\d+)\) at \(([-0-9.]+),([-0-9.]+)\)", block))
        edges = [(String(m[1]), parse(Int, m[2]), parse(Int, m[3]))
                 for m in eachmatch(r"\\draw\[g2tree (internal|pendant)[^\]]*\] \(v(\d+)\)--\(v(\d+)\)", block)]
        function distance(p, (_, a, b))
            (x, y), (u, v) = point[a], point[b]
            t = clamp(((p[1] - x) * (u - x) + (p[2] - y) * (v - y)) / ((u - x)^2 + (v - y)^2), 0.0, 1.0)
            return hypot(p[1] - x - t * (u - x), p[2] - y - t * (v - y))
        end
        weight = Dict{Int,Int}()
        for m in eachmatch(r"\\node\[g2tree weight[^\]]*\] at \(([-0-9.]+),([-0-9.]+)\) \{(\d+)\}", block)
            p = parse.(Float64, (m[1], m[2]))
            e = argmin(k -> distance(p, edges[k]), eachindex(edges))
            haskey(weight, e) && error("$name: two weights nearest to one edge")
            weight[e] = parse(Int, m[3])
        end
        # The taxa reached from a without crossing the edge a-b.
        function side(a, b)
            seen, todo = Set([a]), [a]
            while !isempty(todo)
                v = pop!(todo)
                for (_, x, y) in edges, (s, w) in ((x, y), (y, x))
                    s == v && !(w in seen) && Set([s, w]) != Set([a, b]) && (push!(seen, w); push!(todo, w))
                end
            end
            return sort([v for v in seen if v <= 7])
        end
        split_of(a, b) = (S = side(a, b); 7 in S ? setdiff(1:7, S) : S)
        internal = [(k, a, b) for (k, (kind, a, b)) in enumerate(edges) if kind == "internal"]
        splits = sort([split_of(a, b) for (_, a, b) in internal])
        if isempty(weight)
            push!(trees, name => (; splits, pendant=nothing, lengths=nothing))
            continue
        end
        length(weight) == length(edges) || error("$name: $(length(weight)) weights on $(length(edges)) edges")
        pendant = [weight[only(k for (k, (kind, a, b)) in enumerate(edges) if kind == "pendant" && i in (a, b))] for i in 1:7]
        lengths = Dict(split_of(a, b) => weight[k] for (k, a, b) in internal if weight[k] != 0)
        push!(trees, name => (; splits, pendant, lengths))
    end
    return trees
end

function section_6()
    I = load(g2("grassmannian", "00_ideal.mrdi"))
    P = base_ring(I)
    names = string.(symbols(P))
    p(i, j) = P[findfirst(==("p$i$j"), names)]
    census = Dict(r["orbit"] => r for r in rows(g2("grassmannian", "06_canonical_interior_degenerations.csv")))
    census2 = Dict(r["orbit"] => r for r in rows(g2("grassmannian_gf2", "06_canonical_interior_degenerations.csv")))
    orbit_rows = [r for r in rows(g2("grassmannian_orbits.csv")) if r["kind"] == "cell"]
    qq_orbit_of(type) = only(r["qq_orbit"] for r in orbit_rows if r["orbit_type"] == type)
    # The rows of Table tab:g2-grassmannian-newton-okounkov: type => (orbit size, f-vector, N, Char).
    paper_rows = Dict{String,NTuple{4,String}}()
    for line in readlines(repository("tables", "toric_degeneration_table.tex"))
        m = match(r"oldstylenums\{(\d+)\}\}\((\d+)\)\} & \$\(([\d,]+)\)\$ & .* & (\d+) & (all|\$2\$) \\tabularnewline", line)
        m === nothing || (paper_rows[String(m[1])] = (String(m[2]), String(m[3]), String(m[4]), m[5] == "all" ? "all" : "2"))
    end

    subheading("Prime cones and their Newton-Okounkov polytopes P_C")
    cited("Kaveh and Manon: the toric degenerations correspond to the maximal cones of the tropical variety with prime initial ideal")
    prime_report = studies("g2", "reports", "prime_cone.txt")
    if isfile(prime_report)
        stored("type 1: in_d(G) at four points of relint(C) (supplemental/g2_prime_cone.jl)",
               strip(only(l for l in readlines(prime_report) if occursin("all four initial ideals are equal", l))))
    else
        skipped("type 1: in_d(G) at four points of relint(C) (supplemental/g2_prime_cone.jl)", "report not in this repository")
    end
    if QUICK
        skipped("in_C(G) at a second point of relint(C) for every prime cone", "--quick")
    else
        NR = load(g2("grassmannian", "03_normalized_rays.mrdi"))   # the census column rays indexes its rows
        prime_rows = [r for r in values(census) if r["prime"] == "true"]
        same = map(prime_rows) do r
            J = load(joinpath(g2("grassmannian"), r["initial_ideal"]))
            first_ray = parse(Int, first(split(r["rays"], ';')))
            w = [QQ(parse(Int, x)) for x in split(r["canonical_weight"], ';')] + [QQ(NR[first_ray, j]) for j in 1:21]
            toP = hom(base_ring(J), P, gens(P))
            initial(I, tropical_semiring_map(QQ), w) == ideal(P, toP.(gens(J)))
        end
        claim("in_C(G) = in_d(G) for any d in relint(C): at the canonical weight plus a ray of C it is the stored in_C(G)",
              all(same), "$(count(same)) of $(length(same)) prime cones over QQ")
    end
    for (target, field) in (("grassmannian", "qq"), ("grassmannian_gf2", "gf2"))
        list = toric_orbits(target, field)
        agree = all(list) do o
            B = body(target, o.row)
            dim(B) == 5 && normalized_volume(B) == 18 && string.(f_vector(B)) == split(o.row["body_f_vector"], ';') &&
                length(lattice_points(B)) == parse(Int, o.row["body_lattice_points"])
        end
        claim("$(field == "qq" ? "over QQ" : "over GF(2)"): every P_C has dimension 5 and normalized volume 18 " *
              "(recomputed from the stored polytopes, with the census f-vector and N)", agree, "$(length(list)) orbits")
    end

    subheading("Example: type 1")
    # The trees of the figures, compared with each other and with the census trees of their cells (04_trees.csv).
    type_of_cell = Dict(c => r["orbit_type"] for r in orbit_rows for c in parse.(Int, split(r["qq_members"], ';'; keepempty=false)))
    tree_rows = Dict(parse(Int, r["index"]) => r for r in rows(g2("grassmannian", "04_trees.csv")) if r["kind"] == "cell")
    panel = r"\\oldstylenums\{(\d+)\}\}\{\\large"
    first_prime = first(s6_drawn_trees(read(repository("figures", "toric_degeneration_prime_edge_trees.tex"), String), panel))
    rightmost = only(s6_drawn_trees(read(repository("figures", "triangle_tree.tex"), String), r"^% Coarsest G2 Grassmannian cell (\d+)"))
    cell1 = parse(Int, first(rightmost))
    check("the first tree of Figure atlas_toric_cells_prime (its type) equals the rightmost tree in Figure ex_G2_trees " *
          "(figures/triangle_tree.tex)", ("1", s6_tree_text(last(rightmost))), (first(first_prime), s6_tree_text(last(first_prime))))
    check("it is the tree of cell $cell1 (header of figures/triangle_tree.tex), a cell of type 1 (04_trees.csv)",
          (s6_tree_text(last(first_prime)), "1"), (s6_tree_text(s6_tree(tree_rows[cell1])), type_of_cell[cell1]))
    qq = toric_orbits("grassmannian", "qq")
    type1 = only(o for o in qq if o.type == "1")
    # d of Example ex:fivetrees (max convention); the census computes initial ideals in the min convention.
    d1 = Dict("17" => 7, "56" => 7, "27" => 7, "46" => 7, "37" => 6, "45" => 6, "47" => 7, "23" => 7, "57" => 7,
              "13" => 7, "67" => 6, "12" => 6, "36" => 6, "14" => 8, "25" => 8, "15" => 8, "16" => 1, "24" => 8,
              "26" => 5, "34" => 1, "35" => 5)
    check("the metric d of Example ex:fivetrees is minus the canonical weight of the census cell (QQ orbit $(type1.row["orbit"]))",
          [d1["$i$j"] for (i, j) in PAIRS], [-parse(Int, x) for x in split(type1.row["canonical_weight"], ';')])
    ideal1 = load(joinpath(g2("grassmannian"), type1.row["initial_ideal"]))
    degrees = [total_degree(f) for f in gens(ideal1)]
    check("in_d(G): linear binomials and quadrics", [7, 28], [count(==(1), degrees), count(==(2), degrees)])
    R1 = base_ring(ideal1)
    (p12, p13, p14, p15, p16, p17, p23, p24, p25, p26, p27, p34, p35, p36, p37, p45, p46, p47, p56, p57, p67) = gens(R1)
    paper1 = [p37 - p45, p12 + p67, p27 + p46, p23 + p47, p17 - p56, p13 - p57, p14 + p25,
        p26*p57 + p36*p67, p36*p45 - p35*p46, p45*p56 - p36*p57, p36*p56 + p57*p67,
        p45*p46 - p36*p47, p36*p46 + p47*p67, p25*p35 + p45*p57, p24*p35 + p45*p47,
        p25*p26 - p46*p67, p15*p26 + p56*p67, p47*p56 - p46*p57, p46*p56 - p25*p67,
        p25*p45 + p47*p57, p25*p36 + p46*p57, p24*p36 + p46*p47, p15*p36 - p56*p57,
        p25*p47 - p24*p57, p15*p47 + p25*p57, p25*p46 - p24*p56, p15*p46 + p25*p56,
        p25^2 + p15*p24, p36^2 + p45*p67, p45^2 - p35*p47, p46^2 - p24*p67,
        p47^2 + p24*p45, p56^2 + p15*p67, p57^2 - p15*p45, p67^2 - p26*p56]
    flipped = [f for f in paper1 if !(f in gens(ideal1)) && -f in gens(ideal1)]
    check("the 35 generators printed in the paper are the stored generators of in_d(G), up to sign", (35, true),
          (length(gens(ideal1)), all(f -> f in gens(ideal1) || -f in gens(ideal1), paper1)))
    info("stored with the opposite sign: " * join(string.(flipped), ", "))
    claim("they generate in_d(G)", ideal(R1, paper1) == ideal1)
    P1 = body("grassmannian", type1.row)
    V1 = load(joinpath(g2("grassmannian"), type1.row["degree_one_values"]))   # row k is the point of PAIRS[k]
    point1(k) = [QQ(V1[k, j]) for j in 1:ncols(V1)]
    vertex_set = Set(Vector{QQFieldElem}(collect(v)) for v in vertices(P1))
    squared = [only(k for k in 1:21 if any(e -> e[k] == 2, exponents(f))) for f in paper1[28:35]]
    check("the 7 linear binomials: distinct points among the 21 columns", 14, length(unique(point1.(1:21))))
    check("the last 8 quadrics: distinct points of their squared variables, none a vertex; vertices; dimension (a 5-simplex)",
          (8, true, 6, 5), (length(unique(point1.(squared))), all(k -> !(point1(k) in vertex_set), squared), length(vertex_set), dim(P1)))
    classes = [[k for k in 1:21 if point1(k) == x] for x in unique(point1.(1:21))]
    info("the 21 coordinates on the 14 points: " * join([join(pair_name.(c), " = ") for c in classes if length(c) > 1], ", ") *
         ", and " * join([pair_name(only(c)) for c in classes if length(c) == 1], ", ") * " alone")
    info("the 6 vertices: " * join([pair_name(k) for k in 1:21 if point1(k) in vertex_set], ", "))
    midpoints = map(zip(squared, paper1[28:35])) do (k, f)
        a, b = findall(==(1), only(e for e in exponents(f) if !any(==(2), e)))
        2 * point1(k) == point1(a) + point1(b) ? "2 $(pair_name(k)) = $(pair_name(a)) + $(pair_name(b))" : "$(pair_name(k)) (not a midpoint)"
    end
    info("the squared variables are midpoints: " * join(midpoints, ", "))
    check("type 1 is the Chevalley polytope of Spacek and Wang (census column chevalley_equivalent)", "true", type1.row["chevalley_equivalent"])
    lecture = [1, 3, 2, 3, 1]
    lecture_hall = convex_hull(QQ, [[k > 5 - step ? lecture[k] : 0 for k in 1:5] for step in 0:5])
    claim("P_C is affinely unimodularly equivalent to the lecture hall simplex P_(1,3,2,3,1) (Table, row 1)",
          affine_unimodular_equivalent(P1, lecture_hall))
    stored("Spacek and Wang's G2 polytope (Spacek_Wang/appendix_chevalley.jl)",
           strip(only(l for l in readlines(studies("g2", "degenerations_grassmannian", "Spacek_Wang", "test_summary.txt"))
                      if occursin("lecture-hall simplex", l))))

    subheading("Example ex:twentyone: type 21")
    metrics = Dict(parse(Int, r["number"]) => r for r in rows(g2("maximal_cell_tree_metrics.csv")))
    orbit21 = qq_orbit_of("21")
    cellsqq = fan_cells(fan("grassmannian"))
    triangle_tex = read(repository("figures", "triangle_cell.tex"), String)
    triangle_header = match(r"cell (\d+): rays ([\d,]+), orbit", triangle_tex)
    cell21, rays21 = parse(Int, triangle_header[1]), replace(triangle_header[2], ',' => ';')
    printed = match(r"\n% The cell\n\\node\[g2square cell name[^\n]*\\oldstylenums\{(\d+)\}", triangle_tex)[1]
    check("Figure triangle_cell (figures/triangle_cell.tex: cell $cell21, rays $(triangle_header[2])): the type printed on the triangle; " *
          "its rays in trop(G), those of the census cell (QQ orbit $orbit21), the support of the tree metric of type 21",
          (printed, rays21, rays21, rays21), (type_of_cell[cell21], join(cellsqq[cell21], ';'), census[orbit21]["rays"], metrics[21]["support"]))
    first_nontoric = first(s6_drawn_trees(read(repository("figures", "non_toric_cell_trees.tex"), String), panel))
    on_triangle = only(t for (name, t) in s6_drawn_trees(triangle_tex, r"\n% (The cell)\n") if name == "The cell")
    check("the first tree in Figure atlas_toric_cells_nontoric (its type) equals that on the triangle in Figure triangle_cell " *
          "(drawn without lengths)", ("21", s6_tree_text(on_triangle)), (first(first_nontoric), "splits " * join(join.(last(first_nontoric).splits), " ")))
    check("it is the tree of cell $cell21 (04_trees.csv)", s6_tree_text(last(first_nontoric)), s6_tree_text(s6_tree(tree_rows[cell21])))
    d21 = [parse(Int, metrics[21]["d$(i)$(j)"]) for (i, j) in PAIRS]
    zero_pairs = [PAIRS[k] for k in 1:21 if d21[k] == 0]
    check("its metric: d15 = d26 = d34 = 0 and d_ij = 2 otherwise", ([(1, 5), (2, 6), (3, 4)], [2]),
          (zero_pairs, unique(d for d in d21 if d != 0)))
    check("maximal cells of trop(G) with -d in their relative interior", [cell21], locate(locator(fan("grassmannian")), -QQ.(d21)))
    check("in_d(G) is prime over QQ and not toric (census)", ["true", "false"], [census[orbit21]["prime"], census[orbit21]["toric"]])
    initial21 = initial(I, tropical_semiring_map(QQ, max), d21)
    x, y = p(1, 2)^2 * p(1, 3), p(1, 4)^2 * p(1, 6)
    claim("(p12² p13)² + (p12² p13)(p14² p16) + (p14² p16)² lies in in_d(G)", iszero(normal_form(x^2 + x * y + y^2, initial21)))
    claim("neither factor x - ωy, x - ω²y of x² + xy + y² (x = p12² p13, y = p14² p16) lies in in_d(G) ⊗ C: " *
          "NF(x), NF(y) are not both 0 (a QQ Gröbner basis stays one over Q(ω))",
          !(iszero(normal_form(x, initial21)) && iszero(normal_form(y, initial21))))
    report = replace(join(readlines(studies("g2", "supplemental", "geometric_toric_report.md")), " "), r"\s+" => " ")
    stored("over CC (supplemental/verify_geometric_toric.jl)",
           replace(match(r"Thus the ideal has .*? coordinates\.", report).match, "**" => "", "`" => ""))

    subheading("The maximal cones whose initial ideal is not toric (types 21-35 of Figure atlas_toric_cells_nontoric)")
    census_columns(fields...) = @sprintf("%5s %-6s %-8s %-6s %3s %5s", fields...)
    header = census_columns("orbit", "prime", "binomial", "toric", "GB", "terms")
    println("           type  size   QQ: $header   GF(2): $header")
    for r in orbit_rows
        parse(Int, r["orbit_type"]) >= 21 || continue
        columns = map(((r["qq_orbit"], census), (r["gf2_orbit"], census2))) do (o, c)
            isempty(o) ? census_columns("-", "", "", "", "", "") :
                census_columns(o, c[o]["prime"], c[o]["binomial"], c[o]["toric"], c[o]["groebner_generators"], c[o]["maximum_terms"])
        end
        @printf("           %-5s %4s       %s          %s\n", r["orbit_type"], r["size"], columns[1], columns[2])
    end
    degeneration_count(target, key) = summary(target)["degenerations"][key]
    for (t, l) in (("grassmannian", "G over QQ"), ("grassmannian_gf2", "G over GF(2)"), ("flag", "F over QQ"), ("flag_gf2", "F over GF(2)"))
        info("$l: $(degeneration_count(t, "prime_orbits")) prime and $(degeneration_count(t, "toric_orbits")) toric orbits of maximal cones " *
             "($(degeneration_count(t, "prime_cells")) and $(degeneration_count(t, "toric_cells")) cones)")
    end

    subheading("Theorem thm:NObody: toric degenerations")
    toric(target) = (degeneration_count(target, "toric_cells"), degeneration_count(target, "toric_orbits"))
    istoric(c, o) = !isempty(o) && c[o]["toric"] == "true"
    # The toric orbits over QQ that stay toric orbits of the same size over GF(2), and the further toric orbits over GF(2).
    function in_any_characteristic(target, file)
        qqc = Dict(r["orbit"] => r for r in rows(g2(target, "06_canonical_interior_degenerations.csv")))
        gfc = Dict(r["orbit"] => r for r in rows(g2(target * "_gf2", "06_canonical_interior_degenerations.csv")))
        cells = [r for r in rows(g2(file)) if r["kind"] == "cell"]
        toric_qq = [r for r in cells if istoric(qqc, r["qq_orbit"])]
        stays = all(r -> istoric(gfc, r["gf2_orbit"]) && gfc[r["gf2_orbit"]]["orbit_size"] == qqc[r["qq_orbit"]]["orbit_size"], toric_qq)
        extra = sort([(r["orbit_type"], parse(Int, gfc[r["gf2_orbit"]]["orbit_size"])) for r in cells
                      if istoric(gfc, r["gf2_orbit"]) && !istoric(qqc, r["qq_orbit"])]; by=e -> parse(Int, e[1]))
        return length(toric_qq), stays, extra
    end
    check("G in any characteristic: (degenerations, orbits)", (135, 17), toric("grassmannian"))
    n, stays, extra = in_any_characteristic("grassmannian", "grassmannian_orbits.csv")
    claim("G: the $n toric orbits over QQ are toric orbits of the same size over GF(2)", stays)
    check("G in characteristic 2: further (degenerations, orbits)", (12, 3), toric("grassmannian_gf2") .- toric("grassmannian"))
    check("G in characteristic 2: the further types and their orbit sizes (Table, Char 2)", [("18", 6), ("19", 3), ("20", 3)], extra)
    check("all 20 orbits are listed in Table tab:g2-grassmannian-newton-okounkov: its rows", 20, length(paper_rows))
    check("F in any characteristic: (degenerations, orbits)", (807, 81), toric("flag"))
    n, stays, extra = in_any_characteristic("flag", "flag_orbits.csv")
    claim("F: the $n toric orbits over QQ are toric orbits of the same size over GF(2)", stays)
    check("F in characteristic 2: further (degenerations, orbits)", (18, 2), toric("flag_gf2") .- toric("flag"))
    check("F in characteristic 2: the orbits toric over GF(2) only (degenerations, orbits)", (18, 2),
          (sum(last, extra; init=0), length(extra)))
    info("F in characteristic 2: the further orbits lie over the Grassmannian types " *
         join(("$t ($s cones)" for (t, s) in extra), " and "))

    subheading("Table tab:g2-grassmannian-newton-okounkov, recomputed")
    gf2 = toric_orbits("grassmannian_gf2", "gf2")
    T = run_in_module(repository("julia", "to_tex", "toric_degeneration_table.jl"); name=:ToricTable)
    plain(tex) = replace(replace(tex, r"\^\{([^}]*)\}" => s"^(\1)"), r"[${}]" => "")
    println("           type  size  f-vector            N  volume  normal  char  Newton-Okounkov polytope")
    seen = Dict{String,Any}()
    bodies = Dict{String,Any}()
    chosen = Dict{String,Any}()   # the census row of each type, over QQ when there is one
    for type in string.(1:20)
        over_qq, over_gf2 = [o for o in qq if o.type == type], [o for o in gf2 if o.type == type]
        o, target = isempty(over_qq) ? (first(over_gf2), "grassmannian_gf2") : (first(over_qq), "grassmannian")
        chosen[type] = (o, target)
        bodies[type] = body(target, o.row)
        # The toric variety of the 14 distinct points of the degeneration is normal when these are all the lattice
        # points of P_C and P_C has the integer decomposition property.
        N = parse(Int, o.row["body_lattice_points"])
        normal = N == 14 && has_idp(bodies[type])
        char = isempty(over_qq) ? "2" : isempty(over_gf2) ? "0" : "all"
        f = replace(o.row["body_f_vector"], ';' => ',')
        polytope = call(T, :toric_table_type_tex, o.row["body_type"], parse.(Int, split(o.row["body_f_vector"], ';')),
                        o.row["chevalley_equivalent"] == "true")
        seen[type] = (; N, normal, size=o.row["orbit_size"], f, char)
        @printf("           %-5s %4s  %-18s %3d  %6s  %-6s  %-4s  %s\n", type, o.row["orbit_size"], f, N,
                o.row["body_normalized_volume"], normal, char, plain(polytope))
    end
    computed_rows = Dict(t => (v.size, v.f, string(v.N), v.char) for (t, v) in seen)
    differ = sort([t for t in union(keys(paper_rows), keys(computed_rows)) if get(paper_rows, t, nothing) != get(computed_rows, t, nothing)];
                  by=t -> parse(Int, t))
    claim("its rows agree with the census: type, orbit size, f-vector, N and Char", isempty(differ),
          isempty(differ) ? "$(length(paper_rows)) rows" : "rows that differ: $(join(differ, ", "))")
    check("caption: the orbit sizes of the 20 rows (Char all and 2) sum to 147", 147,
          sum(parse(Int, v[1]) for v in values(paper_rows)))
    info("the 17 rows with Char all sum to $(sum(parse(Int, v[1]) for v in values(paper_rows) if v[4] == "all")), " *
         "the 3 rows with Char 2 to $(sum(parse(Int, v[1]) for v in values(paper_rows) if v[4] == "2"))")
    claim("tables/toric_degeneration_table.tex is the current output of julia/to_tex/toric_degeneration_table.jl (with the polytope column)",
          call(T, :toric_degeneration_table_tex) == read(repository("tables", "toric_degeneration_table.tex"), String))
    if QUICK
        skipped("types 1-17: the polytopes over QQ and over GF(2) are affinely unimodularly equivalent", "--quick")
        line = only(l for l in readlines(studies("g2", "reports", "toric_body_equivalence.txt")) if startswith(l, "equivalent pairs"))
        m = match(r"orbits (\d+) and (\d+)", line)
        type_of_qq = Dict(r["qq_orbit"] => r["orbit_type"] for r in orbit_rows if !isempty(r["qq_orbit"]))
        stored("pairs of the 17 polytopes over QQ that are affinely unimodularly equivalent (supplemental/toric_body_equivalence.jl)",
               line * (m === nothing ? "" : "; QQ orbits $(m[1]) and $(m[2]) are types $(type_of_qq[m[1]]) and $(type_of_qq[m[2]])"))
    else
        gf2_of = Dict(o.type => o for o in gf2)
        claim("types 1-17: the polytopes over QQ and over GF(2) are affinely unimodularly equivalent (one row for all characteristics)",
              all(o -> affine_unimodular_equivalent(body("grassmannian", o.row), body("grassmannian_gf2", gf2_of[o.type].row)), qq))
        equivalent = [(a, b) for a in string.(1:20) for b in string.(parse(Int, a)+1:20) if affine_unimodular_equivalent(bodies[a], bodies[b])]
        check("pairs of the 17 polytopes over QQ that are affinely unimodularly equivalent (of 136)", [("13", "14")],
              [(a, b) for (a, b) in equivalent if parse(Int, b) <= 17])
        info("among all 20 polytopes (types 18-20 over GF(2)), the equivalent pairs: " * join(("$a ~ $b" for (a, b) in equivalent), ", "))
    end

    subheading("Remark: normality and the integral closure")
    claim("a toric degeneration is not normal exactly when |P ∩ Z^5| is 15 or 16",
          all(v -> v.normal == !(v.N in (15, 16)), values(seen)))
    for (target, field) in (("grassmannian", "qq"), ("grassmannian_gf2", "gf2"))
        counts = map(toric_orbits(target, field)) do o
            V = load(joinpath(g2(target), o.row["degree_one_values"]))
            length(unique([[V[i, j] for j in 1:ncols(V)] for i in 1:nrows(V)]))
        end
        check("$(field == "qq" ? "over QQ" : "over GF(2)"): the toric ring has 14 generators (distinct points of the 21 columns) " *
              "for every type", [14], unique(counts))
    end
    # The lattice points of P_C outside the 14 generators, at the census cell of each type that is not normal, written as
    # roots of products of generators: (p_k p_l)^(1/2) when the point is the midpoint of two generators, else a cube root.
    monomial_name(ks) = join([pair_name(k) * ("", "²", "³")[count(==(k), ks)] for k in sort(unique(ks))], " ")
    println("           type  lattice points of P_C outside the 14 generators, as roots of products of generators (census cell)")
    for type in sort([t for (t, v) in seen if !v.normal]; by=t -> parse(Int, t))
        o, target = chosen[type]
        V = load(joinpath(g2(target), o.row["degree_one_values"]))
        point = k -> [QQ(V[k, j]) for j in 1:ncols(V)]
        generators = unique(point.(1:21))
        outside = [q for q in (Vector{QQFieldElem}(collect(x)) for x in lattice_points(bodies[type])) if !(q in generators)]
        roots = map(outside) do q
            two = [(k, l) for k in 1:21 for l in k:21 if point(k) + point(l) == 2q]
            isempty(two) || return join(["($(monomial_name([k, l])))^(1/2)" for (k, l) in two], " = ")
            three = [(k, l, m) for k in 1:21 for l in k:21 for m in l:21 if point(k) + point(l) + point(m) == 3q]
            isempty(three) && return "neither a square nor a cube root"
            return "($(monomial_name(collect(first(three)))))^(1/3)" * (length(three) > 1 ? " and $(length(three) - 1) more" : "")
        end
        @printf("           %-5s %s\n", type, join(roots, "; "))
    end
    # Type 6: the elimination ideal of in_C(G) in p45, p56, p57, p67, at each of the 3 cells of the orbit, in the
    # coordinates p_ij of the paper. Where it is principal, z = p57 p67 / p56 has z² = ±p45 p67; the paper's certificate
    # needs the generator p45 p56² - p57² p67, z² = p45 p67 and z ∉ R, that is p57 p67 ∉ in_C(G) + (p56). (Until
    # 2026-11-11 the paper said type 2, where the generator is p45 p56² + p57² p67 and z ∈ R.)
    raysqq, _ = fan_rays(fan("grassmannian"))
    members_of(type) = parse.(Int, split(only(r["qq_members"] for r in orbit_rows if r["orbit_type"] == type), ';'))
    initials = Dict{Int,Any}()   # in_C(G), at a point of relint(C)
    initial_at(cell) = get!(() -> initial(I, tropical_semiring_map(QQ), sum(raysqq[r] for r in cellsqq[cell])), initials, cell)
    kept = [p(4, 5), p(5, 6), p(5, 7), p(6, 7)]
    others = [v for v in gens(P) if !(v in kept)]
    minus, plus = p(4, 5) * p(5, 6)^2 - p(5, 7)^2 * p(6, 7), p(4, 5) * p(5, 6)^2 + p(5, 7)^2 * p(6, 7)
    generator_name(E) = E == ideal(P, [minus]) ? "p45 p56² - p57² p67" : E == ideal(P, [plus]) ? "p45 p56² + p57² p67" : string(gens(E))
    z_outside(J) = !iszero(normal_form(p(5, 7) * p(6, 7), J + ideal(P, [p(5, 6)])))
    certificate(J) = eliminate(J, others) == ideal(P, [minus]) && z_outside(J) &&
                     iszero(normal_form(p(5, 7)^2 * p(6, 7)^2 - p(4, 5) * p(5, 6)^2 * p(6, 7), J))
    principal, zero_cells, z_inside = Int[], Int[], String[]
    generators6 = String[]
    for cell in members_of("6")
        initial6 = initial_at(cell)
        elimination = eliminate(initial6, others)
        all(iszero, gens(elimination)) && (push!(zero_cells, cell); continue)
        push!(principal, cell)
        push!(generators6, generator_name(elimination))
        z_outside(initial6) && continue
        equal = ["$(s == 1 ? "" : "-")$v" for v in gens(P) for s in (1, -1) if iszero(normal_form(p(5, 7) * p(6, 7) - s * p(5, 6) * v, initial6))]
        square = [s == 1 ? "p45 p67" : "-p45 p67" for s in (1, -1)
                  if iszero(normal_form(p(5, 7)^2 * p(6, 7)^2 - s * p(4, 5) * p(5, 6)^2 * p(6, 7), initial6))]
        push!(z_inside, "at cell $cell z = $(isempty(equal) ? "an element of R" : first(equal)) with z² = $(join(square, ", "))")
    end
    certified = filter(cell -> certificate(initial_at(cell)), members_of("6"))
    claim("type 6, at a cell of its orbit: the elimination ideal of in_C(G) in p45, p56, p57, p67 is generated by " *
          "p45 p56² - p57² p67, and z = p57 p67 / p56 has z² = p45 p67 and lies outside R (the certificate of the paper)",
          !isempty(certified), "at cells $certified")
    info("the elimination ideal is principal at cells $principal (generators $(join(unique(generators6), ", "))) and " *
         "zero at cells $zero_cells; $(join(z_inside, ", "))")
    # Where the certificate of the paper holds as written: at all cells of the toric types over QQ.
    holds = String[]
    for o in sort(qq; by=o -> parse(Int, o.type))
        cells = filter(cell -> certificate(initial_at(cell)), members_of(o.type))
        isempty(cells) || push!(holds, "$(o.type) (cell$(length(cells) == 1 ? "" : "s") $(join(cells, ", ")))")
    end
    info("the certificate of the paper (elimination ideal (p45 p56² - p57² p67), z² = p45 p67, z ∉ R) holds as written in " *
         "the toric types " * join(holds, ", ") * " (all $(length(initials)) cells of the toric types over QQ)")

    subheading("Comparison with the literature (stored reports of degenerations_grassmannian/)")
    D = studies("g2", "degenerations_grassmannian")
    sw = first(l for l in readlines(joinpath(D, "Spacek_Wang", "compare_flag.txt")) if occursin("summand [0, 1]", l))
    stored("type 1 and the Chevalley polytope of Spacek-Wang [Section A.4] (Spacek_Wang/compare_flag.jl)", strip(sw))
    cited("Christophersen and Ilten [Theorem 3.1]: G degenerates to P(K * Δ2) for a 2-sphere K with 11 vertices, one for each term order")
    ci = readlines(joinpath(D, "Christophersen_Ilten", "verify.txt"))
    stored("Christophersen-Ilten: term orders whose 2-sphere K is T_11 (of $(count(l -> occursin(r"^order \d+ signs", l), ci)) term orders; " *
           "Christophersen_Ilten/verify.jl)", count(l -> occursin("isomorphic to T_11", l), ci); paper=2)
    below_t11 = Tuple{String,String}[]
    for l in ci
        m = match(r"^order (\d+)( \(T_11\))?: cells (.*)$", l)
        m === nothing && continue
        types = join(sort(unique([x[1] for x in eachmatch(r"\(type (\d+)\)", m[3])]); by=t -> parse(Int, t)), ", ")
        info("order $(m[1])$(m[2] === nothing ? "" : " (T_11)"): toric types below it $types")
        m[2] === nothing || push!(below_t11, (m[1], types))
    end
    # main.tex: "Each of these degenerations factors through our toric degenerations among types 7, 10, 11, 12 and 17".
    for (order, types) in below_t11
        claim("T_11 order $order: the toric types its degeneration factors through lie among types 7, 10, 11, 12 and 17 " *
              "(Christophersen_Ilten/verify.jl, stored)", issubset(split(types, ", "), ["7", "10", "11", "12", "17"]), types)
    end
    line = only(l for l in ci if startswith(l, "toric types below a T_11 order"))
    stored("the toric types below some T_11 order, the union over the two (Christophersen_Ilten/verify.jl)",
           strip(split(line, ':')[2]); paper="7, 10, 11, 12, 17")
    claim("types 7, 10, 11, 12, 17 (below a T_11 order) are normal", all(t -> seen[t].normal, ("7", "10", "11", "12", "17")))
    cited("Christophersen and Ilten [Corollary 3.3]: V18 degenerates to P(T11 * Δ0) and the K3 surface of genus 10 to P(T11)")
    ffl = readlines(joinpath(D, "Feigin_Fourier_Littelmann", "verify.txt"))
    ffl_g = only(l for l in ffl if startswith(l, "types lattice equivalent to the FFL polytope P(omega_2)"))
    stored("the type whose polytope is lattice equivalent to the FFL polytope P(ω2) [FFL, Definition 11.5] " *
           "(Feigin_Fourier_Littelmann/verify.jl)", strip(split(ffl_g, ':')[2]); paper="7")
    info(only(l for l in ffl if startswith(l, "FFL polytope: ")))
end
