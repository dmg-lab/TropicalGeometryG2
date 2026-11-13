# Section 9 of the paper: run by paper_results.jl, which defines the helpers used here.

function section_9()
    subheading("Hexagons, the variety H and the map psi (eq:hexagon; studies/hexagon/complete_hexagon.jl)")
    cited("the incidence graph of (Q, G, F) has diameter 6 and girth 12, the split Cayley generalized hexagon; hexagons are " *
          "its 12-cycles, the apartments [Aschbacher §5; Van Maldeghem, Def. 2.4.9]")
    m = run_in_module(study_path("hexagon", "complete_hexagon.jl"); name=:ClassicalHexagon)   # includes polynomial_expressions.jl
    e(k) = call(m, :basis_vector, k)
    line(a, b) = call(m, :wedge, e(a), e(b))
    same(u, w) = call(m, :proportional, u, w)
    rotated(H1, H2) = all(same(H2.points[k], H1.points[mod1(k + 1, 6)]) && same(H2.lines[k], H1.lines[mod1(k + 1, 6)]) for k in 1:6)
    # With check=true, complete_hexagon throws an ArgumentError unless the inputs are flags and the result is a
    # non-degenerate hexagon whose twelve pairs (x^(k), p^((k-1)k)), (x^(k), p^(k(k+1))) are flags.
    complete(f1, f4) = try
        call(m, :complete_hexagon, f1, f4)
    catch err
        err isa ArgumentError ? nothing : rethrow()
    end
    order = [1, 5, 3, 4, 2, 6]
    A = complete((e(1), line(1, 5)), (e(4), line(4, 2)))
    claim("psi((e1, e15), (e4, e24)) is the hexagon (eq:hexagon) e1-e15-e5-e35-e3-e34-e4-e24-e2-e26-e6-e16: twelve flags, " *
          "six distinct points and six distinct lines",
          A !== nothing && all(same(A.points[k], e(order[k])) && same(A.lines[k], line(order[k], order[mod1(k + 1, 6)])) for k in 1:6))
    named(v, labels) = (v[findfirst(!iszero, v)] < 0 ? "-e" : "e") * labels[findfirst(!iszero, v)]
    A === nothing || info("psi returns x^(1), ..., x^(6) = " * join(named.(A.points, Ref(string.(1:7))), ", ") *
                          " and p^(12), ..., p^(61) = " * join(named.(A.lines, Ref(["$i$j" for (i, j) in PAIRS])), ", "))
    purity = TOML.parsefile(g2("flag", "00_ideal_purity.toml"))
    check("dim H = 2 dim F, as psi makes H birational to F^2 (dim F = affine dimension - 2, data/g2/flag/00_ideal_purity.toml)",
          12, 2 * (purity["dimension"] - 2))
    P_VARS, X_VARS = latest(m, :P_VARS), latest(m, :X_VARS)
    info("the ideal of H lives in 6 × $(length(X_VARS)) + 6 × $(length(P_VARS)) = $(6 * (length(X_VARS) + length(P_VARS))) " *
         "variables (main.tex: 168); its generators are not computed")
    cited("G2 acts transitively on Q, G, F and on the non-degenerate part of H [Van Maldeghem, Thm. 5.2.9]; every " *
          "non-degenerate hexagon is g·(eq:hexagon), with pointwise stabilizer gTg^(-1) [Aschbacher, (5.8), (5.9)]")
    B = A === nothing ? nothing : complete(A.flags[2], A.flags[5])
    claim("psi of the opposite flags 2 and 5 of this hexagon gives it again (psi inverts the projection H -> F^2)",
          B !== nothing && rotated(A, B))
    Random.seed!(9)
    generic, failed = 0, false
    for _ in 1:10
        H1 = complete(call(m, :random_flag, QQ), call(m, :random_flag, QQ))
        H1 === nothing && continue
        H2 = complete(H1.flags[2], H1.flags[5])
        H2 !== nothing && rotated(H1, H2) ? (generic += 1) : (failed = true)
    end
    claim("random opposite flags over QQ: psi gives a non-degenerate hexagon, and psi of its flags 2 and 5 gives it again",
          generic > 0 && !failed, "$generic of 10 pairs")

    subheading("Theorem thm:completing_hexagons: val(psi(z)) = trop(psi)(val z) at lifts (studies/hexagon/tests/test_tropical_hexagon.jl)")
    # The tests module includes tropical_complete_hexagon.jl; the certificate it loads is used again in the example.
    h = QUICK ? nothing : run_in_module(study_path("hexagon", "tests", "test_tropical_hexagon.jl"); name=:TropicalHexagon)
    if QUICK
        skipped("val(psi(z)) = trop(psi)(val z) for random lifts z over QQ(t) of pairs of flags", "--quick")
    else
        Random.seed!(1)
        lifts_ok = try
            redirect_stdout(devnull) do
                Base.invokelatest(latest(h, :test_tropical_complete_hexagon); lifts=10, fan_pairs=0)
            end
            true
        catch err
            info(sprint(showerror, err))
            false
        end
        claim("val(psi(z)) = trop(psi)(val z), with trop(psi) from the certificate, for 10 random lifts z over QQ(t) of pairs of flags",
              lifts_ok)
    end

    subheading("Lemma lem:complete_hexagon: the polynomial formulas (studies/hexagon/polynomial_expressions.jl)")
    S, xy = polynomial_ring(QQ, vcat(["x$i" for i in 1:7], ["y$i" for i in 1:7]))
    xs, ys = xy[1:7], xy[8:14]
    quad(z) = call(m, :quadric, z)
    polar_form = xs[1] * ys[4] + xs[2] * ys[5] + xs[3] * ys[6] + xs[4] * ys[1] + xs[5] * ys[2] + xs[6] * ys[3] - 2 * xs[7] * ys[7]
    claim("(eq:polarform) q(x+y) - q(x) - q(y) = x1y4 + x2y5 + x3y6 + x4y1 + x5y2 + x6y3 - 2x7y7, the polar form of " *
          "polynomial_expressions.jl and the matrix POLAR",
          quad(xs .+ ys) - quad(xs) - quad(ys) == polar_form && call(m, :polar, xs, ys) == polar_form &&
          sum(POLAR[i, j] * xs[i] * ys[j] for i in 1:7, j in 1:7) == polar_form)
    cited("two points are at distance 6 iff not perpendicular; p^(12) and x^(4) are at distance 5 and have a unique shortest " *
          "path [Van Maldeghem, Rem. 2.4.5, Lemma 1.3.5]")
    H = call(m, :hexagon_expressions)
    bidegree(f) = call(m, :bidegree, f)
    I = call(m, :line_point_ideal)
    in_I(f) = iszero(normal_form(f, I))
    incidence(x, p) = call(m, :incidence_relations, x, p)
    pl(i, j) = call(m, :pl, P_VARS, i, j)            # p_ii = 0, p_ji = -p_ij
    X = X_VARS
    displayed_x2 = [X[4] * pl(1, i) + X[5] * pl(2, i) + X[6] * pl(3, i) + X[1] * pl(4, i) + X[2] * pl(5, i) +
                    X[3] * pl(6, i) - 2 * X[7] * pl(7, i) for i in 1:7]
    claim("x^(2) as displayed: (x^(2))_i = x4 p_1i + x5 p_2i + x6 p_3i + x1 p_4i + x2 p_5i + x3 p_6i - 2 x7 p_7i",
          displayed_x2 == H.x2)
    x2_on_p12 = all(in_I, incidence(H.x2, P_VARS))
    x2_orthogonal = in_I(call(m, :polar, H.x2, X_VARS))
    x4_on_p12 = all(in_I, incidence(X_VARS, P_VARS))
    claim("x^(2) lies on p^(12) and is orthogonal to x^(4); x^(2), x^(3) on Q (modulo the ideal J of Q × G)",
          x2_on_p12 && x2_orthogonal && in_I(call(m, :quadric, H.x2)) && in_I(call(m, :quadric, H.x3)))
    # main.tex:1824: "The intersection $p^{(12)} \cap (x^{(4)})^{\perp}$ is the point $x^{(2)}$": x^(4) is not on p^(12),
    # and x^(2) is the point of p^(12) orthogonal to x^(4).
    check("the point x^(2), as main.tex names it: x^(4) ∉ p^(12) (modulo J), and x^(2) ∈ p^(12) is orthogonal to x^(4)",
          "p^(12) ∩ (x^(4))^⊥", !x4_on_p12 && x2_on_p12 && x2_orthogonal ? "p^(12) ∩ (x^(4))^⊥" : "not determined")
    check("bidegrees (in p, in x) of x^(2), x^(3)", [(1, 1), (1, 2)], [only(unique(bidegree.(v))) for v in (H.x2, H.x3)])
    check("x^(3)_1: number of terms", 24, length(H.x3[1]))
    p_var(i, j) = P_VARS[PAIR_INDEX[(i, j)]]
    # the terms in lex order for p12 > ... > p67 > x1 > ... > x7, as printed (the ring stores them in another order)
    x3_terms = sort(collect(terms(H.x3[1])); by=t -> exponent_vector(t, 1), rev=true)
    check("x^(3)_1: its first three and last two terms (lex order), as displayed",
          [p_var(1, 2) * X[5] * X[7], p_var(1, 3) * X[6] * X[7], p_var(1, 4) * X[1] * X[7], p_var(6, 7) * X[1] * X[3],
           -2 * p_var(6, 7) * X[5] * X[7]],
          x3_terms[[1, 2, 3, end - 1, end]])
    info("x^(3)_1 = $(H.x3[1])")
    times(a, b) = call(m, :cross, a, b)               # <a × b, z> = Φ(a, b, z)
    claim("x^(2) × x^(3) = x^(3) × x^(4) = 0 modulo J", all(in_I, times(H.x2, H.x3)) && all(in_I, times(H.x3, X_VARS)))
    claim("p^(23) = x^(2) ∧ x^(3) and p^(34) = x^(3) ∧ x^(4) in G (linear forms modulo J, Plücker relations); " *
          "x^(2), x^(3) on p^(23); x^(3), x^(4) on p^(34)",
          all(in_I, call(m, :linear_forms, H.p23)) && all(in_I, call(m, :linear_forms, H.p34)) &&
          all(iszero, call(m, :plucker_relations, H.p23)) && all(iszero, call(m, :plucker_relations, H.p34)) &&
          all(iszero, incidence(H.x2, H.p23)) && all(iszero, incidence(H.x3, H.p23)) && all(iszero, incidence(H.x3, H.p34)) &&
          all(iszero, incidence(X_VARS, H.p34)))
    check("bidegrees of p^(23), p^(34)", [(2, 3), (1, 3)], [only(unique(bidegree.(v))) for v in (H.p23, H.p34)])
    info("terms of the coordinates: x^(2) $(length.(H.x2)), x^(3) $(length.(H.x3))")
    info("p^(23) $(length.(H.p23)), p^(34) $(length.(H.p34))")

    subheading("Theorem thm:completing_hexagons: the certificate (studies/hexagon/certificate/; verification/ stored by verify_certificate.jl)")
    E = base_ring(I)
    check("variables of the ring of J: x1, ..., x7, p12, ..., p67", 28, ngens(E))
    G = load(g2("grassmannian", "00_ideal.mrdi"))
    E_names = string.(symbols(E))
    to_E = hom(base_ring(G), E, [E[findfirst(==(s), E_names)] for s in string.(symbols(base_ring(G)))])
    claim("J (line_point_ideal of polynomial_expressions.jl) is the stored ideal of G (data/g2/grassmannian/00_ideal.mrdi) plus q(x)",
          ideal(E, vcat(to_E.(gens(G)), [call(m, :quadric, X_VARS)])) == I)
    R, x = polynomial_ring(QQ, ["x$i" for i in 1:7])
    T = tropical_hypersurface(tropical_polynomial(x[1] * x[4] + x[2] * x[5] + x[3] * x[6] - x[7]^2))
    coarsest = fan("grassmannian")
    check("maximal cones of trop(Q), and of trop(Q) × trop(G) in the coarsest structure (6 × 248 = 1488; main.tex counts the " *
          "6 × 251 = 1506 raw ones, below)", [6, 1488],
          [n_maximal_polyhedra(T), n_maximal_polyhedra(T) * n_maximal_cones(coarsest)])
    directory = study_path("hexagon", "certificate", "verification")
    part_rows(suffix) = vcat([rows(f) for f in sort(readdir(directory; join=true)) if endswith(f, suffix)]...)
    cones = part_rows("_cones.csv")
    representatives = part_rows("_representatives.csv")
    finished = part_rows("_finished.csv")
    maximal = [r for r in cones if r["maximal"] == "true"]
    cells = sum(parse(Int, r["cells"]) for r in cones)
    stored("verify_certificate.jl on the certificate",
           "$(length(cones)) cones ($(length(maximal)) maximal), $cells cells, $(length(representatives)) representatives")
    # main.tex: "subdivides each of the $6\times 251 = 1506$ maximal Gr\"obner cones of $\trop(\pointspace)\times \trop(\linespace)$".
    check("maximal cones of the certificate, on the raw fan of trop(G): its faces, the cones of trop(Q), the pairs " *
          "(main.tex: 6 × 251 = 1506, the 251 raw maximal cones of the proof of Theorem thm:grassmannian)", [251, 6, 6 * 251],
          [length(unique(r["face"] for r in maximal)), length(unique(r["quadric"] for r in maximal)), length(maximal)])
    # Both fans have the same lineality space, and each raw maximal cone (simplicial modulo lineality) lies in the coarsest
    # maximal cone that contains the sum of its rays in its relative interior: its rays lie in that closed cone.
    L = locator(coarsest)
    raw = load(g2("grassmannian", "01_raw_fan.mrdi"))
    raw_rays, raw_lineality = fan_rays(raw)
    coarsest_lineality = last(fan_rays(coarsest))
    same_lineality = rank_of(raw_lineality) == rank_of(coarsest_lineality) == rank_of(vcat(raw_lineality, coarsest_lineality))
    in_closed(k, v) = all(iszero, L[k].equations * v) && all(<=(0), L[k].facets * v)
    containing = map(fan_cells(raw)) do cell
        inner = locate(L, sum(raw_rays[cell]))
        length(inner) == 1 && all(v -> in_closed(only(inner), v), raw_rays[cell]) ? only(inner) : 0
    end
    covered = [count(==(k), containing) for k in eachindex(L)]
    check("the 251 raw maximal cones subdivide the 248 coarsest ones, so the 1506 certified cones subdivide the 1488: each lies " *
          "in one; coarsest cones containing one, two of them (proof of Theorem thm:grassmannian: three pairs have the same span)",
          (true, 245, 3), (same_lineality && all(>(0), containing), count(==(1), covered), count(==(2), covered)))
    # The representatives f_σ = x^(-u) g of the certificate, with x^u f - g computed here: from the certificate that the tests
    # module has loaded (full run), else read here.
    function E_polynomial(p)                       # a polynomial of another ring with the same 28 variables, in E
        B = MPolyBuildCtx(E)
        foreach((c, e) -> push_term!(B, QQ(c), e), AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
        return finish(B)
    end
    certificate = if QUICK
        file = load(study_path("hexagon", "certificate", "hexagon_certificate.mrdi"))
        [(String(n), E_polynomial(g), Int.(u)) for (n, gs, us) in zip(file["coordinate"], file["g"], file["u"]) for (g, u) in zip(gs, us)]
    else
        loaded = call(h, :hexagon_certificate)
        [(n, E_polynomial(r.g), r.u) for (n, rs) in zip(loaded.names, loaded.reps) for r in rs]
    end
    coordinate(n) = (parts = split(n, '_'); getproperty(H, Symbol(parts[1]))[parse(Int, parts[2])])     # "p23_4" -> H.p23[4]
    differences = [prod(gens(E) .^ u) * coordinate(n) - g for (n, g, u) in certificate]
    check("coordinates with representatives in the certificate: x^(2), x^(3), p^(23), p^(34) = 7 + 7 + 21 + 21", 56,
          length(unique(first.(certificate))))
    group(n) = first(split(n, '_'))
    info("representatives by coordinate: " *
         join(["$g: $(count(r -> group(r[1]) == g, certificate))" for g in ("x2", "x3", "p23", "p34")], ", "))
    # main.tex:1871: "a tropically faithful representative $f_{\sigma}\in f+J$ on its relative interior", with J the ideal
    # of Q × G in the 28 variables x_k, p_ij (main.tex:1859).
    claim("main.tex writes f_σ ∈ f+J, J the ideal of Q × G: every representative f_σ lies in f + J (so do the stored rows " *
          "in_f_plus_I of verify_certificate.jl)",
          all(in_I, differences) && length(representatives) == length(certificate) &&
          all(r -> r["in_f_plus_I"] == "true", representatives), "$(length(certificate)) representatives")
    # Until 2026-11-11 main.tex wrote f + I, I the ideal of G: most representatives need the quadric q(x).
    I_G = ideal(E, to_E.(gens(G)))
    in_I_G = count(f -> iszero(normal_form(f, I_G)), differences)
    info("$in_I_G of the $(length(differences)) representatives already lie in f + I, I the ideal of G; the others need q(x)")
    cone_keys = [(r["face"], r["quadric"]) for r in cones]
    claim("the cones of trop(Q) in the certificate: the 11 sets of at least two of the terms x1x4, x2x5, x3x6, x7^2",
          Set(last.(cone_keys)) == Set(join(s) for s in Combinatorics.combinations(1:4) if length(s) >= 2))
    check("cones of the certificate: faces of the raw fan of trop(G) (1 + its f-vector (51, 244, 251) in the proof of Theorem " *
          "thm:grassmannian, the lineality space included), cones of trop(Q), pairs", [1 + 51 + 244 + 251, 11, 6017],
          [length(unique(first.(cone_keys))), length(unique(last.(cone_keys))), length(unique(cone_keys))])
    info("cones by dimension: " * join(["$d: $(count(r -> r["dim"] == d, cones))" for d in sort(unique(r["dim"] for r in cones); by=s -> parse(Int, s))], ", "))
    check("cones on which the saturation of in_v(J) is prime (stored flags J_prime of verify_certificate.jl)", 6017,
          count(r -> r["J_prime"] == "true", cones))
    claim("every check of the certificate passes on every cone (relative interior, maximal cone, Gröbner cone, no failures; " *
          "stored flags of verify_certificate.jl)",
          all(r -> r["point_in_relint"] == "true" && r["in_maximal_cone"] == "true" && r["groebner_cone"] == "true" &&
                   r["failures"] == "0", cones))
    claim("a subdivision tree for each of the 56 coordinates on every maximal cone (stored counts of verify_certificate.jl)",
          all(r -> r["trees"] == "56", maximal))
    claim("all parts of verify_certificate.jl finished, with no trees of cones that are not maximal (stored rows)",
          sort([parse(Int, r["part"]) for r in finished]) == collect(1:parse(Int, first(finished)["parts"])) &&
          all(r -> r["leftover_rows"] == "false", finished))

    subheading("The tropical hexagon of the example (data/g2/hexagon/tropicalized_hexagon.mrdi)")
    hexagon = load(g2("hexagon", "tropicalized_hexagon.mrdi"))
    maxconv(v) = Int.(maximum(v) .- v)
    check("a^(1), ..., a^(6) (max convention)", HEXAGON_POINTS, [maxconv(v) for v in hexagon.points])
    check("d^(12), ..., d^(61) (max convention)", HEXAGON_LINES, [maxconv(v) for v in hexagon.lines])
    if QUICK
        skipped("the completion of the hexagon from the flags (a^(1), d^(12)) and (a^(4), d^(45))", "--quick")
    else
        flag(k) = (Int.(hexagon.points[k]), Int.(hexagon.lines[k]))
        completion = call(h, :tropical_complete_hexagon, flag(1), flag(4))
        check("tropical_complete_hexagon from flags 1 and 4 gives the stored hexagon", (hexagon.points, hexagon.lines),
              (completion.points, completion.lines))
    end
    figure = Meta.parseall(read(project_path("to_tex", "hexagon_tree_ring.jl"), String))
    constant(name) = only(Base.eval(Main, a.args[1].args[2]) for a in figure.args
                          if a isa Expr && a.head == :const && a.args[1] isa Expr && a.args[1].args[1] == name)
    drawn = constant(:HEXAGON)
    claim("Figure fig:tropical_hexagon draws the stored hexagon (const HEXAGON of to_tex/hexagon_tree_ring.jl)",
          isequal((hexagon.points, hexagon.lines), ([QQ.(v) for v in drawn.points], [QQ.(v) for v in drawn.lines])))
    claim("its input flags INPUT_1, INPUT_4 are the flags 1 and 4",
          isequal([(hexagon.points[1], hexagon.lines[1]), (hexagon.points[4], hexagon.lines[4])],
                  [Tuple(QQ.(v) for v in constant(n)) for n in (:INPUT_1, :INPUT_4)]))
    types = orbit_types("grassmannian_orbits.csv", "qq")
    located = [locate(L, QQ.(v)) for v in hexagon.lines]
    claim("all six trees lie in the relative interior of maximal cones of trop(G)", all(c -> length(c) == 1, located))
    paper_types = ["1", "24", "10", "15", "12", "23"]
    check("their orbit types", paper_types, [types[("cell", only(c))] for c in located])
    cell_rows = Dict(r["orbit_type"] => r for r in rows(g2("grassmannian_orbits.csv")) if r["kind"] == "cell" && !isempty(r["qq_orbit"]))
    degeneration = Dict(r["orbit"] => r for r in rows(g2("grassmannian", "06_canonical_interior_degenerations.csv")))
    toric(type) = degeneration[cell_rows[type]["qq_orbit"]]["toric"] == "true"
    line_names = ["12", "23", "34", "45", "56", "61"]
    for (l, line_name) in enumerate(line_names)
        cell_type = types[("cell", only(located[l]))]
        splits = internal_splits(tree_from_min_pluecker(-QQ.(HEXAGON_LINES[l])))
        info("d^($line_name): cell $(only(located[l])), type $cell_type (orbit of $(cell_rows[cell_type]["size"])), " *
             (toric(cell_type) ? "toric" : "not toric") * "; splits " *
             join(["$(join(s.side))|$(join(s.complement)):$(s.length)" for s in splits], " "))
    end
    quadric_terms(a) = [a[1] + a[4], a[2] + a[5], a[3] + a[6], 2 * a[7]]
    minimal_terms(a) = findall(==(minimum(quadric_terms(a))), quadric_terms(a))
    claim("tropicalized hexagons, on the example: a^(1), ..., a^(6) lie in trop(Q), the minimum of x1+x4, x2+x5, x3+x6, 2x7 " *
          "is attained twice (min convention)", all(a -> length(minimal_terms(a)) >= 2, hexagon.points))
    term_names = ["x1x4", "x2x5", "x3x6", "x7^2"]
    info("terms attaining it: " *
         join(["a^($k) " * join(term_names[minimal_terms(a)], ", ") for (k, a) in enumerate(hexagon.points)], "; "))
    projection = Dict(parse(Int, r["flag_cell"]) => parse(Int, r["grassmannian_cell"])
                      for r in rows(g2("flag", "07_projection_cells.csv")))
    flag_cones = maximal_cones(fan("flag"))          # in R^7 × R^21, coordinates (x, p)
    flag_cells(k, l) = sort!([f for (f, c) in projection
                              if c == only(located[l]) && QQ.(vcat(hexagon.points[k], hexagon.lines[l])) in flag_cones[f]])
    incidences = [(k, l) for k in 1:6 for l in (mod1(k - 1, 6), k)]
    incidence_cells = Dict(i => flag_cells(i...) for i in incidences)
    claim("the twelve incidences (a^(k), d^(k-1,k)), (a^(k), d^(k,k+1)) lie in trop(F): in maximal cells of data/g2/flag " *
          "over the cell of the line", all(!isempty, values(incidence_cells)))
    flag_types = orbit_types("flag_orbits.csv", "qq")
    apple = Dict{Int,String}()
    for r in rows(g2("flag_orbits.csv")), member in split(r["qq_members"], ';'; keepempty=false)
        r["kind"] == "cell" && (apple[parse(Int, member)] = r["apple"])
    end
    info("flag cells of the incidences, with the apple on the atlas tree of their type (data/g2/flag_orbits.csv):")
    for (k, l) in incidences
        info("  (a^($k), d^($(line_names[l]))): " *
             join(["cell $f, type $(flag_types[("cell", f)]), $(apple[f])" for f in incidence_cells[(k, l)]], "; "))
    end
    trees_toric = [toric(types[("cell", only(c))]) for c in located]
    check("toric degeneration at each of the six trees (grassmannian/06_canonical_interior_degenerations.csv), as the figure of " *
          "its type says (main.tex: types 1-20 toric, fig:atlas_toric_cells_prime; 21-35 not, fig:atlas_toric_cells_nontoric)",
          [parse(Int, t) <= 20 for t in paper_types], trees_toric)
    claim("so the hexagon appears in both figures", any(trees_toric) && !all(trees_toric),
          "toric: d^(" * join(line_names[trees_toric], "), d^(") * "); not toric: d^(" * join(line_names[.!trees_toric], "), d^(") * ")")
end
