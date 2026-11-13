# Section 4 of the paper: run by paper_results.jl, which defines the helpers used here.

function section_4()
    subheading("Opening of the section: the coefficient 2 in the cubics of eq:threecubics")
    I = load(g2("grassmannian", "00_ideal.mrdi"))
    P = base_ring(I)
    names = string.(symbols(P))
    p(i, j) = P[findfirst(==("p$i$j"), names)]
    cubics = [2p(4, 5)p(4, 6)p(5, 6) + p(2, 3)p(3, 6)p(5, 6) - p(1, 4)p(1, 2)p(4, 5) + p(2, 6)p(1, 3)p(4, 5) - p(3, 5)p(1, 2)p(4, 6),
              2p(2, 3)p(1, 2)p(1, 3) - p(4, 6)p(1, 4)p(1, 3) + p(2, 5)p(5, 6)p(2, 3) + p(1, 6)p(4, 5)p(2, 3) - p(3, 4)p(5, 6)p(1, 2),
              2p(4, 6)p(5, 6)p(4, 5) + p(1, 2)p(2, 5)p(4, 5) - p(3, 6)p(1, 3)p(4, 6) + p(1, 5)p(2, 3)p(4, 6) - p(2, 4)p(1, 3)p(5, 6)]
    check("absolute values of the coefficients of each cubic: the leading 2, then four 1", fill([2, 1, 1, 1, 1], 3),
          [sort([Int(abs(numerator(c))) for c in coefficients(f)]; rev=true) for f in cubics])

    subheading("Theorem thm:char2, odd p: trop_p(G) = trop(G) (sampled here, proved below)")
    for field in ("QQ", "GF2", "GF3", "GF5", "GF7", "GF11", "GF13")
        scan = rows(study_path("g2", "char-p", "supplementary", "char2", "scan-$field.csv"))
        cells = [r for r in scan if startswith(r["point"], "cell")]
        tetrahedra = [r for r in scan if startswith(r["point"], "tetra")]
        stored("scan over $field (char-p/supplementary/char2/scan.jl): canonical weights of the 248 cells in trop",
               "$(count(r -> r["in_trop"] == "true", cells)) of $(length(cells)); interior points of the tetrahedra in trop: " *
               "$(count(r -> r["in_trop"] == "true", tetrahedra)) of $(length(tetrahedra))")
    end

    subheading("Theorem thm:char2, p = 2: trop_2(G)")
    s = summary("grassmannian_gf2")
    check("f-vector modulo lineality (length 3: a 2-dimensional complex)", [54, 253, 254], f_vector_of("grassmannian_gf2"))
    check("2-cells: (triangles, quadrilaterals)", [251, 3], [s["two_dimensional_faces"]["3-gon"], s["two_dimensional_faces"]["4-gon"]])
    info("maximal cells by type: " * join(["$c $t in $o orbits" for ((_, t), (c, o)) in sort(collect(maximal_types("grassmannian_gf2")))], ", "))
    check("reduced integral homology", "H_2 = Z^54", reduced_homology("grassmannian_gf2"))
    check("f-vector modulo D12", [12, 36, 32], s["complex"]["orbit_f_vector"])
    tqq, tgf2 = orbit_types("grassmannian_orbits.csv", "qq"), orbit_types("grassmannian_orbits.csv", "gf2")
    Fgf2 = fan("grassmannian_gf2")
    raysgf2, _ = fan_rays(Fgf2)
    cellsgf2 = fan_cells(Fgf2)
    gf2trees = rows(g2("grassmannian_gf2", "04_trees.csv"))
    check("each cell is encoded by a G2-tree: rays and maximal cells with a tree in grassmannian_gf2/04_trees.csv " *
          "(compared with the fan in section 5)", [54, 254],
          [count(r -> r["kind"] == kind && !isempty(r["pendant_lengths"]), gf2trees) for kind in ("ray", "cell")])
    cited("trop(Gr(3,7)) depends on the characteristic [speyer_sturmfels, Section 7]")
    cited("trop(Gr(2,n)) is characteristic-free: every tree lifts to a line in P^(n-1) [speyer_sturmfels, Theorem 7.1]")

    subheading("The three tetrahedra of the G2-Dressian in characteristic 0 and 2 (Figure fig:tetrahedra)")
    D = fan("dressian")
    drays, dlineality = fan_rays(D)
    dcells = fan_cells(D)
    # The quadrilaterals also have 4 rays; the tetrahedra are the cells of dimension 3 modulo lineality.
    tetrahedra = [cell for cell in dcells if rank_of(vcat(drays[cell], dlineality)) - rank_of(dlineality) == 4]
    check("tetrahedra", 3, length(tetrahedra))
    Fqq = fan("grassmannian")
    Lqq, Lgf2, LD = locator(Fqq), locator(Fgf2), locator(D)
    raysqq, _ = fan_rays(Fqq)
    cellsqq = fan_cells(Fqq)
    centers = [findall(r -> locate(LD, r) == [findfirst(==(cell), dcells)], raysgf2) for cell in tetrahedra]
    check("rays of trop_2(G) in the relative interior of each tetrahedron", [1, 1, 1], length.(centers))
    eps_index = only.(centers)
    check("their types", fill("epsilon", 3), [tgf2[("ray", e)] for e in eps_index])
    gf2_rays = [r for r in rows(g2("grassmannian_orbits.csv")) if r["kind"] == "ray" && isempty(r["qq_orbit"])]
    check("rays of trop_2(G) that are not rays of trop(G): the orbit ε, its size", ["epsilon", "3"],
          [only(gf2_rays)["orbit_type"], only(gf2_rays)["size"]])
    facets_of(cell) = Combinatorics.combinations(cell, 3)
    edges_of(cell) = Combinatorics.combinations(cell, 2)
    claim("p = 2: the 4·3 facets of the tetrahedra are deleted (their barycenters lie in no closed cone of trop_2(G))",
          all(isempty(locate_closed(Lgf2, sum(drays[f]))) for cell in tetrahedra for f in facets_of(cell)))
    gf2_scan = rows(study_path("g2", "char-p", "supplementary", "char2", "scan-GF2.csv"))
    check("cells of trop(G) whose canonical weight is not in trop_2(G) (Bayer test, scan-GF2.csv): the 4·3 facets", 12,
          count(r -> startswith(r["point"], "cell") && r["in_trop"] == "false", gf2_scan))
    qq_cells = Set(Set(raysqq[c]) for c in cellsqq)
    check("maximal cells of trop_2(G) that are maximal cells of trop(G): the 248 of trop(G) but the 4·3 facets", 248 - 4 * 3,
          count(c -> Set(raysgf2[c]) in qq_cells, cellsgf2))
    comparison = TOML.parsefile(g2("grassmannian_gf2", "08_summary.toml"))["comparison"]
    stored("grassmannian_gf2 vs $(comparison["reference"]) (08_summary.toml, census_summary): maximal cells (orbits) in both",
           "$(comparison["common_cells"]) ($(comparison["common_orbits"]))")
    # The cells of trop(G) on the facets, and the cells of trop_2(G) through ε and an edge.
    qq_facets = [reduce(vcat, [locate(Lqq, sum(drays[f])) for f in facets_of(cell)]) for cell in tetrahedra]
    gf2_triangles = [reduce(vcat, [locate(Lgf2, raysgf2[e] + sum(drays[f])) for f in edges_of(cell)]) for (cell, e) in zip(tetrahedra, eps_index)]
    star = [findall(c -> e in c, cellsgf2) for e in eps_index]
    check("maximal cells of trop_2(G) at each ε", [6, 6, 6], length.(star))
    claim("they are the six triangles joining ε to the six edges of the tetrahedron", sort.(gf2_triangles) == star)
    four = ["alpha_2", "alpha_3", "alpha_3", "alpha_4"]
    check("vertices of each tetrahedron: G2-trees of types α2, α4, α3, α3", fill(four, 3),
          [sort([tqq[("ray", v)] for v in union(cellsqq[ks]...)]) for ks in qq_facets])
    check("the other vertices of the six triangles at each ε: the vertices of the tetrahedron", fill(four, 3),
          [sort([tgf2[("ray", v)] for v in setdiff(union(cellsgf2[ks]...), [e])]) for (ks, e) in zip(star, eps_index)])
    by_number(types) = sort(types; by=t -> parse(Int, t))
    check("types of the 4 facets of each tetrahedron, p ≠ 2", fill(["32", "32", "33", "34"], 3),
          [by_number([tqq[("cell", k)] for k in ks]) for ks in qq_facets])
    check("types of the 6 triangles from the center ε to the edges, p = 2", fill(["18", "18", "19", "20", "35", "35"], 3),
          [by_number([tgf2[("cell", k)] for k in ks]) for ks in gf2_triangles])
    face_label(cells, t, k) = "(" * join(sort([t[("ray", v)] for v in cells[k]]), ", ") * "): " * t[("cell", k)]
    facet_tables = [sort([face_label(cellsqq, tqq, k) for k in ks]) for ks in qq_facets]
    triangle_tables = [sort([face_label(cellsgf2, tgf2, k) for k in ks]) for ks in gf2_triangles]
    info("the faces of one tetrahedron by the types of their vertices (not stated in the paper):")
    foreach(l -> info("p ≠ 2: facet " * l), facet_tables[1])
    foreach(l -> info("p = 2: triangle " * l), triangle_tables[1])
    claim("the same face by face for all three tetrahedra", allequal(facet_tables) && allequal(triangle_tables))

    subheading("Proof of thm:char2: the cells off the tetrahedra (Gröbner bases over ZZ: char-p/compute.jl, verify.jl)")
    claim("the 42 generators of I have coefficients ±1, so trop_p(G) lies in the G2-Dressian for every p",
          ngens(I) == 42 && all(g -> all(c -> abs(c) == 1, coefficients(g)), gens(I)))
    verified = Dict{Int,Dict{String,String}}()
    for file in readdir(study_path("g2", "char-p", "results"); join=true)
        startswith(basename(file), "verify-") || continue
        for r in rows(file)
            verified[parse(Int, r["cone"])] = r
        end
    end
    raw = load(g2("grassmannian", "01_raw_fan.mrdi"))
    raw_cones = collect(maximal_cones(raw))
    kinds = Dict(k => count(r -> r["torus_kind"] == k, values(verified)) for k in unique(r["torus_kind"] for r in values(verified)))
    stored("a Gröbner basis over ZZ at an interior point w of each cone of the raw fan, re-checked without Gröbner computations",
           "$(length(verified)) of $(length(raw_cones)) cones; S-pairs checked: $(sum(parse(Int, r["spairs_checked"]) for r in values(verified))); " *
           "torus points: $(kinds)")
    inner_point(C) = sum(QQ.(collect(r)) for r in rays_modulo_lineality(C).rays_modulo_lineality)
    # The boundary triangles of the tetrahedra are the raw cones that leave trop_2(G).
    on_tetrahedra = [k for (k, C) in enumerate(raw_cones) if isempty(locate(Lgf2, inner_point(C)))]
    check("raw cones on the boundary of the tetrahedra (char-p/README.md)", [71, 72, 101, 103, 121, 133, 134, 135, 136, 141, 162, 163],
          on_tetrahedra)
    others = setdiff(eachindex(raw_cones), on_tetrahedra)
    # The maximal cells of D that contain the interior points of the other raw cones.
    halves = Dict{Int,Vector{Int}}()
    for k in others, c in locate(LD, inner_point(raw_cones[k]))
        push!(get!(halves, c, Int[]), k)
    end
    split_cells = sort([c for (c, ks) in halves if length(ks) > 1])
    # main.tex: "(i) the Gr\"obner basis and the $42$ generators generate the same ideal over $\ZZ$, (ii) all leading
    # coefficients are $1$, (iii) the same Gr\"obner basis works for every $w$ in $\relint(C)$".
    stored("(ii) the other 251 − 4·3 cones: all leading coefficients 1 (verify-*.csv, column lcs)",
           count(k -> verified[k]["lcs"] == "1", others); paper=251 - 4 * 3)
    info("leading coefficients 1 and 2 (a power of 2, a unit for odd p) at the cones " *
         join(sort([k for (k, r) in verified if r["lcs"] != "1"]), ", "))
    initial_terms(g, w) = let ws = [sum(w .* e) for e in exponents(g)], m = minimum(ws)   # min convention
        Set((c, e) for (c, e, x) in zip(coefficients(g), exponents(g), ws) if x == m)
    end
    kept = Dict{Int,Any}()   # the stored Gröbner basis G and interior point w of the raw cones in split_cells
    if QUICK
        skipped("largest coefficients of the stored Gröbner bases", "--quick")
    else
        # The certificates were written by Oscar 1.9.0-DEV; Oscar 1.8.2 warns and rethrows when a load fails.
        largest = try
            Dict(k => begin
                     certificate = load(study_path("g2", "char-p", "certificates", "cones", @sprintf("cone_%03d.mrdi", k)))
                     any(c -> k in halves[c], split_cells) && (kept[k] = (certificate["G"], certificate["w"]))
                     maximum(abs(c) for g in certificate["G"] for c in coefficients(g))
                 end for k in others)
        catch e
            sprint(showerror, e)
        end
        if largest isa Dict
            # Only the leading coefficients matter (ii); the paper no longer claims that all coefficients are ±1.
            info("largest |coefficient| of the stored Gröbner basis over ZZ, over the other 251 − 4·3 cones: " *
                 join(["$m for $(count(==(m), values(largest)))" for m in sort(unique(values(largest)))], ", "))
        else
            skipped("largest coefficients of the stored Gröbner bases", "the certificates do not load: " * first(largest, 200))
        end
    end
    stored("(i) the Gröbner basis and the 42 generators generate the same ideal over ZZ, and every S-pair reduces, raw cones " *
           "(verify-*.csv: in_I, I_in_G, spairs_failed)",
           count(r -> r["in_I"] == "true" && r["I_in_G"] == "true" && r["spairs_failed"] == "0", values(verified)); paper=251)
    # main.tex: "for every cell $C$ of the $\Ggroup$-Dressian after subdividing the quadrilaterals": these cells are the 251
    # raw cones of 01_raw_fan, the Gröbner cones.
    stored("(iii) the same Gröbner basis on relint(C), for each of the 251 raw cones of 01_raw_fan (the cells of D after " *
           "subdividing the quadrilaterals): homogeneous for the lineality, initial terms tie at every ray (verify-*.csv: lin_hom, H, I)",
           count(r -> all(r[c] == "true" for c in ("lin_hom", "H", "I")), values(verified)); paper=251)
    info("the quadrilaterals that the paper subdivides are cells $(join(split_cells, ", ")) of D, with " *
         "$(join(unique(length(dcells[c]) for c in split_cells), ", ")) vertices, unions of the raw cones " *
         join([join(halves[c], ";") for c in split_cells], ", "))
    if QUICK
        skipped("the stored Gröbner basis of each raw cone of a quadrilateral fails (iii) on the other one", "--quick")
    elseif length(kept) == sum(c -> length(halves[c]), split_cells; init=0)
        moves(a, b) = any(g -> initial_terms(g, kept[a][2]) != initial_terms(g, kept[b][2]), kept[a][1])
        claim("the stored Gröbner basis of each raw cone of a quadrilateral has other initial forms at the interior point of the " *
              "other raw cone (char-p/certificates/cones)",
              all(c -> all(moves(a, b) for a in halves[c] for b in halves[c] if a != b), split_cells))
    else
        skipped("the stored Gröbner basis of each raw cone of a quadrilateral fails (iii) on the other one", "the certificates do not load")
    end
    good = [r for r in values(verified) if r["spairs_failed"] == "0" &&
            all(r[c] == "true" for c in ("in_I", "I_in_G", "lcs_powers_of_2", "lin_hom", "H", "I", "torus"))]
    stored("raw cones certified for every odd p: (ii), (iii), leading coefficients powers of 2, a torus point of in_w(G) (verify-*.csv)",
           length(good); paper=251)

    subheading("Proof of thm:char2: the open tetrahedra for p ≥ 3 (char-p/witness.jl)")
    witness = filter(l -> startswith(l, "tetrahedron "), readlines(study_path("g2", "char-p", "results", "witness.txt")))
    stored("p odd, open tetrahedra: a Weyl image of the cubic f in I over ZZ whose term ±2 p17 p27 p37 is strictly minimal",
           join(witness, "\n           "))
    stored("tetrahedra with such a witness (witness.txt)",
           count(l -> occursin("2-term strictly minimal on the open cell; in I over ZZ with cofactors: true", l), witness); paper=3)

    subheading("Proof of thm:char2: the cubics in characteristic 2")
    check("terms of each cubic modulo 2", [4, 4, 4], [count(c -> isodd(numerator(c)), coefficients(f)) for f in cubics])
    term_weights(f, v) = [(isodd(numerator(c)), sum(v[k] * a[k] for k in eachindex(a))) for (c, a) in zip(coefficients(f), exponents(f))]
    odd_weights(f, v) = [x for (o, x) in term_weights(f, v) if o]
    minimal_terms(f, v) = (ws = odd_weights(f, v); findall(==(minimum(ws)), ws))
    function ties_at(f, v)
        t = term_weights(f, v)
        odd, even = [x for (o, x) in t if o], [x for (o, x) in t if !o]
        return allequal(odd) && all(<(first(odd)), even)
    end
    eps_vectors = raysgf2[eps_index]
    for (c, f) in enumerate(cubics)
        info("cubic $c, weights of its terms (the 2-term first) at ε = " *
             join(["ray $e $(show_value([x for (_, x) in sort(term_weights(f, v); by=first)]))" for (e, v) in zip(eps_index, eps_vectors)], "; "))
    end
    matches = [findall(v -> ties_at(f, v), eps_vectors) for f in cubics]
    check("rays ε at which the four odd terms of a cubic tie and its term 2m is strictly smaller (min convention)", [1, 1, 1], length.(matches))
    if all(==(1), length.(matches))
        claim("the three cubics belong to the three tetrahedra", allunique(only.(matches)))
        edge_minima, facet_minima, missing_terms = Int[], Int[], Vector{Vector{Int}}[]
        for (f, m) in zip(cubics, matches)
            e = eps_vectors[only(m)]
            cell = dcells[only(locate(LD, e))]
            append!(edge_minima, [length(minimal_terms(f, e + sum(drays[g]))) for g in edges_of(cell)])
            append!(facet_minima, [length(minimal_terms(f, sum(drays[g]))) for g in facets_of(cell)])
            push!(missing_terms, sort([setdiff(1:4, minimal_terms(f, drays[a])) for a in cell]))
        end
        check("at each vertex of the tetrahedron three of the four terms are minimal, a different one missing at each: " *
              "near ε the tropical surface of the cubic mod 2 is the cone over the six edges", fill([[1], [2], [3], [4]], 3), missing_terms)
        check("the cubic mod 2 attains its minimum twice on the 6 triangles at ε", fill(2, 18), edge_minima)
        check("and once on the 4 facets of its tetrahedron", fill(1, 12), facet_minima)
    end

    subheading("Theorem thm:flagchar2: trop_2(F)")
    s = summary("flag_gf2")
    check("f-vector modulo lineality (length 4: a 3-dimensional complex)", [156, 1445, 3758, 2792], f_vector_of("flag_gf2"))
    mt = maximal_types("flag_gf2")
    cell_kinds = [(3, "3-simplex"), (3, "Pyramid(4-gon)"), (3, "Δ1 × Δ2")]
    check("3-cells by type: (3-simplex, Pyramid(4-gon), Δ1 × Δ2)", [2630, 144, 18], [get(mt, k, (0, 0))[1] for k in cell_kinds])
    info("their orbits: $(join([get(mt, k, (0, 0))[2] for k in cell_kinds], ", ")); types of maximal cells: $(length(mt))")
    check("reduced integral homology", "H_3 = Z^324", reduced_homology("flag_gf2"))
    check("f-vector modulo D12", [31, 192, 416, 281], s["complex"]["orbit_f_vector"])
    certificate2 = TOML.parsefile(g2("flag_gf2", "07_projection_certificate.toml"))
    claim("the exact projection certificate flag_gf2/07_projection_certificate.toml is valid (ray images, cells, D12-transport)",
          projection_certificate_valid("flag_gf2"))
    images2 = [c["target_cell"] for c in certificate2["cells"]]
    check("each maximal cell of trop_2(F) maps onto a maximal cell of trop_2(G); the images (onto)", [2792, 254],
          [length(images2), length(unique(images2))])
    claim("flag_gf2/07_projection_cells.csv lists the same images",
          [parse(Int, r["grassmannian_cell"]) for r in rows(g2("flag_gf2", "07_projection_cells.csv"))] == images2)
    ftrees2 = rows(g2("flag_gf2", "04_trees.csv"))
    check("rays with a G2-tree, vertical rays (the 7 star trees, empty row), maximal cells with a G2-tree (flag_gf2/04_trees.csv)",
          [156 - 7, 7, 2792],
          [count(r -> r["kind"] == "ray" && r["vertical"] == "false" && !isempty(r["pendant_lengths"]), ftrees2),
           count(r -> r["kind"] == "ray" && r["vertical"] == "true", ftrees2),
           count(r -> r["kind"] == "cell" && !isempty(r["pendant_lengths"]), ftrees2)])
    info("the decorations (apples on the edges of the trees of trop_2(G)) are checked over GF(2) in section 7")

    subheading("Proof of thm:flagchar2")
    fiber = reduce(vcat, [rows(f) for f in readdir(study_path("g2", "char-p", "supplementary", "flag"); join=true)
                          if startswith(basename(f), "fiber-")])
    for field in ("QQ", "GF3", "GF5", "GF2")
        sample = [r for r in fiber if r["field"] == field]
        stored("sampled over $field (char-p/supplementary/flag/fiber.jl): points (a, b) with b inside a raw cone of trop(G) and a on the " *
               "line L_b in trop; points with a off L_b in trop",
               "$(count(r -> r["in_trop"] == "true", sample)) of $(length(sample)); " *
               "$(count(r -> r["perturbed_on_line"] == "false" && r["perturbed_in_trop"] == "true", sample)) of " *
               "$(count(r -> r["perturbed_on_line"] == "false", sample))")
    end
    skipped("odd p: trop_p(F) = trop(F), a certificate like verify.jl",
            "the stored check char-p/supplementary/flag/fstep2.jl covers an older raw flag fan of 3879 cones")
    support = TOML.parsefile(g2("flag_gf2", "01_support_certificate.toml"))
    stored("p = 2: Algorithm alg:balanced_support over GF(2) (01_support_certificate.toml)",
           "$(support["certification_strategy"]): $(support["proof"]); orbits checked $(support["checked_orbits"]), kept $(support["selected_orbits"])")
    basis = filter(!isempty, readlines(g2("grassmannian_gf2", "01_tropical_basis.txt")))
    monomial_degree(m) = sum(occursin('^', x) ? parse(Int, split(x, '^')[2]) : 1 for x in split(strip(m), '*'))
    shape(line) = (monomial_degree(first(split(line, " + "))), length(split(line, " + ")))
    shapes = shape.(basis)
    check("tropical basis of I over GF(2) (grassmannian_gf2/01_tropical_basis.txt): the 42 generators (linear forms and Plücker " *
          "quadrics) and 12 additional polynomials (4 terms)", [42, 12],
          [count(s -> s[1] == 1 || s == (2, 3), shapes), count(s -> s[2] == 4, shapes)])
    info("the 42: $(count(s -> s[1] == 1, shapes)) linear forms and $(count(==((2, 3)), shapes)) quadrics with 3 terms; the 12: " *
         "$(count(==((2, 4)), shapes)) quadrics and $(count(==((3, 4)), shapes)) cubics with 4 terms")
    foreach(l -> info(l), [l for (l, s) in zip(basis, shapes) if s[2] == 4])
    candidates = filter(!isempty, readlines(g2("flag_gf2", "01_candidate_basis.txt")))
    J2 = load(g2("flag_gf2", "00_ideal.mrdi"))
    R2 = base_ring(J2)
    variables2 = Dict(string(v) => v for v in gens(R2))
    factor2(x) = occursin('^', x) ? variables2[split(x, '^')[1]]^parse(Int, split(x, '^')[2]) : variables2[x]
    # Over GF(2) every coefficient is 1: a line is a sum of products of variables.
    as_polynomial2(line) = sum(prod(factor2, split(strip(t), '*')) for t in split(strip(line), " + "))
    generators = ngens(J2)
    check("the set S over GF(2) (01_candidate_basis.txt): generators of the flag ideal, further polynomials", [78, 12],
          [generators, length(candidates) - generators])
    S2, J2gens = Set(as_polynomial2.(candidates)), Set(gens(J2))
    twelve = Set(as_polynomial2.([l for (l, s) in zip(basis, shapes) if s[2] == 4]))
    claim("S = the 78 generators of flag_gf2/00_ideal.mrdi and the 12 additional polynomials of the GF(2) tropical basis of I",
          length(S2) == length(candidates) && issubset(J2gens, S2) && setdiff(S2, J2gens) == twelve)
    prevariety = load(g2("flag_gf2", "01_prevariety_fan.mrdi"))
    check("cones of the prevariety of S; maximal cones after merging, their dimension, their orbits", [3038, 2792, 8, 281],
          [n_maximal_cones(prevariety), s["complex"]["maximal_cells"], dim(fan("flag_gf2")), s["complex"]["maximal_cell_orbits"]])
    check("prevariety of S: dimension, lineality dimension, pure (so none has dimension 9)", (8, 4, true),
          (dim(prevariety), lineality_dim(prevariety), is_pure(prevariety)))
    # f_vector(prevariety) takes about 3 minutes; census.md stores it.
    census_line = only(filter(l -> startswith(l, "| Candidate prevariety"), readlines(g2("flag_gf2", "census.md"))))
    fp = parse.(Int, split(strip(last(filter(!isempty, strip.(split(census_line, '|')))), ['(', ')']), ", "))[2:end]
    euler(f) = sum((-1)^k * f[k+1] for k in 0:length(f)-1) - 1
    stored("its f-vector modulo lineality (census.md, run_census_gf2.jl)",
           "$(show_value(fp)), reduced Euler characteristic $(euler(fp)); that of trop_2(F): $(euler(f_vector_of("flag_gf2")))")
    membership = rows(g2("flag_gf2", "01_support_membership.csv"))
    check("second loop: orbits tested, of dimension 8, kept; their cells", [281, 281, 281, 2792],
          [length(membership), count(r -> r["dimension"] == "8", membership), count(r -> r["tropical"] == "true", membership),
           sum(parse(Int, r["orbit_size"]) for r in membership)])

    subheading("Proof of thm:flagchar2: the output in the flag Dressian (thm:flag_dressian), and Figure fig:tetrahedra decorated")
    flag_rows = rows(g2("flag_orbits.csv"))
    rays_f = [r for r in flag_rows if r["kind"] == "ray" && !isempty(r["gf2_orbit"])]
    new_rays = [r for r in rays_f if isempty(r["qq_orbit"])]
    check("rays of trop_2(F): orbits and rays of trop(F) (the rays of the flag Dressian), new orbits and rays (31 − 27, 156 − 141)",
          [27, 141, 31 - 27, 156 - 141],
          [length(rays_f) - length(new_rays), sum(parse(Int, r["size"]) for r in rays_f if !isempty(r["qq_orbit"])),
           length(new_rays), sum(parse(Int, r["size"]) for r in new_rays)])
    check("type of the new rays", ["epsilon"], unique(r["orbit_type"] for r in new_rays))
    info("sizes of their orbits: " * join([r["size"] for r in new_rays], ", "))
    # The new cells of trop_2(F) lie in the 4-simplices of the flag Dressian.
    simplices = [C for C in maximal_cones(fan("flag_dressian")) if dim(C) == 9]   # 5 rays and the 4-dimensional lineality
    Lsimplices = [locator_entry(C) for C in simplices]
    rays2f, _ = fan_rays(fan("flag_gf2"))
    cells2f = fan_cells(fan("flag_gf2"))
    new_cells = [parse(Int, m) for r in flag_rows if r["kind"] == "cell" && isempty(r["qq_orbit"]) for m in split(r["gf2_members"], ';')]
    hosts = [locate(Lsimplices, sum(rays2f[cells2f[k]])) for k in new_cells]
    per_simplex = [count(==([s]), hosts) for s in eachindex(simplices)]
    claim("the maximal cells of trop_2(F) not in trop(F) lie in the relative interiors of the 33 4-simplices of the flag Dressian, " *
          "the same number in each", length(simplices) == 33 && all(h -> length(h) == 1, hosts) && allequal(per_simplex),
          "$(length(new_cells)) cells, $(join(unique(per_simplex), " or ")) in each; not stated in the paper")
    comparison = TOML.parsefile(g2("flag_gf2", "08_summary.toml"))["comparison"]
    stored("flag_gf2 vs $(comparison["reference"]) (08_summary.toml, census_summary): maximal cells (orbits) in both",
           "$(comparison["common_cells"]) ($(comparison["common_orbits"]))")
    above = Dict{String,Vector{Int}}()
    for (target, t, types) in (("flag", tqq, ("32", "33", "34")), ("flag_gf2", tgf2, ("18", "19", "20", "35")))
        per_cell = Dict{Int,Int}()
        for r in rows(g2(target, "07_projection_cells.csv"))
            k = parse(Int, r["grassmannian_cell"])
            t[("cell", k)] in types && (per_cell[k] = get(per_cell, k, 0) + 1)
        end
        info("$target, maximal cells over the triangles of type " *
             join(["$ty: $(sum(n for (k, n) in per_cell if t[("cell", k)] == ty))" for ty in types], ", "))
        above[target] = [length(per_cell); sort(unique(values(per_cell)))]
    end
    check("triangles of the tetrahedra and flag cells above each: char 0 (the 4·3 facets), char 2 (the 6·3 triangles at ε) " *
          "(10 and 11: the maximal cells of the fiber complexes 𝔄 and 𝔈 of Table tab:local_flag_complex)",
          [[4 * 3, 10], [6 * 3, 11]], [above["flag"], above["flag_gf2"]])
end
