# Section 3 of the paper: run by paper_results.jl, which defines the helpers used here.

"The basis of the lineality space of a fan, as a matrix over QQ."
s3_lineality(F) = matrix(QQ, rays_modulo_lineality(F).lineality_basis)

"Whether two fans have the same lineality space."
function s3_same_lineality(A, B)
    LA, LB = s3_lineality(A), s3_lineality(B)
    return rank(LA) == rank(LB) == rank(vcat(LA, LB))
end

"The linear span of each maximal cone of a fan, as the reduced row echelon form of its equations."
s3_spans(F) = [Matrix{QQFieldElem}(rref(matrix(QQ, Matrix{QQFieldElem}(linear_equation_matrix(linear_span(C)))))[2])
               for C in maximal_cones(F)]

"""The image of a ray of the flag Dressian (in R^7 x R^21) in trop(G): the signature of its Plücker part modulo the
lineality space of the Dressian (columns of `quotient`), or nothing when that part lies in the lineality space."""
function s3_plucker_image(ray, quotient)
    plucker = QQ.(collect(ray))[8:28]
    return all(iszero, vec(transpose(plucker) * quotient)) ? nothing : quotient_ray_signature(plucker, quotient)
end

"The census members of the rows of an orbit table, as integers."
s3_members(row, field) = isempty(row["$(field)_members"]) ? Int[] : parse.(Int, split(row["$(field)_members"], ';'))

function section_3()
    G, D = fan("grassmannian"), fan("dressian")
    raw = load(g2("grassmannian", "01_raw_fan.mrdi"))
    Dcones = collect(maximal_cones(D))

    subheading("The opening paragraph: trop(G) inside the G2-Dressian")
    check("dimensions in R^21/1 of the cones over the three tetrahedra of the Dressian (larger than dim G = 5)", [6, 6, 6],
          [dim(C) - 1 for C in Dcones if dim(C) == dim(D)])
    check("dimension of trop(G) in R^21/1", 5, dim(G) - 1)
    claim("trop(G) is pure", is_pure(G))
    info("that the tetrahedra are the only cells in the difference is checked in the proof of Theorem thm:grassmannian")

    subheading("Theorem thm:grassmannian: the tropical G2-Grassmannian trop(G)")
    s = summary("grassmannian")
    check("f-vector modulo lineality", [51, 241, 248], f_vector_of("grassmannian"))
    check("2-cells: (triangles, quadrilaterals)", [245, 3], [s["two_dimensional_faces"]["3-gon"], s["two_dimensional_faces"]["4-gon"]])
    check("reduced integral homology", "H_2 = Z^57", reduced_homology("grassmannian"))
    check("f-vector modulo D12", [11, 33, 31], s["complex"]["orbit_f_vector"])
    gtrees = rows(g2("grassmannian", "04_trees.csv"))
    check("each cell is encoded by a G2-tree: rays and maximal cells with a tree in grassmannian/04_trees.csv", [51, 248],
          [count(r -> r["kind"] == kind && !isempty(r["pendant_lengths"]), gtrees) for kind in ("ray", "cell")])

    subheading("Proof of thm:grassmannian: from the Groebner fan structure to the coarsest fan")
    claim("the raw tropical variety (01_raw_fan) has the lineality space of the Dressian", s3_same_lineality(raw, D),
          "dimension $(rank(s3_lineality(D))) in R^21")
    raw_f = Int.(f_vector(raw))[lineality_dim(raw)+1:end]
    check("the fan of Oscar's tropical_variety modulo lineality: f-vector", [51, 244, 251], raw_f)
    claim("it is simplicial", all(C -> length(rays_modulo_lineality(C).rays_modulo_lineality) == 3, maximal_cones(raw)))
    spans = s3_spans(raw)
    distinct = unique(spans)
    pairs = count(sp -> count(==(sp), spans) == 2, distinct)
    check("distinct linear spans of the 251 maximal cones, and pairs with the same span", [248, 3], [length(distinct), pairs])
    quotient = Matrix{QQFieldElem}(kernel(s3_lineality(D); side=:right))
    Gcells = Set(TropicalGeometryG2.cell_signatures(G, quotient))
    Dcells = TropicalGeometryG2.cell_signatures(D, quotient)          # in the order of maximal_cones(D)
    raw_cells = TropicalGeometryG2.cell_signatures(raw, quotient)     # in the order of spans
    quadrilaterals = Set(c for c in Gcells if length(c) == 4)
    merged = Set(union(raw_cells[findall(==(sp), spans)]...) for sp in distinct if count(==(sp), spans) == 2)
    claim("each pair of triangles with the same span merges into one of the three quadrilaterals of trop(G)",
          length(merged) == 3 && merged == quadrilaterals)
    coarse = rows(g2("grassmannian", "02_coarsest_cells.csv"))
    stored("the merged pairs of raw cells -> rays of the quadrilateral (stage coarsen, grassmannian/02_coarsest_cells.csv)",
           join(["$(r["original_cells"]) -> $(r["extremal_original_rays"])" for r in coarse if r["original_cell_count"] == "2"], ", "))
    check("distinct linear spans of the 248 maximal cones after merging (so the fan structure is coarsest)", 248, length(unique(s3_spans(G))))
    tetrahedra = [c for (c, C) in zip(Dcells, Dcones) if dim(C) == dim(D)]
    polygons = Set(c for (c, C) in zip(Dcells, Dcones) if dim(C) < dim(D))
    tetrahedron_facets = Set(Set(facet) for tet in tetrahedra for facet in Combinatorics.combinations(collect(tet), 3))
    claim("the three quadrilaterals are those of the G2-Dressian", quadrilaterals == Set(c for c in polygons if length(c) == 4))
    check("maximal cells of trop(G): maximal polygons of the Dressian, facets of its tetrahedra", [236, 12],
          [length(intersect(Gcells, polygons)), length(intersect(Gcells, tetrahedron_facets))])
    claim("trop(G) is the Dressian without the three open tetrahedra", Gcells == union(polygons, tetrahedron_facets))

    subheading("The atlas of trop(G): types of rays (fig:atlas_rays) and of maximal cells (fig:atlas_toric_cells_*)")
    check("types of rays and of maximal cells (orbits under D12)", [11, 31], [s["complex"]["orbit_f_vector"][1], s["complex"]["maximal_cell_orbits"]])
    orbit_rows = rows(g2("grassmannian_orbits.csv"))
    tree_of = Dict(parse(Int, r["index"]) => r for r in gtrees if r["kind"] == "ray")
    info("ray type   rays  census rays                         G2-tree of the first (pendant lengths; splits:length)")
    for r in orbit_rows
        r["kind"] == "ray" && !isempty(r["qq_members"]) || continue
        members = s3_members(r, "qq")
        info(rpad(r["orbit_type"], 10) * lpad(length(members), 5) * "  " * rpad(join(members, ";"), 36) *
             tree_of[first(members)]["pendant_lengths"] * "; " * tree_of[first(members)]["internal_splits"])
    end

    subheading("Remark on multiplicities (stored: studies/g2/weights/weights.jl --absolute, weights/weights.csv)")
    weights = rows(study_path("g2", "weights", "weights.csv"))
    double = [r for r in weights if r["mult_QQbar"] != "1"]
    L = locator(G)
    types = orbit_types("grassmannian_orbits.csv", "qq")
    double_cells = sort([only(locate(L, QQ.(parse.(Int, split(r["weight"], ';'))))) for r in double])
    stored("multiplicities of the 251 cones of the raw fan over QQbar (cones)",
           join(["$m ($(count(r -> r["mult_QQbar"] == m, weights)))" for m in sort(unique(r["mult_QQbar"] for r in weights))], ", "))
    # main.tex: "One of the $31$ maximal cones has multiplicity $2$. This is type \celltype{21}" (the 31 types).
    check("types of maximal cells of trop(G) of multiplicity 2 (\"one of the 31\", \"type 21\")", ["21"],
          unique(types[("cell", c)] for c in double_cells))
    cells21 = sort(only(s3_members(r, "qq") for r in orbit_rows if r["kind"] == "cell" && r["orbit_type"] == "21"))
    check("the maximal cells of multiplicity 2 are the whole orbit of type 21", cells21, double_cells)
    info("that orbit has $(length(double_cells)) cells ($(join(double_cells, ", "))): \"maximal cones\" in the remark means types")
    balancing = run_in_module(study_path("g2", "weights", "balancing.jl"); name=:Balancing)
    failures = redirect_stdout(() -> call(balancing, :check_balancing), devnull)
    check("walls of the raw fan where trop(G) with the QQbar multiplicities is not balanced (weights/balancing.jl)", 0,
          length(failures.QQbar))
    info("with all weights 1 (the multiplicities over QQ) balancing fails at $(length(failures.QQ)) walls")

    # The table of the 31 types of maximal cells.
    certificate = TOML.parsefile(g2("flag", "07_projection_certificate.toml"))
    above = Dict{Int,Int}()
    for c in certificate["cells"]
        above[c["target_cell"]] = get(above, c["target_cell"], 0) + 1
    end
    degenerations = Dict(r["orbit"] => r for r in rows(g2("grassmannian", "06_canonical_interior_degenerations.csv")))
    Gcell_rays = fan_cells(G)
    cell_rows = sort([r for r in orbit_rows if r["kind"] == "cell" && !isempty(r["qq_members"])]; by=r -> parse(Int, r["orbit_type"]))
    info("the 31 types of maximal cells of trop(G):")
    info("type  cells  rays  multiplicity  prime  toric  cells of trop(F) above each")
    for r in cell_rows
        members = s3_members(r, "qq")
        degeneration = degenerations[r["qq_orbit"]]
        info(lpad(r["orbit_type"], 4) * lpad(length(members), 7) * lpad(length(Gcell_rays[first(members)]), 6) *
             lpad(isempty(intersect(members, double_cells)) ? 1 : 2, 14) * lpad(degeneration["prime"], 7) *
             lpad(degeneration["toric"], 7) * lpad(get(above, first(members), 0), 8))
    end

    subheading("Theorem thm:flag: the tropical G2 flag variety trop(F)")
    s = summary("flag")
    orbits = face_orbits("flag")
    check("f-vector modulo lineality", [141, 1352, 3584, 2714], f_vector_of("flag"))
    check("3-cells by number of vertices: (tetrahedra, square pyramids, triangular prisms)", [2588, 108, 18],
          [face_count(orbits, 3, 4), face_count(orbits, 3, 5), face_count(orbits, 3, 6)])
    flag_types = maximal_types("flag")
    shapes = [(3, "3-simplex"), (3, "Pyramid(4-gon)"), (3, "Δ1 × Δ2")]
    check("3-cells by combinatorial type: 3-simplex, Pyramid(4-gon), Δ1 × Δ2", [2588, 108, 18], [first(flag_types[t]) for t in shapes])
    info("their orbits under D12: " * join([last(flag_types[t]) for t in shapes], ", "))
    check("reduced integral homology", "H_3 = Z^342", reduced_homology("flag"))
    check("f-vector modulo D12", [27, 173, 388, 270], s["complex"]["orbit_f_vector"])
    ftrees = rows(g2("flag", "04_trees.csv"))
    check("rays with a G2-tree, vertical rays (star tree, empty row), maximal cells with a G2-tree (flag/04_trees.csv)",
          [134, 7, 2714],
          [count(r -> r["kind"] == "ray" && r["vertical"] == "false" && !isempty(r["pendant_lengths"]), ftrees),
           count(r -> r["kind"] == "ray" && r["vertical"] == "true", ftrees),
           count(r -> r["kind"] == "cell" && !isempty(r["pendant_lengths"]), ftrees)])
    info("the decorations (the point 8) are in the column apple of flag_orbits.csv, checked in the proof below")

    subheading("Proof of thm:flag: Algorithm alg:balanced_support")
    J = load(g2("flag", "00_ideal.mrdi"))
    check("S: the generators of the ideal of F", 78, ngens(J))
    check("Gamma = D12: order of the symmetry group of trop(F)", 12, s["complex"]["group_order"])
    purity = TOML.parsefile(g2("flag", "00_ideal_purity.toml"))
    stored("the initial ideal in_{w0}(J) (FlagPurityVerifier.jl, data/g2/flag/00_ideal_purity.toml)",
           "binomial $(purity["initial_binomial"]), prime $(purity["initial_geometrically_prime"]), $(purity["initial_generators"]) generators, " *
           "lattice index $(purity["initial_lattice_index"]), dim in_{w0}(J) = $(purity["initial_dimension"]), d = dim J = $(purity["dimension"])")
    if QUICK
        skipped("in_{w0}(J) recomputed: 83 binomials, lattice index 1, saturated, d = 8")
    else
        J0 = initial(J, tropical_semiring_map(QQ, min), QQ.(-purity["weight_max"]))
        reduced = collect(gens(groebner_basis(J0; complete_reduction=true)))
        check("in_{w0}(J) recomputed: elements of its reduced Groebner basis", 83, length(reduced))
        claim("they are binomials, and their exponent differences span a saturated lattice (index 1)",
              all(g -> length(g) == 2, reduced) && TropicalGeometryG2.binomial_lattice_index(reduced) == 1)
        RJ = base_ring(J)
        claim("in_{w0}(J) is saturated by the product of the variables, hence, with the saturated lattice, prime over QQbar",
              saturation(J0, ideal(RJ, prod(gens(RJ)))) == J0)
        # The algorithm requires J prime of dimension d; main.tex no longer states d = 8 here.
        check("d = dim J (the algorithm requires J prime of dimension d), and dim in_{w0}(J)", [8, 8], [Int(dim(J)), Int(dim(J0))])
        stored_J0 = load(g2("flag", "00_verified_prime_initial_ideal.mrdi"))
        to_J = hom(base_ring(stored_J0), RJ, gens(RJ))
        claim("it equals data/g2/flag/00_verified_prime_initial_ideal.mrdi", ideal(RJ, to_J.(gens(stored_J0))) == J0)
    end
    fd = summary("flag_dressian")
    check("the prevariety of S (the flag Dressian): maximal cells, and their dimensions modulo lineality", [2627, [3, 4]],
          [fd["complex"]["maximal_cells"], sort(unique(t["dimension"] for t in fd["maximal_cell_types"]))])
    simplices = maximal_types("flag_dressian")[(4, "4-simplex")]
    check("first loop: 4-simplices removed, and their orbits", [33, 6], collect(simplices))
    check("3-cells of the flag Dressian, and their orbits", [2729, 274], [fd["complex"]["f_vector"][4], fd["complex"]["orbit_f_vector"][4]])
    check("second loop: orbits kept, and their cells", [270, 2714], [s["complex"]["maximal_cell_orbits"], s["complex"]["maximal_cells"]])
    run_log = TOML.parsefile(study_path("g2", "supplemental", "full28_benchmark.toml"))
    stored("the run of the algorithm (studies/g2/supplemental/full28_benchmark.toml)",
           "Gfan prevariety $(run_log["initial_prevariety_rays"]) rays, $(run_log["initial_prevariety_maximal_cones"]) cones; " *
           "first loop excludes $(run_log["oversized_orbits_excluded"]) orbits; $(run_log["candidate_skeleton_maximal_cones"]) 3-cells " *
           "in $(run_log["candidate_skeleton_orbits"]) orbits; second loop drops $(run_log["excluded_dimension_eight_orbits"]) orbits, " *
           "keeps $(run_log["selected_orbits"]) orbits, $(run_log["selected_maximal_cones"]) cells, $(run_log["selected_rays"]) rays")
    info(@sprintf("its timings: prevariety %.0f s, certification %.0f s, membership tests %.1f s in all (at most %.1f s each)",
                  run_log["initial_intersection_seconds"], run_log["certification_seconds"],
                  run_log["membership_seconds_cumulative"], run_log["membership_max_seconds"]))
    claim("the log agrees with the counts computed from the fans above",
          [run_log[k] for k in ("oversized_orbits_excluded", "candidate_skeleton_maximal_cones", "candidate_skeleton_orbits",
                                "excluded_dimension_eight_orbits", "selected_orbits", "selected_maximal_cones", "selected_rays")] ==
          [last(simplices), fd["complex"]["f_vector"][4], fd["complex"]["orbit_f_vector"][4],
           fd["complex"]["orbit_f_vector"][4] - s["complex"]["maximal_cell_orbits"], s["complex"]["maximal_cell_orbits"],
           s["complex"]["maximal_cells"], f_vector_of("flag")[1]])
    check("the output (01_raw_fan) is already coarsest: rays, maximal cones", [141, 2714],
          [s["stage1_fan"]["rays"], s["stage1_fan"]["maximal_cones"]])
    info("its f-vector and 3-cells are those of the theorem, checked above")

    subheading("Proof of thm:flag: the projection to trop(G), the decorations and the rays")
    claim("the exact projection certificate 07_projection_certificate.toml is valid (ray images, cells, D12-transport)",
          projection_certificate_valid("flag"))
    images = [c["target_cell"] for c in certificate["cells"]]
    check("each of the 2714 maximal cells maps onto a maximal cell of trop(G); the images (onto)", [2714, 248],
          [length(images), length(unique(images))])
    claim("07_projection_cells.csv lists the same images", [parse(Int, r["grassmannian_cell"]) for r in rows(g2("flag", "07_projection_cells.csv"))] == images)
    flag_orbit_rows = rows(g2("flag_orbits.csv"))
    flag_cell_rows = [r for r in flag_orbit_rows if r["kind"] == "cell" && !isempty(r["qq_orbit"])]
    decoration(r) = unique(first.(split.(split(r["apple"], '|'))))
    check("cell orbits of trop(F), and their cells", [270, 2714], [length(flag_cell_rows), sum(parse(Int, r["size"]) for r in flag_cell_rows)])
    claim("each is a maximal G2-tree decorated by a point on an edge (pendant or inner; column apple of flag_orbits.csv)",
          all(r -> all(in(("pendant", "inner")), decoration(r)), flag_cell_rows))
    info("cells with the point on a pendant edge: $(sum(parse(Int, r["size"]) for r in flag_cell_rows if decoration(r) == ["pendant"])), " *
         "on an inner edge: $(sum(parse(Int, r["size"]) for r in flag_cell_rows if decoration(r) == ["inner"]))")
    F, FD = fan("flag"), fan("flag_dressian")
    claim("trop(F) and the flag Dressian have the same lineality space", s3_same_lineality(F, FD))
    q = Matrix{QQFieldElem}(kernel(s3_lineality(FD); side=:right))
    FD_rays = [quotient_ray_signature(r, q) for r in rays_modulo_lineality(FD).rays_modulo_lineality]
    Fcells = Set(TropicalGeometryG2.cell_signatures(F, q))
    FDcones = collect(maximal_cones(FD))
    FDcells = TropicalGeometryG2.cell_signatures(FD, q)
    maximal3 = Set(c for (c, C) in zip(FDcells, FDcones) if dim(C) < dim(FD))
    simplex_facets = Set(Set(facet) for (c, C) in zip(FDcells, FDcones) if dim(C) == dim(FD) for facet in Combinatorics.combinations(collect(c), 4))
    removed = setdiff(union(maximal3, simplex_facets), Fcells)
    check("3-cells of the flag Dressian; those in trop(F); the others; tetrahedra among them", [2729, 2714, 15, 15],
          [length(union(maximal3, simplex_facets)), length(intersect(union(maximal3, simplex_facets), Fcells)), length(removed),
           count(c -> length(c) == 4, removed)])
    ray_index = Dict(signature => k for (k, signature) in enumerate(FD_rays))
    removed_rays = sort([sort([ray_index[r] for r in c]) for c in removed])
    info("of these, facets of the 4-simplices: $(length(intersect(removed, simplex_facets))); their rays (flag Dressian census): " *
         join(join.(removed_rays, ";"), ", "))
    FD_vectors = rays_modulo_lineality(FD).rays_modulo_lineality
    over = [Set(filter(!isnothing, [s3_plucker_image(FD_vectors[k], quotient) for k in c])) for c in removed_rays]
    check("the removed tetrahedra over each of the three tetrahedra of the Dressian (node-decorated tetrahedral G2-trees)",
          [5, 5, 5], [count(==(tet), over) for tet in tetrahedra])
    claim("the 141 rays of trop(F) are those of the flag Dressian, as vectors modulo lineality",
          Set(quotient_ray_signature(r, q) for r in rays_modulo_lineality(F).rays_modulo_lineality) == Set(FD_rays))
    flag_rays = sort([(r["orbit_type"], r["apple"], r["size"]) for r in flag_orbit_rows if r["kind"] == "ray" && !isempty(r["qq_orbit"])])
    flag_dressian_rays = sort([(r["orbit_type"], r["apple"], r["size"]) for r in rows(g2("flag_dressian", "flag_dressian_orbits.csv")) if r["kind"] == "ray"])
    claim("and the orbit tables agree (orbit types, apples and sizes)", flag_rays == flag_dressian_rays)
    targets = [r["target_ray"] for r in certificate["rays"]]
    check("rays of trop(F): over rays of trop(G), over its lineality space", [134, 7], [count(>(0), targets), count(==(0), targets)])
    check("rays of trop(G) that occur below a ray of trop(F)", 51, length(unique(filter(>(0), targets))))
    flag_type = Dict(m => r["orbit_type"] for r in flag_orbit_rows if r["kind"] == "ray" for m in s3_members(r, "qq"))
    claim("each ray over a ray of trop(G) has the G2-tree type of its image (a minimal G2-tree with a marked point 8)",
          all(flag_type[k] == types[("ray", t)] for (k, t) in enumerate(targets) if t > 0))
    FR = Matrix{QQFieldElem}(load(g2("flag", "03_normalized_rays.mrdi")))
    vertical = [k for k in 1:size(FR, 1) if all(iszero, FR[k, 8:28])]
    check("rays of trop(F) with zero Pluecker part", 7, length(vertical))
    claim("they are the rays over the lineality space, of type star, and the vertical rays of flag/04_trees.csv",
          vertical == findall(==(0), targets) && vertical == sort([k for (k, t) in flag_type if t == "star"]) &&
          vertical == sort([parse(Int, r["index"]) for r in ftrees if r["kind"] == "ray" && r["vertical"] == "true"]))
    unit(i) = [QQ(i == j ? 1 : 0) for j in 1:7]
    claim("their point parts are e_1, ..., e_7 modulo (1,...,1) (min convention): the point 8 on the pendant edge of each leaf",
          Set(FR[k, 1:7] .- minimum(FR[k, 1:7]) for k in vertical) == Set(unit.(1:7)))
    info("the 27 ray orbits of trop(F) and their images in trop(G):")
    for r in unique(r -> r["qq_orbit"], [r for r in flag_orbit_rows if r["kind"] == "ray" && !isempty(r["qq_orbit"])])
        members = s3_members(r, "qq")
        image = unique(t == 0 ? "the lineality space" : "rays of type " * types[("ray", t)] for t in targets[members])
        info(rpad(r["orbit_type"], 8) * "orbit of size " * rpad(length(members), 3) * rpad(join(members, ";"), 46) * "-> " * join(image, ", "))
    end

    subheading("Corollary cor:tropicalbasis: the three cubics")
    I = load(g2("grassmannian", "00_ideal.mrdi"))
    P = base_ring(I)
    names = string.(symbols(P))
    p(i, j) = P[findfirst(==("p$i$j"), names)]
    cubics = [2p(4, 5)p(4, 6)p(5, 6) + p(2, 3)p(3, 6)p(5, 6) - p(1, 4)p(1, 2)p(4, 5) + p(2, 6)p(1, 3)p(4, 5) - p(3, 5)p(1, 2)p(4, 6),
              2p(2, 3)p(1, 2)p(1, 3) - p(4, 6)p(1, 4)p(1, 3) + p(2, 5)p(5, 6)p(2, 3) + p(1, 6)p(4, 5)p(2, 3) - p(3, 4)p(5, 6)p(1, 2),
              2p(4, 6)p(5, 6)p(4, 5) + p(1, 2)p(2, 5)p(4, 5) - p(3, 6)p(1, 3)p(4, 6) + p(1, 5)p(2, 3)p(4, 6) - p(2, 4)p(1, 3)p(5, 6)]
    claim("the three cubics of eq:threecubics lie in I", all(c -> iszero(normal_form(c, I)), cubics))
    check("the candidates of the tropical basis: the 42 = 35 + 7 generators of I and the 3 cubics", 45, ngens(I) + length(cubics))

    subheading("Figure fig:tropbasis and the proof of cor:tropicalbasis: the G2-trees of the 5 + 5 + 5 monomials")
    R = Matrix{QQFieldElem}(load(g2("grassmannian", "03_normalized_rays.mrdi")))
    E = Matrix{QQFieldElem}(load(g2("grassmannian_gf2", "03_normalized_rays.mrdi")))
    gf2_types = orbit_types("grassmannian_orbits.csv", "gf2")
    monomials = [[[(4, 5), (4, 6), (5, 6)], [(2, 3), (3, 6), (5, 6)], [(1, 4), (1, 2), (4, 5)], [(2, 6), (1, 3), (4, 5)], [(3, 5), (1, 2), (4, 6)]],
                 [[(2, 3), (1, 2), (1, 3)], [(4, 6), (1, 4), (1, 3)], [(2, 5), (5, 6), (2, 3)], [(1, 6), (4, 5), (2, 3)], [(3, 4), (5, 6), (1, 2)]],
                 [[(4, 6), (5, 6), (4, 5)], [(1, 2), (2, 5), (4, 5)], [(3, 6), (1, 3), (4, 6)], [(1, 5), (2, 3), (4, 6)], [(2, 4), (1, 3), (5, 6)]]]
    vertices = [[34, 9, 40, 6], [33, 10, 41, 7], [8, 30, 12, 11]]   # QQ census rays = atlas rays 1:51 (to_tex/cubic_tree_table.sh)
    centres = [3, 39, 36]                                           # GF(2) census rays of type epsilon (atlas rays 52, 54, 53)
    tropical_value(v, m) = sum(v[PAIR_INDEX[pair]] for pair in m)
    lin = s3_lineality(G)
    G_rays = Set(quotient_ray_signature(r, quotient) for r in rays_modulo_lineality(G).rays_modulo_lineality)
    tree_types, at_vertices, at_centres, interior, rows_are_tetrahedra, witnesses = String[], true, true, true, true, true
    rng = MersenneTwister(3)
    for c in 1:3
        push!(tree_types, gf2_types[("ray", centres[c])], (types[("ray", r)] for r in vertices[c])...)
        rows_are_tetrahedra &= Set(quotient_ray_signature(R[r, :], quotient) for r in vertices[c]) in tetrahedra
        info("row $c: values of the 5 monomials of the cubic at the leftmost tree (GF(2) ray $(centres[c])) and the rays " *
             join(vertices[c], ", ") * ":")
        vals = [tropical_value(E[centres[c], :], m) for m in monomials[c]]
        at_centres &= findall(==(minimum(vals)), vals) == [1]
        info("    " * show_value(vals))
        for (k, r) in enumerate(vertices[c])
            vals = [tropical_value(R[r, :], m) for m in monomials[c]]
            at_vertices &= findall(>(minimum(vals)), vals) == [k + 1]
            info("    " * show_value(vals))
        end
        solvable, x = can_solve_with_solution(vcat(matrix(QQ, R[vertices[c], :]), lin), matrix(QQ, 1, 21, E[centres[c], :]); side=:left)
        interior &= solvable && all(k -> x[1, k] > 0, 1:4)
        solvable && info("    the leftmost tree = " * join(["$(x[1, k]) ray $(vertices[c][k])" for k in 1:4], " + ") * " modulo lineality")
        for _ in 1:100
            a = rand(rng, 1:20, 4)
            w = sum(a[k] .* R[vertices[c][k], :] for k in 1:4)
            vals = [tropical_value(w, m) for m in monomials[c]]
            witnesses &= findall(==(minimum(vals)), vals) == [1]
        end
    end
    check("orbit types of the 15 trees, row by row", repeat(["epsilon", "alpha_2", "alpha_4", "alpha_3", "alpha_3"], 3), tree_types)
    # main.tex, caption of fig:tropbasis, of all 15 trees: "These are rays of $\trop_2(\linespace)$".
    figure_trees = vcat([quotient_ray_signature(R[r, :], quotient) for v in vertices for r in v],
                        [quotient_ray_signature(E[r, :], quotient) for r in centres])
    G2_rays = Set(quotient_ray_signature(E[k, :], quotient) for k in axes(E, 1))
    check("trees of the figure that are rays of trop_2(G) (caption of fig:tropbasis: all 15)", 15, count(in(G2_rays), figure_trees))
    info("$(count(in(G_rays), figure_trees)) of them are rays of trop(G) as well: all but the leftmost, of type epsilon")
    claim("type epsilon has no ray over QQ", !("epsilon" in values(types)))
    qtrees = Dict(parse(Int, r["index"]) => r for r in gtrees if r["kind"] == "ray")
    etrees = Dict(parse(Int, r["index"]) => r for r in rows(g2("grassmannian_gf2", "04_trees.csv")) if r["kind"] == "ray")
    edge_lengths(r) = (t = tree_row(r); Set(filter(!iszero, vcat(t.pendant, collect(values(t.lengths))))))
    claim("all edges of the 15 trees have length 1 or 2",
          all(issubset(edge_lengths(t), Set([1, 2])) for t in vcat([qtrees[r] for v in vertices for r in v], [etrees[r] for r in centres])))
    claim("second column: splits", all(r -> qtrees[r]["pendant_lengths"] == "0;0;0;0;0;0;0" && !occursin(';', qtrees[r]["internal_splits"]), (34, 33, 8)))
    check("last tree of the first row (ray 6): d13, d23, d12, d35", [3, 2, 1, 0], [-R[6, PAIR_INDEX[pair]] for pair in ((1, 3), (2, 3), (1, 2), (3, 5))])
    claim("each row is one of the three tetrahedra of the Dressian: its last four trees are the vertices", rows_are_tetrahedra)
    claim("the leftmost tree lies in the relative interior of that tetrahedron", interior)
    claim("at each vertex tree d, in_d(f) is the sum of the four other monomials (min convention: only its own monomial is above the minimum)", at_vertices)
    claim("at the leftmost tree d, in_d(f) is the monomial above it (the unique minimum)", at_centres)
    claim("at 100 random positive combinations of the vertices of each row, in_d(f) is that monomial", witnesses)
    lineality_homogeneous(L) = (M = Matrix{QQFieldElem}(L);
        all(allequal([tropical_value(M[k, :], m) for m in ms]) for ms in monomials for k in 1:size(M, 1)))
    claim("the cubics are homogeneous for the lineality spaces of trop(G) and trop_2(G), so in_d(f) does not depend on the representative of d",
          lineality_homogeneous(lin) && lineality_homogeneous(s3_lineality(fan("grassmannian_gf2"))))
end
