# Section 2 of the paper: run by paper_results.jl, which defines the helpers used here.

"The star tree d_ij = s_i + s_j of the pendant lengths s, in the order of PAIRS."
s2_star_pluecker_vector(s) = [s[i] + s[j] for (i, j) in PAIRS]


"The rows of eq:lineality."
const s2_lineality_rows = [1 -1 0 -1 1 0 0; 0 1 -1 0 -1 1 0; -1 0 1 1 0 -1 0]

"Whether an inner node of the tree meets three inner edges (three pairwise disjoint sides of splits); else a caterpillar."
function s2_branching(r)
    sides = [Set(S) for S in keys(tree_splits(r))]
    sides = vcat(sides, [setdiff(Set(1:7), A) for A in sides])
    return any(isempty(A ∩ B) && isempty(A ∩ C) && isempty(B ∩ C) for (A, B, C) in Combinatorics.combinations(sides, 3))
end

"Whether the tree with pendant lengths `s` and internal splits `splits` (side => length) is normalized as in prop:normalized."
s2_normalized(s, splits) = minimum(s) == 0 && iszero(s2_lineality_rows * s) && gcd(vcat(s, Int[l for (_, l) in splits])) == 1

"Whether a row of a 04_trees.csv is normalized as in prop:normalized."
s2_normalized_tree(r) = s2_normalized(tree_row(r)...)


"The reduced row echelon basis of the span of `vectors`."
function s2_row_space(vectors)
    r, R = rref(matrix(QQ, reduce(vcat, permutedims.(vectors))))
    return Matrix{QQFieldElem}(R)[1:r, :]
end

"The maximal cones of a fan grouped by their linear span: a dictionary from the span to the cones."
function s2_span_classes(F)
    rs, lin = fan_rays(F)
    classes = Dict{Matrix{QQFieldElem},Vector{Int}}()
    for (k, cell) in enumerate(fan_cells(F))
        push!(get!(classes, s2_row_space(vcat(lin, rs[cell])), Int[]), k)
    end
    return classes
end

"A key of the ray of `v` modulo the space with reduced row echelon basis `R`: equal keys mean the same ray."
function s2_ray_key(v, R)
    w = copy(v)
    for k in 1:size(R, 1)
        p = findfirst(!iszero, R[k, :])
        w = w .- w[p] .* R[k, :]
    end
    return w .// abs(w[findfirst(!iszero, w)])
end

function section_2()
    D = fan("dressian")
    rays, lineality = fan_rays(D)
    LD = locator(D)
    dorbits = rows(g2("dressian", "dressian_orbits.csv"))

    subheading("eq:dressian: the tropicalizations of the seven linear forms")
    form_supports = Set(Set(PAIRS[k] for k in 1:21 if !iszero(f[k])) for f in phi_linear_forms())
    paper_supports = Set([Set([(1, 7), (5, 6)]), Set([(2, 7), (4, 6)]), Set([(3, 7), (4, 5)]), Set([(4, 7), (2, 3)]),
                          Set([(5, 7), (1, 3)]), Set([(6, 7), (1, 2)]), Set([(1, 4), (2, 5), (3, 6)])])
    claim("the supports of the 7 linear forms: d17=d56, d27=d46, d37=d45, d47=d23, d57=d13, d67=d12 and {d14, d25, d36}",
          form_supports == paper_supports)
    # The leaves of Figure fig:shape_hexagon: the colours {1,2,3} and {4,5,6}, the shapes {1,4}, {2,5}, {3,6}, the heart 7.
    taxon_colour(i) = i <= 3 ? 1 : 2
    taxon_shape(i) = mod1(i, 3)
    other_shapes(i) = Tuple(sort([j for j in 1:6 if taxon_colour(j) != taxon_colour(i) && taxon_shape(j) != taxon_shape(i)]))
    claim("with the colours {1,2,3}, {4,5,6} and the shapes {1,4}, {2,5}, {3,6} of Figure fig:shape_hexagon, (i) \"the " *
          "distance from any shape to the heart is the same as the distance between the two other shapes of the opposite " *
          "color\" and (ii) \"the maximum of the three distances between equal shapes is attained at least twice\" are eq:dressian",
          union(Set(Set([(i, 7), other_shapes(i)]) for i in 1:6),
                Set([Set((i, j) for i in 1:6 for j in i+1:6 if taxon_shape(i) == taxon_shape(j))])) == paper_supports)

    subheading("Figure fig:ex_G2_trees: the five trees")
    gtrees = Dict((r["kind"], parse(Int, r["index"])) => r for r in rows(g2("grassmannian", "04_trees.csv")))
    qq_ray_count = count(key -> first(key) == "ray", keys(gtrees))
    gf2_trees = Dict(parse(Int, r["index"]) => r for r in rows(g2("grassmannian_gf2", "04_trees.csv")) if r["kind"] == "ray")
    epsilon_rays = sort([k for ((kind, k), t) in orbit_types("grassmannian_orbits.csv", "gf2") if kind == "ray" && t == "epsilon"])
    gf2_only_types = [r["orbit_type"] for r in rows(g2("grassmannian_orbits.csv")) if r["kind"] == "ray" && isempty(r["qq_orbit"])]
    claim("the fourth tree exists only over GF(2): the rays of trop_2(G) that are not rays of trop(G) are the 3 rays of " *
          "type epsilon (grassmannian_orbits.csv), and their trees are snowflakes (a node meets three inner edges)",
          gf2_only_types == ["epsilon"] && length(epsilon_rays) == length(gf2_trees) - qq_ray_count &&
          all(k -> s2_branching(gf2_trees[k]), epsilon_rays), "GF(2) census rays $(join(epsilon_rays, ", "))")

    subheading("Example ex:fivetrees: the caterpillar G2-tree of that triangle")
    metrics = Dict(parse(Int, r["number"]) => r for r in rows(g2("maximal_cell_tree_metrics.csv")))
    paper = Dict("17" => 7, "56" => 7, "27" => 7, "46" => 7, "37" => 6, "45" => 6, "47" => 7, "23" => 7, "57" => 7,
                 "13" => 7, "67" => 6, "12" => 6, "36" => 6, "14" => 8, "25" => 8, "15" => 8, "16" => 1, "24" => 8,
                 "26" => 5, "34" => 1, "35" => 5)
    paper_d = [paper["$i$j"] for (i, j) in PAIRS]
    row = metrics[1]
    triangle = parse.(Int, split(row["support"], ';'))
    cell_number = findfirst(==(triangle), fan_cells(fan("grassmannian")))
    check("its metric d (d_ij in lexicographic order)", paper_d, [parse(Int, row["d$i$j"]) for (i, j) in PAIRS])
    d = QQ.(paper_d)
    dist(i, j) = d[PAIR_INDEX[minmax(i, j)]]
    twice(v) = count(==(maximum(v)), v) >= 2
    equalities = [((1, 7), (5, 6)), ((2, 7), (4, 6)), ((3, 7), (4, 5)), ((4, 7), (2, 3)), ((5, 7), (1, 3)), ((6, 7), (1, 2))]
    claim("d is a G2-tree: the six equations and max(d14, d25, d36) attained twice (eq:dressian)",
          all(dist(a...) == dist(b...) for (a, b) in equalities) && twice([dist(1, 4), dist(2, 5), dist(3, 6)]))
    claim("d satisfies the four-point condition eq:4point for all 35 quadruples",
          all(twice([dist(i, j) + dist(k, l), dist(i, k) + dist(j, l), dist(i, l) + dist(j, k)])
              for (i, j, k, l) in Combinatorics.combinations(1:7, 4)))
    hits = locate(LD, [-x for x in d])   # the stored fans use Oscar's min convention, so they contain -d
    dcell = Dict(parse(Int, m) => r["orbit_type"] for r in dorbits for m in split(r["members"], ';') if r["kind"] == "cell")
    claim("-d lies in the relative interior of exactly one maximal cell of D", length(hits) == 1,
          length(hits) == 1 ? "cell $(only(hits)), $(length(fan_cells(D)[only(hits)])) rays, type $(dcell[only(hits)])" :
                              "$(length(hits)) cells")
    if length(hits) == 1
        tree = only(r for r in rows(g2("dressian", "04_trees.csv")) if r["kind"] == "cell" && parse(Int, r["index"]) == only(hits))
        check("the metric of the normalized tree of that cell (dressian/04_trees.csv)", paper_d, tree_metric(tree))
        caterpillar = gtrees[("cell", cell_number)]
        claim("that tree is the tree of cell $cell_number of trop(G), with the rays of row 1 (same pendant lengths and splits), " *
              "and a caterpillar (no inner node meets three inner edges)",
              tree["pendant_lengths"] == caterpillar["pendant_lengths"] && tree_splits(tree) == tree_splits(caterpillar) &&
              !s2_branching(caterpillar))
    end

    subheading("The lineality space of trop(Gr(2,7)): the star trees")
    star_space = [QQ.(s2_star_pluecker_vector([Int(m == i) for m in 1:7])) for i in 1:7]
    ones21 = fill(QQ(1), 21)
    check("star trees d_ij = s_i + s_j: dimension in R^21/1 (they contain 1)", 6, rank_of(vcat(star_space, [ones21])) - 1)

    subheading("Lemma lem:lineality-space: the lineality space of the G2-Dressian")
    check("dimension of the lineality space in R^21/1", 2, lineality_dim(D) - 1)
    H = [[QQ(Int(k == PAIR_INDEX[a]) - Int(k == PAIR_INDEX[b])) for k in 1:21] for (a, b) in equalities]
    S = matrix(QQ, reduce(hcat, star_space))                 # 21 x 7: s -> the star tree of s
    _, N = nullspace(matrix(QQ, reduce(vcat, permutedims.(H))) * S)
    T = S * N
    meet = [[T[r, c] for r in 1:21] for c in 1:ncols(T)]
    check("star trees satisfying the six equations of eq:dressian: dimension in R^21/1", 2, rank_of(meet) - 1)
    claim("they form exactly the lineality space of D", rank_of(vcat(lineality, meet)) == rank_of(lineality) == rank_of(meet))
    claim("on them also d14 = d25 = d36, the lineality of max(d14, d25, d36): so they are the star trees in the lineality of trop(L)",
          all(v -> v[PAIR_INDEX[(1, 4)]] == v[PAIR_INDEX[(2, 5)]] == v[PAIR_INDEX[(3, 6)]], meet))
    stars = [s2_lineality_rows[k, :] for k in 1:3]
    vectors = [QQ.(s2_star_pluecker_vector(s)) for s in stars]
    claim("the three star trees of eq:lineality lie in it, and any two of them with 1 span it",
          all(v -> rank_of(vcat(lineality, [v])) == rank_of(lineality), vectors) &&
          all(rank_of([vectors[a], vectors[b], ones21]) == lineality_dim(D) for a in 1:3 for b in a+1:3))
    check("dimension of the linear span of the Dressian (six independent equations eq:dressian)", 15,
          rank_of(vcat(rays, lineality)))

    subheading("Figure fig:root_lattice: the weights of the 28 coordinates (columns of eq:lineality)")
    wx(i) = [stars[1][i], stars[2][i], stars[3][i]]
    variable_weights = vcat([("x$i", wx(i)) for i in 1:7], [("p$i$j", wx(i) + wx(j)) for (i, j) in PAIRS])
    root_weight(name) = only(w for (n, w) in variable_weights if n == name)
    classes = Dict{Vector{Int},Vector{String}}()
    for (name, w) in variable_weights
        push!(get!(classes, w, String[]), name)
    end
    for (w, names) in sort(collect(classes); by=p -> (sum(abs2, p[1]), p[1]))
        info(rpad(string(w), 12) * "squared length $(sum(abs2, w)): " * join(names, ", "))
    end
    check("the variables of weight 0", ["x7", "p14", "p25", "p36"], classes[[0, 0, 0]])
    check("alpha = (1,0,-1) is the weight of x1, p17, p56; beta = (-1,-1,2) is the weight of", (["x1", "p17", "p56"], ["p34"]),
          (classes[[1, 0, -1]], classes[[-1, -1, 2]]))
    check("nonzero weights: short roots (squared length 2), long roots (squared length 6)", (6, 6),
          (count(w -> sum(abs2, w) == 2, keys(classes)), count(w -> sum(abs2, w) == 6, keys(classes))))
    claim("the leaves x1..x6 carry the six short roots, the heart x7 the origin; the long roots are the hexagon edges",
          Set(wx(i) for i in 1:6) == Set(w for w in keys(classes) if sum(abs2, w) == 2) && iszero(wx(7)) &&
          Set(only(classes[w]) for w in keys(classes) if sum(abs2, w) == 6) == Set(["p15", "p16", "p24", "p26", "p34", "p35"]))
    figure = [["x1", "p17", "p56"], ["x5", "p13", "p57"], ["x3", "p37", "p45"], ["x4", "p23", "p47"], ["x2", "p27", "p46"],
              ["x6", "p12", "p67"], ["p15"], ["p35"], ["p34"], ["p24"], ["p26"], ["p16"], ["x7", "p14", "p25", "p36"]]
    check("the labels of figures/root_lattice_V2.tex (drawn by hand) are the classes of equal weight",
          sort(sort.(figure)), sort(sort.(collect(values(classes)))))
    monomial_weight(f, e) = sum(e[k] * root_weight(string(symbols(parent(f))[k])) for k in eachindex(e))
    homogeneous(f) = length(unique(monomial_weight(f, e) for e in AbstractAlgebra.exponent_vectors(f))) == 1
    claim("the torus acts with these weights: the 42 generators of I and the 78 of J are homogeneous",
          all(homogeneous, gens(load(g2("grassmannian", "00_ideal.mrdi")))) &&
          all(homogeneous, gens(load(g2("flag", "00_ideal.mrdi")))))

    subheading("Proposition prop:normalized: the rows of eq:lineality")
    claim("\"the plane spanned by the three rows\" (rank 2), and \"Their entries sum to zero\"",
          rank_of(stars) == 2 && all(iszero, sum(s2_lineality_rows; dims=2)))

    subheading("The Weyl group D12 acting on the coordinates")
    rotation, reflection = [5, 6, 4, 2, 3, 1, 7], [1, 3, 2, 4, 6, 5, 7]
    generated = Set([collect(1:7)])
    frontier = [collect(1:7)]
    while !isempty(frontier)
        frontier = unique([s[g] for g in frontier for s in (rotation, reflection) if !(s[g] in generated)])
        union!(generated, frontier)
    end
    hexagon = Set([Set([1, 5]), Set([5, 3]), Set([3, 4]), Set([4, 2]), Set([2, 6]), Set([6, 1])])
    claim("D12 = <(1 5 3 4 2 6), (2 3)(5 6)> preserves the hexagon 1-5-3-4-2-6 of Figure fig:root_lattice, " *
          "so it is its full symmetry group (order 12)",
          length(generated) == 12 && all(g -> Set(Set(g[collect(e)]) for e in hexagon) == hexagon, generated))
    for target in ("dressian", "flag_dressian")
        census_group = Set(load_symmetry_artifacts(joinpath(TARGET_DIRECTORIES[target], "05_symmetry")).permutations)
        claim("the symmetry group of the $target census is this D12", census_group == generated,
              "order $(length(census_group))")
    end
    check("order of the symmetry group of the census of D", 12, summary("dressian")["complex"]["group_order"])
    check("order of the symmetry group of the census of the flag Dressian", 12,
          summary("flag_dressian")["complex"]["group_order"])
    witness = readlines(study_path("g2", "char-p", "results", "witness.txt"))
    stored("the signed action (char-p/witness.jl)", only(l for l in witness if startswith(l, "all 12 taxon permutations")))

    subheading("Theorem thm:dressian: the G2-Dressian")
    s = summary("dressian")
    check("f-vector modulo lineality", [51, 241, 248, 3], f_vector_of("dressian"))
    check("2-dimensional cells: (triangles, quadrilaterals)", [245, 3],
          [s["two_dimensional_faces"]["3-gon"], s["two_dimensional_faces"]["4-gon"]])
    types = maximal_types("dressian")
    check("the three 3-cells are tetrahedra", 3, first(types[(3, "3-simplex")]))
    check("f-vector modulo D12", [11, 33, 31, 1], s["complex"]["orbit_f_vector"])
    check("reduced integral homology", "H_2 = Z^54", reduced_homology("dressian"))
    prevariety = load(g2("dressian", "01_prevariety_fan.mrdi"))
    check("the prevariety of the 42 generators: f-vector of the fan D modulo lineality", [51, 241, 248, 3],
          Int.(f_vector(prevariety))[lineality_dim(prevariety)+1:end])
    stored("generators of the prevariety (data/g2/dressian/provenance/tropicalize.toml)",
           TOML.parsefile(g2("dressian", "provenance", "tropicalize.toml"))["metadata"]["defining_generators"])
    claim("the rays and the lineality space of D satisfy the six equations of eq:dressian, which are linearly independent",
          rank_of(H) == 6 && all(v -> all(h -> iszero(sum(h .* v)), H), vcat(rays, lineality)))
    check("maximal cells: (triangles, quadrilaterals, tetrahedra)", [233, 3, 3],
          [first(types[(2, "3-gon")]), first(types[(2, "4-gon")]), first(types[(3, "3-simplex")])])
    Dcells = fan_cells(D)
    tetrahedra = [k for (k, c) in enumerate(Dcells) if rank_of(vcat(lineality, rays[c])) - rank_of(lineality) == 4]
    tetrahedron_facets = length(unique(setdiff(Dcells[k], [r]) for k in tetrahedra for r in Dcells[k]))
    check("polygons that are facets of the tetrahedra (distinct), and polygons that are maximal cells", [12, 236],
          [tetrahedron_facets, first(types[(2, "3-gon")]) + first(types[(2, "4-gon")])])
    claim("together they are the 248 polygons",
          tetrahedron_facets + first(types[(2, "3-gon")]) + first(types[(2, "4-gon")]) == f_vector_of("dressian")[3])
    check("maximal cells, and their orbits under D12 (28 of polygons, 1 of tetrahedra)", [239, 29],
          [s["complex"]["maximal_cells"], s["complex"]["maximal_cell_orbits"]])
    check("distinct linear spans of the 239 maximal cones (so the fan structure is the coarsest)", 239,
          length(s2_span_classes(D)))
    maximal_polygon_orbits = last(types[(2, "3-gon")]) + last(types[(2, "4-gon")])
    check("orbits of vertices, and orbits of maximal polygons", [11, 28],
          [s["complex"]["orbit_f_vector"][1], maximal_polygon_orbits])
    polygon_orbits = s["complex"]["orbit_f_vector"][3]
    info("\"The maximal triangles and quadrilaterals form 28 orbits\"; all 248 polygons form " *
         "$polygon_orbits orbits, the 12 facets of the tetrahedra $(polygon_orbits - maximal_polygon_orbits)")
    # figures/ray_orbits_styled_catalog.tex, without epsilon (3)
    figure_rays = [("alpha_1", 6), ("alpha_2", 3), ("alpha_3", 6), ("alpha_4", 3), ("alpha_5", 2), ("beta_1", 12), ("beta_2", 6),
                   ("gamma_1", 6), ("gamma_2", 3), ("gamma_3", 3), ("delta", 1)]
    check("the 11 ray orbits of D: the types of Figure fig:atlas_rays without epsilon, with their orbit sizes", sort(figure_rays),
          sort([(r["orbit_type"], parse(Int, r["size"])) for r in dorbits if r["kind"] == "ray"]))
    by_number(t) = parse(Int, t)
    polygons = [r for r in dorbits if r["kind"] == "cell" && r["orbit_type"] != "epsilon"]
    figure_cells = string.(setdiff(1:35, [18, 19, 20, 35, 32, 33, 34]))   # the panels "char 2" (18-20, 35), "char ≠ 2" (32-34)
    check("the 28 polygon orbits of D: the cell types of Figures fig:atlas_toric_cells_prime/nontoric without " *
          "'char 2' and 'char != 2'", figure_cells, sort([r["orbit_type"] for r in polygons]; by=by_number))
    grassmannian_cells = [r for r in rows(g2("grassmannian_orbits.csv")) if r["kind"] == "cell"]
    claim("these are the cell types of trop(G) over both QQ and GF(2) (grassmannian_orbits.csv)",
          figure_cells == sort([r["orbit_type"] for r in grassmannian_cells if !isempty(r["qq_orbit"]) && !isempty(r["gf2_orbit"])];
                               by=by_number))
    grassmannian_size = Dict(r["orbit_type"] => r["size"] for r in grassmannian_cells)
    claim("with the orbit sizes of grassmannian_orbits.csv, which the two figures print",
          all(r -> grassmannian_size[r["orbit_type"]] == r["size"], polygons))
    epsilon = only(r for r in dorbits if r["kind"] == "cell" && r["orbit_type"] == "epsilon")
    check("the orbit of the 3 tetrahedra is labelled epsilon (dressian_orbits.csv)", tetrahedra,
          sort(parse.(Int, split(epsilon["members"], ';'))))
    gf2_rays, _ = fan_rays(fan("grassmannian_gf2"))
    epsilon_cells = [locate(LD, gf2_rays[k]) for k in epsilon_rays]
    check("the maximal cells of D whose relative interiors contain the 3 epsilon rays of trop_2(G)", tetrahedra,
          sort(reduce(vcat, epsilon_cells)))
    dressian_trees = Dict(parse(Int, r["index"]) => r for r in rows(g2("dressian", "04_trees.csv")) if r["kind"] == "cell")
    topology(r) = Set(keys(tree_splits(r)))
    claim("the tetrahedral trees have the splits of the epsilon trees (the same G2-tree up to edge lengths)",
          Set(topology(dressian_trees[k]) for k in tetrahedra) == Set(topology(gf2_trees[k]) for k in epsilon_rays))
    for (k, cells) in zip(epsilon_rays, epsilon_cells), c in cells
        info("tetrahedron $c: pendant $(dressian_trees[c]["pendant_lengths"]), splits $(dressian_trees[c]["internal_splits"]); " *
             "epsilon ray $k over GF(2): pendant $(gf2_trees[k]["pendant_lengths"]), splits $(gf2_trees[k]["internal_splits"])")
    end
    for target in ("dressian", "grassmannian", "grassmannian_gf2")
        trees = rows(joinpath(TARGET_DIRECTORIES[target], "04_trees.csv"))
        claim("the $(length(trees)) trees of the rays and maximal cells of the $target census are normalized (prop:normalized)",
              all(s2_normalized_tree, trees))
    end

    subheading("The G2 flag Dressian")
    FD = fan("flag_dressian")
    fp = TOML.parsefile(g2("flag_dressian", "provenance", "tropicalize.toml"))["metadata"]
    stored("the Gfan prevariety of the flag Dressian (flag_dressian/provenance/tropicalize.toml)",
           "$(fp["original_generators"]) generators, reduced to $(fp["reduced_generators"]); reduction verified: " *
           "$(fp["generator_reduction_verified"]), every cone checked against the original generators: " *
           "$(fp["all_cones_verified_against_original_generators"])")
    check("dimension of its lineality space modulo the two all-ones vectors", 2, lineality_dim(FD) - 2)
    FDrays, FDlin = fan_rays(FD)
    concatenations = [QQ.(vcat(t, s2_star_pluecker_vector(t))) for t in stars]
    all_ones = [vcat(fill(QQ(1), 7), fill(QQ(0), 21)), vcat(fill(QQ(0), 7), fill(QQ(1), 21))]
    claim("its lineality space is spanned by the concatenations (s, star tree of s) of the rows of eq:lineality and " *
          "the two all-ones vectors",
          rank_of(vcat(concatenations, all_ones)) == rank_of(FDlin) == rank_of(vcat(FDlin, concatenations, all_ones)))
    # main.tex: "The former with the vector (1,1,1,-1,-1,-1,0) span the lineality of the 3-dimensional tropical hypersurface
    # trop(Q) and the latter span the lineality of the G2-Dressian."
    _, xq = polynomial_ring(QQ, ["x$i" for i in 1:7])
    TQ = tropical_hypersurface(xq[1] * xq[4] + xq[2] * xq[5] + xq[3] * xq[6] - xq[7]^2, tropical_semiring_map(QQ, min))
    Qlin = [QQ.(collect(v)) for v in lineality_space(TQ)]
    ones7 = fill(QQ(1), 7)
    check("trop(Q) ⊂ R^7/1, the tropical hypersurface of q: the dimension of its lineality space (w1 + w4 = w2 + w5 = w3 + w6 = 2 w7)",
          3, lineality_dim(TQ) - 1)
    info("trop(Q) itself has dimension $(dim(TQ) - 1) in R^7/1: the \"3-dimensional\" of the paper is its lineality space")
    extra = QQ.([1, 1, 1, -1, -1, -1, 0])
    claim("the rows of eq:lineality with (1,1,1,-1,-1,-1,0) span the lineality of trop(Q) (modulo 1)",
          rank_of(Qlin) == rank_of(vcat(Qlin, [QQ.(t) for t in stars], [extra])) ==
          rank_of(vcat([QQ.(t) for t in stars], [extra], [ones7])))
    claim("the star trees of the rows of eq:lineality span the lineality of the G2-Dressian (modulo 1)",
          rank_of(lineality) == rank_of(vcat(lineality, vectors)) == rank_of(vcat(vectors, [ones21])))

    subheading("Theorem thm:flag_dressian: the G2 flag Dressian")
    s = summary("flag_dressian")
    orbits = face_orbits("flag_dressian")
    flag_types = maximal_types("flag_dressian")
    check("f-vector modulo lineality", [141, 1352, 3584, 2729, 33], f_vector_of("flag_dressian"))
    pyramids, prisms = first(flag_types[(3, "Pyramid(4-gon)")]), first(flag_types[(3, "Δ1 × Δ2")])
    check("3-cells: (tetrahedra, square pyramids, triangular prisms)", [2603, 108, 18],
          [face_count(orbits, 3, 4), pyramids, prisms])
    claim("these are all the 3-cells: the square pyramids and prisms (all maximal) are the 3-cells with 5 and 6 vertices",
          face_count(orbits, 3, 5) == pyramids && face_count(orbits, 3, 6) == prisms &&
          face_count(orbits, 3, 4) + pyramids + prisms == f_vector_of("flag_dressian")[4])
    check("4-cells: 4-simplices", 33, first(flag_types[(4, "4-simplex")]))
    check("reduced integral homology", "H_3 = Z^324", reduced_homology("flag_dressian"))
    check("f-vector modulo D12", [27, 173, 388, 274, 6], s["complex"]["orbit_f_vector"])

    subheading("Proof of thm:flag_dressian: from the Gfan prevariety to the coarsest fan")
    prevariety = load(g2("flag_dressian", "01_prevariety_fan.mrdi"))
    check("maximal cones of the prevariety of the 78 generators", 2726, n_maximal_cones(prevariety))
    pre_classes = s2_span_classes(prevariety)
    fd_classes = s2_span_classes(FD)
    class_sizes = [length(g) for g in values(pre_classes)]
    check("maximal cones of the prevariety by linear span: (alone, in pairs, in larger groups)", [2528, 99, 0],
          [count(==(1), class_sizes), count(==(2), class_sizes), count(>(2), class_sizes)])
    Prays, Plin = fan_rays(prevariety)
    Pcells = fan_cells(prevariety)
    fd_cones = collect(maximal_cones(FD))
    as_matrix(vs) = matrix(QQ, reduce(vcat, permutedims.(vs)))
    merged_pairs = [g for g in values(pre_classes) if length(g) == 2]
    claim("the positive hull of each pair is the maximal cone of the coarsest fan with the same span (a convex union)",
          all(haskey(fd_classes, sp) &&
              positive_hull(as_matrix(Prays[union(Pcells[g[1]], Pcells[g[2]])]), as_matrix(Plin)) == fd_cones[only(fd_classes[sp])]
              for (sp, g) in pre_classes if length(g) == 2))
    pair_dimensions = [rank_of(vcat(Plin, Prays[Pcells[first(g)]])) - rank_of(Plin) - 1 for g in merged_pairs]
    info("the merged pairs: $(count(==(3), pair_dimensions)) of 3-cells and $(count(==(4), pair_dimensions)) of 4-cells")
    R = s2_row_space(FDlin)
    fd_keys = Set(s2_ray_key(r, R) for r in FDrays)
    removed = [k for k in eachindex(Prays) if !(s2_ray_key(Prays[k], R) in fd_keys)]
    paired = Set(vcat(merged_pairs...))
    containing(k) = [c for c in eachindex(Pcells) if k in Pcells[c]]
    check("rays removed by merging; maximal cones left", [6, 2627], [length(removed), n_maximal_cones(FD)])
    claim("each removed ray lies only in cones of the merged pairs", all(k -> all(in(paired), containing(k)), removed),
          "in $(join(length.(containing.(removed)), ", ")) cones")
    check("pairwise distinct linear spans of the maximal cones of the coarsest fan", 2627, length(fd_classes))
    stored("convexity and support (flag_dressian/02_coarsest_report.md)",
           only(l for l in readlines(g2("flag_dressian", "02_coarsest_report.md")) if startswith(l, "Maximal cones are grouped")))

    subheading("Proof of thm:flag_dressian: the rays are decorated G2-trees")
    base = rank_of(lineality)
    ray_image = map(FDrays) do r                # coordinates x1..x7, then the p_ij in lexicographic order
        p = r[8:28]
        rank_of(vcat(lineality, [p])) == base && return 0
        return something(findfirst(q -> rank_of(vcat(lineality, [p, q])) == base + 1, rays), -1)
    end
    check("rays of the flag Dressian over rays of D, and over its lineality space (star trees)", [134, 7],
          [count(>(0), ray_image), count(==(0), ray_image)])
    ray_trees = Dict(parse(Int, r["index"]) => r for r in rows(g2("dressian", "04_trees.csv")) if r["kind"] == "ray")
    nodes = [length(tree_splits(ray_trees[j])) + 1 for j in eachindex(rays)]
    check("rays of D by inner nodes of their trees (2, 3, 4), from dressian/04_trees.csv: 2*20 + 3*30 + 4*1", [20, 30, 1],
          [count(==(k), nodes) for k in 2:4])
    claim("over each ray of D lie as many rays of the flag Dressian as its tree has inner nodes",
          all(j -> count(==(j), ray_image) == nodes[j], eachindex(rays)))
    for r in dorbits
        r["kind"] == "ray" || continue
        members = parse.(Int, split(r["members"], ';'))
        info(rpad(r["orbit_type"], 8) * lpad(r["size"], 2) * (r["size"] == "1" ? " ray" : " rays") * " of D, with " *
             join(unique(nodes[members]), "/") * " inner nodes and " *
             join(unique(count(==(j), ray_image) for j in members), "/") * " rays of the flag Dressian above each")
    end
    flag_dressian_rays = [r for r in rows(g2("flag_dressian", "flag_dressian_orbits.csv")) if r["kind"] == "ray"]
    ray_orbit_types = unique(r["orbit_type"] for r in flag_dressian_rays)
    info("the $(length(flag_dressian_rays)) ray orbits of the flag Dressian, by type of the tree (flag_dressian_orbits.csv): " *
         join(["$t " * join([r["size"] for r in flag_dressian_rays if r["orbit_type"] == t], "+") for t in ray_orbit_types], ", "))
    ones7 = fill(QQ(1), 7)
    units = [[QQ(Int(m == i)) for m in 1:7] for i in 1:7]
    check("the 7 star rays: the leaf whose pendant edge carries point 8 (x-part e_i modulo 1)", [[i] for i in 1:7],
          sort([findall(i -> rank_of([r[1:7], units[i], ones7]) == 2, 1:7) for (r, j) in zip(FDrays, ray_image) if j == 0]))

    subheading("Proof of thm:flag_dressian: decorated trees, the faces over the cells of D (by orbits)")
    cone_dim_D(o) = rank_of(vcat(lineality, rays[o])) - base - 1
    over(face) = sort!(unique(filter(>(0), ray_image[face])))
    decorations = Dict{Tuple{Int,Int},Int}()     # (dimension of the face, dimension of its image in D) => faces
    for (face, (k, size)) in orbits
        k >= 3 || continue
        key = (k, cone_dim_D(over(face)))
        decorations[key] = get(decorations, key, 0) + size
    end
    check("4-cells over tetrahedra; 3-cells over 2-cells (apple on an edge) and over tetrahedra (apple on a node)",
          [33, 2714, 15], [get(decorations, (4, 3), 0), get(decorations, (3, 2), 0), get(decorations, (3, 3), 0)])
    inner_edges(k) = length(tree_splits(dressian_trees[k]))
    check("the 3 tetrahedral trees: (edges, inner nodes); hence 33 = 3*11 and 15 = 3*5", [(11, 5)],
          unique([(7 + inner_edges(k), inner_edges(k) + 1) for k in tetrahedra]))
    Dsym = load_symmetry_artifacts(g2("dressian", "05_symmetry"))
    Dorbit = Dict(c => o for (o, orbit) in enumerate(Dsym.cell_orbits) for c in orbit)
    Dindex = Dict(c => k for (k, c) in enumerate(Dcells))
    above, on_facets = Dict{Int,Int}(), 0
    for (face, (k, size)) in orbits
        (k == 3 && cone_dim_D(over(face)) == 2) || continue
        o = over(face)
        if haskey(Dindex, o)
            above[Dorbit[Dindex[o]]] = get(above, Dorbit[Dindex[o]], 0) + size
        else
            on_facets += size
        end
    end
    claim("over each of the 28 orbits of maximal polygons of D lie as many 3-cells as its G2-tree has edges (10 or 11)",
          length(above) == 28 &&
          all(n == length(Dsym.cell_orbits[o]) * (7 + inner_edges(first(Dsym.cell_orbits[o]))) for (o, n) in above))
    check("3-cells over the 12 facets of the tetrahedra (12 x 10 edges)", 12 * 10, on_facets)
    labelled = Dict{String,Int}()
    for r in rows(g2("flag_dressian", "flag_dressian_orbits.csv"))
        r["kind"] == "cell" && (labelled[r["orbit_type"]] = get(labelled, r["orbit_type"], 0) + parse(Int, r["size"]))
    end
    table = String[]
    for r in sort([r for r in dorbits if r["kind"] == "cell"]; by=r -> something(tryparse(Int, r["orbit_type"]), 99))
        representative = first(parse.(Int, split(r["members"], ';')))
        shape = representative in tetrahedra ? "tetrahedra" : length(Dcells[representative]) == 3 ? "triangles" : "quadrilaterals"
        cells_above = representative in tetrahedra ? get(decorations, (4, 3), 0) : get(above, Dorbit[representative], 0)
        push!(table, "$(r["orbit_type"]): $(r["size"]) $shape x $(7 + inner_edges(representative)) edges = $cells_above")
    end
    info("maximal cells of the flag Dressian over each orbit of maximal cells of D (type: orbit, edges of the tree):")
    foreach(k -> info("  " * join(table[k:min(k + 3, end)], ";  ")), 1:4:length(table))
    claim("flag_dressian_orbits.csv labels the same number of maximal cells with each type, $(sum(values(labelled))) in all",
          sum(values(labelled)) == n_maximal_cones(FD) &&
          all(r -> r["kind"] != "cell" || get(labelled, r["orbit_type"], 0) ==
                   parse(Int, r["size"]) * (7 + inner_edges(first(parse.(Int, split(r["members"], ';'))))), dorbits))
end
