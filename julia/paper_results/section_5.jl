# Section 5 of the paper: run by julia/paper_results.jl, which defines the helpers used here.

"The internal splits of the tree of a ray (the side of each split without taxon 7)."
function s5_ray_splits(ray)
    d = Dict{Tuple{Int,Int},QQFieldElem}()
    for (k, (i, j)) in enumerate(PAIRS)
        d[(i, j)] = -ray[k]     # the tree metric is minus the min-convention point
    end
    # Four-point condition: {i,j}|{k,l} is a split of the quartet when d_ij + d_kl < the other two sums.
    dist(i, j) = i == j ? QQ(0) : d[minmax(i, j)]
    splits = Set{Vector{Int}}()
    for S in (s for n in 2:5 for s in Combinatorics.combinations(1:6, n))
        complement = setdiff(1:7, S)
        ok = all(dist(i, j) + dist(k, l) < min(dist(i, k) + dist(j, l), dist(i, l) + dist(j, k))
                 for i in S for j in S if i < j for k in complement for l in complement if k < l)
        ok && push!(splits, sort(S))
    end
    return collect(splits)
end

"The internal splits of a row of 04_trees.csv: the side without taxon 7 (sorted) => length."
function s5_split_lengths(r)
    lengths = Dict{Vector{Int},Int}()
    for s in split(r["internal_splits"], ';'; keepempty=false)
        side, len = split(s, ':')
        A = Set(c - '0' for c in first(split(side, '|')))
        lengths[sort(collect(7 in A ? setdiff(Set(1:7), A) : A))] = parse(Int, len)
    end
    return lengths
end

"The tree metric d_ij (in the order of PAIRS) of the pendant lengths of 1-7 and the internal splits => length."
s5_metric(pendant, lengths) = [pendant[i] + pendant[j] + sum((l for (S, l) in lengths if (i in S) != (j in S)); init=0) for (i, j) in PAIRS]

"The tree metric d_ij (in the order of PAIRS) of a row of 04_trees.csv."
s5_tree_metric(r) = s5_metric(parse.(Int, split(r["pendant_lengths"], ';')), s5_split_lengths(r))

"The rows of 04_trees.csv of a census target, by (kind, index)."
s5_tree_rows(target) = Dict((r["kind"], parse(Int, r["index"])) => r
                            for r in rows(joinpath(TARGET_DIRECTORIES[target], "04_trees.csv")))

"The metric d' with d'_{g(i)g(j)} = d_ij, for a permutation g of the taxa."
function s5_permuted_metric(d, g)
    image = similar(d)
    for (k, (i, j)) in enumerate(PAIRS)
        image[PAIR_INDEX[minmax(g[i], g[j])]] = d[k]
    end
    return image
end

"Two splits, given by their sides without taxon 7, are compatible."
s5_compatible(A, B) = issubset(A, B) || issubset(B, A) || isempty(intersect(A, B))

"""
The edges of a census fan modulo lineality (the edges of its maximal cells), the pairs of rays whose trees have
compatible splits, and whether the rays of every maximal cell are pairwise compatible.
"""
function s5_compatibility(target)
    rays, lineality = fan_rays(fan(target))
    L = matrix(QQ, reduce(vcat, permutedims.(lineality)))
    cells = fan_cells(fan(target))
    splits_of = [s5_ray_splits(r) for r in rays]
    compatible_rays(cell) = all(s5_compatible(A, B) for a in cell for b in cell if a < b for A in splits_of[a] for B in splits_of[b])
    cell_edges = Set{Vector{Int}}()
    for cell in cells   # a maximal cell with three rays is a triangle
        union!(cell_edges, length(cell) == 3 ? [collect(e) for e in Combinatorics.combinations(cell, 2)] :
                           first(cell_faces(cell, rays, L))[2])
    end
    compatible_pairs = Set([a, b] for a in eachindex(rays) for b in a+1:length(rays) if compatible_rays([a, b]))
    return (cells=length(cells), edges=cell_edges, pairs=compatible_pairs, cells_compatible=all(compatible_rays, cells))
end

"""
The tree drawn in one block of a figure of G2-trees, or nothing if it draws none: the pendant lengths of the taxa 1-7
(drawn at the vertices v1-v7) and the internal splits (side without 7) => length. Each weight is read on its edge (near
its middle, not at its ends), or with `nearest` on its nearest edge, since the atlas figures draw some weights beside it.
"""
function s5_picture_tree(block; nearest=false)
    point = Dict(parse(Int, m[1]) => parse.(Float64, (m[2], m[3]))
                 for m in eachmatch(r"\\coordinate \(v(\d+)\) at \(([-0-9.]+),([-0-9.]+)\)", block))
    edges = [(String(m[1]), parse(Int, m[2]), parse(Int, m[3]))
             for m in eachmatch(r"\\draw\[g2tree (internal|pendant)[^\]]*\] \(v(\d+)\)--\(v(\d+)\)", block)]
    isempty(edges) && return nothing
    # The position t of the foot of p along the edge a-b (0 at a, 1 at b), and the distance from p to the segment.
    function foot(p, (_, a, b))
        (x, y), (u, v) = point[a], point[b]
        t = ((p[1] - x) * (u - x) + (p[2] - y) * (v - y)) / ((u - x)^2 + (v - y)^2)
        s = clamp(t, 0.0, 1.0)
        return t, hypot(p[1] - x - s * (u - x), p[2] - y - s * (v - y))
    end
    on = [Int[] for _ in edges]
    for m in eachmatch(r"\\node\[g2tree weight[^\]]*\] at \(([-0-9.]+),([-0-9.]+)\) \{(\d+)\}", block)
        p, w = parse.(Float64, (m[1], m[2])), parse(Int, m[3])
        if nearest
            push!(on[argmin(k -> last(foot(p, edges[k])), eachindex(edges))], w)
        else
            for k in eachindex(edges)
                t, distance = foot(p, edges[k])
                0.1 < t < 0.9 && distance < 1e-3 && push!(on[k], w)
            end
        end
    end
    all(ws -> length(ws) == 1, on) || error("a drawn tree with an edge without exactly one weight")
    # The taxa reached from a without crossing the edge a-b.
    function side(a, b)
        seen, todo = Set([a]), [a]
        while !isempty(todo)
            v = pop!(todo)
            for (_, x, y) in edges, (u, w) in ((x, y), (y, x))
                if u == v && !(w in seen) && Set([u, w]) != Set([a, b])
                    push!(seen, w)
                    push!(todo, w)
                end
            end
        end
        return sort([v for v in seen if v <= 7])
    end
    pendant = [only(on[only(k for (k, (kind, a, b)) in enumerate(edges) if kind == "pendant" && i in (a, b))]) for i in 1:7]
    lengths = Dict{Vector{Int},Int}()
    for (k, (kind, a, b)) in enumerate(edges)
        kind == "internal" || continue
        S = side(a, b)
        lengths[7 in S ? setdiff(1:7, S) : S] = only(on[k])
    end
    return (pendant, lengths)
end

"The trees drawn in figures/square_cell.tex, by block (\"Ray 13\", ..., \"The cell\"), as `s5_picture_tree` reads them."
function s5_drawn_trees(tex)
    trees = Dict{String,Tuple{Vector{Int},Dict{Vector{Int},Int}}}()
    for block in split(tex, "\n% ")[2:end]
        tree = s5_picture_tree(block)
        tree === nothing || (trees[String(first(split(block, '\n')))] = tree)
    end
    return trees
end

"The trees drawn in an atlas figure, by the panel label (the first capture of `label`): those of its tikzpictures."
function s5_panel_trees(tex, label)
    marks = [(m.offset, String(m[1])) for m in eachmatch(label, tex)]
    trees = Dict{String,Vector{Tuple{Vector{Int},Dict{Vector{Int},Int}}}}()
    for (n, (start, name)) in enumerate(marks)
        stop = n < length(marks) ? prevind(tex, first(marks[n+1])) : lastindex(tex)
        panel = get!(trees, name, Tuple{Vector{Int},Dict{Vector{Int},Int}}[])
        for picture in eachmatch(r"\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}"s, SubString(tex, start, stop))
            tree = s5_picture_tree(picture.match; nearest=true)
            tree === nothing || push!(panel, tree)
        end
    end
    return trees
end

function section_5()
    orbit_rows = rows(g2("grassmannian_orbits.csv"))
    ray_rows = [r for r in orbit_rows if r["kind"] == "ray"]
    cell_rows = [r for r in orbit_rows if r["kind"] == "cell"]
    orbit_row(t) = only(r for r in orbit_rows if r["orbit_type"] == t)
    members(r, field) = isempty(r["$(field)_members"]) ? Int[] : parse.(Int, split(r["$(field)_members"], ';'))
    field_name(field) = field == "qq" ? "QQ" : "GF(2)"
    targets = (("grassmannian", "qq"), ("grassmannian_gf2", "gf2"))
    Fqq, Fgf2 = fan("grassmannian"), fan("grassmannian_gf2")
    raysqq, linqq = fan_rays(Fqq)
    raysgf2, lingf2 = fan_rays(Fgf2)
    cellsqq, cellsgf2 = fan_cells(Fqq), fan_cells(Fgf2)
    tqq = orbit_types("grassmannian_orbits.csv", "qq")
    tgf2 = orbit_types("grassmannian_orbits.csv", "gf2")
    trees = Dict("qq" => s5_tree_rows("grassmannian"), "gf2" => s5_tree_rows("grassmannian_gf2"))
    tex(name) = read(repository("figures", name * ".tex"), String)
    hexagon_tex = tex("taxon_hexagon")

    subheading("The D12-orbits of the G2-trees (Figures fig:atlas_rays, fig:atlas_toric_cells_prime, fig:atlas_toric_cells_nontoric)")
    for (field, rays, lineality) in (("qq", raysqq, linqq), ("gf2", raysgf2, lingf2))
        same = all(eachindex(rays)) do k
            haskey(trees[field], ("ray", k)) || return false
            m = QQ.(s5_tree_metric(trees[field][("ray", k)]))
            rank_of(vcat([m, rays[k]], lineality)) == rank_of(lineality) + 1 && sort(s5_ray_splits(-m)) == sort(s5_ray_splits(rays[k]))
        end
        claim("$(field_name(field)): the trees of 04_trees.csv are the rays of the fan (the same line modulo lineality, the same splits)",
              same, "$(length(rays)) rays")
    end
    ray_label = Dict(String(chop(m[1]; head=1, tail=0)) => parse(Int, m[2]) for m in
                     eachmatch(r"\\boldmath \\\((\\[a-z]+(?:_\d)?)\\\)\}\{\\large\((\d+)\)\}", tex("ray_orbits_styled_catalog")))
    cell_label(name) = Dict(String(m[1]) => parse(Int, m[2]) for m in
                            eachmatch(r"\\oldstylenums\{(\d+)\}\}\{\\large\((\d+)\)\}", tex(name)))
    toric_label, nontoric_label = cell_label("toric_degeneration_prime_edge_trees"), cell_label("non_toric_cell_trees")
    check("types drawn in fig:atlas_rays (12 rays), fig:atlas_toric_cells_prime (1-20) and fig:atlas_toric_cells_nontoric (21-35)",
          [12, collect(1:20), collect(21:35)],
          [length(ray_label), sort(parse.(Int, collect(keys(toric_label)))), sort(parse.(Int, collect(keys(nontoric_label))))])
    figure_sizes = merge(ray_label, toric_label, nontoric_label)
    check("types whose label (orbit size) in the figures differs from the size in grassmannian_orbits.csv", String[],
          sort([r["orbit_type"] for r in orbit_rows if get(figure_sizes, r["orbit_type"], 0) != parse(Int, r["size"])]))
    # The trees themselves, against trop(G) and trop_2(G) (Section 2 checks them against the Dressian, and that the
    # vertices v1-v7 of these panels carry the shapes of the taxa 1-7).
    drawn_atlas = merge(s5_panel_trees(tex("ray_orbits_styled_catalog"), r"\\boldmath \\\(\\([a-z]+(?:_\d)?)\\\)\}\{\\large"),
                        s5_panel_trees(tex("toric_degeneration_prime_edge_trees"), r"\\oldstylenums\{(\d+)\}\}\{\\large"),
                        s5_panel_trees(tex("non_toric_cell_trees"), r"\\oldstylenums\{(\d+)\}\}\{\\large"))
    member_metrics(r) = Set(s5_tree_metric(trees[field][(r["kind"], k)]) for (_, field) in targets for k in members(r, field))
    check("types whose panel in these figures does not draw one tree, the tree of a member of its orbit in trop(G) or trop_2(G) " *
          "(grassmannian(_gf2)/04_trees.csv)", String[],
          sort([r["orbit_type"] for r in orbit_rows if length(get(drawn_atlas, r["orbit_type"], [])) != 1 ||
                !(s5_metric(only(drawn_atlas[r["orbit_type"]])...) in member_metrics(r))]))
    # The hexagon of Figure fig:shape_hexagon: its six sides are the solid lines, the dotted ones are the diagonals.
    hexagon_edges = Set(Set(parse.(Int, (m[1], m[2]))) for m in
                        eachmatch(r"\\draw\[draw=black,line width=\.6pt\] \(v(\d)\)--\(v(\d)\);", hexagon_tex))
    size_from_trees = Dict{String,Set{Int}}()
    images_are_orbits = true
    for (target, field) in targets
        perms = load_symmetry_artifacts(joinpath(TARGET_DIRECTORIES[target], "05_symmetry")).permutations
        claim("$(field_name(field)): the census acts by 12 permutations of the taxa, the symmetries of the hexagon of Figure fig:shape_hexagon fixing the heart 7",
              length(unique(perms)) == 12 &&
              all(g -> g[7] == 7 && Set(Set(g[a] for a in e) for e in hexagon_edges) == hexagon_edges, perms),
              "$(length(hexagon_edges)) sides of the hexagon")
        metric = Dict(key => s5_tree_metric(r) for (key, r) in trees[field])
        key_of = Dict(d => key for (key, d) in metric)
        for r in orbit_rows
            m = sort(members(r, field))
            isempty(m) && continue
            images = unique(s5_permuted_metric(metric[(r["kind"], first(m))], g) for g in perms)
            images_are_orbits &= sort([get(key_of, d, ("none", 0)) for d in images]) == [(r["kind"], k) for k in m]
            push!(get!(size_from_trees, r["orbit_type"], Set{Int}()), length(images))
        end
    end
    claim("the 12 permutations map the tree of one member of each orbit onto the trees of its members (QQ and GF(2))", images_are_orbits)
    check("types whose orbit size, from the D12-action on the trees, differs from the label in the figures", String[],
          sort([t for t in keys(figure_sizes) if get(size_from_trees, t, Set{Int}()) != Set([figure_sizes[t]])]))
    info("orbit sizes from the trees: " * join(["$(r["orbit_type"]) $(join(sort(collect(size_from_trees[r["orbit_type"]])), "/"))"
                                                for r in orbit_rows], ", "))

    subheading("The 12 types of rays (Figure fig:atlas_rays)")
    ray_trees(target) = Set(s5_tree_metric(r) for ((kind, _), r) in s5_tree_rows(target) if kind == "ray")
    qq_ray_trees = ray_trees("grassmannian")
    claim("the rays of the Dressian are the rays of trop(G) (the same normalized trees)", ray_trees("dressian") == qq_ray_trees,
          "$(length(qq_ray_trees)) trees")
    check("the rays of trop_2(G) that are not rays of trop(G): the orbit ε", sort(members(orbit_row("epsilon"), "gf2")),
          sort([k for ((kind, k), r) in trees["gf2"] if kind == "ray" && !(s5_tree_metric(r) in qq_ray_trees)]))
    check("types of rays: in characteristic ≠ 2, and only in characteristic 2", [11, 1],
          [count(r -> !isempty(r["qq_orbit"]), ray_rows), count(r -> isempty(r["qq_orbit"]), ray_rows)])
    println("           type      size  inner edges  shape")
    shapes = Dict{String,String}()
    # Three disjoint cherries without the heart, and a fourth split joining two of them: the heart hangs on the arm of the third.
    function heart_on_arm(splits)
        cherries = [s for s in splits if length(s) == 2]
        return length(splits) == 4 && length(cherries) == 3 && length(union(cherries...)) == 6 &&
               any(s -> s in [sort(union(a, b)) for a in cherries for b in cherries if a < b], splits)
    end
    for r in ray_rows
        ray = isempty(r["qq_orbit"]) ? raysgf2[parse(Int, r["gf2_orbit"])] : raysqq[parse(Int, r["qq_orbit"])]
        splits = s5_ray_splits(ray)
        # The splits are given by their sides without taxon 7. Two of them are nested when the heart 7 sits at an end
        # of the caterpillar, and disjoint when it sits in the middle.
        shape = length(splits) == 1 ? "split" :
                length(splits) == 2 ? (isempty(intersect(splits...)) ? "caterpillar, heart in the middle" : "caterpillar, heart at an end") :
                length(splits) == 3 && all(isempty(intersect(a, b)) for a in splits for b in splits if a != b) ? "snowflake, heart at the center" :
                heart_on_arm(splits) ? "snowflake with the heart on an arm (4 inner edges)" : "other ($(length(splits)) inner edges)"
        shapes[r["orbit_type"]] = shape
        @printf("           %-9s %4s  %11d  %s\n", r["orbit_type"], r["size"], length(splits), shape)
    end
    expected_shape(t) = startswith(t, "alpha") ? "split" : startswith(t, "beta") ? "caterpillar, heart at an end" :
                        startswith(t, "gamma") ? "caterpillar, heart in the middle" : t == "delta" ? "snowflake, heart at the center" :
                        "snowflake with the heart on an arm (4 inner edges)"
    claim("α are splits, β caterpillars with the heart at an end, γ with the heart in the middle, δ and ε snowflakes",
          all(shapes[t] == expected_shape(t) for t in keys(shapes)))
    info("the trees of one ray per type (04_trees.csv): pendant lengths of the taxa 1-7, internal splits:length")
    for r in ray_rows
        field = isempty(r["qq_orbit"]) ? "gf2" : "qq"
        k = parse(Int, r["$(field)_orbit"])
        tree = trees[field][("ray", k)]
        info(@sprintf("%-9s %-11s %-15s %s", r["orbit_type"], "$(field_name(field)) ray $k", tree["pendant_lengths"],
                      replace(tree["internal_splits"], ';' => ", ")))
    end
    check("δ is fixed under D12: its orbit size, from the D12-action on its tree", [1], sort(collect(size_from_trees["delta"])))
    check("the ray type that appears only in characteristic 2", "epsilon", only(r["orbit_type"] for r in ray_rows if isempty(r["qq_orbit"])))

    subheading("The 35 types of maximal cells (Figures fig:atlas_toric_cells_prime and fig:atlas_toric_cells_nontoric)")
    check("types of maximal cells", collect(1:35), sort(parse.(Int, [r["orbit_type"] for r in cell_rows])))
    toric_of = Dict{Tuple{String,Int},Bool}()     # (field, type) => toric
    for (target, field) in targets
        census = rows(g2(target, "06_canonical_interior_degenerations.csv"))
        orbit_type = Dict(r["$(field)_orbit"] => parse(Int, r["orbit_type"]) for r in cell_rows if !isempty(r["$(field)_orbit"]))
        for r in census
            toric_of[(field, orbit_type[r["orbit"]])] = r["toric"] == "true"
        end
    end
    toric_types(field) = sort([t for ((f, t), toric) in toric_of if f == field && toric])
    info("type size  fields      toric      tree of one cell (04_trees.csv): pendant lengths, internal splits:length")
    for r in sort(cell_rows; by=r -> parse(Int, r["orbit_type"]))
        t = parse(Int, r["orbit_type"])
        fields = [field for (_, field) in targets if !isempty(r["$(field)_orbit"])]
        status = unique(get(toric_of, (field, t), missing) for field in fields)
        tree = trees[first(fields)][("cell", first(members(r, first(fields))))]
        info(@sprintf("%-4s %4s  %-11s %-10s %-15s %s", r["orbit_type"], r["size"], join(field_name.(fields), ", "),
                      length(status) == 1 ? (only(status) === true ? "yes" : only(status) === false ? "no" : "?") : "differs",
                      tree["pendant_lengths"], replace(tree["internal_splits"], ';' => ", ")))
    end
    check("types that give a toric degeneration (in some characteristic)", collect(1:20), sort(union(toric_types("qq"), toric_types("gf2"))))
    check("types that give a toric degeneration over QQ (18-20 exist only in characteristic 2)", collect(1:17), toric_types("qq"))
    check("types that give a toric degeneration over GF(2)", collect(1:20), toric_types("gf2"))
    check("types that give no toric degeneration (in any characteristic)", collect(21:35),
          [t for t in 1:35 if !any(get(toric_of, (field, t), false) for (_, field) in targets)])
    only_gf2 = sort([parse(Int, r["orbit_type"]) for r in cell_rows if isempty(r["qq_orbit"])])
    only_qq = sort([parse(Int, r["orbit_type"]) for r in cell_rows if isempty(r["gf2_orbit"])])
    check("types only in characteristic 2", [18, 19, 20, 35], only_gf2)
    check("types only in characteristic ≠ 2", [32, 33, 34], only_qq)
    quads = [[k for (k, cell) in enumerate(cells) if length(cell) == 4] for cells in (cellsqq, cellsgf2)]
    type30 = sort(members(orbit_row("30"), "qq"))
    check("type 30 is the orbit of the three quadrilaterals: over QQ, over GF(2)", quads,
          [type30, sort(members(orbit_row("30"), "gf2"))])

    subheading("The split type α2 (Figure fig:shape_hexagon)")
    check("α2 has an orbit of size 3 (D12-action on its trees)", [3], sort(collect(size_from_trees["alpha_2"])))
    # The colors and forms of the six shapes, read from the figure.
    taxon = Dict{Int,String}()
    for m in eachmatch(r"g2tree (dark|light),(rectangle|circle),[^\n]*at \(v(\d)\)", hexagon_tex)
        taxon[parse(Int, m[3])] = (m[1] == "dark" ? "blue" : "yellow") * " " * (m[2] == "rectangle" ? "square" : "circle")
    end
    for m in eachmatch(r"shift=\{\(v(\d)\)\}[^\n]*\n\\path\[fill=g2Taxon(Blue|Yellow)", hexagon_tex)
        taxon[parse(Int, m[1])] = lowercase(m[2]) * " triangle"
    end
    info("the shapes of Figure fig:shape_hexagon: " * join(["$i $(get(taxon, i, "?"))" for i in 1:6], ", ") * ", 7 the heart")
    alpha2 = sort(members(orbit_row("alpha_2"), "qq"))
    alpha2_splits = [s5_ray_splits(raysqq[k]) for k in alpha2]
    claim("α2: each of its 3 trees has one split, grouping the heart with two shapes against the other four leaves",
          all(s -> length(s) == 1 && length(only(s)) == 4, alpha2_splits), "QQ rays $(join(alpha2, ", "))")
    groups = sort([setdiff(1:6, only(s)) for s in alpha2_splits if length(s) == 1])
    form(i) = last(split(get(taxon, i, "?")))
    claim("α2: \"grouping the heart with two same shapes\" (main.tex), of the two colors",
          !isempty(groups) && all(g -> length(unique(form.(g))) == 1 && length(unique(first.(split.(get.(Ref(taxon), g, "?"))))) == 2, groups),
          "the heart is grouped with " * join(["$(g[1]) and $(g[2]) ($(get(taxon, g[1], "?")), $(get(taxon, g[2], "?")))" for g in groups], "; "))

    subheading("Table tab:cell_ray_orbits: the ray types of the maximal cells")
    table_tex = read(repository("tables", "cell_ray_orbits.tex"), String)
    table_body = table_tex[last(findfirst("\\toprule", table_tex)):first(findfirst("\\bottomrule", table_tex))]
    table = Dict{Int,Vector{String}}()
    for chunk in split(table_body, r"&|\\tabularnewline")
        m = match(r"oldstylenums\{(\d+)\}", chunk)
        m === nothing && continue
        table[parse(Int, m[1])] = [String(x[1]) for x in eachmatch(r"\$\\([a-z]+(?:_\d)?)\$", chunk)]
    end
    check("cell types listed in tables/cell_ray_orbits.tex", collect(1:35), sort(collect(keys(table))))
    differ = Int[]
    for r in sort(cell_rows; by=r -> parse(Int, r["orbit_type"]))
        t = parse(Int, r["orbit_type"])
        types_by_field = String[]
        for (field, cells, types) in (("qq", cellsqq, tqq), ("gf2", cellsgf2, tgf2))
            isempty(r["$(field)_orbit"]) && continue
            ray_types = sort([types[("ray", ray)] for ray in cells[first(members(r, field))]])
            isempty(types_by_field) && append!(types_by_field, ray_types)
            ray_types == sort(get(table, t, String[])) || push!(differ, t)
        end
        entry = get(table, t, String[])
        @printf("           type %-3s size %-3s rays %s%s\n", r["orbit_type"], r["size"], join(types_by_field, ", "),
                sort(entry) != types_by_field ? "   (table: " * join(entry, ", ") * ")" :
                entry != types_by_field ? "   (table, in cyclic order: " * join(entry, ", ") * ")" : "")
    end
    check("cell types whose ray types (over QQ or GF(2)) differ from the table", Int[], sort(unique(differ)))
    # The order of the entry of type 30: along the boundary of the quadrilateral 29.
    square = cellsqq[first(type30)]
    square_faces, _ = cell_faces(square, raysqq, matrix(QQ, reduce(vcat, permutedims.(linqq))))
    neighbours(v) = sort([only(setdiff(e, [v])) for e in square_faces[2] if v in e])
    boundary = [first(square)]
    while length(boundary) < length(square)
        push!(boundary, first(setdiff(neighbours(boundary[end]), boundary)))
    end
    cycle_types = [tqq[("ray", k)] for k in boundary]
    turns = [circshift(c, s) for c in (cycle_types, reverse(cycle_types)) for s in 0:length(c)-1]
    claim("type 30 in the table lists the ray types in the cyclic order of the quadrilateral (up to rotation and reflection)",
          get(table, 30, String[]) in turns, "boundary of cell $(first(type30)): " * join(["$k $t" for (k, t) in zip(boundary, cycle_types)], ", "))

    compatibility = Dict(field => s5_compatibility(target) for (target, field) in targets)
    for ((target, field), f_vector) in zip(targets, ((51, 241, 248), (54, 253, 254)))
        c = compatibility[field]
        types = field == "qq" ? tqq : tgf2
        claim("$(field_name(field)): the rays of each of the $(c.cells) maximal cells have pairwise compatible splits", c.cells_compatible)
        check("$(field_name(field)): edges (f-vector $f_vector)", f_vector[2], length(c.edges))
        claim("$(field_name(field)): compatible splits do not suffice: the pairs of rays with compatible splits contain the edges, and are more",
              issubset(c.edges, c.pairs) && length(c.pairs) > length(c.edges), "$(length(c.pairs)) compatible pairs")
        info("compatible pairs that are not edges: " *
             join(["$a-$b ($(types[("ray", a)]), $(types[("ray", b)]))" for (a, b) in sort(collect(setdiff(c.pairs, c.edges)))], ", "))
    end
    alpha4 = sort(members(orbit_row("alpha_4"), "qq"))
    delta = only(members(orbit_row("delta"), "qq"))
    claim("δ and α4 are compatible: the splits of δ (ray $delta) are compatible with those of each ray α4 ($(join(alpha4, ", "))), and each pair is an edge",
          all(s5_compatible(A, B) for a in alpha4 for A in s5_ray_splits(raysqq[delta]) for B in s5_ray_splits(raysqq[a])) &&
          all(sort([delta, a]) in compatibility["qq"].edges for a in alpha4))
    for t in ("13", "14")
        cell = cellsqq[first(members(orbit_row(t), "qq"))]
        cell_types = [tqq[("ray", k)] for k in cell]
        claim("type $t is a triangle with one ray δ and one ray α4, so δ and α4 span one of its edges",
              length(cell) == 3 && count(==("delta"), cell_types) == 1 && count(==("alpha_4"), cell_types) == 1,
              "cell $(first(members(orbit_row(t), "qq"))): " * join(["$k $s" for (k, s) in zip(cell, cell_types)], ", "))
    end

    subheading("The quadrilateral of type 30 (Figure fig:square_cell)")
    square_tex = tex("square_cell")
    figure_cell = parse(Int, match(r"cell (\d+): rays", square_tex)[1])
    check("the cell drawn in figures/square_cell.tex: its type and its rays (main.tex: rays 13, 14, 17, 18)",
          ["30", [13, 14, 17, 18]], [get(tqq, ("cell", figure_cell), "?"), cellsqq[figure_cell]])
    # The rays of the figure, clockwise from the top left: by the angle of the centers of their circles.
    centers = Dict(parse(Int, m[1]) => atand(parse(Float64, m[3]), parse(Float64, m[2]))
                   for m in eachmatch(r"% Ray (\d+)\n\\draw\[g2square circle\] \(([-0-9.]+),([-0-9.]+)\)", square_tex))
    clockwise = sort(collect(keys(centers)); by=k -> mod(135 - centers[k], 360))
    check("the rays of the figure, clockwise from the top left (main.tex: 13, 14, 17, 18)", [13, 14, 17, 18], clockwise)
    check("the edges of the quadrilateral join consecutive rays of this cyclic order",
          sort([sort([clockwise[i], clockwise[mod1(i + 1, end)]]) for i in eachindex(clockwise)]),
          figure_cell == first(type30) ? square_faces[2] : Vector{Int}[])
    check("G2-trees on the quadrilateral: vertices + edges + interior", [4, 4, 1], length.(square_faces))
    label_of = Dict(parse(Int, m[1]) => String(m[2]) for m in
                    eachmatch(r"% Ray (\d+)\n[^\n]*\n[^\n]*\\boldmath \\\(\\([a-z]+_\d'?)\\\)", square_tex))
    labels = [get(label_of, k, "?") for k in clockwise]
    check("the labels of the rays in the figure are their types: two vertices in the orbit γ1", replace.(labels, "'" => ""),
          [tqq[("ray", k)] for k in clockwise])
    ray_of = Dict(zip(labels, clockwise))
    paper_d = Dict{Tuple{Int,Int},Int}()
    for (value, list) in ((2, "12 16 26 34 35 36 37 45 67"), (4, "13 17 23 27 46 47 56 57"), (6, "14 15 24 25"))
        for s in split(list)
            paper_d[(parse(Int, s[1]), parse(Int, s[2]))] = value
        end
    end
    d = [paper_d[p] for p in PAIRS]
    support = sort(s5_ray_splits(-QQ.(d)))
    side_names = join(join.(support), ", ")
    drawn = s5_drawn_trees(square_tex)
    drawn_names = vcat(["Ray $k" for k in clockwise], ["The cell"])
    check("the inner edges of the 5 trees drawn: all those of d, $side_names (sides without 7), some of length zero",
          [support], unique([sort(collect(keys(get(drawn, n, (Int[], Dict{Vector{Int},Int}()))[2]))) for n in drawn_names]))
    check("the union of the splits of the trees of the four rays: the splits of d", support,
          sort(unique(vcat([s5_ray_splits(raysqq[k]) for k in clockwise]...))))
    stored_rows = vcat([trees["qq"][("ray", k)] for k in clockwise], [trees["qq"][("cell", figure_cell)]])
    stored_trees = [(parse.(Int, split(r["pendant_lengths"], ';')), merge(Dict(S => 0 for S in support), s5_split_lengths(r)))
                    for r in stored_rows]
    check("trees drawn whose weights differ from 04_trees.csv (rays 13, 14, 17, 18 and the cell)", String[],
          [n for (n, tree) in zip(drawn_names, stored_trees) if get(drawn, n, nothing) != tree])
    info("tree               pendant lengths of 1-7   lengths of the inner edges $side_names")
    for (n, label, (pendant, lengths)) in zip(drawn_names, vcat(labels, ["d"]), stored_trees)
        info(@sprintf("%-18s %-24s %s", "$label ($(lowercase(n)))", join(pendant, ";"), join([lengths[S] for S in support], ", ")))
    end
    tree_of(k) = QQ.(s5_tree_metric(trees["qq"][("ray", k)]))
    s = fill(QQ(2), length(PAIRS))     # the star tree with all pendant edges of length 1
    check("d = (γ1 + γ1' - 2s)/6 = γ2 + γ3 (the trees of the rays of the figure)", [QQ.(d), QQ.(d)],
          [(tree_of(ray_of["gamma_1"]) + tree_of(ray_of["gamma_1'"]) - 2 * s) .* QQ(1, 6), tree_of(ray_of["gamma_2"]) + tree_of(ray_of["gamma_3"])])
    row = Dict(parse(Int, r["number"]) => r for r in rows(g2("maximal_cell_tree_metrics.csv")))[30]
    check("d_ij (lexicographic order): row 30 of maximal_cell_tree_metrics.csv", d, [parse(Int, row["d$(i)$(j)"]) for (i, j) in PAIRS])
    check("d_ij: the normalized tree of the cell (grassmannian/04_trees.csv)", d, s5_tree_metric(trees["qq"][("cell", figure_cell)]))
    check("maximal cells of trop(G) with -d in their relative interior", [figure_cell], locate(locator(Fqq), -QQ.(d)))
end
