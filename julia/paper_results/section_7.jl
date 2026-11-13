# Section 7 of the paper: run by paper_results.jl, which defines the helpers used here.

"The local flag cell atlas: one row per 3-cell (location_kind edge) and per node 2-cell (node) of each fiber complex."
const s7_ATLAS = study_path("g2", "atlas", "flag", "local_flag_cell_atlas.csv")

"""
The letters of the vertices R_g^X of F_{A_1} (type 21) and F_D (type 5) in Figures fig:triangle_cell and
fig:prism_projection, as listed in figures/triangle_incidence.tex and figures/triangle_incidence_5.tex.
"""
const s7_LETTERS = Dict(
    "21" => Dict((7, 'A') => 'A', (7, 'B') => 'B', (12, 'A') => 'C', (12, 'B') => 'D', (40, 'A') => 'F', (40, 'B') => 'E'),
    "5" => Dict((2, 'A') => 'B', (2, 'B') => 'A', (2, 'C') => 'C', (13, 'A') => 'E', (13, 'B') => 'F', (13, 'C') => 'D',
                (14, 'A') => 'H', (14, 'B') => 'I', (14, 'C') => 'G'))

s7_set(s) = Set(parse.(Int, split(s, ';')))
"The inner edge [a, b] of the tree of C under an atlas 3-cell, or nothing for a pendant edge."
s7_edge(r) = (hit = match(r"^\(:edge, (\d+), (\d+)\)$", r["location"]); hit === nothing ? nothing : sort(parse.(Int, hit.captures)))
"The leaf k of an atlas 3-cell over the pendant edge k, or nothing for an inner edge."
s7_leaf(r) = (hit = match(r"^\(:ray, (\d+)\)$", r["location"]); hit === nothing ? nothing : parse(Int, hit[1]))
"An atlas 3-cell by its edge of the tree of C: the leaf k of a pendant edge, a-b for the inner edge ab."
s7_place(r) = (edge = s7_edge(r); edge === nothing ? string(s7_leaf(r)) : "$(edge[1])-$(edge[2])")
"The 3-cells of a fiber complex, the pendant edges by leaf first, then the inner edges."
s7_sorted_cells(F) = sort(F.cells; by=r -> (edge = s7_edge(r); edge === nothing ? (0, s7_leaf(r), 0) : (1, edge[1], edge[2])))
"The vertices of an atlas row: (Grassmannian ray g, node letter X) for each R_g^X, and the leaves k of the shapes e_k."
s7_symbols(r) = ([(parse(Int, x[1]), only(x[2])) for x in eachmatch(r"R_\{(\d+)\}\^\{(\w)\}", r["symbolic_vertices"])],
                 [parse(Int, x[1]) for x in eachmatch(r"e_\{(\d+)\}", r["symbolic_vertices"])])
"The name of an atlas row in the letters of the paper: its letters and the leaves of its shapes, sorted (? for a vertex without a letter)."
s7_name(r, letters) = join(sort(vcat(string.(last(s7_symbols(r))), [string(get(letters, x, '?')) for x in first(s7_symbols(r))])))

"""
The fiber complex of an atlas label: its 3-cells, its node 2-cells by node, the node 2-cells inside each 3-cell, the
leaves at each node of the tree of C, and the rays of C (support, in the atlas numbering: that of QQ, and 52-54 for the
rays of characteristic 2 only).
"""
function s7_fiber(atlas, label)
    entries = [r for r in atlas if r["label"] == label]
    fields = split(first(entries)["fields"], '+')
    vertex_set(r) = s7_set(r["$(first(fields))_flag_rays"])
    cells = [r for r in entries if r["location_kind"] == "edge"]
    nodes = Dict(parse(Int, split(r["location"])[2]) => vertex_set(r) for r in entries if r["location_kind"] == "node")
    inside = Dict(r["location"] => sort([n for (n, s) in nodes if s ⊆ vertex_set(r)]) for r in cells)
    leaves = Dict(n => sort(Int[s7_leaf(r) for r in cells if s7_leaf(r) !== nothing && inside[r["location"]] == [n]])
                  for n in keys(nodes))
    support = parse.(Int, split(first(entries)["support"], ';'))
    return (; entries, fields, support, vertices=vertex_set, cells, nodes, inside, leaves)
end

"The split of the leaves made by the inner edge ab of the tree of a fiber complex, as an unordered pair of leaf sets."
function s7_split(F, a, b)
    adjacent = Dict(n => Int[] for n in keys(F.nodes))
    for r in F.cells
        edge = s7_edge(r)
        edge === nothing || (push!(adjacent[edge[1]], edge[2]); push!(adjacent[edge[2]], edge[1]))
    end
    side, seen, stack = Set{Int}(), Set([a]), [a]
    while !isempty(stack)
        x = pop!(stack)
        union!(side, F.leaves[x])
        for y in adjacent[x]
            (x == a && y == b) || y in seen || (push!(seen, y); push!(stack, y))
        end
    end
    return Set([side, setdiff(Set(1:7), side)])
end

"The inner splits of the trees of a Grassmannian census (04_trees.csv), as unordered pairs of leaf sets, by (kind, index)."
function s7_splits(target)
    pair(S) = Set([Set(S), setdiff(Set(1:7), S)])
    return Dict((r["kind"], parse(Int, r["index"])) => Set{Set{Set{Int}}}(pair(S) for S in keys(tree_splits(r)))
                for r in rows(g2(target, "04_trees.csv")))
end

"The cells of a stage-2 coarsest fan by their rays (the numbering of 07_projection_cells.csv and of the census)."
s7_cell_rays(target) = Dict(parse(Int, r["coarsest_cell"]) => s7_set(r["extremal_original_rays"]) for r in rows(g2(target, "02_coarsest_cells.csv")))

"The rays of a flag fan over each ray of its Grassmannian (0: over the lineality space), from the projection certificate."
function s7_over_ray(target)
    over = Dict{Int,Set{Int}}()
    for r in TOML.parsefile(g2(target, "07_projection_certificate.toml"))["rays"]
        push!(get!(over, r["target_ray"], Set{Int}()), r["source_ray"])
    end
    return over
end

"""
The edges of the polygon C, a maximal cell of a Grassmannian census given by its rays: every pair of rays of a triangle;
for a quadrilateral, the pairs that lie in another maximal cell too (its walls; a diagonal lies in no other cell).
"""
function s7_polygon_edges(cell_rays, c)
    R = sort(collect(cell_rays[c]))
    pairs = [Set([g, h]) for (i, g) in enumerate(R) for h in R[i+1:end]]
    length(R) == 3 && return pairs
    return [p for p in pairs if any(d != c && p ⊆ other for (d, other) in cell_rays)]
end

"""
The f-vector of F_C counted as the decorations of the trees of the faces of C: an apple on one of the 7 pendant edges of
the star tree (a vertex); on a node or an edge of the tree of each vertex of C (a vertex or an edge), of each edge of C
(an edge or a 2-cell) and of C itself (a 2-cell or a 3-cell). The tree of an edge of C has the splits of its two ends.
"""
function s7_decorations(splits, rays_of_C, edges_of_C, cell)
    node_count(s) = length(s) + 1
    edge_count(s) = 7 + length(s)
    vertex_trees = [splits[("ray", g)] for g in rays_of_C]
    edge_trees = [reduce(union, [splits[("ray", g)] for g in edge]) for edge in edges_of_C]
    tree = splits[("cell", cell)]
    return (7 + sum(node_count, vertex_trees), sum(edge_count, vertex_trees) + sum(node_count, edge_trees),
            sum(edge_count, edge_trees) + node_count(tree), edge_count(tree))
end

"""
Whether every 3-cell of the fiber complex F has as vertices the images of the endpoints of its edge e of the tree of C
under the contractions to the trees of the vertices of C: two vertices over a vertex whose tree keeps the split of e,
one over the others (one over each vertex of C for a pendant edge).
"""
function s7_contractions(F, splits)
    return all(F.cells) do r
        over_vertex = [count(x -> first(x) == g, first(s7_symbols(r))) for g in F.support]
        edge = s7_edge(r)
        edge === nothing && return over_vertex == ones(Int, length(F.support))
        return over_vertex == [s7_split(F, edge...) in splits[("ray", g)] ? 2 : 1 for g in F.support]
    end
end

"The node triangles of a fiber complex, in the letters of the paper, each with the 3-cells that contain it."
s7_incidence(F, name_of) = join(["$(name_of["node $n"]) (node $n): " *
                                 join(sort([name_of[r["location"]] for r in F.cells if n in F.inside[r["location"]]]), ", ")
                                 for n in sort(collect(keys(F.nodes)))], "; ")

function section_7()
    flag_rows = rows(g2("flag_orbits.csv"))
    grass_orbit_rows = rows(g2("grassmannian_orbits.csv"))
    fields = (("qq", "QQ", "flag", "grassmannian"), ("gf2", "GF(2)", "flag_gf2", "grassmannian_gf2"))
    over = Dict(f => s7_over_ray(target) for (f, _, target, _) in fields)
    splits = Dict(f => s7_splits(grass) for (f, _, _, grass) in fields)
    grass_types = Dict(f => orbit_types("grassmannian_orbits.csv", f) for (f, _, _, _) in fields)
    stars = Dict(f => Set(parse(Int, x) for r in flag_rows if r["kind"] == "ray" && r["orbit_type"] == "star"
                          for x in split(r["$(f)_members"], ';') if !isempty(x)) for (f, _, _, _) in fields)
    atlas = rows(s7_ATLAS)
    labels = sort(unique(r["label"] for r in atlas if r["location_kind"] in ("edge", "node")); by=l -> parse(Int, l))
    fibers = Dict(l => s7_fiber(atlas, l) for l in labels)

    subheading("The vertices and maximal cells of trop(F) as decorated G2-trees")
    rays = [r for r in flag_rows if r["kind"] == "ray" && !isempty(r["qq_orbit"])]
    check("vertices: star trees with an apple on an edge, decorated trees with an apple on a node", [7, 134],
          [sum(parse(Int, r["size"]) for r in rays if r["orbit_type"] == "star"),
           sum(parse(Int, r["size"]) for r in rays if r["orbit_type"] != "star")])
    inner_nodes = Dict(k => length(s) + 1 for ((kind, k), s) in splits["qq"] if kind == "ray")
    check("over each ray of trop(G): one ray of trop(F) per inner node of its tree; their total", (true, 134),
          (all(length(get(over["qq"], k, Set{Int}())) == n for (k, n) in inner_nodes), sum(values(inner_nodes))))
    grass_ray_types = Set(r["orbit_type"] for r in grass_orbit_rows if r["kind"] == "ray" && !isempty(r["qq_orbit"]))
    node_types = Set(r["orbit_type"] for r in rays if r["orbit_type"] != "star")
    check("the trees of the 134 node-decorated rays: the ray types of trop(G) in Figure fig:atlas_rays (without ε)", (11, true),
          (length(node_types), node_types == grass_ray_types))
    info("rays of trop(G) by type × rays of trop(F) over each (one per inner node of its tree):")
    for (f, name, _, _) in fields
        terms = String[]
        for r in grass_orbit_rows
            r["kind"] == "ray" && !isempty(r["$(f)_members"]) || continue
            members = parse.(Int, split(r["$(f)_members"], ';'))
            per_tree = unique(length(get(over[f], k, Set{Int}())) for k in members)
            push!(terms, "$(r["orbit_type"]) $(length(members))×$(join(per_tree, '/'))")
        end
        node_rays = sum(length(s) for (k, s) in over[f] if k != 0)
        info("$name: " * join(terms, " + ") * " = $node_rays; with the 7 star trees $(7 + node_rays)")
    end
    for (f, name, _, _) in fields
        edges_of = Dict(k => 7 + length(s) for ((kind, k), s) in splits[f] if kind == "cell")
        present = [r for r in grass_orbit_rows if r["kind"] == "cell" && !isempty(r["$(f)_orbit"])]
        flag_cells = sum(parse(Int, r["size"]) for r in flag_rows if r["kind"] == "cell" && !isempty(r["$(f)_orbit"]))
        paper_count = f == "qq" ? 2714 : 2792
        check("maximal cells of $(f == "qq" ? "trop(F)" : "trop_2(F)"), and Σ (orbit size) × (edges of the tree of C) " *
              "over the maximal cells C of $(f == "qq" ? "trop(G)" : "trop_2(G)")", (paper_count, paper_count),
              (flag_cells, sum(parse(Int, r["size"]) * edges_of[parse(Int, first(split(r["$(f)_members"], ';')))] for r in present)))
    end
    cell_types(f) = sort(unique(parse(Int, r["orbit_type"]) for r in flag_rows if r["kind"] == "cell" && !isempty(r["$(f)_orbit"])))
    # Known to the authors and left as is: a warning, not a failure.
    warned("the types of the trees of the 2714 maximal cells of trop(F) (main.tex: \"of type 1-35\"; but 18-20 and 35 " *
          "\"appear only in characteristic 2\")", collect(1:35), cell_types("qq"))
    claim("the types 1-35 are those of trop(F) and trop_2(F) together", sort(union(cell_types("qq"), cell_types("gf2"))) == 1:35)

    subheading("The fiber complexes F_C over the maximal cells C of trop(G)")
    claim("the 7 rays of trop(F) over the lineality space of trop(G) are the star trees with an apple on an edge (QQ and GF(2))",
          all(length(stars[f]) == 7 && over[f][0] == stars[f] for (f, _, _, _) in fields))
    fiber = Dict{String,Any}()
    for (f, _, target, grass) in fields
        flag_cell_rays, grass_cell_rays = s7_cell_rays(target), s7_cell_rays(grass)
        vertices_over, cells_over = Dict{Int,Set{Int}}(), Dict{Int,Int}()
        for r in rows(g2(target, "07_projection_cells.csv"))
            c = parse(Int, r["grassmannian_cell"])
            union!(get!(vertices_over, c, Set{Int}()), flag_cell_rays[parse(Int, r["flag_cell"])])
            cells_over[c] = get(cells_over, c, 0) + 1
        end
        fiber[f] = (; flag_cell_rays, grass_cell_rays, vertices_over, cells_over)
    end
    claim("the 7 star vertices lie in the fiber complex F_C of every maximal cell C of trop(G)",
          all(stars[f] ⊆ V for (f, _, _, _) in fields for V in values(fiber[f].vertices_over)),
          join(["$(length(fiber[f].vertices_over)) cells over $name" for (f, name, _, _) in fields], ", "))
    fA, fB, fC, fD, fE, fF = (13, 33, 31, 10), (14, 36, 34, 11), (15, 38, 35, 11), (16, 39, 35, 11), (16, 40, 36, 11), (19, 52, 45, 11)
    paper_table = [("A_{1}", fA, [21, 34]), ("A_{2}", fA, [32]), ("A_{3}", fA, [33]),
                   ("B_{1}", fB, [1, 2]), ("B_{2}", fB, [6]), ("B_{3}", fB, [7, 12, 17]), ("B_{4}", fB, [9, 26]), ("B_{5}", fB, [11, 24]),
                   ("B_{6}", fB, [15, 23, 25, 27, 31]),
                   ("C_{1}", fC, [3]), ("C_{2}", fC, [4]), ("C_{3}", fC, [8]), ("C_{4}", fC, [10]), ("C_{5}", fC, [13]), ("C_{6}", fC, [14]),
                   ("C_{7}", fC, [16, 22, 29]), ("C_{8}", fC, [28]),
                   ("D", fD, [5]), ("E_{1}", fE, [18]), ("E_{2}", fE, [19]), ("E_{3}", fE, [20]), ("E_{4}", fE, [35]), ("F", fF, [30])]
    paper_rows = [(fv, sort(reduce(vcat, [types for (_, other, types) in paper_table if other == fv]))) for fv in (fA, fB, fC, fD, fE, fF)]
    decorated = Dict{Int,Set{NTuple{4,Int}}}()
    for (f, _, _, _) in fields
        cell_rays = fiber[f].grass_cell_rays
        for (c, R) in cell_rays
            push!(get!(decorated, parse(Int, grass_types[f][("cell", c)]), Set{NTuple{4,Int}}()),
                  s7_decorations(splits[f], R, s7_polygon_edges(cell_rays, c), c))
        end
    end
    check("the cells of F_C decorate the trees of the faces of C: their f-vectors, with the types of C (as in Table tab:local_flag_complex)",
          paper_rows, [(fv, sort([t for (t, s) in decorated if fv in s])) for fv in sort(collect(reduce(union, values(decorated))))])
    for (f, name, _, _) in fields
        edges_of = Dict(k => 7 + length(s) for ((kind, k), s) in splits[f] if kind == "cell")
        check("over $name: numbers of maximal cells of the fiber complexes", [10, 11], sort(unique(values(fiber[f].cells_over))))
        claim("over $name: the 3-cells of each F_C are the apples on the edges of the tree of C (7 pendant + inner)",
              fiber[f].cells_over == edges_of, "$(length(edges_of)) cells C")
    end
    claim("every F_C: a 3-cell over the inner edge ab contains the node 2-cells of a and b, a pendant one exactly one " *
          "(the 3-cells are glued along them; the incidence graph is the tree of C with its edges subdivided)",
          all(F -> all(r -> let edge = s7_edge(r)
                               edge === nothing ? length(F.inside[r["location"]]) == 1 : F.inside[r["location"]] == edge
                           end, F.cells), values(fibers)), "$(length(labels)) types")
    claim("every node 2-cell has the shape of C and one vertex over each vertex of C (it bijects to C)",
          all(r -> r["polytope"] == r["kind"] && sort(first.(first(s7_symbols(r)))) == sort(fibers[r["label"]].support) &&
                   length(fibers[r["label"]].vertices(r)) == length(fibers[r["label"]].support),
              [r for r in atlas if r["location_kind"] == "node"]))
    claim("the other vertices of F_C decorate the nodes of the trees of the vertices of C (QQ and GF(2))",
          all(setdiff(V, stars[f]) == reduce(union, [get(over[f], g, Set{Int}()) for g in fiber[f].grass_cell_rays[c]])
              for (f, _, _, _) in fields for (c, V) in fiber[f].vertices_over))

    subheading("The 23 combinatorial types of fiber complexes (Proposition and Table tab:local_flag_complex)")
    words = Dict("tetrahedron" => "tetrahedra", "square pyramid" => "square pyramids", "triangular prism" => "triangular prisms")
    paper_cells = Dict(fA => [("tetrahedra", 10)], fB => [("tetrahedra", 11)], fC => [("square pyramids", 1), ("tetrahedra", 10)],
                       fD => [("tetrahedra", 10), ("triangular prisms", 1)], fE => [("square pyramids", 2), ("tetrahedra", 9)],
                       fF => [("square pyramids", 7), ("triangular prisms", 4)])
    labelled = nothing
    if QUICK
        skipped("23 combinatorial types of fiber complexes, 6 distinct f-vectors, Table tab:local_flag_complex from the face lattices", "--quick")
    else
        m = run_in_module(project_path("to_tex", "local_flag_complex_table.jl"); name=:LocalFlagComplexes)
        types = call(m, :local_flag_types)
        check("combinatorial types of fiber complexes, and distinct f-vectors", [23, 6],
              [length(types), length(unique(t.complex.f_vector for t in types))])
        labelled = call(m, :local_flag_labelled_types, types)
        check("Table tab:local_flag_complex: label, f-vector and orbits of every type (A_1 over 21, D over 5)", paper_table,
              [(replace(t.label, "\\mathfrak " => ""), Tuple(t.complex.f_vector), sort(unique(x.number for x in t.members))) for t in labelled])
        census_words = Dict("3-simplex" => "tetrahedra", "Pyramid(4-gon)" => "square pyramids", "Δ1 × Δ2" => "triangular prisms")
        check("Table tab:local_flag_complex: the maximal cells of each f-vector", [(fv, paper_cells[fv]) for fv in (fA, fB, fC, fD, fE, fF)],
              unique((Tuple(t.complex.f_vector), sort([(census_words[k], v) for (k, v) in t.complex.types])) for t in labelled))
    end
    atlas_vertices = Dict(l => length(reduce(union, [fibers[l].vertices(r) for r in fibers[l].cells])) for l in labels)
    atlas_cells = Dict{String,Vector{Tuple{String,Int}}}()
    for l in labels
        counts = Dict{String,Int}()
        foreach(r -> counts[words[r["polytope"]]] = get(counts, words[r["polytope"]], 0) + 1, fibers[l].cells)
        atlas_cells[l] = sort([(k, v) for (k, v) in counts])
    end
    by_row = Dict{Any,Vector{Int}}()
    foreach(l -> push!(get!(by_row, (atlas_vertices[l], atlas_cells[l]), Int[]), parse(Int, l)), labels)
    check("the atlas: vertices and maximal cells of each F_C, with the types of C, as in the rows of Table tab:local_flag_complex",
          sort([(fv[1], paper_cells[fv], types) for (fv, types) in paper_rows]), sort([(k..., sort(v)) for (k, v) in by_row]))
    label_of = labelled === nothing ? Dict{Int,String}() :
               Dict(x.number => replace(t.label, "\\mathfrak " => "") for t in labelled for x in t.members)
    size_of = Dict(r["orbit_type"] => r["size"] for r in grass_orbit_rows if r["kind"] == "cell")
    singular = Dict("tetrahedra" => "tetrahedron", "square pyramids" => "square pyramid", "triangular prisms" => "triangular prism")
    println("           type  fields  rays of C     orbit  3-cells  node 2-cells  f-vector (decorations)  maximal cells (atlas)                  label")
    for l in labels
        F = fibers[l]
        @printf("           %4s  %-6s  %-12s  %5s  %7d  %12d  %-22s  %-37s  %s\n", l, join(F.fields, '+'), join(F.support, ','),
                size_of[l], length(F.cells), length(F.nodes), show_value(only(decorated[parse(Int, l)])),
                join(["$v $(v == 1 ? singular[k] : k)" for (k, v) in atlas_cells[l]], ", "), get(label_of, parse(Int, l), ""))
    end

    subheading("The fiber complexes F_{A_1} over type 21 and F_D over type 5 (Figures fig:triangle_cell and fig:prism_projection)")
    A1, letters = fibers["21"], s7_LETTERS["21"]
    name_of = Dict(r["location"] => s7_name(r, letters) for r in A1.entries if r["location_kind"] in ("edge", "node"))
    check("F_{A_1} over type 21: tetrahedra (one per edge of the tree); node triangles (one per inner node); " *
          "vertices besides the 7 shapes; the letters A-F of figures/triangle_incidence.tex name exactly the vertices R_g^X of the atlas",
          (10, 4, 6, true),
          (count(r -> r["polytope"] == "tetrahedron", A1.cells), length(A1.nodes),
           length(setdiff(reduce(union, [A1.vertices(r) for r in A1.cells]), stars["qq"])),
           Set(x for r in A1.cells for x in first(s7_symbols(r))) == Set(keys(letters))))
    BDE = only(n for n in keys(A1.nodes) if name_of["node $n"] == "BDE")
    inner_tetrahedra = sort([name_of[r["location"]] for r in A1.cells if s7_edge(r) !== nothing])
    check("F_{A_1}: the tetrahedra over the inner edges of the tree, and those of them on the node triangle BDE",
          (["ABDE", "BCDE", "BDEF"], ["ABDE", "BCDE", "BDEF"]),
          (inner_tetrahedra, sort([name_of[r["location"]] for r in A1.cells if s7_edge(r) !== nothing && BDE in A1.inside[r["location"]]])))
    info("F_{A_1}, each node triangle with the tetrahedra on it: " * s7_incidence(A1, name_of))
    qq_labels = [l for l in labels if "qq" in fibers[l].fields]
    claim("every 3-cell over an edge e of the tree of C: over each vertex of C, the endpoints of e in its tree (two vertices), " *
          "or the node that e is contracted to (one)", all(l -> s7_contractions(fibers[l], splits["qq"]), qq_labels),
          "the $(length(qq_labels)) types over QQ")
    # figures/triangle_incidence.tex draws the tetrahedron ABDE at the middle of the inner edge v10-v11 (ie10).
    ABDE = only(r for r in A1.cells if r["location"] == "(:edge, 10, 11)")
    abde_edge = s7_edge(ABDE)
    check("the tetrahedron over the inner edge e = v10-v11 of figures/triangle_incidence.tex: its name, and the rays of C " *
          "whose trees keep e (e is AB there, contracted to D and E in the others)", ("ABDE", [7]),
          (name_of[ABDE["location"]], [g for g in A1.support if s7_split(A1, abde_edge...) in splits["qq"][("ray", g)]]))
    claim("the 3-cell over the pendant edge k has the shape e_k as a vertex, a 3-cell over an inner edge no shape",
          all(r -> last(s7_symbols(r)) == (s7_leaf(r) === nothing ? Int[] : [s7_leaf(r)]), [r for r in atlas if r["location_kind"] == "edge"]))
    D5, letters5 = fibers["5"], s7_LETTERS["5"]
    name5 = Dict(r["location"] => s7_name(r, letters5) for r in D5.entries if r["location_kind"] in ("edge", "node"))
    prisms = [r for r in D5.cells if r["polytope"] == "triangular prism"]
    check("F_D over type 5: its triangular prisms; the prism; its triangles that are node triangles of F_D",
          (1, "BCDEGH", ["BEH", "CDG"]),
          (length(prisms), name5[only(prisms)["location"]], sort([name5["node $n"] for n in D5.inside[only(prisms)["location"]]])))
    prism_edge = s7_edge(only(prisms))
    cherry = Set([Set([1, 2]), Set(3:7)])
    check("the prism lies over an inner edge with the cherry {1,2} (blue square, blue triangle) at one end; that cherry is in " *
          "the trees of the 3 vertices of C; the edge (v8-v10 in figures/triangle_incidence_5.tex)", (true, true, "(:edge, 8, 10)"),
          (prism_edge !== nothing && D5.inside[only(prisms)["location"]] == prism_edge && any(n -> D5.leaves[n] == [1, 2], prism_edge),
           all(g -> cherry in splits["qq"][("ray", g)], D5.support), only(prisms)["location"]))
    info("F_D, each node triangle with the cells on it: " * s7_incidence(D5, name5))

    subheading("Simplicial fiber complexes, pyramids and prisms")
    labelled === nothing || check("the simplicial fiber complexes: the f-vectors of the first two rows", [fA, fB],
                                  sort(unique(Tuple(t.complex.f_vector) for t in labelled if keys(t.complex.types) == Set(["3-simplex"]))))
    cell_shapes(l) = sort(unique(r["polytope"] for r in fibers[l].cells))
    check("fiber complexes with only tetrahedra (rows A and B of Table tab:local_flag_complex)",
          [1, 2, 6, 7, 9, 11, 12, 15, 17, 21, 23, 24, 25, 26, 27, 31, 32, 33, 34],
          sort([parse(Int, l) for l in labels if cell_shapes(l) == ["tetrahedron"]]))
    check("types with a triangular prism (rows D, F), types with a square pyramid (rows C, E, F)",
          ([5, 30], [3, 4, 8, 10, 13, 14, 16, 18, 19, 20, 22, 28, 29, 30, 35]),
          (sort([parse(Int, l) for l in labels if "triangular prism" in cell_shapes(l)]),
           sort([parse(Int, l) for l in labels if "square pyramid" in cell_shapes(l)])))

    subheading("Prime cones of trop(F) over prime cones of trop(G)")
    grass_cell_type = Dict(k => t for ((kind, k), t) in grass_types["qq"] if kind == "cell")
    gf2_type = Dict(k => t for ((kind, k), t) in grass_types["gf2"] if kind == "cell")
    flag_census = rows(g2("flag", "06_canonical_interior_degenerations.csv"))
    toric_flag = Set(parse(Int, c) for r in flag_census if r["toric"] == "true" for c in split(r["cells"], ';'))
    grass_census = rows(g2("grassmannian", "06_canonical_interior_degenerations.csv"))
    toric_grass = Set(parse(Int, c) for r in grass_census if r["toric"] == "true" for c in split(r["cells"], ';'))
    projection = Dict(parse(Int, r["flag_cell"]) => parse(Int, r["grassmannian_cell"]) for r in rows(g2("flag", "07_projection_cells.csv")))
    images = Set(projection[c] for c in toric_flag)
    check("prime flag cones, and their images (prime cones of trop(G))", [807, 135], [length(toric_flag), length(images)])
    claim("the images are exactly the prime cones of trop(G)", images == toric_grass)
    flag_orbit_types = Set(grass_cell_type[projection[c]] for c in toric_flag)
    check("the types of toric degenerations of F (orbits of prime cones), and the types of toric degenerations of G below them",
          (81, 17), (count(r -> r["toric"] == "true", flag_census), length(flag_orbit_types)))
    gf2_census = rows(g2("flag_gf2", "06_canonical_interior_degenerations.csv"))
    toric_gf2 = Set(parse(Int, c) for r in gf2_census if r["toric"] == "true" for c in split(r["cells"], ';'))
    gf2_projection = Dict(parse(Int, r["flag_cell"]) => parse(Int, r["grassmannian_cell"]) for r in rows(g2("flag_gf2", "07_projection_cells.csv")))
    grass_gf2_census = rows(g2("grassmannian_gf2", "06_canonical_interior_degenerations.csv"))
    toric_grass_gf2 = Set(parse(Int, c) for r in grass_gf2_census if r["toric"] == "true" for c in split(r["cells"], ';'))
    images_gf2 = Set(gf2_projection[c] for c in toric_gf2)
    claim("characteristic 2 (\"an analogous correspondence\"): every prime cone of trop_2(F) lies over a prime cone of trop_2(G)",
          images_gf2 ⊆ toric_grass_gf2, "they lie over $(length(images_gf2)) of the $(length(toric_grass_gf2))")
    claim("a 3-cell of a fiber complex is prime exactly when it is toric (census columns prime and toric, QQ and GF(2))",
          all(r -> r["prime"] == r["toric"], flag_census) && all(r -> r["prime"] == r["toric"], gf2_census))
    green, sizes = Dict{String,Int}(), Dict{String,Int}()
    for r in grass_census
        r["toric"] == "true" || continue
        cell = parse(Int, r["representative_cell"])
        green[grass_cell_type[cell]] = count(c -> projection[c] == cell, toric_flag)
        sizes[grass_cell_type[cell]] = parse(Int, r["orbit_size"])
    end
    green2, sizes2 = Dict{String,Int}(), Dict{String,Int}()
    for r in grass_gf2_census
        r["toric"] == "true" || continue
        cell = parse(Int, r["representative_cell"])
        green2[gf2_type[cell]] = count(c -> gf2_projection[c] == cell, toric_gf2)
        sizes2[gf2_type[cell]] = parse(Int, r["orbit_size"])
    end
    order = sort(collect(keys(green)); by=t -> parse(Int, t))
    claim("the green numbers of the types 1-17 are the same in characteristic 2", all(green2[t] == green[t] && sizes2[t] == sizes[t] for t in order))
    check("Figure fig:atlas_toric_cells_prime: (orbit size, green number) of the types 1-20 (18-20 in characteristic 2)",
          [(6, 4), (6, 4), (12, 6), (12, 6), (6, 8), (3, 6), (3, 11), (12, 2), (12, 2), (6, 9), (6, 10), (6, 10), (6, 6), (3, 6),
           (12, 6), (12, 6), (12, 8), (6, 2), (3, 0), (3, 2)],
          [(sizes2[string(t)], green2[string(t)]) for t in 1:20])
    toric_cells = Dict("qq" => toric_flag, "gf2" => toric_gf2)
    flag_cell_of = Dict(f => Dict(v => k for (k, v) in fiber[f].flag_cell_rays) for (f, _, _, _) in fields)
    green_edges = Dict{Int,Set{String}}()
    fields_agree = true
    for t in 1:20
        F = fibers[string(t)]
        by_field = [Set(r["location"] for r in F.cells if flag_cell_of[f][s7_set(r["$(f)_flag_rays"])] in toric_cells[f]) for f in F.fields]
        fields_agree &= length(unique(by_field)) == 1
        green_edges[t] = first(by_field)
    end
    claim("the green edges from the atlas: the same in both fields, as many as the green numbers of the census",
          fields_agree && all(length(green_edges[t]) == green2[string(t)] for t in 1:20))
    if drawings("Figure fig:atlas_toric_cells_prime: its green and red edges, against the prime cones of the atlas",
                "figures/toric_degeneration_prime_edge_trees.tex")
        figure = figure_text("toric_degeneration_prime_edge_trees")
        heads = collect(eachmatch(r"\\oldstylenums\{(\d+)\}\}\{\\large\((\d+)\)\}", figure))
        figure_edges = Dict{Int,String}()
        for (n, h) in enumerate(heads)
            panel = SubString(figure, h.offset, n < length(heads) ? prevind(figure, heads[n+1].offset) : lastindex(figure))
            drawn = [(d[2], minmax(parse(Int, d[3]), parse(Int, d[4])), d[1] == "pendant")
                     for d in eachmatch(r"\\draw\[g2tree (internal|pendant),draw=(g2Prime|g2Barren)\] \(v(\d+)\)--\(v(\d+)\);", panel)]
            places = colour -> [pendant ? string(a) : "$a-$b" for (c, (a, b), pendant) in sort(drawn; by=d -> (!d[3], d[2])) if c == colour]
            figure_edges[parse(Int, h[1])] = join(places("g2Prime"), " ") * " | " * join(places("g2Barren"), " ")
        end
        info("the green | red edges (leaf k, inner a-b) of figures/toric_degeneration_prime_edge_trees.tex = the prime | other 3-cells:")
        for t in 1:20
            cells = s7_sorted_cells(fibers[string(t)])
            check(@sprintf("Figure fig:atlas_toric_cells_prime, type %2d (%2d)", t, sizes2[string(t)]), get(figure_edges, t, "no panel"),
                  join([s7_place(r) for r in cells if r["location"] in green_edges[t]], " ") * " | " *
                  join([s7_place(r) for r in cells if !(r["location"] in green_edges[t])], " "))
        end
    end
    total = sum(green[t] * sizes[t] for t in order)
    # main.tex: "by adding up all 17 products, we arrive at the total count: 807".
    check("807 = Σ (orbit size) × (prime cones over one tree), all 17 products, over the types 1-17: " *
          join(["$(sizes[t])×$(green[t])" for t in order], " + "), 807, total)
    check("the displayed products 6×4, 6×4, 12×6 (types 1-3), 3×11 (type 7), 12×8 (type 17)",
          [(6, 4), (6, 4), (12, 6), (3, 11), (12, 8)], [(sizes[t], green[t]) for t in ("1", "2", "3", "7", "17")])
    info("the 20 products, types 18-20 included, give $(sum(sizes2[string(t)] * green2[string(t)] for t in 1:20)); types 18-20 occur only in characteristic 2 (primality over GF(2)): their products 6×2 + 3×0 + 3×2 = 18 give " *
         "825 = 807 + 18, the prime cones of trop_2(F); the 17 products of the types 1-17 give 807")
    over19 = count(c -> gf2_type[gf2_projection[c]] == "19", toric_gf2)
    check("characteristic 2 (\"type 19 has no prime cone in its fiber complex\"): prime cones of trop_2(F) over the cells of " *
          "type 19; the types of prime cones of trop_2(G) with none over them", (0, ["19"]),
          (over19, sort(unique(gf2_type[c] for c in setdiff(toric_grass_gf2, images_gf2)))))

    subheading("The toric degenerations of the flag variety (in words in the paper; here as a table)")
    D = study_path("g2", "degenerations_grassmannian")
    equivalent(file) = strip(last(split(first(l for l in readlines(joinpath(D, file)) if startswith(l, "QQ flag:")), "lattice equivalent to it:")))
    spacek_wang = equivalent(joinpath("Spacek_Wang", "compare_flag.txt"))
    ffl = equivalent(joinpath("Feigin_Fourier_Littelmann", "verify.txt"))
    remark = Dict(spacek_wang => "Spacek-Wang", ffl => "Feigin-Fourier-Littelmann")
    println("           field  orbit  size  type below  apple                      f-vector                    N   volume  Remark")
    for (target, field) in (("flag", "qq"), ("flag_gf2", "gf2"))
        for o in sort(toric_orbits(target, field); by=o -> (parse(Int, o.type), parse(Int, o.row["orbit"])))
            @printf("           %-5s  %5s  %4s  %-10s  %-25s  %-26s %3s  %6s  %s\n", field, o.row["orbit"], o.row["orbit_size"], o.type,
                    o.apple, replace(o.row["body_f_vector"], ';' => ','), o.row["body_lattice_points"], o.row["body_normalized_volume"],
                    field == "qq" ? get(remark, o.row["orbit"], "") : "")
        end
    end
    flag_toric = toric_orbits("flag", "qq")
    claim("every toric degeneration of F has a Newton-Okounkov polytope of dimension 6",
          all(o -> length(split(o.row["body_f_vector"], ';')) == 6, flag_toric))
    stored("Spacek-Wang flag degeneration (Spacek_Wang/compare_flag.jl)",
           strip(first(l for l in readlines(joinpath(D, "Spacek_Wang", "compare_flag.txt")) if startswith(l, "QQ flag:"))))
    flag_orbit_of = Dict(parse(Int, c) => r["orbit"] for r in flag_census for c in split(r["cells"], ';'))
    orbit_over(label, location) =
        flag_orbit_of[flag_cell_of["qq"][s7_set(only(r for r in fibers[label].cells if r["location"] == location)["qq_flag_rays"])]]
    check("Spacek-Wang: their flag orbit (compare_flag.txt) = the apples on the circle leaves 3 and 6 of type 1",
          [spacek_wang], unique([orbit_over("1", "(:ray, 3)"), orbit_over("1", "(:ray, 6)")]))
    T7 = fibers["7"]
    cherries = [n for n in keys(T7.nodes) if length(T7.leaves[n]) == 2]
    supporting = [r["location"] for r in T7.cells if let edge = s7_edge(r)
                      edge !== nothing && !isempty(intersect(edge, cherries))
                  end]
    check("Feigin-Fourier-Littelmann: the inner edges at the cherries of type 7, and their flag orbit = that of P(ω1+ω2) (verify.txt)",
          (2, [ffl]), (length(supporting), unique([orbit_over("7", l) for l in supporting])))
end
