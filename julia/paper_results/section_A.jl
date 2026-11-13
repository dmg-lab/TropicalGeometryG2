# Appendix A of the paper: run by paper_results.jl, which defines the helpers used here.

import SHA
import Tar

"The Zenodo tarball: at the repository root (..), or else in julia/."
const sA_TARBALL = let name = "TropicalGeometryG2_data.tar.gz", paths = [repository(name), project_path(name)]
    paths[something(findfirst(isfile, paths), 1)]
end
const sA_PIECES = data_path("mukai_curves", "genus_10_method_1_pieces.mrdi")
const sA_IDEALS = ("pointspace", "pointspace_gf2", "grassmannian", "grassmannian_gf2", "flag", "flag_gf2")
const sA_FANS = ("dressian", "flag_dressian", "grassmannian", "grassmannian_gf2", "flag", "flag_gf2")
const sA_TORIC = (("grassmannian", "qq"), ("grassmannian_gf2", "gf2"), ("flag", "qq"), ("flag_gf2", "gf2"))
"The name in the tarball of the fan of a census target."
const sA_FAN_NAME = Dict("dressian" => "01_dressians/dressian", "flag_dressian" => "01_dressians/flag_dressian",
                         "grassmannian" => "02_tropical_varieties/grassmannian_qq",
                         "grassmannian_gf2" => "02_tropical_varieties/grassmannian_gf2",
                         "flag" => "02_tropical_varieties/flag_qq", "flag_gf2" => "02_tropical_varieties/flag_gf2")

"The Oscar version in the first bytes of an .mrdi file."
function sA_oscar_version(head)
    m = match(r"Oscar\.jl\",\"([0-9.]+)\"", head)
    return m === nothing ? "no Oscar header" : "Oscar " * m[1]
end

"""
The first `bytes` bytes and the SHA-256 digest of each file of a .tar.gz archive, by name below the top folder: the
names and sizes from the Tar standard library (not from the layout of `tar -tv`, which differs between GNU tar and
bsdtar), the contents from one stream of `tar -xzOf`, which writes the files in the same order.
"""
function sA_tarball_heads(tarball; bytes=300)
    headers = open(`gzip -dc $tarball`) do io
        list = Tar.list(io; strict=false)
        read(io)   # the padding after the end of the archive
        list
    end
    files = [(h.path, h.size) for h in headers if h.type == :file]
    buffer = Vector{UInt8}(undef, 1 << 20)
    return open(`tar -xzOf $tarball`) do io
        map(files) do (name, size)
            context = SHA.SHA256_CTX()
            head = read(io, min(size, bytes))
            length(head) == min(size, bytes) || error("$tarball ends inside $name")
            SHA.update!(context, head)
            rest = size - length(head)
            while rest > 0
                n = readbytes!(io, buffer, min(rest, length(buffer)))
                n > 0 || error("$tarball ends inside $name")
                SHA.update!(context, buffer, n)
                rest -= n
            end
            replace(name, r"^[^/]*/" => "") => (head=String(head), digest=bytes2hex(SHA.digest!(context)))
        end
    end
end

"The source of each file of the tarball but README.md, as ../production/make_zenodo_tarball.sh copies it: name => path."
function sA_sources()
    folder(target, field) = replace(target, "_gf2" => "") * "_" * field
    ideals = ["00_ideals/$(s)_$(f).mrdi" for s in ("pointspace", "linespace", "flagspace") for f in ("qq", "gf2")]
    return vcat([name => g2(t, "00_ideal.mrdi") for (name, t) in zip(ideals, sA_IDEALS)],
                [sA_FAN_NAME[t] * ".mrdi" => g2(t, "03_normalized_coarsest_fan.mrdi") for t in sA_FANS],
                ["03_no_polytopes/$(folder(t, f))/orbit_$(lpad(o.row["orbit"], 3, '0')).mrdi" =>
                     joinpath(g2(t), o.row["newton_okounkov_body"]) for (t, f) in sA_TORIC for o in toric_orbits(t, f)],
                ["04_hexagons/tropicalized_hexagon.mrdi" => g2("hexagon", "tropicalized_hexagon.mrdi"),
                 "05_curves/genus_10_method_1.mrdi" => data_path("mukai_curves", "genus_10_method_1.mrdi"),
                 "05_curves/genus_10_method_2.mrdi" => data_path("mukai_curves", "genus_10_method_2.mrdi"),
                 "05_curves/genus_10_method_1_pieces.mrdi" => sA_PIECES,
                 "grassmannian_orbits.csv" => g2("grassmannian_orbits.csv"), "flag_orbits.csv" => g2("flag_orbits.csv"),
                 "dressian_orbits.csv" => g2("dressian", "dressian_orbits.csv"),
                 "flag_dressian_orbits.csv" => g2("flag_dressian", "flag_dressian_orbits.csv")])
end

"The file names of Table tab:data, with orbit_NNN for the toric orbits NNN of the four censuses."
function sA_table_names()
    folder(target, field) = replace(target, "_gf2" => "") * "_" * field
    return vcat(["00_ideals/$(s)_$(f).mrdi" for s in ("pointspace", "linespace", "flagspace") for f in ("qq", "gf2")],
                ["01_dressians/$d.mrdi" for d in ("dressian", "flag_dressian")],
                ["02_tropical_varieties/$(s)_$(f).mrdi" for s in ("grassmannian", "flag") for f in ("qq", "gf2")],
                ["03_no_polytopes/$(folder(t, f))/orbit_$(lpad(o.row["orbit"], 3, '0')).mrdi"
                 for (t, f) in sA_TORIC for o in toric_orbits(t, f)],
                ["04_hexagons/tropicalized_hexagon.mrdi"],
                ["05_curves/genus_10_method_$m.mrdi" for m in ("1", "2", "1_pieces")],
                ["$(t)_orbits.csv" for t in ("grassmannian", "flag", "dressian", "flag_dressian")], ["README.md"])
end

"Integers as runs, such as 1-17, 21-31."
function sA_runs(numbers)
    numbers = sort(numbers)
    runs = String[]
    i = 1
    while i <= length(numbers)
        j = i
        while j < length(numbers) && numbers[j+1] == numbers[j] + 1
            j += 1
        end
        push!(runs, i == j ? string(numbers[i]) : "$(numbers[i])-$(numbers[j])")
        i = j + 1
    end
    return join(runs, ", ")
end

"The orbit types present in the lines of an orbit table: ray types, then cell types (numbers as runs)."
function sA_types(lines)
    ray = unique(String(r["orbit_type"]) for r in lines if r["kind"] == "ray")
    cell = unique(String(r["orbit_type"]) for r in lines if r["kind"] == "cell")
    numbers = [parse(Int, t) for t in cell if all(isdigit, t)]
    return "rays " * join(ray, ", ") * "; cells " * join(vcat(isempty(numbers) ? String[] : [sA_runs(numbers)],
                                                            [t for t in cell if !all(isdigit, t)]), ", ")
end

"The normalized rays of a census target (03_normalized_rays.mrdi), numbered as the rays of its fan."
function sA_normalized_rays(target)
    M = Matrix{QQFieldElem}(load(joinpath(TARGET_DIRECTORIES[target], "03_normalized_rays.mrdi")))
    return [M[i, :] for i in 1:size(M, 1)]
end

"The number of connected components of the graph on the vertices of `edges`."
function sA_components(edges)
    parent = Dict{Int,Int}()
    root(v) = (p = get!(parent, v, v); p == v ? v : (parent[v] = root(p)))
    for (a, b) in edges
        parent[root(a)] = root(b)
    end
    return length(unique(root(v) for v in collect(keys(parent))))
end

"The rational number of a polymake string such as -36387/7."
function sA_rational(s)
    parts = split(s, '/')
    return length(parts) == 1 ? QQ(parse(BigInt, s)) : QQ(parse(BigInt, parts[1]), parse(BigInt, parts[2]))
end

"""
The polyhedra of a stored vector of polyhedra (the pieces of Method 1), read from the text of the .mrdi file rather
than loaded by Oscar (which takes minutes): their number, the number of empty ones, and the vertices and directions of
the one-dimensional ones, in R^21.
"""
function sA_pieces(text)
    bytes = codeunits(text)
    objects = UnitRange{Int}[]   # the entries of the top-level array "data"
    depth, quoted, escaped, opening = 0, false, false, 0
    for i in (last(findfirst("\"data\":[", text)) + 1):length(bytes)
        b = bytes[i]
        if quoted
            escaped ? (escaped = false) : b == UInt8('\\') ? (escaped = true) : b == UInt8('"') && (quoted = false)
        elseif b == UInt8('"')
            quoted = true
        elseif b == UInt8('{') || b == UInt8('[')
            depth == 0 && (opening = i)
            depth += 1
        elseif b == UInt8('}') || b == UInt8(']')
            depth -= 1
            depth < 0 && break
            depth == 0 && push!(objects, opening:i)
        end
    end
    entries = [SubString(text, first(r), last(r)) for r in objects]
    one_dimensional = filter(e -> occursin(r"\"CONE_DIM\":2[,}]", e), entries)
    pieces = map(one_dimensional) do e
        generators = [[sA_rational(strip(x, '"')) for x in split(row, ',')]
                      for row in split(match(r"\"VERTICES\":\[\[(.*?)\]\]", e)[1], "],[")]
        (points=[g[2:end] .// g[1] for g in generators if !iszero(g[1])],
         directions=[g[2:end] for g in generators if iszero(g[1])])
    end
    return (total=length(entries), empty=count(e -> occursin("\"FEASIBLE\":false", e), entries), pieces=pieces)
end

"t with x = t d, or nothing if x is not a multiple of the nonzero vector d."
function sA_multiple(x, d)
    k = findfirst(!iszero, d)
    t = x[k] // d[k]
    return x == t .* d ? t : nothing
end

"A direction scaled so that its first nonzero coordinate is 1 or -1."
sA_direction(r) = r .// abs(r[findfirst(!iszero, r)])

"Whether a piece (points and directions) lies in a maximal cell (its vertices, and its ray if unbounded) of a curve."
function sA_in_cell(piece, cell)
    u = cell.vertices[1]
    d = cell.bounded ? cell.vertices[2] - u : cell.ray
    on(p) = (t = sA_multiple(p - u, d); t !== nothing && t >= 0 && (!cell.bounded || t <= 1))
    along(r) = !cell.bounded && (t = sA_multiple(r, d); t !== nothing && t > 0)
    return all(on, piece.points) && all(along, piece.directions)
end

"The ends of a piece or cell: its vertices, and the direction of its ray."
sA_ends(points, directions) = vcat([("point", p) for p in points], [("ray", sA_direction(r)) for r in directions])

function sA_files()
    subheading("The 222 files of the Zenodo data, and its copy of data/")
    polytopes = [summary(t)["degenerations"]["toric_orbits"] for (t, _) in sA_TORIC]
    files = 6 + 4 + 2 + 4 + sum(polytopes) + 1 + 3 + 1
    check("files of Table tab:data: 6 ideals, 4 orbit tables, 2 Dressians, 4 tropical varieties, the polytopes, the hexagon, 3 curve files, README",
          222, files)
    sources = sA_sources()
    heads = nothing
    if isfile(sA_TARBALL)
        all_heads = sA_tarball_heads(sA_TARBALL)
        # The top level holds the files of Table tab:data; its data/ a copy of data/
        # (../production/make_zenodo_tarball.sh).
        heads = [(name, h) for (name, h) in all_heads if !startswith(name, "data/")]
        copies = Dict(name[length("data/")+1:end] => h.digest for (name, h) in all_heads if startswith(name, "data/"))
        names = first.(heads)
        check("files at the top level of TropicalGeometryG2_data.tar.gz (../production/make_zenodo_tarball.sh)", 222,
              length(names))
        expected = sA_table_names()
        claim("the files of the tarball are named as in Table tab:data, with orbit_NNN for the toric orbits NNN of the censuses",
              sort(expected) == sort(names),
              sort(expected) == sort(names) ? "" : "missing $(setdiff(expected, names)), extra $(setdiff(names, expected))")
        folder(name) = startswith(name, "03_") ? join(split(name, '/')[1:2], '/') : occursin('/', name) ? first(split(name, '/')) : "top level"
        counts = Dict{String,Int}()
        foreach(name -> counts[folder(name)] = get(counts, folder(name), 0) + 1, names)
        info("files per folder: " * join(("$f $(counts[f])" for f in sort(collect(keys(counts)))), ", ") *
             "; and $(length(copies)) files in data/")
        # The numbers of this appendix are computed from data/: the files of the table must be copies of its files.
        source = Dict(sources)
        data_name(path) = relpath(path, data_path())
        mismatched = [name for (name, h) in heads if haskey(source, name) && get(copies, data_name(source[name]), "") != h.digest]
        claim("each file of the tarball but README.md is a byte-identical (SHA-256) copy of its source in the folder data/ of the tarball",
              isempty(mismatched), "$(count(name -> haskey(source, name), names)) compared" *
              (isempty(mismatched) ? "" : "; different or absent from data/: " * join(mismatched, ", ")))
        # A Markdown or TOML file with private blocks is left out: the tarball's copy lacks their lines.
        private_blocks(path) = (endswith(path, ".md") || endswith(path, ".toml")) &&
                               any(l -> strip(l) in ("<!-- private -->", "# private"), eachline(path))
        present = [(name, data_path(name)) for name in keys(copies) if isfile(data_path(name))]
        local_files = [(name, path) for (name, path) in present if !private_blocks(path)]
        differ = sort!([name for (name, path) in local_files if bytes2hex(open(SHA.sha256, path)) != copies[name]])
        kinds = Dict{String,Int}()
        foreach(name -> kinds[splitext(name)[2]] = get(kinds, splitext(name)[2], 0) + 1, differ)
        info("data/ of the tarball against data/: $(length(local_files) - length(differ)) of the " *
             "$(length(local_files)) files identical, $(length(differ)) different, " *
             "$(length(copies) - length(present)) absent from data/" *
             (length(present) == length(local_files) ? "" : ", $(length(present) - length(local_files)) with private blocks") *
             (isempty(differ) ? "" : length(differ) <= 5 ? ": " * join(differ, ", ") :
                                     ", by extension " * join(("$e $n" for (e, n) in sort!(collect(kinds))), ", ") *
                                     " (a computed data/ has new identifiers in its .mrdi files and new times and " *
                                     "hashes in its provenance)"))
    else
        skipped("files in TropicalGeometryG2_data.tar.gz", "not at the repository root (..) nor in julia/: " *
                "download it from Zenodo, DOI 10.5281/zenodo.24104003, to the repository root")
    end
    # main.tex: "All files except the four tables *_orbits.csv are written by OSCAR 1.8.2 in its format .mrdi".
    if heads === nothing   # the sources of the tarball in data/
        paths = [path for (name, path) in sources if !endswith(name, "_orbits.csv")]
        versions = [relpath(p, project_path()) => sA_oscar_version(String(read(p, 300))) for p in paths]
        scope = "of the data (no tarball)"
    else
        versions = [name => sA_oscar_version(h.head) for (name, h) in heads if !endswith(name, "_orbits.csv")]
        scope = "of the tarball"
    end
    # Left as it is by the authors: a warning, not a failure.
    warned("files $scope other than the four *_orbits.csv tables that are not .mrdi files of Oscar 1.8.2 (the paper: none, " *
          "of $(length(versions)))", String[],["$name ($version)" for (name, version) in versions if version != "Oscar 1.8.2"])
end

function sA_conventions()
    subheading("Conventions")
    xs, ps = ["x$i" for i in 1:7], [pair_name(k) for k in 1:21]
    rings = [base_ring(load(g2(t, "00_ideal.mrdi"))) for t in sA_IDEALS]
    check("the suffixes qq and gf2 name the field of coefficients: characteristics of the coefficient fields of the rings " *
          "of Q, G, F over QQ and GF(2)", [0, 2, 0, 2, 0, 2], [Int(characteristic(coefficient_ring(R))) for R in rings])
    claim("the rings of Q, G, F over QQ and GF(2) have the coordinates x1..x7, p12..p67 in lexicographic order, x1..x7,p12..p67",
          [string.(symbols(R)) for R in rings] == [xs, xs, ps, ps, vcat(xs, ps), vcat(xs, ps)])
    check("ambient dimensions of the fans (Dressian, flag Dressian, trop(G) over QQ, GF(2), trop(F) over QQ, GF(2))",
          [21, 28, 21, 21, 28, 28], [ambient_dim(fan(t)) for t in sA_FANS])
    claim("the six fans are PolyhedralFans", all(t -> fan(t) isa PolyhedralFan, sA_FANS))
end

function sA_table()
    subheading("Table tab:data: the contents of the files")
    ideals = Dict(t => load(g2(t, "00_ideal.mrdi")) for t in sA_IDEALS)
    check("generators of the ideals of Q, G, F over QQ and GF(2)", [1, 1, 42, 42, 78, 78], [ngens(ideals[t]) for t in sA_IDEALS])
    degree_counts(t) = (d = total_degree.(gens(ideals[t])); [count(==(1), d), count(==(2), d)])
    check("linear forms and quadrics among them: Q one quadric, G 7 linear forms and 35 quadrics (QQ, GF(2))",
          [[0, 1], [0, 1], [7, 35], [7, 35]], [degree_counts(t) for t in sA_IDEALS[1:4]])
    for t in ("flag", "flag_gf2")
        bidegree(f) = (e = exponent_vector(f, 1); (sum(e[1:7]), sum(e[8:28])))
        degrees = bidegree.(gens(ideals[t]))
        info("bidegrees in (x, p) of the generators of F over $(t == "flag" ? "QQ" : "GF(2)"): " *
             join(("$d x$(count(==(d), degrees))" for d in sort(unique(degrees))), ", "))
    end
    rays_cells(target) = [f_vector_of(target)[1], summary(target)["complex"]["maximal_cells"]]
    check("G2-Dressian: rays, maximal cones", [51, 239], rays_cells("dressian"))
    check("G2 flag Dressian: rays, maximal cones", [141, 2627], rays_cells("flag_dressian"))
    check("trop(G) over QQ, GF(2): rays, maximal cones", [51, 248, 54, 254], vcat(rays_cells("grassmannian"), rays_cells("grassmannian_gf2")))
    check("trop(F) over QQ, GF(2): rays, maximal cones", [141, 2714, 156, 2792], vcat(rays_cells("flag"), rays_cells("flag_gf2")))
    check("Newton-Okounkov polytopes: G over QQ, GF(2); F over QQ, GF(2)", [17, 20, 81, 83],
          [summary(t)["degenerations"]["toric_orbits"] for (t, _) in sA_TORIC])
end

function sA_orbit_tables()
    subheading("Orbit types: the four tables *_orbits.csv")
    symmetry = Dict(t => load_symmetry_artifacts(joinpath(TARGET_DIRECTORIES[t], "05_symmetry")) for t in sA_FANS)
    groups = map(sA_FANS) do t
        s = symmetry[t]
        S = symmetric_group(length(first(s.permutations)))
        G, _ = sub(S, [perm(S, p) for p in s.permutations])
        cells = fan_cells(fan(t))   # maximal_cones, with the rays numbered as in rays_modulo_lineality
        order(G) == 12 && is_isomorphic(G, dihedral_group(12)) && length(s.ray_actions) == length(s.cell_actions) == 12 &&
            all(sort(s.ray_actions[g][cells[c]]) == cells[s.cell_actions[g][c]] for g in 1:12, c in eachindex(cells))
    end
    claim("the census numbers the rays and maximal cones of each of the six fans as rays_modulo_lineality and maximal_cones, " *
          "and D12 (dihedral of order 12) permutes them", all(groups))
    # The actions are those of the permutations of the taxa on the coordinates (x1..x7 first in R^28); the orbits of
    # load_symmetry_artifacts are those of the actions.
    geometric = map(sA_FANS) do t
        try
            data = TropicalGeometryG2.projection_fan_data(fan(t), load(joinpath(TARGET_DIRECTORIES[t], "03_normalized_rays.mrdi")))
            TropicalGeometryG2.verify_projection_symmetry(data, symmetry[t], p -> tropical_coordinate_permutation(p, ambient_dim(fan(t))))
        catch error
            info("$t: " * sprint(showerror, error))
            false
        end
    end
    claim("each of the 12 permutations of the taxa, acting on the coordinates, maps the normalized rays and maximal cones " *
          "of each of the six fans as the stored actions say (modulo lineality)", all(geometric))

    ray_labels = Set(vcat(["alpha_$k" for k in 1:5], ["beta_1", "beta_2", "gamma_1", "gamma_2", "gamma_3", "delta", "epsilon", "star"]))
    cell_labels = Set(vcat(string.(1:35), ["epsilon"]))
    for (file, target, prefix) in (("grassmannian_orbits.csv", "grassmannian", "qq_"), ("grassmannian_orbits.csv", "grassmannian_gf2", "gf2_"),
                                   ("flag_orbits.csv", "flag", "qq_"), ("flag_orbits.csv", "flag_gf2", "gf2_"),
                                   (joinpath("dressian", "dressian_orbits.csv"), "dressian", ""),
                                   (joinpath("flag_dressian", "flag_dressian_orbits.csv"), "flag_dressian", ""))
        table = filter(r -> !isempty(r[prefix * "members"]), rows(g2(file)))
        members(r) = sort(parse.(Int, split(r[prefix * "members"], ';')))
        of_kind(kind) = [r for r in table if r["kind"] == kind]
        s = symmetry[target]
        ray_ok = startswith(target, "flag") ? ray_labels : setdiff(ray_labels, ["star"])   # star: over the lineality space
        ok = sort(members.(of_kind("ray"))) == s.ray_orbits && sort(members.(of_kind("cell"))) == s.cell_orbits &&
             all(r -> length(members(r)) == parse(Int, r["size"]), table) &&
             all(r -> parse(Int, r[prefix * "orbit"]) == minimum(members(r)), of_kind("ray")) &&
             all(r -> s.cell_orbits[parse(Int, r[prefix * "orbit"])] == members(r), of_kind("cell")) &&
             all(r -> r["orbit_type"] in ray_ok, of_kind("ray")) && all(r -> r["orbit_type"] in cell_labels, of_kind("cell"))
        claim("$(basename(file)) for $(sA_FAN_NAME[target]): the lines are the D12 orbits of its rays and maximal cones " *
              "(members, size, orbit number; types among the labels alpha_1..epsilon" * (startswith(target, "flag") ? ", star" : "") *
              ", 1..35, epsilon)", ok, "$(length(of_kind("ray"))) ray orbits, $(length(of_kind("cell"))) cell orbits")
    end
    for file in ("grassmannian_orbits.csv", "flag_orbits.csv")
        table = rows(g2(file))
        lines_over(f) = count(r -> !isempty(r["$(f)_members"]), table)
        both = count(r -> !isempty(r["qq_members"]) && !isempty(r["gf2_members"]), table)
        info("$file: $(length(table)) lines, $both over both fields, $(lines_over("qq") - both) over QQ only, " *
             "$(lines_over("gf2") - both) over GF(2) only; " * sA_types(table))
    end
    for file in (joinpath("dressian", "dressian_orbits.csv"), joinpath("flag_dressian", "flag_dressian_orbits.csv"))
        table = rows(g2(file))
        info("$(basename(file)): $(length(table)) lines; " * sA_types(table))
    end

    for (file, a, b) in (("grassmannian_orbits.csv", "grassmannian", "grassmannian_gf2"), ("flag_orbits.csv", "flag", "flag_gf2"))
        normalized = Dict("qq" => (sA_normalized_rays(a), fan_cells(fan(a))), "gf2" => (sA_normalized_rays(b), fan_cells(fan(b))))
        ids(r, f) = parse.(Int, split(r["$(f)_members"], ';'))
        function content(r, f)   # the orbit as a set of normalized rays, or of cones as sets of normalized rays
            R, C = normalized[f]
            return (r["kind"], r["kind"] == "ray" ? Set(R[ids(r, f)]) : Set(Set(R[C[c]]) for c in ids(r, f)))
        end
        table = rows(g2(file))
        on_line(r, f) = !isempty(r["$(f)_members"])
        shared = [r for r in table if on_line(r, "qq") && on_line(r, "gf2")]
        only_qq = Set(content(r, "qq") for r in table if on_line(r, "qq") && !on_line(r, "gf2"))
        only_gf2 = Set(content(r, "gf2") for r in table if !on_line(r, "qq") && on_line(r, "gf2"))
        claim("$file: an orbit over QQ and an orbit over GF(2) share a line exactly when they consist of the same rays or cones " *
              "(normalized rays)", all(r -> content(r, "qq") == content(r, "gf2"), shared) && isempty(intersect(only_qq, only_gf2)),
              "$(length(shared)) shared lines")
    end

    # The types as in Figures fig:atlas_rays, fig:atlas_toric_cells_prime and fig:atlas_toric_cells_nontoric: the rays and
    # maximal cones of the G2-Dressian other than the tetrahedra are those of trop(G).
    D, FD = fan("dressian"), fan("flag_dressian")
    dressian_file, flag_dressian_file = joinpath("dressian", "dressian_orbits.csv"), joinpath("flag_dressian", "flag_dressian_orbits.csv")
    types_of(file) = Dict((r["kind"], parse(Int, m)) => r["orbit_type"] for r in rows(g2(file)) for m in split(r["members"], ';'))
    gtype = orbit_types("grassmannian_orbits.csv", "qq")
    RD, RG = sA_normalized_rays("dressian"), sA_normalized_rays("grassmannian")
    CD, CG = fan_cells(D), fan_cells(fan("grassmannian"))
    to_G = [something(findfirst(==(r), RG), 0) for r in RD]   # the ray of trop(G) equal to a ray of D
    gcell = Dict(c => k for (k, c) in enumerate(CG))
    type_below(rays_of_G) = get(gtype, ("cell", get(gcell, sort(unique(rays_of_G)), 0)), "none")
    dtype = types_of(dressian_file)
    claim("dressian_orbits.csv: each ray, and each maximal cone not of type epsilon, is a ray or maximal cone of trop(G) " *
          "over QQ of the same type (normalized rays)",
          all(>(0), to_G) && allunique(to_G) && all(k -> dtype[("ray", k)] == gtype[("ray", to_G[k])], eachindex(RD)) &&
          all(c -> dtype[("cell", c)] == "epsilon" || dtype[("cell", c)] == type_below(to_G[CD[c]]), eachindex(CD)))

    for (target, field) in (("flag", "qq"), ("flag_gf2", "gf2"))
        flag_types, grass_types = orbit_types("flag_orbits.csv", field), orbit_types("grassmannian_orbits.csv", field)
        certificate = TOML.parsefile(g2(target, "07_projection_certificate.toml"))
        images = [(parse(Int, r["flag_cell"]), parse(Int, r["grassmannian_cell"])) for r in rows(g2(target, "07_projection_cells.csv"))]
        if target == "flag_gf2"   # Section 3 verifies the certificate over QQ
            claim("flag_gf2/07_projection_certificate.toml is valid (ray images, cells, D12-transport), and " *
                  "07_projection_cells.csv lists its images",
                  projection_certificate_valid(target) && images == [(c["source_cell"], c["target_cell"]) for c in certificate["cells"]])
        end
        rays_ok = all(r -> flag_types[("ray", r["source_ray"])] ==
                           (r["target_ray"] == 0 ? "star" : grass_types[("ray", r["target_ray"])]), certificate["rays"])
        cells_ok = all(((f, g),) -> flag_types[("cell", f)] == grass_types[("cell", g)], images)
        stars = count(r -> r["target_ray"] == 0, certificate["rays"])
        claim("$(sA_FAN_NAME[target]): each ray and maximal cone has the type of the ray or cone of trop(G) below it, " *
              "star over the lineality space", rays_ok && cells_ok, "$stars rays of type star")
    end
    # The flag Dressian over the Dressian: the image of the p-part of a ray modulo the lineality space of D is a ray of D
    # (then of trop(G)), or 0.
    Dprojection = TropicalGeometryG2.projection_fan_data(D, load(g2("dressian", "03_normalized_rays.mrdi")))
    Dray = Dict(s => k for (k, s) in enumerate(Dprojection.signatures))
    FDnormal = Matrix{QQFieldElem}(load(g2("flag_dressian", "03_normalized_rays.mrdi")))
    below = map(1:size(FDnormal, 1)) do k   # the ray of trop(G) below, 0 over the lineality space, -1 if none
        p = FDnormal[k, 8:28]
        all(iszero, permutedims(p) * Dprojection.quotient) && return 0
        j = get(Dray, TropicalGeometryG2.quotient_ray_signature(p, Dprojection.quotient), 0)
        return j > 0 ? to_G[j] : -1
    end
    fdtype = types_of(flag_dressian_file)
    CFD = fan_cells(FD)
    claim("01_dressians/flag_dressian: each ray, and each maximal cone not of type epsilon, has the type of the ray or cone " *
          "of trop(G) below it, star over the lineality space",
          all(>=(0), below) && all(k -> fdtype[("ray", k)] == (below[k] == 0 ? "star" : gtype[("ray", below[k])]), eachindex(below)) &&
          all(c -> fdtype[("cell", c)] == "epsilon" || fdtype[("cell", c)] == type_below(filter(>(0), below[CFD[c]])), eachindex(CFD)),
          "$(count(==(0), below)) rays of type star")

    epsilon_lines(file) = [r for r in rows(g2(file)) if r["kind"] == "cell" && r["orbit_type"] == "epsilon"]
    epsilon_cells(file) = sort([parse(Int, m) for r in epsilon_lines(file) for m in split(r["members"], ';')])
    dressian_cones = collect(maximal_cones(D))
    tetrahedra = [i for (i, c) in enumerate(fan_cells(D)) if length(c) == 4 && dim(dressian_cones[i]) - lineality_dim(D) == 4]
    claim("dressian_orbits.csv: the maximal cones of type epsilon are the tetrahedra of the G2-Dressian",
          epsilon_cells(dressian_file) == tetrahedra, "$(length(tetrahedra)) cones")
    flag_rays, _ = fan_rays(FD)
    _, dressian_lineality = fan_rays(D)
    over = [i for (i, c) in enumerate(fan_cells(FD))
            if rank_of(vcat([r[8:28] for r in flag_rays[c]], dressian_lineality)) - length(dressian_lineality) == 4]
    claim("flag_dressian_orbits.csv: the maximal cones of type epsilon are those over the tetrahedra (an image in R^21 of dimension 4)",
          epsilon_cells(flag_dressian_file) == over, "$(length(over)) cones")
    for (file, F) in ((dressian_file, D), (flag_dressian_file, FD))
        ray_type = Dict(parse(Int, m) => r["orbit_type"] for r in rows(g2(file)) if r["kind"] == "ray" for m in split(r["members"], ';'))
        for r in epsilon_lines(file)
            cell = parse(Int, first(split(r["members"], ';')))
            info("$(basename(file)) epsilon orbit $(r["orbit"]) of size $(r["size"]): cell $cell with rays " *
                 join(("$k ($(ray_type[k]))" for k in fan_cells(F)[cell]), ", "))
        end
    end

    for file in ("flag_orbits.csv", joinpath("flag_dressian", "flag_dressian_orbits.csv"))
        table = rows(g2(file))
        apples = [a for r in table if r["kind"] == "cell" && !isempty(r["apple"]) for a in split(r["apple"], '|')]
        valid(a) = occursin(r"^pendant [1-7]$", a) || occursin(r"^inner (8|9|10|11|12)-(8|9|10|11|12)$", a)
        claim("$(basename(file)): the apple entries are 'pendant k' (k = 1..7) or 'inner a-b' (a, b = 8..12), separated by |",
              all(valid, apples), "$(length(apples)) entries")
        without = unique(r["orbit_type"] for r in table if r["kind"] == "cell" && isempty(r["apple"]))
        info("apple empty on $(count(r -> r["kind"] == "ray" && isempty(r["apple"]), table)) of " *
             "$(count(r -> r["kind"] == "ray", table)) ray lines; cell lines without apple: " *
             (isempty(without) ? "none" : join(without, ", ")))
    end

    line = "ray,alpha_2,3,8,1,8;33;34,1;19;22"
    claim("the example line of grassmannian_orbits.csv", line in readlines(g2("grassmannian_orbits.csv")), line)
    qq, gf2 = sA_normalized_rays("grassmannian"), sA_normalized_rays("grassmannian_gf2")
    claim("the rays 8, 33, 34 of grassmannian_qq are the rays 1, 19, 22 of grassmannian_gf2 (normalized rays)",
          Set(qq[[8, 33, 34]]) == Set(gf2[[1, 19, 22]]))
    for i in [8, 33, 34]
        info("ray $i over QQ = ray $(something(findfirst(==(qq[i]), gf2), "none")) over GF(2): ($(join(string.(qq[i]), ", ")))")
    end
end

function sA_polytopes()
    subheading("Newton-Okounkov polytopes")
    for (target, field) in sA_TORIC
        file = startswith(target, "grassmannian") ? "grassmannian_orbits.csv" : "flag_orbits.csv"
        table = Dict(r["$(field)_orbit"] => r for r in rows(g2(file)) if r["kind"] == "cell" && !isempty(r["$(field)_orbit"]))
        census = rows(g2(target, "06_canonical_interior_degenerations.csv"))
        cells(s) = sort(parse.(Int, split(s, ';')))
        claim("$(sA_FAN_NAME[target]): the census orbit NNN of orbit_NNN has the cells of the line with $(field)_orbit = NNN",
              length(census) == length(table) && all(r -> cells(r["cells"]) == cells(table[r["orbit"]]["$(field)_members"]), census),
              "$(count(r -> r["toric"] == "true", census)) toric of $(length(census)) orbits")
    end
    for (target, field) in sA_TORIC
        n = startswith(target, "grassmannian") ? 5 : 6
        bodies = [body(target, o.row) for o in toric_orbits(target, field)]
        claim("$(sA_FAN_NAME[target]): the $(length(bodies)) polytopes are lattice polytopes of dimension $n in R^$n",
              all(P -> ambient_dim(P) == n && dim(P) == n && is_lattice_polytope(P), bodies))
    end
    for (target, field) in sA_TORIC[3:4]
        minkowski = all(toric_orbits(target, field)) do o
            V = Matrix{QQFieldElem}(load(joinpath(g2(target), o.row["degree_one_values"])))
            value(i, k) = V[i+7*(k-1), :]   # the degree-one value of x_i p_k
            all(value(i, k) - value(1, k) == value(i, 1) - value(1, 1) for i in 1:7, k in 1:21) && body(target, o.row) == convex_hull(V)
        end
        claim("$(sA_FAN_NAME[target]): P_C is the convex hull of the 7 x 21 sums of the values of x_i and p_k, i.e. conv(x) + conv(p)",
              minkowski)
    end
end

function sA_hexagon_and_curves()
    subheading("Hexagon and curves")
    hexagon = load(g2("hexagon", "tropicalized_hexagon.mrdi"))
    claim("tropicalized_hexagon: fields points (6 vectors in QQ^7) and lines (6 vectors in QQ^21), each with minimum 0",
          keys(hexagon) == (:points, :lines) && length(hexagon.points) == 6 && length(hexagon.lines) == 6 &&
          all(v -> length(v) == 7 && iszero(minimum(v)), hexagon.points) && all(v -> length(v) == 21 && iszero(minimum(v)), hexagon.lines))
    maxconv(v) = Int.(maximum(v) .- v)
    claim("max(v)1 - v gives the a^(1), ..., a^(6) and d^(12), ..., d^(61) displayed in Section 9",
          [maxconv(v) for v in hexagon.points] == HEXAGON_POINTS && [maxconv(v) for v in hexagon.lines] == HEXAGON_LINES)
    for (name, vectors) in (("point", hexagon.points), ("line", hexagon.lines)), (k, v) in enumerate(vectors)
        info("stored $name $k: ($(join(string.(v), ", ")))")
    end

    method_1 = nothing
    for (file, expected) in (("genus_10_method_1.mrdi", [449, 239, 201, 10]), ("genus_10_method_2.mrdi", [447, 210, 228, 10]))
        curve = load(data_path("mukai_curves", file))
        C = curve.complex
        file == "genus_10_method_1.mrdi" && (method_1 = C)
        graph = curve_graph(C)
        used = unique(vcat(first.(graph.edges), last.(graph.edges)))
        genus = length(graph.edges) - length(used) + sA_components(graph.edges)
        check("$file: maximal cells, unbounded ones, vertices of the graph of the bounded ones, its genus", expected,
              [n_maximal_polyhedra(C), length(graph.rays), length(used), genus])
        directions = vcat([p - graph.points[1] for p in graph.points[2:end]], [r for (_, r) in graph.rays])
        claim("$file: a one-dimensional PolyhedralComplex in R^21 with one point over each point of R^21/R1",
              C isa PolyhedralComplex && dim(C) == 1 && ambient_dim(C) == 21 && lineality_dim(C) == 0 &&
              rank_of(directions) < rank_of(vcat(directions, [fill(QQ(1), 21)])))
        if haskey(curve, :multiplicities)
            check("$file: fields; one multiplicity per maximal cell; their values", ((:complex, :multiplicities), n_maximal_polyhedra(C), [1, 2, 4, 5]),
                  (keys(curve), length(curve.multiplicities), sort(unique(Int.(curve.multiplicities)))))
        else
            check("$file: fields", (:complex,), keys(curve))
        end
        valence = zeros(Int, length(graph.points))
        foreach(e -> valence[collect(e)] .+= 1, graph.edges)
        foreach(r -> valence[first(r)] += 1, graph.rays)
        sums, p12 = unique(sum.(graph.points)), unique(first.(graph.points))
        info("$(length(graph.points)) vertices, $(length(graph.edges)) bounded cells, " *
             "$(length(unique(sA_direction(r) for (_, r) in graph.rays))) ray directions, valences " *
             join(("$v x$(count(==(v), valence))" for v in sort(unique(valence))), ", ") * "; " *
             (length(sums) == 1 ? "every vertex has coordinate sum $(only(sums))" :
              length(p12) == 1 ? "every vertex has p12 = $(only(p12))" : "no common normalization"))
        if haskey(curve, :multiplicities)
            m = Int.(curve.multiplicities)
            unbounded = [m[k] for (k, P) in enumerate(maximal_polyhedra(C)) if !is_bounded(P)]
            tally(list) = join(("$v x$(count(==(v), list))" for v in sort(unique(list))), ", ")
            info("multiplicities $(tally(m)); on the unbounded cells $(tally(unbounded))")
        end
    end

    # The pieces, from data/ or else from the tarball.
    pieces_text = isfile(sA_PIECES) ? read(sA_PIECES, String) :
                  isfile(sA_TARBALL) ? read(`tar -xzOf $sA_TARBALL TropicalGeometryG2_data/05_curves/genus_10_method_1_pieces.mrdi`, String) :
                  nothing
    if pieces_text !== nothing
        stored_pieces = sA_pieces(pieces_text)
        pieces_text = nothing
        pieces = stored_pieces.pieces
        cells = map(collect(maximal_polyhedra(method_1))) do P
            (bounded=is_bounded(P), vertices=[QQ.(collect(v)) for v in vertices(P)],
             ray=is_bounded(P) ? QQFieldElem[] : QQ.(collect(first(Oscar.rays(P)))))
        end
        containing = [findall(cell -> sA_in_cell(piece, cell), cells) for piece in pieces]
        inside = [findall(==([c]), containing) for c in eachindex(cells)]
        made_up = all(eachindex(cells)) do c
            ends = Dict{Any,Int}()
            for k in inside[c], e in sA_ends(pieces[k].points, pieces[k].directions)
                ends[e] = get(ends, e, 0) + 1
            end
            cell = cells[c]
            all(<=(2), values(ends)) &&
                Set(e for (e, n) in ends if n == 1) == Set(sA_ends(cell.vertices, cell.bounded ? [] : [cell.ray]))
        end
        check("genus_10_method_1_pieces: one-dimensional polyhedra, maximal cells of genus_10_method_1, cells of two collinear pieces",
              [454, 449, 5], [length(pieces), length(cells), count(k -> length(k) == 2, inside)])
        claim("each one-dimensional piece lies in exactly one maximal cell, and the pieces in a cell make it up",
              all(c -> length(c) == 1, containing) && made_up)
        info("$(isfile(sA_PIECES) ? relpath(sA_PIECES, project_path()) : "from the tarball"): " *
             "$(stored_pieces.total) polyhedra, $(stored_pieces.empty) of them empty; the one-dimensional ones: " *
             "$(count(p -> isempty(p.directions), pieces)) bounded, $(count(p -> !isempty(p.directions), pieces)) unbounded; " *
             "cells of two pieces: " * join(findall(k -> length(k) == 2, inside), ", "))
    else
        skipped("the pieces of Method 1", "data/mukai_curves/genus_10_method_1_pieces.mrdi is missing: run " *
                "make_paper_data.jl or unpack the Zenodo tarball")
    end
end

function section_A()
    sA_files()
    sA_conventions()
    sA_table()
    sA_orbit_tables()
    sA_polytopes()
    sA_hexagon_and_curves()
end
