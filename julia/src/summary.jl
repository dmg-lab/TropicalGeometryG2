# Stage `summarize`: the statements of the paper about one census target,
# computed from its stored artifacts. Modulo lineality the fan is a polyhedral
# complex: its faces of dimension k are the cones of dimension k + 1 modulo
# lineality, and a vertex is a ray. Faces and their types are computed on one
# maximal cell per orbit of the taxon group and counted by orbit sizes.

"The file of the summary of the census target stored in `directory`."
summary_path(directory::AbstractString) = joinpath(directory, "08_summary.toml")

"""
The faces of the maximal cell `cell` (a sorted list of ray indices), from the
facets of its cone alone: every proper face of a polytope is an intersection
of facets. Returns the faces as sorted ray lists, entry `k + 1` listing those
of dimension `k` modulo lineality, and the facet-by-ray incidence matrix.
"""
function cell_faces(cell, rays, lineality)
    generators = matrix(QQ, [rays[r][j] for r in cell, j in eachindex(first(rays))])
    cone = positive_hull(generators, lineality)
    values = linear_inequality_matrix(facets(cone)) * transpose(generators)
    facet_sets = unique!([[cell[j] for j in eachindex(cell) if iszero(values[i, j])] for i in 1:nrows(values)])
    filter!(!isempty, facet_sets)
    faces = Set{Vector{Int}}(vcat(facet_sets, [cell]))
    while true
        meets = Set(sort!(intersect(a, b)) for a in faces for b in faces)
        new = setdiff(filter(!isempty, meets), faces)
        isempty(new) && break
        union!(faces, new)
    end
    # The face lattice is graded: a face without smaller faces is a vertex.
    dimension = Dict{Vector{Int}, Int}()
    for face in sort!(collect(faces); by=length)
        below = [dimension[small] for small in keys(dimension) if length(small) < length(face) && issubset(small, face)]
        dimension[face] = isempty(below) ? 0 : maximum(below) + 1
    end
    top = dimension[cell]
    top == dim(cone) - lineality_dim(cone) - 1 || error("the face lattice of cell $cell has the wrong rank")
    incidence = BitMatrix([ray in facet for facet in facet_sets, ray in cell])
    return [sort!([face for face in faces if dimension[face] == k]) for k in 0:top], incidence
end

"""
Reduced integral homology of the polyhedral complex with maximal cells `cells`
(sorted vertex lists). It is that of the simplicial complex with these facets:
in both, a family of maximal cells meets exactly when their vertex sets meet,
and a nonempty intersection is contractible, so the nerve lemma applies.
"""
function reduced_integral_homology(cells)
    complex = simplicial_complex(cells)
    return map(0:dim(complex)) do k
        divisors = elementary_divisors(homology(complex, k))
        (dimension=k, rank=count(iszero, divisors) - (k == 0),
         torsion=[Int(d) for d in divisors if !iszero(d) && !isone(d)])
    end
end

function describe_reduced_homology(groups)
    parts = String[]
    for group in groups
        summands = vcat(iszero(group.rank) ? String[] : ["Z^$(group.rank)"], ["Z/$t" for t in group.torsion])
        isempty(summands) || push!(parts, "H_$(group.dimension) = " * join(summands, " + "))
    end
    return isempty(parts) ? "zero" : join(parts, ", ")
end

"Split one CSV line, honouring double quotes."
function csv_fields(line)
    fields, buffer, quoted = String[], IOBuffer(), false
    for character in line
        if character == '"'
            quoted = !quoted
        elseif character == ',' && !quoted
            push!(fields, String(take!(buffer)))
        else
            print(buffer, character)
        end
    end
    quoted && error("unterminated CSV quote")
    push!(fields, String(take!(buffer)))
    return fields
end

function csv_rows(path)
    lines = filter(!isempty, readlines(path))
    header = csv_fields(first(lines))
    return map(lines[2:end]) do line
        fields = csv_fields(line)
        length(fields) == length(header) || error("unexpected column count in $path")
        Dict(zip(header, fields))
    end
end

"The maximal cells of the orbits of the stage-6 census with `column` true."
degeneration_cells(rows, column) =
    sort!([parse(Int, cell) for row in rows if row[column] == "true" for cell in split(row["cells"], ';')])

function degeneration_counts(rows)
    counts = Dict{String, Any}("orbits" => length(rows), "cells" => sum(parse(Int, row["orbit_size"]) for row in rows))
    for column in ("toric", "prime", "binomial", "monomial_free")
        selected = filter(row -> row[column] == "true", rows)
        counts["$(column)_orbits"] = length(selected)
        counts["$(column)_cells"] = sum((parse(Int, row["orbit_size"]) for row in selected); init=0)
    end
    return counts
end

"Each maximal cell as the set of the signatures of its rays modulo lineality, by `quotient`."
function cell_signatures(fan, quotient)
    signatures = [quotient_ray_signature(ray, quotient) for ray in rays_modulo_lineality(fan).rays_modulo_lineality]
    return [Set(signatures[cell]) for cell in fan_cells(fan)]
end

"""
The maximal cells of `fan` that are maximal cells of the fan `reference_fan`,
and the toric ones among them, from the stage-6 census `rows` of each. The
two groups agree, so one representative per orbit decides for its orbit.
"""
function compare_with_reference(fan, rows, cell_orbits, reference_fan, reference_rows)
    lineality = matrix(QQ, rays_modulo_lineality(fan).lineality_basis)
    reference_lineality = matrix(QQ, rays_modulo_lineality(reference_fan).lineality_basis)
    rank(lineality) == rank(reference_lineality) == rank(vcat(lineality, reference_lineality)) ||
        error("the fan and its reference have different lineality spaces")
    quotient = Matrix{QQFieldElem}(kernel(lineality; side=:right))
    cells, reference_cells = cell_signatures(fan, quotient), cell_signatures(reference_fan, quotient)
    here = Dict(cell => index for (index, cell) in enumerate(cells))
    there = Dict(cell => index for (index, cell) in enumerate(reference_cells))
    toric, reference_toric = Set(degeneration_cells(rows, "toric")), Set(degeneration_cells(reference_rows, "toric"))
    counts = Dict{String, Any}(name * suffix => 0
        for name in ("common", "toric_in_both", "toric_here_only", "reference_toric_not_toric_here")
        for suffix in ("_cells", "_orbits"))
    add!(name, size) = (counts[name * "_cells"] += size; counts[name * "_orbits"] += 1)
    for orbit in cell_orbits
        representative = first(orbit)
        common = haskey(there, cells[representative])
        common && add!("common", length(orbit))
        representative in toric || continue
        common && there[cells[representative]] in reference_toric ?
            add!("toric_in_both", length(orbit)) : add!("toric_here_only", length(orbit))
    end
    for row in reference_rows
        row["toric"] == "true" || continue
        cell = get(here, reference_cells[parse(Int, row["representative_cell"])], 0)
        cell in toric || add!("reference_toric_not_toric_here", parse(Int, row["orbit_size"]))
    end
    return counts
end

"""
    census_summary(directory; stage1_fan, reference=nothing)

The statements of the paper about the census target stored in `directory`,
from its normalized coarsest fan, its symmetry artifacts and, when present,
its stage-6 census. Only one maximal cell per orbit is examined. `stage1_fan`
is the fan of stage 1, whose rays and maximal cones are counted. `reference`
is the directory of a target over QQ to compare the maximal cells and toric
degenerations with.
"""
function census_summary(directory; stage1_fan, reference=nothing)
    fan = load(joinpath(directory, "03_normalized_coarsest_fan.mrdi"))
    symmetry = load_symmetry_artifacts(joinpath(directory, "05_symmetry"))
    data = rays_modulo_lineality(fan)
    rays = [QQ.(collect(ray)) for ray in data.rays_modulo_lineality]
    lineality = matrix(QQ, data.lineality_basis)
    cells = fan_cells(fan)
    images(face) = unique(sort!([action[ray] for ray in face]) for action in symmetry.ray_actions)
    # Every face is a face of a maximal cell, hence in the orbit of a face of a
    # representative: the orbits of faces are those of the faces of the representatives.
    face_orbits = Dict{Vector{Int}, Tuple{Int, Int}}()
    types = Dict{Tuple{Int, String}, Vector{Int}}()
    for orbit in symmetry.cell_orbits
        faces, incidence = cell_faces(cells[first(orbit)], rays, lineality)
        for (level, list) in enumerate(faces), face in list
            orbit_images = images(face)
            get!(face_orbits, minimum(orbit_images), (level - 1, length(orbit_images)))
        end
        length(images(cells[first(orbit)])) == length(orbit) ||
            error("the ray actions and the stored orbit of cell $(first(orbit)) disagree")
        counts = get!(types, (length(faces) - 1, describe_type(determine_polytope(incidence, length(faces) - 1))), [0, 0])
        counts .+= (length(orbit), 1)
    end
    top = maximum(first, values(face_orbits))
    f = [sum(size for (k, size) in values(face_orbits) if k == dimension) for dimension in 0:top]
    orbit_f = [count(record -> first(record) == dimension, values(face_orbits)) for dimension in 0:top]
    first(f) == length(rays) || error("the faces of the representatives miss rays")
    sum(first, values(types)) == length(cells) || error("the cell types do not cover the maximal cells")
    homology = reduced_integral_homology(cells)
    # The reduced Euler characteristic, once from the faces and once from the homology of the nerve.
    sum((-1)^k * f[k + 1] for k in 0:top) - 1 == sum((-1)^g.dimension * g.rank for g in homology) ||
        error("the f-vector and the homology have different Euler characteristics")
    stage1 = load(stage1_fan)
    summary = Dict{String, Any}(
        "complex" => Dict{String, Any}(
            "dimension" => top,
            "lineality_dimension" => lineality_dim(fan),
            "pure" => allequal(first.(keys(types))),
            "f_vector" => f,
            "orbit_f_vector" => orbit_f,
            "maximal_cells" => length(cells),
            "maximal_cell_orbits" => length(symmetry.cell_orbits),
            "group_order" => length(symmetry.permutations),
        ),
        "maximal_cell_types" => [
            Dict{String, Any}("dimension" => d, "type" => type, "cells" => counts[1], "orbits" => counts[2])
            for ((d, type), counts) in sort!(collect(types); by=item -> (item[1][1], -item[2][1]))
        ],
        "reduced_homology" => Dict{String, Any}(
            "description" => describe_reduced_homology(homology),
            "ranks" => [group.rank for group in homology],
            "torsion" => [join(group.torsion, ";") for group in homology],
        ),
        "stage1_fan" => Dict{String, Any}(
            "path" => relpath(stage1_fan, directory),
            "rays" => length(rays_modulo_lineality(stage1).rays_modulo_lineality),
            "maximal_cones" => n_maximal_cones(stage1),
        ),
    )
    if top >= 2
        # Every polygon is determined by its number of vertices.
        polygons = Dict{String, Int}()
        for (face, (k, size)) in face_orbits
            k == 2 && (polygons["$(length(face))-gon"] = get(polygons, "$(length(face))-gon", 0) + size)
        end
        summary["two_dimensional_faces"] = polygons
    end
    census = joinpath(directory, "06_canonical_interior_degenerations.csv")
    if isfile(census)
        rows = csv_rows(census)
        summary["degenerations"] = degeneration_counts(rows)
        if !isnothing(reference)
            summary["comparison"] = compare_with_reference(
                fan, rows, symmetry.cell_orbits,
                load(joinpath(reference, "03_normalized_coarsest_fan.mrdi")),
                csv_rows(joinpath(reference, "06_canonical_interior_degenerations.csv")),
            )
            summary["comparison"]["reference"] = relpath(reference, directory)
        end
    end
    return summary
end

"Write `census_summary(directory; ...)` to `08_summary.toml` and print its main lines."
function write_census_summary(directory; label, kwargs...)
    summary = census_summary(directory; kwargs...)
    summary["label"] = label
    open(io -> TOML.print(io, summary; sorted=true), summary_path(directory), "w")
    complex = summary["complex"]
    println("  $label: dimension $(complex["dimension"]) modulo lineality, f-vector $(Tuple(complex["f_vector"])), ",
            "orbit f-vector $(Tuple(complex["orbit_f_vector"]))")
    println("  maximal cells: ", join(("$(t["cells"]) $(t["type"]) ($(t["orbits"]) orbits)" for t in summary["maximal_cell_types"]), ", "))
    println("  reduced integral homology: ", summary["reduced_homology"]["description"])
    if haskey(summary, "degenerations")
        d = summary["degenerations"]
        println("  toric degenerations: $(d["toric_cells"]) in $(d["toric_orbits"]) orbits")
    end
    return summary
end
