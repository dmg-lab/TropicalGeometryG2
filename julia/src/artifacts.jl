fan_cells(fan) = begin
    incidence = maximal_cones(IncidenceMatrix, fan)
    [sort!(collect(row(incidence, index))) for index in 1:nrows(incidence)]
end

function rows_matrix(rows; ring=ZZ)
    isempty(rows) && return zero_matrix(ring, 0, 0)
    return matrix(ring, reduce(vcat, permutedims.(rows)))
end

"The four files of the symmetry artifacts with path prefix `prefix`, in this order."
symmetry_artifact_paths(prefix::AbstractString) = (
    permutations=prefix * "_taxon_permutations.mrdi",
    rays=prefix * "_ray_actions.mrdi",
    cells=prefix * "_cell_actions.mrdi",
    orbits=prefix * "_cell_orbits.csv",
)

function save_symmetry_artifacts(prefix::AbstractString, symmetry)
    mkpath(dirname(prefix))
    paths = symmetry_artifact_paths(prefix)
    save(paths.permutations, rows_matrix(symmetry.permutations))
    save(paths.rays, rows_matrix(symmetry.ray_actions))
    save(paths.cells, rows_matrix(symmetry.cell_actions))
    open(paths.orbits, "w") do io
        println(io, "orbit,size,representative,cells")
        for (index, orbit) in enumerate(symmetry.cell_orbits)
            println(io, "$index,$(length(orbit)),$(first(orbit)),\"$(join(orbit, ';'))\"")
        end
    end
    return paths
end

matrix_rows(value) = [
    [Int(value[i, j]) for j in 1:ncols(value)] for i in 1:nrows(value)
]

function load_symmetry_artifacts(prefix::AbstractString)
    paths = symmetry_artifact_paths(prefix)
    permutations = matrix_rows(load(paths.permutations))
    ray_actions = matrix_rows(load(paths.rays))
    cell_actions = matrix_rows(load(paths.cells))
    return (
        permutations=permutations,
        ray_actions=ray_actions,
        cell_actions=cell_actions,
        ray_orbits=action_orbits(ray_actions, isempty(ray_actions) ? 0 : length(first(ray_actions))),
        cell_orbits=action_orbits(cell_actions, isempty(cell_actions) ? 0 : length(first(cell_actions))),
    )
end

"Save `value` unless `path` holds it: each save gives the file a new parent identifier, hence a new hash."
function save_unless_current(path, value)
    isfile(path) && load(path) == value && return path
    save(path, value)
    return path
end

function write_tree_census(path, normalized_rays, weights)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "kind,index,pendant_lengths,internal_splits")
        for (kind, vectors) in (("ray", eachrow(normalized_rays)), ("cell", weights))
            for (index, vector) in enumerate(vectors)
                tree = tree_from_min_pluecker(collect(vector))
                pendant = join(rational_key.(pendant_lengths(tree)), ';')
                splits = join(
                    (join(split.side, "") * "|" * join(split.complement, "") * ":" * rational_key(split.length)
                     for split in internal_splits(tree)),
                    ';',
                )
                println(io, "$kind,$index,\"$pendant\",\"$splits\"")
            end
        end
    end
    return path
end

csv_escape(value) = begin
    text = string(value)
    (occursin(',', text) || occursin('"', text) || occursin('\n', text)) ?
        "\"" * replace(text, "\"" => "\"\"") * "\"" : text
end

"""
Split one CSV line, honouring double quotes: a quoted field may contain commas,
and a doubled quote in it is one quote, as `csv_escape` writes them.
"""
function csv_fields(line)
    characters = collect(line)
    fields, buffer, quoted = String[], IOBuffer(), false
    index = 1
    while index <= length(characters)
        character = characters[index]
        if character == '"'
            if quoted && index < length(characters) && characters[index + 1] == '"'
                print(buffer, '"')
                index += 1
            else
                quoted = !quoted
            end
        elseif character == ',' && !quoted
            push!(fields, String(take!(buffer)))
        else
            print(buffer, character)
        end
        index += 1
    end
    quoted && error("unterminated CSV quote in $line")
    push!(fields, String(take!(buffer)))
    return fields
end

"""
    csv_rows(path; header=nothing, comments=false)

The rows of the CSV file `path`, each a dictionary by column name. Empty lines
are skipped, and with `comments` also lines starting with #. With `header`, the
column names must be these.
"""
function csv_rows(path; header=nothing, comments::Bool=false)
    lines = filter(line -> !isempty(line) && !(comments && startswith(line, '#')), readlines(path))
    isempty(lines) && error("empty CSV file $path")
    names = csv_fields(first(lines))
    isnothing(header) || names == header || error("$path: unexpected columns $(join(names, ','))")
    return map(lines[2:end]) do line
        fields = csv_fields(line)
        length(fields) == length(names) ||
            error("$path: a line with $(length(fields)) fields, the header has $(length(names)): $line")
        Dict(zip(names, fields))
    end
end

function degeneration_object_directory(path::AbstractString)
    stem, _ = splitext(path)
    return stem * "_objects"
end

function degeneration_artifact_paths(path::AbstractString)
    directory = degeneration_object_directory(path)
    objects = isdir(directory) ? sort!(readdir(directory; join=true)) : String[]
    return vcat([String(path)], objects)
end

const DEGENERATION_COLUMNS = [
        "orbit", "orbit_size", "representative_cell", "cells", "rays",
        "canonical_weight", "initial_ideal", "torus_saturation", "degree_one_values",
        "newton_okounkov_body", "prime", "binomial", "monomial_free",
        "torus_saturated", "toric", "dimension", "degree",
        "groebner_generators", "maximum_terms", "body_f_vector",
        "body_lattice_points", "body_normalized_volume", "body_type", "body_slice_types",
        "chevalley_equivalent", "binomial_lattice_index", "geometrically_prime",
]

write_degenerations_header(io) = println(io, join(DEGENERATION_COLUMNS, ','))

function write_degenerations_record(io, path::AbstractString, record)
    mkpath(dirname(path))
    object_directory = degeneration_object_directory(path)
    mkpath(object_directory)
    algebra = record.algebra
    body = record.body
    prefix = joinpath(object_directory, @sprintf("orbit_%03d", record.orbit))
    ideal_path = prefix * "_initial_ideal.mrdi"
    saturation_path = isnothing(algebra.saturation) ? "" : prefix * "_torus_saturation.mrdi"
    save(ideal_path, algebra.ideal)
    isempty(saturation_path) || save(saturation_path, algebra.saturation)
    values_path = ""
    body_path = ""
    if !isnothing(body)
        values_path = prefix * "_degree_one_values.mrdi"
        body_path = prefix * "_newton_okounkov_body.poly"
        save(values_path, body.points)
        save(body_path, body.polytope)
        for (index, slice) in enumerate(hasproperty(body, :slices) ? body.slices : ())
            save(prefix * "_newton_okounkov_body_slice_$index.poly", slice.polytope)
        end
    end
    values = Any[
        record.orbit,
        length(record.cells),
        record.representative,
        join(record.cells, ';'),
        join(record.rays, ';'),
        join(rational_key.(record.canonical_weight), ';'),
        relpath(ideal_path, dirname(path)),
        isempty(saturation_path) ? "" : relpath(saturation_path, dirname(path)),
        isempty(values_path) ? "" : relpath(values_path, dirname(path)),
        isempty(body_path) ? "" : relpath(body_path, dirname(path)),
        algebra.prime,
        algebra.binomial,
        algebra.monomial_free,
        ismissing(algebra.torus_saturated) ? "" : algebra.torus_saturated,
        algebra.toric,
        algebra.dimension,
        ismissing(algebra.degree) ? "" : algebra.degree,
        algebra.groebner_generators,
        algebra.maximum_terms,
        isnothing(body) ? "" : join(body.f_vector, ';'),
        isnothing(body) ? "" : body.lattice_points,
        isnothing(body) ? "" : body.normalized_volume,
        !hasproperty(record, :body_type) || ismissing(record.body_type) ? "" : record.body_type,
        !hasproperty(record, :body_slice_types) || ismissing(record.body_slice_types) ?
            "" : record.body_slice_types,
        ismissing(record.chevalley_equivalent) ? "" : record.chevalley_equivalent,
        !hasproperty(algebra, :binomial_lattice_index) || ismissing(algebra.binomial_lattice_index) ?
            "" : algebra.binomial_lattice_index,
        !hasproperty(algebra, :geometrically_prime) || ismissing(algebra.geometrically_prime) ?
            "" : algebra.geometrically_prime,
    ]
    println(io, join(csv_escape.(values), ','))
    flush(io)
end

function write_degenerations(path::AbstractString, result)
    mkpath(dirname(path))
    open(path, "w") do io
        write_degenerations_header(io)
        for record in result.representatives
            write_degenerations_record(io, path, record)
        end
    end
    return degeneration_artifact_paths(path)
end
