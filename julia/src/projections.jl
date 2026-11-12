"Exact witnesses for projections whose nonvertical rays map to target rays."
function projection_lineality(fan)
    return Matrix{QQFieldElem}(matrix(QQ, rays_modulo_lineality(fan).lineality_basis))
end

function projection_rational(value)
    parts = split(string(value), '/')
    length(parts) in (1, 2) || error("invalid rational witness")
    return QQ(parse(BigInt, parts[1]), length(parts) == 1 ? BigInt(1) : parse(BigInt, parts[2]))
end
projection_rationals(values) = projection_rational.(values)
projection_strings(values) = rational_key.(collect(values))

function projection_lineality_coefficients(rows, lineality)
    if size(lineality, 1) == 0
        all(iszero, rows) || error("projection residual is outside zero target lineality")
        return Matrix{QQFieldElem}(undef, size(rows, 1), 0)
    end
    solvable, coefficients = can_solve_with_solution(
        matrix(QQ, lineality), matrix(QQ, rows); side=:left,
    )
    solvable || error("projection residual is outside target lineality")
    return Matrix{QQFieldElem}(coefficients)
end

function projection_fan_data(fan, normalized_rays)
    rays = Matrix{QQFieldElem}(normalized_rays)
    data = rays_modulo_lineality(fan)
    raw = Matrix{QQFieldElem}(matrix(QQ, data.rays_modulo_lineality))
    size(rays) == size(raw) || error("projection ray matrix has the wrong dimensions")
    lineality = projection_lineality(fan)
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    signatures = [quotient_ray_signature(ray, quotient) for ray in eachrow(rays)]
    length(unique(signatures)) == length(signatures) || error("duplicate target or source quotient rays")
    for index in axes(rays, 1)
        signatures[index] == quotient_ray_signature(raw[index, :], quotient) ||
            error("normalized ray $index disagrees with its fan ray")
    end
    return (; rays, lineality, quotient, signatures, cells=fan_cells(fan))
end

"Verify the stored ray and cell actions against the actual rational fan geometry."
function verify_projection_symmetry(data, symmetry, coordinate_indices)
    permutations = symmetry.permutations
    isempty(permutations) && error("projection symmetry group is empty")
    keys = Set(Tuple.(permutations))
    length(keys) == length(permutations) || error("duplicate symmetry permutations")
    ntaxa = length(first(permutations))
    all(p -> sort(p) == collect(1:ntaxa), permutations) || error("invalid taxon permutation")
    Tuple(1:ntaxa) in keys || error("projection symmetry lacks the identity")
    for left in permutations, right in permutations
        Tuple(left[right]) in keys || error("projection permutations are not a group")
    end
    length(symmetry.ray_actions) == length(symmetry.cell_actions) == length(permutations) ||
        error("projection symmetry action counts disagree")
    cell_index = Dict(Tuple(cell) => index for (index, cell) in enumerate(data.cells))
    for (group, permutation) in enumerate(permutations)
        indices = coordinate_indices(permutation)
        sort(indices) == collect(axes(data.rays, 2)) || error("invalid coordinate action")
        preserves_lineality(data.lineality, indices) || error("symmetry changes fan lineality")
        action = symmetry.ray_actions[group]
        sort(action) == collect(axes(data.rays, 1)) || error("invalid symmetry ray action")
        for ray in axes(data.rays, 1)
            quotient_ray_signature(data.rays[ray, indices], data.quotient) ==
                data.signatures[action[ray]] || error("incorrect symmetry image of ray $ray")
        end
        induced_cell_action(action, data.cells, cell_index) == symmetry.cell_actions[group] ||
            error("incorrect symmetry cell action")
    end
    action_orbits(symmetry.cell_actions, length(data.cells)) == symmetry.cell_orbits ||
        error("incorrect symmetry orbit partition")
    return true
end

"""
    certify_ray_projection(source_fan, source_rays, target_fan, target_rays,
                           source_symmetry, target_symmetry; ...)

Certify equality of projected maximal cones using positive ray multiples and
exact lineality residuals. Cone equality is checked once per source orbit;
verified coordinate actions transport the result to every other source cell.
Throws if a projected ray is neither vertical nor a positive target ray, or if
any orbit representative fails to project onto a target maximal cell.
"""
function certify_ray_projection(
    source_fan, source_rays, target_fan, target_rays, source_symmetry, target_symmetry;
    projection_columns,
    source_coordinate_indices,
    target_coordinate_indices=coordinate_permutation,
)
    source = projection_fan_data(source_fan, source_rays)
    target = projection_fan_data(target_fan, target_rays)
    columns = collect(projection_columns)
    projected = source.rays[:, columns]
    projected_q = projected * target.quotient
    target_q = target.rays * target.quotient
    lookup = Dict(signature => index for (index, signature) in enumerate(target.signatures))
    ray_records = Dict{String, Any}[]
    for ray in axes(projected, 1)
        image = projected_q[ray, :]
        if all(iszero, image)
            target_ray, scale = 0, QQ(0)
            residual = projected[ray:ray, :]
        else
            signature = quotient_ray_signature(projected[ray, :], target.quotient)
            target_ray = get(lookup, signature, 0)
            target_ray > 0 || error("source ray $ray does not project to a target ray")
            pivot = findfirst(!iszero, target_q[target_ray, :])
            scale = image[pivot] / target_q[target_ray, pivot]
            residual = projected[ray:ray, :] - scale .* target.rays[target_ray:target_ray, :]
        end
        coefficients = projection_lineality_coefficients(residual, target.lineality)
        push!(ray_records, Dict(
            "source_ray" => ray, "target_ray" => target_ray,
            "scale" => rational_key(scale),
            "lineality_coefficients" => projection_strings(vec(coefficients)),
        ))
    end
    lineality_coefficients = projection_lineality_coefficients(
        source.lineality[:, columns], target.lineality,
    )
    target_lookup = Dict(Tuple(cell) => index for (index, cell) in enumerate(target.cells))
    target_orbits = zeros(Int, length(target.cells))
    for (orbit, members) in enumerate(target_symmetry.cell_orbits)
        target_orbits[members] .= orbit
    end
    target_group = Dict(Tuple(p) => g for (g, p) in enumerate(target_symmetry.permutations))
    orbit_records, cell_records = Dict{String, Any}[], Dict{String, Any}[]
    ray_images = [record["target_ray"] for record in ray_records]
    for (orbit, members) in enumerate(source_symmetry.cell_orbits)
        representative = first(members)
        source_cell = source.cells[representative]
        image_rays = sort!(unique(filter(>(0), ray_images[source_cell])))
        target_cell = get(target_lookup, Tuple(image_rays), 0)
        target_cell > 0 || error("source orbit $orbit projects to $image_rays, not a target maximal cell")
        push!(orbit_records, Dict(
            "orbit" => orbit, "orbit_size" => length(members),
            "source_cell" => representative, "source_cells" => members,
            "source_rays" => source_cell, "target_cell" => target_cell,
            "target_orbit" => target_orbits[target_cell], "target_rays" => image_rays,
        ))
        for cell in members
            group = findfirst(action -> action[representative] == cell, source_symmetry.cell_actions)
            isnothing(group) && error("source cell $cell has no orbit transport")
            target_action = target_symmetry.cell_actions[target_group[Tuple(source_symmetry.permutations[group])]]
            push!(cell_records, Dict(
                "source_cell" => cell, "source_orbit" => orbit,
                "group_element" => group, "target_cell" => target_action[target_cell],
            ))
        end
    end
    sort!(cell_records; by=record -> record["source_cell"])
    certificate = Dict{String, Any}(
        "schema_version" => 1, "projection_columns" => columns,
        "source_maximal_cells" => length(source.cells), "target_maximal_cells" => length(target.cells),
        "group_permutations" => source_symmetry.permutations,
        "source_lineality_dimension" => size(source.lineality, 1),
        "target_lineality_dimension" => size(target.lineality, 1),
        "lineality_coefficients" => [projection_strings(row) for row in eachrow(lineality_coefficients)],
        "rays" => ray_records, "orbits" => orbit_records, "cells" => cell_records,
    )
    verify_ray_projection_certificate(
        certificate, source_fan, source_rays, target_fan, target_rays, source_symmetry, target_symmetry;
        projection_columns=columns, source_coordinate_indices, target_coordinate_indices,
    )
    return certificate
end

"Verify a saved certificate using only rational identities, incidence, and group actions."
function verify_ray_projection_certificate(
    certificate, source_fan, source_rays, target_fan, target_rays, source_symmetry, target_symmetry;
    projection_columns,
    source_coordinate_indices,
    target_coordinate_indices=coordinate_permutation,
)
    certificate["schema_version"] == 1 || error("unknown projection certificate schema")
    columns = collect(projection_columns)
    certificate["projection_columns"] == columns || error("incorrect certified projection")
    source = projection_fan_data(source_fan, source_rays)
    target = projection_fan_data(target_fan, target_rays)
    certificate["source_maximal_cells"] == length(source.cells) || error("source cell count mismatch")
    certificate["target_maximal_cells"] == length(target.cells) || error("target cell count mismatch")
    certificate["source_lineality_dimension"] == size(source.lineality, 1) || error("source lineality count mismatch")
    certificate["target_lineality_dimension"] == size(target.lineality, 1) || error("target lineality count mismatch")
    coefficients = isempty(certificate["lineality_coefficients"]) ?
        Matrix{QQFieldElem}(undef, 0, size(target.lineality, 1)) :
        reduce(vcat, permutedims.(projection_rationals.(certificate["lineality_coefficients"])))
    coefficients * target.lineality == source.lineality[:, columns] || error("incorrect lineality witness")
    rank(matrix(QQ, source.lineality[:, columns])) == size(target.lineality, 1) ||
        error("source lineality does not project onto target lineality")
    records = certificate["rays"]
    length(records) == size(source.rays, 1) || error("ray witnesses are incomplete")
    ray_images = zeros(Int, length(records))
    for (ray, record) in enumerate(records)
        record["source_ray"] == ray || error("incorrect source ray index")
        image = record["target_ray"]
        0 <= image <= size(target.rays, 1) || error("incorrect target ray index")
        scale = projection_rational(record["scale"])
        (image == 0 ? iszero(scale) : scale > 0) || error("ray scale must be positive, or zero for a vertical ray")
        residual = vec(transpose(projection_rationals(record["lineality_coefficients"])) * target.lineality)
        expected = image == 0 ? residual : scale .* target.rays[image, :] + residual
        source.rays[ray, columns] == expected || error("incorrect exact projection witness for source ray $ray")
        ray_images[ray] = image
    end
    verify_projection_symmetry(source, source_symmetry, source_coordinate_indices)
    verify_projection_symmetry(target, target_symmetry, target_coordinate_indices)
    certificate["group_permutations"] == source_symmetry.permutations || error("incorrect certificate group")
    Set(Tuple.(source_symmetry.permutations)) == Set(Tuple.(target_symmetry.permutations)) ||
        error("source and target symmetry groups disagree")
    target_group = Dict(Tuple(p) => g for (g, p) in enumerate(target_symmetry.permutations))
    for permutation in source_symmetry.permutations
        source_coordinate_indices(permutation)[columns] == columns[target_coordinate_indices(permutation)] ||
            error("projection does not commute with the coordinate actions")
    end
    orbits = certificate["orbits"]
    length(orbits) == length(source_symmetry.cell_orbits) || error("orbit witnesses are incomplete")
    for (orbit, record) in enumerate(orbits)
        members = source_symmetry.cell_orbits[orbit]
        record["orbit"] == orbit && record["orbit_size"] == length(members) &&
            record["source_cells"] == members && record["source_cell"] == first(members) ||
            error("incorrect source orbit witness")
        record["source_rays"] == source.cells[first(members)] || error("incorrect source cone generators")
        cell = record["target_cell"]
        1 <= cell <= length(target.cells) || error("invalid target maximal cell")
        cell in target_symmetry.cell_orbits[record["target_orbit"]] || error("incorrect target orbit")
        image = sort!(unique(filter(>(0), ray_images[record["source_rays"]])))
        image == record["target_rays"] == target.cells[cell] || error("projected cone is not equal to target maximal cone")
    end
    cells = certificate["cells"]
    length(cells) == length(source.cells) || error("cell transports are incomplete")
    for (cell, record) in enumerate(cells)
        record["source_cell"] == cell || error("incorrect transported source cell")
        orbit, group = record["source_orbit"], record["group_element"]
        representative = orbits[orbit]["source_cell"]
        source_symmetry.cell_actions[group][representative] == cell || error("incorrect source transport")
        target_action = target_symmetry.cell_actions[target_group[Tuple(source_symmetry.permutations[group])]]
        target_action[orbits[orbit]["target_cell"]] == record["target_cell"] || error("incorrect target transport")
    end
    return true
end

function write_projection_tables(orbit_io, cell_io, certificate)
    println(orbit_io, "flag_orbit,orbit_size,flag_cell,flag_rays,grassmannian_cell,grassmannian_orbit,grassmannian_rays")
    for record in certificate["orbits"]
        println(orbit_io, join([
            record["orbit"], record["orbit_size"], record["source_cell"], join(record["source_rays"], ';'),
            record["target_cell"], record["target_orbit"], join(record["target_rays"], ';'),
        ], ','))
    end
    println(cell_io, "flag_cell,flag_orbit,group_element,grassmannian_cell")
    for record in certificate["cells"]
        println(cell_io, join([record["source_cell"], record["source_orbit"], record["group_element"], record["target_cell"]], ','))
    end
end

function write_ray_projection_certificate(prefix, certificate)
    mkpath(dirname(prefix))
    open(prefix * "_certificate.toml", "w") do io
        TOML.print(io, certificate; sorted=true)
    end
    open(prefix * "_orbits.csv", "w") do orbit_io
        open(prefix * "_cells.csv", "w") do cell_io
            write_projection_tables(orbit_io, cell_io, certificate)
        end
    end
    return prefix
end
