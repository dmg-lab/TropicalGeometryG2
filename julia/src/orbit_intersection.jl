module OrbitIntersection

using SHA
using TOML

export expand_restricted_prevariety, canonical_all_ones_ray

"Canonical positive ray modulo all-ones lineality; opposite rays stay distinct."
function canonical_all_ones_ray(row)
    values = BigInt.(collect(row))
    isempty(values) && throw(ArgumentError("a ray has no coordinates"))
    shifted = values .- last(values)
    divisor = foldl(gcd, abs.(shifted); init=big(0))
    iszero(divisor) && throw(ArgumentError("a purported ray lies in all-ones lineality"))
    return Tuple(div.(shifted, divisor))
end

"""
Reduced row echelon basis of a lineality space, with pivots on the last
possible coordinates, as `(pivot, row)` pairs with `row[pivot] == 1`.
"""
function _lineality_reduction(rows, ambient)
    work = [Rational{BigInt}.(collect(row)) for row in rows]
    reduction = Tuple{Int,Vector{Rational{BigInt}}}[]
    for column in ambient:-1:1
        index = findfirst(row -> !iszero(row[column]), work)
        index === nothing && continue
        row = popat!(work, index)
        row = row ./ row[column]
        for other in work
            iszero(other[column]) || (other .-= other[column] .* row)
        end
        for (_, previous) in reduction
            iszero(previous[column]) || (previous .-= previous[column] .* row)
        end
        push!(reduction, (column, row))
    end
    all(iszero, Iterators.flatten(work)) || error("lineality reduction left a nonzero row")
    return reduction
end

function _reduce_modulo(row, reduction)
    values = Rational{BigInt}.(collect(row))
    for (pivot, basis_row) in reduction
        value = values[pivot]
        iszero(value) || (values .-= value .* basis_row)
    end
    return values
end

function _primitive(values)
    scale = foldl(lcm, denominator.(values); init=big(1))
    integers = numerator.(values .* scale)
    divisor = foldl(gcd, abs.(integers); init=big(0))
    iszero(divisor) && throw(ArgumentError("a purported ray lies in the lineality space"))
    return Tuple(div.(integers, divisor))
end

"""
Canonical positive ray modulo the lineality space: its pivot coordinates are
cleared and the result divided by a positive gcd. For the all-ones line this is
`canonical_all_ones_ray`.
"""
_canonical_ray(row, reduction) = _primitive(_reduce_modulo(row, reduction))

_in_span(row, reduction) = all(iszero, _reduce_modulo(row, reduction))

function _read_fan(path::AbstractString)
    text = read(path, String)
    sections = Dict{String,Vector{String}}()
    current = nothing
    for raw in split(text, '\n')
        line = strip(raw)
        isempty(line) && continue
        if occursin(r"^[A-Z][A-Z_]*$", line)
            current = String(line)
            sections[current] = String[]
        elseif current !== nothing
            push!(sections[current], String(line))
        end
    end
    scalar(name) = parse(Int, only(sections[name]))
    ambient = scalar("AMBIENT_DIM")
    rows(name) = [parse.(BigInt, split(strip(first(split(line, '#'; limit=2)))))
                  for line in get(sections, name, String[])
                  if !isempty(strip(first(split(line, '#'; limit=2))))]
    rays = rows("RAYS")
    all(length(ray) == ambient for ray in rays) || throw(ArgumentError("wrong ray ambient dimension"))
    lineality = rows("LINEALITY_SPACE")
    all(length(row) == ambient for row in lineality) || throw(ArgumentError("wrong lineality ambient dimension"))
    reduction = _lineality_reduction(lineality, ambient)
    length(reduction) == length(lineality) == scalar("LINEALITY_DIM") ||
        throw(ArgumentError("lineality rows are not a basis of LINEALITY_DIM dimensions"))
    name = haskey(sections, "MAXIMAL_CONES") ? "MAXIMAL_CONES" : "MAXIMAL_CONES_OF_CLOSURE"
    haskey(sections, name) || throw(ArgumentError("missing maximal cones"))
    cones = Vector{Vector{Int}}()
    for line in sections[name]
        data = strip(first(split(line, '#'; limit=2)))
        isempty(data) && continue
        match_value = match(r"^\{([^}]*)\}$", data)
        match_value === nothing && throw(ArgumentError("malformed maximal cone incidence"))
        inside = strip(only(match_value.captures))
        cone = isempty(inside) ? Int[] : parse.(Int, split(inside)) .+ 1
        all(1 <= index <= length(rays) for index in cone) || throw(ArgumentError("cone references a missing ray"))
        length(unique(cone)) == length(cone) || throw(ArgumentError("a cone repeats a ray"))
        push!(cones, sort!(cone))
    end
    isempty(cones) && throw(ArgumentError("fan has no maximal cones"))
    return (; text, ambient, rays, cones, lineality, reduction)
end

function _integer_rank(rows, ambient)
    isempty(rows) && return 0
    work = [BigInt.(collect(row)) for row in rows]
    rank = 0
    for column in 1:ambient
        pivot = findfirst(index -> !iszero(work[index][column]), (rank + 1):length(work))
        pivot === nothing && continue
        pivot_index = rank + pivot
        rank += 1
        work[rank], work[pivot_index] = work[pivot_index], work[rank]
        for index in (rank + 1):length(work)
            iszero(work[index][column]) && continue
            left, right = work[rank][column], work[index][column]
            work[index] = left .* work[index] .- right .* work[rank]
            divisor = foldl(gcd, abs.(work[index]); init=big(0))
            !iszero(divisor) && (work[index] = div.(work[index], divisor))
        end
        rank == length(work) && break
    end
    return rank
end

function _closed_group(permutations, ambient)
    group = sort!(unique(Tuple(Int.(collect(permutation))) for permutation in permutations))
    isempty(group) && throw(ArgumentError("coordinate group is empty"))
    all(length(permutation) == ambient && sort(collect(permutation)) == collect(1:ambient)
        for permutation in group) || throw(ArgumentError("invalid coordinate permutation"))
    group_set = Set(group)
    Tuple(1:ambient) in group_set || throw(ArgumentError("coordinate group omits the identity"))
    all(Tuple(left[index] for index in right) in group_set for left in group for right in group) ||
        throw(ArgumentError("supply the full closed coordinate group, not just generators"))
    return group
end

function _previous_symmetry(previous, group, reduction)
    rays = [_canonical_ray(ray, reduction) for ray in previous.rays]
    length(unique(rays)) == length(rays) || throw(ArgumentError("previous quotient rays are not distinct"))
    ray_index = Dict(ray => index for (index, ray) in enumerate(rays))
    cone_set = Set(Tuple(cone) for cone in previous.cones)
    length(cone_set) == length(previous.cones) || throw(ArgumentError("previous maximal cones repeat"))
    for permutation in group
        action = Int[]
        for ray in previous.rays
            image = _canonical_ray(ray[collect(permutation)], reduction)
            haskey(ray_index, image) || throw(ArgumentError("coordinate group fails to preserve previous rays"))
            push!(action, ray_index[image])
        end
        all(Tuple(sort!(action[cone])) in cone_set for cone in previous.cones) ||
            throw(ArgumentError("coordinate group fails to preserve previous maximal cones"))
    end
    return true
end

_coordinate_group_digest(group) = bytes2hex(sha256(join((join(permutation, ",") for permutation in group), '\n') * "\n"))

function _restriction_proof(path, previous_path, restricted_path, group, previous;
                            restriction_input_path, basis_input_path)
    path === nothing && return nothing
    previous_path === nothing && throw(ArgumentError("a restriction certificate requires the previous fan"))
    restriction_input_path === nothing && throw(ArgumentError("a restriction certificate requires its restriction input"))
    basis_input_path === nothing && throw(ArgumentError("a restriction certificate requires its basis input"))
    proof = TOML.parsefile(path)
    required = ("support_preserved", "all_previous_cells_covered", "previous_basis_subset",
                "invariant_basis_closed", "closed_restrictions", "group_closed")
    all(get(proof, flag, false) === true for flag in required) ||
        throw(ArgumentError("restriction coverage/common-refinement certificate is incomplete"))
    for (key, input) in (("source_fan_sha256", previous_path),
                         ("restricted_fan_sha256", restricted_path),
                         ("restriction_input_sha256", restriction_input_path),
                         ("basis_input_sha256", basis_input_path))
        get(proof, key, "") == bytes2hex(sha256(read(input))) ||
            throw(ArgumentError("restriction certificate digest mismatch: $key"))
    end
    get(proof, "group_order", 0) == length(group) || throw(ArgumentError("restriction group order mismatch"))
    get(proof, "coordinate_group_sha256", "") == _coordinate_group_digest(group) ||
        throw(ArgumentError("restriction coordinate group digest mismatch"))
    get(proof, "number_of_previous_cells", 0) == length(previous.cones) ||
        throw(ArgumentError("restriction certificate omits previous cells"))
    get(proof, "number_of_representatives", 0) > 0 || throw(ArgumentError("restriction certificate has no representatives"))
    return proof
end

function _atomic_write(path, text)
    mkpath(dirname(path))
    temporary, io = mktemp(dirname(path))
    try
        write(io, text)
        close(io)
        mv(temporary, path; force=true)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
end

function _fan_text(ambient, rays, cones, dimensions, reduction)
    lines = String["_application PolyhedralFan", "_version 2.2", "_type PolyhedralFan", "",
                   "AMBIENT_DIM", string(ambient), "", "DIM", string(maximum(dimensions)), "",
                   "LINEALITY_DIM", string(length(reduction)), "", "RAYS"]
    append!(lines, [join(ray, " ") for ray in rays])
    append!(lines, ["", "N_RAYS", string(length(rays)), "", "LINEALITY_SPACE"])
    append!(lines, [join(_primitive(row), " ") for (_, row) in reverse(reduction)])
    append!(lines, ["", "ORTH_LINEALITY_SPACE"])
    # One kernel vector per free coordinate; for the all-ones line, e_i - e_n.
    pivots = Set(first.(reduction))
    for free in 1:ambient
        free in pivots && continue
        row = zeros(Rational{BigInt}, ambient)
        row[free] = 1
        for (pivot, basis_row) in reduction
            row[pivot] = -basis_row[free]
        end
        push!(lines, join(_primitive(row), " "))
    end
    append!(lines, ["", "MAXIMAL_CONES"])
    for (cone, dimension) in zip(cones, dimensions)
        push!(lines, "{" * join(cone .- 1, " ") * "} # Dimension $dimension")
    end
    append!(lines, ["", "SIMPLICIAL", all(length(cone) == dimension - length(reduction) for (cone, dimension) in zip(cones, dimensions)) ? "1" : "0",
                    "", "PURE", length(unique(dimensions)) == 1 ? "1" : "0", ""])
    return join(lines, '\n')
end

"""
    expand_restricted_prevariety(restricted_fan_path, output_path, permutations; ...)

Expand a max-convention Gfan fan by a complete group of coordinate *pullback*
permutations (`row[permutation]`), which must preserve the common lineality
space.  Quotient rays are normalized by clearing the pivot coordinates of that
space, chosen last (for the all-ones line: subtracting the final coordinate),
and dividing by a positive gcd.  Deduplication and removal of contained faces
preserve support exactly, including shared closed boundaries.

A `restriction_certificate_path`, together with the previous fan and both
input paths, authorizes `support_checked=true` only after checking the hashes
and the complete closed-restriction/common-invariant-basis proof.  Under those
conditions all pieces refine the same normal fan, so their overlaps agree.
Without that certificate the exact expansion is still available for fixtures,
but compatibility/coverage is explicitly left uncertified.

`progress(description, total)` returns a function called with the number of
steps done since its last call, for reporting the two long loops.
"""
function expand_restricted_prevariety(
    restricted_fan_path::AbstractString,
    output_path::AbstractString,
    permutations;
    previous_fan_path=nothing,
    provenance_path=String(output_path) * ".expansion.toml",
    restriction_certificate_path=nothing,
    restriction_input_path=nothing,
    basis_input_path=nothing,
    progress=(description, total) -> Returns(nothing),
)
    abspath(restricted_fan_path) == abspath(output_path) && throw(ArgumentError("output would overwrite the restricted fan"))
    previous_fan_path !== nothing && abspath(previous_fan_path) == abspath(output_path) &&
        throw(ArgumentError("output would overwrite the previous fan"))
    restricted = _read_fan(restricted_fan_path)
    group = _closed_group(permutations, restricted.ambient)
    previous = previous_fan_path === nothing ? nothing : _read_fan(previous_fan_path)
    reduction = restricted.reduction
    all(_in_span(row[collect(permutation)], reduction) for row in restricted.lineality for permutation in group) ||
        throw(ArgumentError("coordinate group fails to preserve the lineality space"))
    if previous !== nothing
        previous.ambient == restricted.ambient || throw(ArgumentError("previous fan ambient dimension differs"))
        length(previous.lineality) == length(reduction) && all(row -> _in_span(row, reduction), previous.lineality) ||
            throw(ArgumentError("previous fan lineality space differs"))
        _previous_symmetry(previous, group, reduction)
    end
    proof = _restriction_proof(restriction_certificate_path, previous_fan_path, restricted_fan_path, group, previous;
                               restriction_input_path, basis_input_path)
    normalized = [_canonical_ray(ray, reduction) for ray in restricted.rays]
    source_dimensions = [_integer_rank([normalized[index] for index in cone], restricted.ambient) + length(reduction)
                         for cone in restricted.cones]
    rays = Tuple[]
    ray_index = Dict{Tuple,Int}()
    cone_dimensions = Dict{Tuple,Int}()
    expanded = progress("orbit expansion", length(group) * length(restricted.cones))
    for permutation in group
        action = Int[]
        for ray in restricted.rays
            image = _canonical_ray(ray[collect(permutation)], reduction)
            target = get!(ray_index, image) do
                push!(rays, image)
                length(rays)
            end
            push!(action, target)
        end
        for (cone, dimension) in zip(restricted.cones, source_dimensions)
            image = Tuple(sort!(unique(action[cone])))
            get(cone_dimensions, image, dimension) == dimension || error("orbit cone dimension changed")
            cone_dimensions[image] = dimension
        end
        expanded(length(restricted.cones))
    end
    before_pruning = length(cone_dimensions)
    # Search only previously retained larger-dimensional cells containing the
    # rarest ray.  Incidence inclusion is an exact containment certificate.
    candidates = sort!(collect(keys(cone_dimensions)); by=cone -> (-cone_dimensions[cone], -length(cone), cone))
    retained = Tuple[]
    containing = Dict{Int,Vector{Int}}()
    removed = 0
    pruned = progress("removal of contained cones", length(candidates))
    for cone in candidates
        pruned(1)
        possible = if isempty(cone)
            eachindex(retained)
        else
            ray = argmin(index -> length(get(containing, index, Int[])), cone)
            get(containing, ray, Int[])
        end
        contained = any(index -> cone_dimensions[retained[index]] > cone_dimensions[cone] &&
                                 all(ray in retained[index] for ray in cone), possible)
        if contained
            removed += 1
            continue
        end
        push!(retained, cone)
        for ray in cone
            push!(get!(containing, ray, Int[]), length(retained))
        end
    end
    used = sort!(unique(index for cone in retained for index in cone); by=index -> rays[index])
    renumber = Dict(old => new for (new, old) in enumerate(used))
    output_rays = [rays[index] for index in used]
    output_pairs = [(Tuple(sort!([renumber[index] for index in cone])), cone_dimensions[cone]) for cone in retained]
    sort!(output_pairs; by=pair -> (last(pair), first(pair)))
    output_cones = [collect(first(pair)) for pair in output_pairs]
    dimensions = last.(output_pairs)
    output_text = _fan_text(restricted.ambient, output_rays, output_cones, dimensions, reduction)
    _atomic_write(output_path, output_text)
    checked = proof !== nothing
    metadata = Dict{String,Any}(
        "method" => "complete coordinate-orbit expansion of closed common-normal-fan refinements",
        "ambient_dimension" => restricted.ambient, "lineality_dimension" => length(reduction),
        "dimension" => maximum(dimensions), "input_rays" => length(restricted.rays),
        "input_cones" => length(restricted.cones), "output_rays" => length(output_rays),
        "output_cones" => length(output_cones), "group_order" => length(group), "group_closed" => true,
        "expanded_cones_before_dedup" => length(group) * length(restricted.cones),
        "unique_cones_before_pruning" => before_pruning, "removed_contained_cones" => removed,
        "support_checked" => checked, "fan_compatibility_checked" => checked,
        "input_sha256" => bytes2hex(sha256(restricted.text)), "output_sha256" => bytes2hex(sha256(output_text)),
        "previous_sha256" => previous_fan_path === nothing ? "" : bytes2hex(sha256(read(previous_fan_path))),
        "previous_group_symmetry_checked" => previous !== nothing,
        "coordinate_group_sha256" => _coordinate_group_digest(group),
        "restriction_certificate_sha256" => restriction_certificate_path === nothing ? "" : bytes2hex(sha256(read(restriction_certificate_path))),
        "maximal_cones_by_dimension" => Dict(string(dimension) => count(==(dimension), dimensions) for dimension in unique(dimensions)),
    )
    io = IOBuffer()
    TOML.print(io, metadata; sorted=true)
    _atomic_write(provenance_path, String(take!(io)))
    return (; fan_path=String(output_path), provenance_path=String(provenance_path),
             input_rays=length(restricted.rays), input_cones=length(restricted.cones),
             output_rays=length(output_rays), output_cones=length(output_cones),
             expanded_cones_before_dedup=length(group) * length(restricted.cones),
             unique_cones_before_pruning=before_pruning, removed_contained_cones=removed,
             ambient_dimension=restricted.ambient, lineality_dimension=length(reduction), dimension=maximum(dimensions),
             group_order=length(group), group_closed=true, support_checked=checked,
             fan_compatibility_checked=checked, previous_group_symmetry_checked=previous !== nothing,
             input_sha256=metadata["input_sha256"], output_sha256=metadata["output_sha256"],
             previous_sha256=metadata["previous_sha256"], coordinate_group_sha256=metadata["coordinate_group_sha256"])
end

end
