# Geometry and labels used by the local G2 flag cell atlas.
# Keep these calculations separate from the removed historical renderers.

include(joinpath(@__DIR__, "..", "grassmannian", "atlas_common.jl"))

const CLEAN_ATLAS_FAN = TropicalCensus.data_path("g2", "flag", "02_coarsest_fan.mrdi")
const CLEAN_ATLAS_GRASS_FAN =
    TropicalCensus.data_path("g2", "grassmannian", "01_raw_fan.mrdi")

clean_vector(M, i::Int) = [QQ(M[i, j]) for j in 1:ncols(M)]
clean_tex_set(values) = "\\{" * join(values, ",\\,") * "\\}"
clean_quoted(values) = "\"" * join(values, ';') * "\""

function clean_positive_ray_match(v, rays)
    pivot = findfirst(!iszero, v)
    pivot === nothing && return 0
    for i in 1:nrows(rays)
        iszero(rays[i, pivot]) && continue
        scale = v[pivot] / rays[i, pivot]
        scale > 0 || continue
        v == scale .* clean_vector(rays, i) && return i
    end
    return -1
end

function clean_graph_distance(tree::TreeMetric, source::Int, target::Int)
    source == target && return QQ(0)
    stack = [(source, 0, QQ(0))]
    while !isempty(stack)
        vertex, parent, distance = pop!(stack)
        for (neighbor, edge_length) in tree.adjacency[vertex]
            neighbor == parent && continue
            neighbor == target && return distance + edge_length
            push!(stack, (neighbor, vertex, distance + edge_length))
        end
    end
    error("tree vertices $source and $target are disconnected")
end

function clean_internal_vertex(tree::TreeMetric, vertex::Int)
    x = [-clean_graph_distance(tree, vertex, taxon) for taxon in 1:tree.ntaxa]
    return x .- x[7]
end

function clean_taxon_direction(taxon::Int)
    direction = zeros(QQFieldElem, 7)
    if taxon <= 6
        direction[taxon] = 1
    else
        direction[1:6] .= -1
    end
    return direction
end

function clean_direction_taxon(direction)
    iszero(direction[7]) || error("fiber direction is not in the x_7=0 gauge")
    for taxon in 1:6
        if direction[taxon] > 0 &&
           all(iszero(direction[index]) for index in 1:6 if index != taxon)
            return taxon
        end
    end
    direction[1] < 0 && all(direction[index] == direction[1] for index in 2:6) &&
        return 7
    error("unrecognized tropical-line direction $direction")
end

function clean_segment_parameter(point, left, right)
    difference = right .- left
    pivot = findfirst(!iszero, difference)
    pivot === nothing && return nothing
    parameter = (point[pivot] - left[pivot]) / difference[pivot]
    left .+ parameter .* difference == point || return nothing
    0 <= parameter <= 1 || return nothing
    return parameter
end

function clean_ray_parameter(point, origin, direction)
    pivot = findfirst(!iszero, direction)
    parameter = (point[pivot] - origin[pivot]) / direction[pivot]
    origin .+ parameter .* direction == point || return nothing
    parameter >= 0 || return nothing
    return parameter
end

function clean_locate_point(point, tree, vertex_points)
    internal = sort([vertex for vertex in keys(tree.adjacency) if vertex > tree.ntaxa])
    for vertex in internal
        point == vertex_points[vertex] && return (kind=:vertex, vertex=vertex)
    end
    for left in internal, (right, _) in tree.adjacency[left]
        left < right && right > tree.ntaxa || continue
        parameter = clean_segment_parameter(point, vertex_points[left], vertex_points[right])
        parameter === nothing || return (
            kind=:edge, left=left, right=right, parameter=QQ(parameter),
        )
    end
    for taxon in 1:tree.ntaxa
        attachment = only(tree.adjacency[taxon])[1]
        parameter = clean_ray_parameter(
            point, vertex_points[attachment], clean_taxon_direction(taxon),
        )
        parameter === nothing || return (
            kind=:ray, taxon=taxon, attachment=attachment, parameter=QQ(parameter),
        )
    end
    error("flag point is not on its reconstructed tropical line: $point")
end

clean_update_lower(current, candidate) = current === nothing ? candidate : max(current, candidate)
clean_update_upper(current, candidate) = current === nothing ? candidate : min(current, candidate)

"""Slice a simplicial four-ray subcone by p=w and x_7=0."""
function clean_slice_four(flag_rays, flag_lineality, weight, cell::Int, ray_ids)
    generators = vcat(flag_rays[ray_ids, :], flag_lineality)
    equations = matrix(QQ, vcat(
        transpose(generators[:, 8:28]), transpose(generators[:, 7:7]),
    ))
    target = matrix(QQ, reshape(vcat(weight, [QQ(0)]), :, 1))
    solvable, solution = can_solve_with_solution(equations, target; side=:right)
    solvable || return nothing
    kernel_basis = kernel(equations; side=:right)
    size(kernel_basis) == (8, 1) || return nothing
    z0 = vec(Matrix{QQFieldElem}(solution))
    dz = vec(Matrix{QQFieldElem}(kernel_basis))
    lower = nothing
    upper = nothing
    for index in 1:4
        if dz[index] > 0
            lower = clean_update_lower(lower, -z0[index] / dz[index])
        elseif dz[index] < 0
            upper = clean_update_upper(upper, -z0[index] / dz[index])
        elseif z0[index] < 0
            return nothing
        end
    end
    lower !== nothing && upper !== nothing && lower > upper && return nothing
    xmap = matrix(QQ, transpose(generators[:, 1:7]))
    x0 = vec(Matrix{QQFieldElem}(xmap * solution))
    dx = vec(Matrix{QQFieldElem}(xmap * kernel_basis))
    point(parameter) = x0 .+ parameter .* dx
    if lower !== nothing && upper !== nothing
        return (
            cell=cell, flag_rays=ray_ids, kind=:bounded, left=point(lower),
            right=point(upper), taxon=0, degenerate=lower == upper,
        )
    end
    xor(lower === nothing, upper === nothing) || return nothing
    endpoint = lower === nothing ? upper : lower
    outward = lower === nothing ? -dx : dx
    return (
        cell=cell, flag_rays=ray_ids, kind=:unbounded, left=point(endpoint),
        right=QQFieldElem[], taxon=clean_direction_taxon(outward), degenerate=false,
    )
end

function clean_four_subsets(ray_ids)
    return [ray_ids[[a, b, c, d]]
            for a in 1:(length(ray_ids) - 3)
            for b in (a + 1):(length(ray_ids) - 2)
            for c in (b + 1):(length(ray_ids) - 1)
            for d in (c + 1):length(ray_ids)]
end

function clean_merge_slices(parts, cell::Int, ray_ids)
    isempty(parts) && return nothing
    finite = Vector{Vector{QQFieldElem}}()
    for part in parts
        push!(finite, part.left)
        part.kind == :bounded && push!(finite, part.right)
    end
    points = collect(values(Dict(Tuple(point) => point for point in finite)))
    origin = first(points)
    direction = nothing
    for point in points
        difference = point .- origin
        any(!iszero, difference) && (direction = difference; break)
    end
    if direction === nothing
        unbounded = findfirst(part -> part.kind == :unbounded, parts)
        direction = unbounded === nothing ? zeros(QQFieldElem, 7) :
                    clean_taxon_direction(parts[unbounded].taxon)
    end
    pivot = findfirst(!iszero, direction)
    pivot === nothing && return (
        cell=cell, flag_rays=ray_ids, kind=:bounded, left=origin, right=origin,
        taxon=0, degenerate=true,
    )
    function scalar(point)
        parameter = (point[pivot] - origin[pivot]) / direction[pivot]
        origin .+ parameter .* direction == point || error(
            "four-ray slices of clean cell $cell are not collinear",
        )
        return parameter
    end
    parameters = Dict(scalar(point) => point for point in points)
    lower, upper = minimum(keys(parameters)), maximum(keys(parameters))
    lower_infinite = false
    upper_infinite = false
    lower_taxon = 0
    upper_taxon = 0
    for part in parts
        part.kind == :unbounded || continue
        ray_direction = clean_taxon_direction(part.taxon)
        scale = ray_direction[pivot] / direction[pivot]
        if scale < 0
            lower_infinite, lower_taxon = true, part.taxon
        else
            upper_infinite, upper_taxon = true, part.taxon
        end
    end
    lower_infinite && upper_infinite && error("clean cell $cell contains a full fiber line")
    if upper_infinite
        return (cell=cell, flag_rays=ray_ids, kind=:unbounded,
                left=parameters[lower], right=QQFieldElem[], taxon=upper_taxon,
                degenerate=false)
    elseif lower_infinite
        return (cell=cell, flag_rays=ray_ids, kind=:unbounded,
                left=parameters[upper], right=QQFieldElem[], taxon=lower_taxon,
                degenerate=false)
    end
    return (cell=cell, flag_rays=ray_ids, kind=:bounded,
            left=parameters[lower], right=parameters[upper], taxon=0,
            degenerate=lower == upper)
end

function clean_slice_cell(flag_rays, flag_lineality, incidence, weight, cell::Int)
    ray_ids = findall(incidence[cell, :])
    length(ray_ids) >= 4 || error("clean cell $cell has fewer than four rays")
    parts = Any[]
    for subset in clean_four_subsets(ray_ids)
        part = clean_slice_four(flag_rays, flag_lineality, weight, cell, subset)
        part === nothing || push!(parts, part)
    end
    return clean_merge_slices(parts, cell, ray_ids)
end

function clean_facet_data(flag_rays, flag_lineality, incidence)
    owners = Dict{Tuple{Vararg{Int}}, Set{Int}}()
    for cell in axes(incidence, 1)
        ray_ids = findall(incidence[cell, :])
        cone = positive_hull(flag_rays[ray_ids, :], flag_lineality; non_redundant=true)
        for halfspace in facets(Halfspace, cone)
            normal = QQ.(collect(normal_vector(halfspace)))
            key = Tuple(sort!([
                ray for ray in ray_ids if iszero(sum(
                    normal[j] * flag_rays[ray, j] for j in 1:ncols(flag_rays)
                ))
            ]))
            rank(matrix(QQ, vcat(flag_lineality, flag_rays[collect(key), :]))) == 7 ||
                error("a clean maximal-cone facet has the wrong dimension")
            push!(get!(owners, key, Set{Int}()), cell)
        end
    end
    keys_sorted = sort!(collect(keys(owners)))
    return [(wall=index, rays=collect(key), owners=sort!(collect(owners[key])))
            for (index, key) in enumerate(keys_sorted)]
end

"""
The point where the fiber over `weight` meets the relative interior of a wall,
or `nothing`.  The coefficient vectors of the wall's rays over `weight` form a
point or, for a square wall, a line; the wall is met when one of them is
positive on every ray.  Testing the whole wall, rather than its three-ray
subcones, also finds square walls over the crossing of the square's diagonals.
"""
function clean_slice_wall(flag_rays, flag_lineality, weight, wall)
    count = length(wall.rays)
    generators = vcat(flag_rays[wall.rays, :], flag_lineality)
    equations = matrix(QQ, vcat(
        transpose(generators[:, 8:28]), transpose(generators[:, 7:7]),
    ))
    target = matrix(QQ, reshape(vcat(weight, [QQ(0)]), :, 1))
    solvable, solution = can_solve_with_solution(equations, target; side=:right)
    solvable || return nothing
    kernel_basis = kernel(equations; side=:right)
    xmap = matrix(QQ, transpose(generators[:, 1:7]))
    iszero(xmap * kernel_basis) || error("clean wall $(wall.wall) has an x-fiber")
    z0 = vec(Matrix{QQFieldElem}(solution))[1:count]
    relations = Matrix{QQFieldElem}(kernel_basis)[1:count, :]
    rank(matrix(QQ, relations)) <= 1 ||
        error("clean wall $(wall.wall) has several relations among its rays")
    column = findfirst(index -> any(!iszero, relations[:, index]), axes(relations, 2))
    lower = nothing
    upper = nothing
    for index in 1:count
        slope = column === nothing ? QQ(0) : relations[index, column]
        if slope > 0
            lower = clean_update_lower(lower, -z0[index] / slope)
        elseif slope < 0
            upper = clean_update_upper(upper, -z0[index] / slope)
        elseif z0[index] <= 0
            return nothing
        end
    end
    lower !== nothing && upper !== nothing && lower >= upper && return nothing
    return vec(Matrix{QQFieldElem}(xmap * solution))
end

function clean_projection_data(flag, balanced_rays, support_specs)
    data = rays_modulo_lineality(flag)
    flag_rays = Matrix{QQFieldElem}(matrix(QQ, data.rays_modulo_lineality))
    flag_lineality = Matrix{QQFieldElem}(matrix(QQ, data.lineality_basis))
    incidence = Matrix{Bool}(maximal_cones(IncidenceMatrix, flag))
    grass = load(CLEAN_ATLAS_GRASS_FAN)
    grass_lineality = matrix(QQ, rays_modulo_lineality(grass).lineality_basis)
    quotient = Matrix{QQFieldElem}(kernel(grass_lineality; side=:right))
    balanced = Matrix{QQFieldElem}(balanced_rays)
    balanced_q = balanced * quotient
    flag_q = flag_rays[:, 8:28] * quotient
    matches = [clean_positive_ray_match(clean_vector(flag_q, ray), balanced_q)
               for ray in 1:nrows(flag_q)]
    all(match -> match >= 0, matches) || error("a surviving flag ray projects inside a Grassmannian cell")
    support_cells = [spec.representative for spec in support_specs]
    append!(support_cells, [member for spec in support_specs for member in spec.members[2:end]])
    support_lookup = Dict(Tuple(sort!(collect(cell))) => cell for cell in support_cells)
    by_support = Dict{Tuple{Vararg{Int}}, Vector{Int}}()
    projected = Vector{Vector{Int}}(undef, nrows(incidence))
    for cell in axes(incidence, 1)
        ray_ids = findall(incidence[cell, :])
        image = sort!(unique(filter(>(0), matches[ray_ids])))
        key = Tuple(image)
        haskey(support_lookup, key) || error(
            "clean flag cell $cell projects to rays $image, not a support maximal cell",
        )
        projected[cell] = image
        push!(get!(by_support, key, Int[]), cell)
    end
    length(by_support) == length(support_cells) ||
        error("clean flag cones do not cover all $(length(support_cells)) support cells")
    return (
        rays=flag_rays, lineality=flag_lineality, incidence=incidence,
        ray_matches=matches, projected=projected, by_support=by_support,
        grass_quotient=quotient, balanced_q=balanced_q,
    )
end

"""Put one flag ray in the exact p-gauge of its balanced Grassmannian ray."""
function clean_align_flag_ray(flag_data, balanced_rays, flag_ray::Int, grass_ray::Int)
    vector = clean_vector(flag_data.rays, flag_ray)
    vector_q = vec(reshape(vector[8:28], 1, :) * flag_data.grass_quotient)
    target_q = clean_vector(flag_data.balanced_q, grass_ray)
    pivot = findfirst(!iszero, target_q)
    pivot === nothing && error("cannot align a nonvertical flag ray to the zero ray")
    scale = vector_q[pivot] / target_q[pivot]
    scale > 0 && vector_q == scale .* target_q || error("flag-ray quotient mismatch")
    aligned = vector ./ scale
    target_p = clean_vector(balanced_rays, grass_ray)
    right_hand_side = vcat(target_p .- aligned[8:28], [-aligned[7]])
    equations = matrix(QQ, vcat(
        transpose(flag_data.lineality[:, 8:28]),
        transpose(flag_data.lineality[:, 7:7]),
    ))
    solvable, coefficients = can_solve_with_solution(
        equations,
        matrix(QQ, reshape(right_hand_side, :, 1));
        side=:right,
    )
    solvable || error("could not lift the balanced p-gauge through flag lineality")
    coefficient_matrix = Matrix{QQFieldElem}(coefficients)
    aligned .+= vec(transpose(flag_data.lineality) * coefficient_matrix)
    aligned[8:28] == target_p || error("balanced p-gauge lift failed")
    iszero(aligned[7]) || error("balanced x_7 gauge lift failed")
    return aligned
end

function clean_support_trials(number_of_rays::Int)
    trials = [Tuple(entries) for entries in Iterators.product(fill(1:4, number_of_rays)...)
              if foldl(gcd, entries) == 1]
    sort!(trials; by=entries -> (sum(entries), maximum(entries), entries))
    return trials
end

function clean_support_weight(balanced_rays, support_rays, coefficients)
    weight = zeros(QQFieldElem, 21)
    for (coefficient, ray) in zip(coefficients, support_rays)
        weight .+= QQ(coefficient) .* clean_vector(balanced_rays, ray)
    end
    return weight
end

function clean_decorate_sample(sample, walls, flag_data)
    tree = tree_from_min_pluecker(sample.weight)
    internal = sort([vertex for vertex in keys(tree.adjacency) if vertex > tree.ntaxa])
    vertex_points = Dict(vertex => clean_internal_vertex(tree, vertex) for vertex in internal)
    decorated = Any[]
    actual_portions = Set{Tuple}()
    for slice in sample.slices
        left_location = clean_locate_point(slice.left, tree, vertex_points)
        left_location.kind == :vertex || error(
            "clean cell $(slice.cell) still starts inside a tree edge/ray",
        )
        if slice.kind == :bounded
            right_location = clean_locate_point(slice.right, tree, vertex_points)
            right_location.kind == :vertex || error(
                "clean cell $(slice.cell) still ends inside a bounded tree edge",
            )
            key = (:edge, min(left_location.vertex, right_location.vertex),
                   max(left_location.vertex, right_location.vertex))
        else
            right_location = nothing
            key = (:ray, slice.taxon)
        end
        key in actual_portions && error("two clean cells still cover the same tree portion $key")
        push!(actual_portions, key)
        push!(decorated, merge(slice, (
            left_location=left_location, right_location=right_location, portion=key,
        )))
    end
    expected = Set{Tuple}((:ray, taxon) for taxon in 1:tree.ntaxa)
    for left in internal, (right, _) in tree.adjacency[left]
        left < right && right > tree.ntaxa && push!(expected, (:edge, left, right))
    end
    actual_portions == expected || error(
        "fiber portions $(sort!(collect(actual_portions))) do not equal tree portions " *
        "$(sort!(collect(expected)))",
    )

    active = Set(slice.cell for slice in decorated)
    wall_targets = Any[]
    for wall in walls
        count(owner -> owner in active, wall.owners) >= 2 || continue
        point = clean_slice_wall(flag_data.rays, flag_data.lineality, sample.weight, wall)
        point === nothing && continue
        location = clean_locate_point(point, tree, vertex_points)
        location.kind == :vertex || error(
            "surviving clean wall $(wall.wall) lands inside a tree edge",
        )
        push!(wall_targets, merge(wall, (location=location,)))
    end
    return merge(sample, (
        tree=tree, vertex_points=vertex_points, slices=sort!(decorated; by=s -> s.cell),
        walls=sort!(wall_targets; by=w -> w.wall),
    ))
end

function clean_generic_sample(spec, balanced_rays, walls, flag_data)
    support_rays = collect(spec.representative)
    candidates = flag_data.by_support[Tuple(support_rays)]
    for coefficients in clean_support_trials(length(support_rays))
        weight = clean_support_weight(balanced_rays, support_rays, coefficients)
        slices = Any[]
        valid = true
        for cell in candidates
            slice = clean_slice_cell(
                flag_data.rays, flag_data.lineality, flag_data.incidence, weight, cell,
            )
            if slice === nothing || slice.degenerate
                valid = false
                break
            end
            push!(slices, slice)
        end
        valid || continue
        sample = (weight=weight, coefficients=coefficients, support=spec, slices=slices)
        try
            return clean_decorate_sample(sample, walls, flag_data)
        catch exception
            exception isa ErrorException || rethrow()
        end
    end
    error("no unsplit generic fiber sample for support representative $(spec.representative)")
end

"""
The taxon towards which a flag ray with zero support projection moves the
point: its direction in the gauge with zero Pluecker part and x_7 = 0.
"""
function clean_vertical_taxon(flag_data, flag_ray::Int)
    vector = clean_vector(flag_data.rays, flag_ray)
    equations = matrix(QQ, vcat(
        transpose(flag_data.lineality[:, 8:28]),
        transpose(flag_data.lineality[:, 7:7]),
    ))
    solvable, coefficients = can_solve_with_solution(
        equations,
        matrix(QQ, reshape(vcat(-vector[8:28], [-vector[7]]), :, 1));
        side=:right,
    )
    solvable || error("flag ray $flag_ray has a nonzero support projection")
    vector .+= vec(transpose(flag_data.lineality) * Matrix{QQFieldElem}(coefficients))
    return clean_direction_taxon(vector[1:7])
end

const FLAG_TREE_LETTERS = [string(Char(Int('A') + index - 1)) for index in 1:12]

flag_tree_tex_tuple(values) = "(" * join(values, ",\\,") * ")"

const FLAG_TREE_TAXON_MACROS = (
    "taxonone", "taxontwo", "taxonthree", "taxonfour",
    "taxonfive", "taxonsix", "taxonseven",
)

function flag_tree_edges(tree)
    result = Tuple{Int,Int,QQFieldElem}[]
    for left in sort!(collect(keys(tree.adjacency)))
        for (right, edge_length) in tree.adjacency[left]
            left < right && push!(result, (left, right, QQ(edge_length)))
        end
    end
    sort!(result; by=edge -> (
        edge[1] <= tree.ntaxa || edge[2] <= tree.ntaxa,
        edge[1], edge[2],
    ))
    return result
end

function flag_tree_exception_endpoints(tree, portion)
    if portion[1] == :edge
        return portion[2], portion[3]
    elseif portion[1] == :ray
        taxon = portion[2]
        return taxon, only(tree.adjacency[taxon])[1]
    end
    error("unknown tree portion $portion")
end

function flag_tree_ray_lifts(flag_data, balanced_rays)
    ray_records = [balance_outer_edges(clean_vector(balanced_rays, ray))
                   for ray in 1:nrows(balanced_rays)]
    by_grass = Dict{Int, Dict{Int,Int}}()
    flag_to_vertex = Dict{Int,Int}()
    for flag_ray in 1:nrows(flag_data.rays)
        grass_ray = flag_data.ray_matches[flag_ray]
        grass_ray == 0 && continue
        record = ray_records[grass_ray]
        internal = sort([vertex for vertex in keys(record.tree.adjacency)
                         if vertex > record.tree.ntaxa])
        vertex_points = Dict(vertex => clean_internal_vertex(record.tree, vertex)
                             for vertex in internal)
        aligned = clean_align_flag_ray(flag_data, balanced_rays, flag_ray, grass_ray)
        location = clean_locate_point(aligned[1:7], record.tree, vertex_points)
        location.kind == :vertex || error("flag ray $flag_ray is not at an internal node")
        node_map = get!(by_grass, grass_ray, Dict{Int,Int}())
        haskey(node_map, location.vertex) && error(
            "two flag rays project to node $(location.vertex) of G2 ray $grass_ray",
        )
        node_map[location.vertex] = flag_ray
        flag_to_vertex[flag_ray] = location.vertex
    end

    letters = Dict{Int, Dict{Int,String}}()
    for grass_ray in 1:length(ray_records)
        internal = sort([vertex for vertex in keys(ray_records[grass_ray].tree.adjacency)
                         if vertex > ray_records[grass_ray].tree.ntaxa])
        node_map = get(by_grass, grass_ray, Dict{Int,Int}())
        Set(keys(node_map)) == Set(internal) || error(
            "G2 ray $grass_ray does not have exactly one flag ray at every internal node",
        )
        letters[grass_ray] = Dict(
            vertex => FLAG_TREE_LETTERS[position]
            for (position, vertex) in enumerate(internal)
        )
    end
    count(>(0), flag_data.ray_matches) == length(flag_to_vertex) ||
        error("a nonvertical flag ray is missing from the node atlas")
    vertical = findall(==(0), flag_data.ray_matches)
    vertical_taxon = Dict(flag_ray => clean_vertical_taxon(flag_data, flag_ray)
                          for flag_ray in vertical)
    sort!(collect(values(vertical_taxon))) == collect(1:7) ||
        error("the vertical flag rays are not the seven taxon directions")
    return (
        records=ray_records,
        by_grass=by_grass,
        flag_to_vertex=flag_to_vertex,
        letters=letters,
        vertical=vertical,
        vertical_taxon=vertical_taxon,
    )
end
