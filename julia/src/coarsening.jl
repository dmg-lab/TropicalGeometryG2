module CoarsestFan

using Oscar
import ..rational_key

export LinkCoarseningResult,
       MergeRound,
       merge_cells,
       merge_cells_file,
       coarsen_fan,
       coarsen_fan_file,
       open_equal_span_walls,
       open_wall_certificate

struct IterativeCell
    original_cells::Vector{Int}
    original_extreme_rays::Vector{Int}
    volume_sum::Union{Nothing, QQFieldElem}
    hull_volume::Union{Nothing, QQFieldElem}
    span_block::Int
    boundary_facets::Vector{Vector{Int}}
end

"""One fixed affine cross-section for an invariant equal-span block."""
struct SpanBlockGeometry
    affine_coordinates::Matrix{QQFieldElem}
end

struct MergeRound
    codimension::Int
    cells_before::Int
    candidate_faces::Int
    dimension_candidates::Int
    spherical_links::Int
    convex_candidates::Int
    merged_stars::Int
    cells_after::Int
end

struct LinkCoarseningResult
    fan
    strategy::Symbol
    input_cones::Int
    input_rays::Int
    output_cones::Int
    output_rays::Int
    ambient_dimension::Int
    fan_dimension::Int
    lineality_dimension::Int
    rounds::Vector{MergeRound}
    original_cells::Vector{Vector{Int}}
    open_equal_span_walls::Int
    # Connected sets of cones across open walls whose union was checked not convex.
    checked_open_wall_unions::Int
    coarsest_certified::Bool
end

mutable struct DisjointSets
    parent::Vector{Int}
end

DisjointSets(n::Int) = DisjointSets(collect(1:n))

matrix_row(M, i::Int) = [M[i, j] for j in 1:ncols(M)]

function root!(sets::DisjointSets, x::Int)
    while sets.parent[x] != x
        sets.parent[x] = sets.parent[sets.parent[x]]
        x = sets.parent[x]
    end
    return x
end

function unite!(sets::DisjointSets, a::Int, b::Int)
    a, b = root!(sets, a), root!(sets, b)
    a == b || (sets.parent[b] = a)
end

function canonical_ray_signature(row, quotient_basis)
    q = vec(reshape(row, 1, :) * quotient_basis)
    pivot = findfirst(!iszero, q)
    pivot === nothing && error("a purported ray lies in lineality")
    return Tuple(q ./ abs(q[pivot]))
end

function signature_map(rays, lineality)
    _, quotient_basis_raw = nullspace(matrix(QQ, lineality))
    quotient_basis = Matrix{QQFieldElem}(quotient_basis_raw)
    result = Dict{Tuple{Vararg{QQFieldElem}}, Int}()
    for i in 1:nrows(rays)
        key = canonical_ray_signature(matrix_row(rays, i), quotient_basis)
        haskey(result, key) && error("duplicate quotient rays $(result[key]) and $i")
        result[key] = i
    end
    return result, quotient_basis
end

function add_combinations!(owners, ray_set, cell::Int, size::Int)
    length(ray_set) >= size || return
    positions = Vector{Int}(undef, size)
    function visit(depth::Int, first_position::Int)
        if depth > size
            key = Tuple(ray_set[positions])
            push!(get!(owners, key, Int[]), cell)
            return
        end
        final_position = length(ray_set) - (size - depth)
        for position in first_position:final_position
            positions[depth] = position
            visit(depth + 1, position + 1)
        end
    end
    visit(1, 1)
end

"""Candidate adjacent pairs obtained from the rays of a quotient facet."""
function facet_candidate_pairs(ray_sets, quotient_dimension::Int)
    facet_ray_minimum = quotient_dimension - 1
    facet_ray_minimum >= 1 || error("the pointed quotient must have dimension at least two")
    owners = Dict{Tuple{Vararg{Int}}, Vector{Int}}()
    for (cell, ray_set) in enumerate(ray_sets)
        add_combinations!(owners, ray_set, cell, facet_ray_minimum)
    end
    pairs = Set{Tuple{Int, Int}}()
    for cells in values(owners)
        sort!(unique!(cells))
        for a in 1:(length(cells) - 1), b in (a + 1):length(cells)
            i, j = cells[a], cells[b]
            push!(pairs, i < j ? (i, j) : (j, i))
        end
    end
    return sort!(collect(pairs))
end

function facet_adjacencies(ray_sets, rays, lineality, quotient_dimension)
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    projected = rays * quotient
    facets_by_rays = Dict{Tuple{Vararg{Int}}, Set{Int}}()
    for (i, j) in facet_candidate_pairs(ray_sets, quotient_dimension)
        common = sort!(intersect(ray_sets[i], ray_sets[j]))
        rank(matrix(QQ, projected[common, :])) == quotient_dimension - 1 || continue
        owners = get!(facets_by_rays, Tuple(common), Set{Int}())
        push!(owners, i, j)
    end
    return Dict(key => sort!(collect(owners)) for (key, owners) in facets_by_rays)
end

"""Interior codimension-two faces obtained from pairs of shared facets."""
function codimension_two_faces(maximal_ray_sets, facets, rays, lineality, quotient_dimension)
    quotient_dimension == 2 && return [Int[]]
    quotient_dimension == 3 && return [[ray] for ray in 1:nrows(rays)]
    quotient_dimension > 4 && return lower_shared_faces(
        maximal_ray_sets, facets, rays, lineality, quotient_dimension, 2,
    )
    facet_keys = sort!(collect(keys(facets)))
    facet_ray_sets = [collect(facet) for facet in facet_keys]
    facets_by_ray = ray_to_cells(facet_ray_sets, nrows(rays))

    # Every pointed two-dimensional cone has exactly two extreme rays.  It is
    # therefore much cheaper to index ray pairs and verify that the shared
    # facets containing them intersect in exactly that pair than to enumerate
    # all pairs of facets in every maximal cone.
    candidate_pairs = Set{Tuple{Int, Int}}()
    for ray_set in maximal_ray_sets
        for left in 1:(length(ray_set) - 1), right in (left + 1):length(ray_set)
            push!(candidate_pairs, (ray_set[left], ray_set[right]))
        end
    end
    result = Set{Tuple{Vararg{Int}}}()
    for pair in candidate_pairs
        containing_facets = intersect(facets_by_ray[pair[1]], facets_by_ray[pair[2]])
        length(containing_facets) >= 2 || continue
        common_rays = copy(facet_ray_sets[first(containing_facets)])
        for facet_number in Iterators.drop(containing_facets, 1)
            common_rays = intersect(common_rays, facet_ray_sets[facet_number])
            length(common_rays) < 2 && break
        end
        Tuple(common_rays) == pair && push!(result, pair)
    end
    return [collect(face) for face in sort!(collect(result))]
end

"""Index lower interior faces by exact intersections of shared facets in a cell."""
function lower_shared_faces(maximal_ray_sets, facets, rays, lineality, quotient_dimension, codimension)
    codimension >= 1 || throw(ArgumentError("face codimension must be positive"))
    codimension < quotient_dimension || throw(ArgumentError("face must contain a quotient ray"))
    facet_keys = sort!(collect(keys(facets)))
    codimension == 1 && return collect.(facet_keys)
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    projected = rays * quotient
    facets_by_cell = [Vector{Vector{Int}}() for _ in maximal_ray_sets]
    for facet in facet_keys, cell in facets[facet]
        push!(facets_by_cell[cell], collect(facet))
    end
    dimensions = Dict{Tuple{Vararg{Int}}, Int}()
    result = Set{Tuple{Vararg{Int}}}()
    for cell_facets in facets_by_cell
        current = Set(Tuple.(cell_facets))
        for level in 2:codimension
            next = Set{Tuple{Vararg{Int}}}()
            for face in current, facet in cell_facets
                common = Tuple(sort!(intersect(collect(face), facet)))
                length(common) >= quotient_dimension - level || continue
                dimension = get!(dimensions, common) do
                    rank(matrix(QQ, projected[collect(common), :]))
                end
                dimension == quotient_dimension - level && push!(next, common)
            end
            current = next
            isempty(current) && break
        end
        union!(result, current)
    end
    return collect.(sort!(collect(result)))
end

function equal_span_components(ray_sets, facets, rays, lineality, quotient_dimension)
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    projected = rays * quotient
    span_keys = map(ray_sets) do ray_set
        r, reduced = rref(matrix(QQ, projected[ray_set, :]))
        (r, Tuple(vec(Matrix{QQFieldElem}(reduced)[1:r, :])))
    end
    all(key -> first(key) == quotient_dimension, span_keys) ||
        error("input fan is not pure modulo its common lineality space")
    sets = DisjointSets(length(ray_sets))
    for owners in values(facets)
        for a in 1:(length(owners) - 1), b in (a + 1):length(owners)
            i, j = owners[a], owners[b]
            span_keys[i] == span_keys[j] || continue
            unite!(sets, i, j)
        end
    end
    groups = Dict{Int, Vector{Int}}()
    for i in eachindex(ray_sets)
        push!(get!(groups, root!(sets, i), Int[]), i)
    end
    return sort!(collect(values(groups)); by=first)
end

function independent_coordinate_columns(rows, dimension::Int)
    pivots = Int[]
    for column in 1:ncols(rows)
        rank(matrix(QQ, rows[:, vcat(pivots, column)])) > length(pivots) &&
            push!(pivots, column)
        length(pivots) == dimension && return pivots
    end
    error("could not find $dimension independent quotient coordinates")
end

cross2(a, b, c) =
    (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])

function convex_polygon_order(coordinates)
    order = sort!(collect(1:nrows(coordinates)); by=i -> (
        coordinates[i, 1],
        coordinates[i, 2],
    ))
    length(order) <= 2 && return order
    lower = Int[]
    for point in order
        while length(lower) >= 2 && cross2(
            matrix_row(coordinates, lower[end - 1]),
            matrix_row(coordinates, lower[end]),
            matrix_row(coordinates, point),
        ) <= 0
            pop!(lower)
        end
        push!(lower, point)
    end
    upper = Int[]
    for point in reverse(order)
        while length(upper) >= 2 && cross2(
            matrix_row(coordinates, upper[end - 1]),
            matrix_row(coordinates, upper[end]),
            matrix_row(coordinates, point),
        ) <= 0
            pop!(upper)
        end
        push!(upper, point)
    end
    return vcat(lower[1:(end - 1)], upper[1:(end - 1)])
end

"""Exact hull data, using the lightweight implementation through dimension three."""
function exact_affine_hull(coordinates)
    point_count = nrows(coordinates)
    dimension = ncols(coordinates)
    point_count >= dimension + 1 || return (
        volume=QQ(0), vertices=Int[], facets=Vector{Vector{Int}}(),
    )
    if dimension == 1
        minimum_index = argmin(coordinates[:, 1])
        maximum_index = argmax(coordinates[:, 1])
        return (
            volume=coordinates[maximum_index, 1] - coordinates[minimum_index, 1],
            vertices=sort!(unique([minimum_index, maximum_index])),
            facets=[[minimum_index], [maximum_index]],
        )
    end
    if dimension == 2
        polygon = convex_polygon_order(coordinates)
        twice_area = sum((
            coordinates[polygon[i], 1] * coordinates[polygon[mod1(i + 1, length(polygon))], 2] -
            coordinates[polygon[i], 2] * coordinates[polygon[mod1(i + 1, length(polygon))], 1]
            for i in eachindex(polygon)
        ); init=QQ(0))
        polygon_facets = [
            sort!([polygon[i], polygon[mod1(i + 1, length(polygon))]])
            for i in eachindex(polygon)
        ]
        return (
            volume=abs(twice_area) / 2,
            vertices=sort!(copy(polygon)),
            facets=polygon_facets,
        )
    end
    if dimension > 3
        if point_count == dimension + 1
            differences = matrix(QQ, [
                coordinates[i, j] - coordinates[1, j]
                for i in 2:point_count, j in 1:dimension
            ])
            simplex_volume = abs(det(differences)) / factorial(big(dimension))
            iszero(simplex_volume) && return (
                volume=QQ(0), vertices=Int[], facets=Vector{Vector{Int}}(),
            )
            return (
                volume=simplex_volume,
                vertices=collect(1:point_count),
                facets=[setdiff(1:point_count, [i]) for i in 1:point_count],
            )
        end
        hull = convex_hull(matrix(QQ, coordinates))
        input_index = Dict(Tuple(coordinates[i, :]) => i for i in 1:point_count)
        hull_vertices = sort!([
            input_index[Tuple(vertex)] for vertex in vertices(hull)
        ])
        hull_facets = Vector{Vector{Int}}()
        for (normal, offset) in facets(Pair, hull)
            push!(hull_facets, [
                i for i in 1:point_count if
                sum(normal[j] * coordinates[i, j] for j in 1:dimension; init=QQ(0)) == offset
            ])
        end
        return (volume=QQ(volume(hull)), vertices=hull_vertices, facets=hull_facets)
    end

    facets_by_vertices = Dict{Tuple{Vararg{Int}}, Vector{QQFieldElem}}()
    for i in 1:(point_count - 2), j in (i + 1):(point_count - 1), k in (j + 1):point_count
        u = [coordinates[j, column] - coordinates[i, column] for column in 1:3]
        v = [coordinates[k, column] - coordinates[i, column] for column in 1:3]
        normal = [
            u[2] * v[3] - u[3] * v[2],
            u[3] * v[1] - u[1] * v[3],
            u[1] * v[2] - u[2] * v[1],
        ]
        all(iszero, normal) && continue
        values = [
            sum(normal[column] * (coordinates[point, column] - coordinates[i, column])
                for column in 1:3; init=QQ(0))
            for point in 1:point_count
        ]
        any(!iszero, values) || continue
        (all(value -> value >= 0, values) || all(value -> value <= 0, values)) ||
            continue
        get!(facets_by_vertices, Tuple(findall(iszero, values)), normal)
    end

    center = [
        sum(coordinates[:, column]; init=QQ(0)) / point_count for column in 1:3
    ]
    result = QQ(0)
    hull_vertices = Set{Int}()
    for (facet, normal) in facets_by_vertices
        facet_vertices = collect(facet)
        drop::Int = findfirst(!iszero, normal)
        kept = setdiff(1:3, [drop])
        projected = coordinates[facet_vertices, kept]
        boundary = facet_vertices[convex_polygon_order(projected)]
        union!(hull_vertices, boundary)
        anchor = boundary[1]
        for position in 2:(length(boundary) - 1)
            tetrahedron = matrix(QQ, [
                coordinates[vertex, column] - center[column]
                for vertex in [anchor, boundary[position], boundary[position + 1]],
                    column in 1:3
            ])
            result += abs(det(tetrahedron)) / 6
        end
    end
    return (
        volume=result,
        vertices=sort!(collect(hull_vertices)),
        facets=[collect(facet) for facet in keys(facets_by_vertices)],
    )
end

exact_affine_volume(coordinates) = exact_affine_hull(coordinates).volume

"""Find `g` with every row of `coordinates` evaluating positively on `g`."""
function positive_grading(coordinates)
    dimension = ncols(coordinates)
    values = coordinates * reshape(vec(sum(coordinates; dims=1)), :, 1)
    all(>(0), values) && return vec(sum(coordinates; dims=1))

    # Searching combinations of tight rows grows rapidly in higher dimensions.
    # The negative sum of the exact facet normals grades every nonzero ray of
    # a pointed full-dimensional cone positively.
    if dimension > 4
        cone = positive_hull(matrix(QQ, coordinates))
        lineality_dim(cone) == 0 || return nothing
        grading = zeros(QQFieldElem, dimension)
        for halfspace in facets(Halfspace, cone)
            grading .-= collect(normal_vector(halfspace))
        end
        heights = coordinates * reshape(grading, :, 1)
        return all(>(0), heights) ? grading : nothing
    end

    # The nonempty polyhedron A*g >= 1 has a vertex because A has full column
    # rank.  A vertex is determined by `dimension` independent tight rows.
    chosen = Vector{Int}(undef, dimension)
    result = Ref{Union{Nothing, Vector{QQFieldElem}}}(nothing)
    right_hand_side = matrix(QQ, ones(Int, dimension, 1))
    function visit(depth::Int, first_row::Int)
        if depth > dimension
            basis = matrix(QQ, coordinates[chosen, :])
            rank(basis) == dimension || return false
            grading = vec(Matrix{QQFieldElem}(inv(basis) * right_hand_side))
            heights = coordinates * reshape(grading, :, 1)
            all(height -> height >= 1, heights) || return false
            result[] = grading
            return true
        end
        final_row = nrows(coordinates) - (dimension - depth)
        for row_number in first_row:final_row
            chosen[depth] = row_number
            visit(depth + 1, row_number + 1) && return true
        end
        return false
    end
    visit(1, 1)
    return result[]
end

"""Certify exactly that a component union equals its convex hull."""
function convex_component_certificate(
    component,
    ray_sets,
    rays,
    lineality,
    quotient_dimension::Int,
)
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    union_rays = sort!(unique(vcat((ray_sets[cell] for cell in component)...)))
    projected = rays[union_rays, :] * quotient
    pointed_dimension = rank(matrix(QQ, projected))
    pointed_dimension == quotient_dimension || error(
        "component $(first(component)) has quotient dimension $pointed_dimension " *
        "instead of $quotient_dimension",
    )
    pivots = independent_coordinate_columns(projected, pointed_dimension)
    coordinates = projected[:, pivots]
    hull_cone = positive_hull(coordinates)
    lineality_dim(hull_cone) == 0 || error("component convex hull is not pointed")

    grading = zeros(QQFieldElem, pointed_dimension)
    for halfspace in facets(Halfspace, hull_cone)
        grading .-= collect(normal_vector(halfspace))
    end
    heights = coordinates * reshape(grading, :, 1)
    all(>(0), heights) || error("failed to find a positive cross-section grading")

    normalized = Matrix{QQFieldElem}(undef, length(union_rays), pointed_dimension)
    for i in eachindex(union_rays)
        normalized[i, :] = coordinates[i, :] ./ heights[i, 1]
    end
    drop = findfirst(!iszero, grading)
    affine_coordinates = normalized[:, setdiff(1:pointed_dimension, [drop])]
    local_index = Dict(ray => i for (i, ray) in enumerate(union_rays))
    hull_volume = exact_affine_volume(affine_coordinates)
    volume_sum = sum((
        exact_affine_volume(affine_coordinates[
            [local_index[ray] for ray in ray_sets[cell]],
            :,
        ]) for cell in component
    ); init=QQ(0))
    return (
        convex=volume_sum == hull_volume,
        volume_sum=volume_sum,
        hull_volume=hull_volume,
        union_rays=union_rays,
    )
end

function extreme_original_rays(cone, allowed, signature_to_ray, quotient_basis)
    result = Int[]
    cone_rays = matrix(QQ, rays_modulo_lineality(cone).rays_modulo_lineality)
    allowed_set = Set(allowed)
    for i in 1:nrows(cone_rays)
        key = canonical_ray_signature(matrix_row(cone_rays, i), quotient_basis)
        ray = get(signature_to_ray, key, 0)
        ray == 0 && error("coarsening created an unknown ray")
        ray in allowed_set || error("coarsening created a ray outside its component")
        push!(result, ray)
    end
    return sort!(unique!(result))
end

function validated_output_fan(cells, rays, lineality)
    used_original_rays = sort!(unique(vcat(
        (cell.original_extreme_rays for cell in cells)...,
    )))
    compact_index = Dict(ray => i for (i, ray) in enumerate(used_original_rays))
    incidence = IncidenceMatrix([
        sort!([compact_index[ray] for ray in cell.original_extreme_rays]) for cell in cells
    ])
    fan = polyhedral_fan(
        incidence,
        matrix(QQ, rays[used_original_rays, :]),
        matrix(QQ, lineality);
        non_redundant=true,
    )
    n_maximal_cones(fan) == length(cells) || error("fan construction changed cell count")
    return fan
end

function output_order(fan, cells, signature_to_ray, quotient_basis)
    output_rays = matrix(QQ, rays_modulo_lineality(fan).rays_modulo_lineality)
    output_to_original = Int[]
    for i in 1:nrows(output_rays)
        key = canonical_ray_signature(matrix_row(output_rays, i), quotient_basis)
        push!(output_to_original, get(signature_to_ray, key, 0))
        last(output_to_original) > 0 || error("unrecognized output ray")
    end
    by_ray_set = Dict(
        Tuple(cell.original_extreme_rays) => index for (index, cell) in enumerate(cells)
    )
    incidence = maximal_cones(IncidenceMatrix, fan)
    order = Int[]
    for i in 1:nrows(incidence)
        original_rays = sort!([output_to_original[j] for j in row(incidence, i)])
        state = get(by_ray_set, Tuple(original_rays), 0)
        state > 0 || error("could not match output cone $i to a span component")
        push!(order, state)
    end
    length(unique(order)) == length(cells) || error("output order is not bijective")
    return order, output_to_original
end

rational_string(x) = isnothing(x) ? "" : rational_key(x)

quoted_list(xs) = "\"" * join(xs, ';') * "\""

function write_provenance(path, order, cells)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "coarsest_cell,original_cells,original_cell_count,extremal_original_rays," *
                    "ray_count,cross_section_volume_sum,cross_section_hull_volume")
        for (output_cell, state) in enumerate(order)
            cell = cells[state]
            println(io, join([
                string(output_cell),
                quoted_list(cell.original_cells),
                string(length(cell.original_cells)),
                quoted_list(cell.original_extreme_rays),
                string(length(cell.original_extreme_rays)),
                rational_string(cell.volume_sum),
                rational_string(cell.hull_volume),
            ], ','))
        end
    end
end

face_contained_in(small, large) = all(ray -> ray in large, small)

function connected_graph(adjacency)
    isempty(adjacency) && return false
    seen = Set{Int}()
    stack = [first(eachindex(adjacency))]
    while !isempty(stack)
        vertex = pop!(stack)
        vertex in seen && continue
        push!(seen, vertex)
        append!(stack, adjacency[vertex])
    end
    return length(seen) == length(adjacency)
end

function ray_to_cells(ray_sets, number_of_rays::Int)
    result = [Int[] for _ in 1:number_of_rays]
    for (cell, rays) in enumerate(ray_sets), ray in rays
        push!(result[ray], cell)
    end
    return result
end

function owners_of_face(face, ray_owners, number_of_cells::Int)
    isempty(face) && return collect(1:number_of_cells)
    owners = copy(ray_owners[first(face)])
    for ray in Iterators.drop(face, 1)
        owners = intersect(owners, ray_owners[ray])
        isempty(owners) && break
    end
    return sort!(owners)
end

function cofaces_containing(face, cofaces, ray_to_cofaces)
    indices = owners_of_face(face, ray_to_cofaces, length(cofaces))
    return cofaces[indices]
end

function circle_link(
    face,
    owners,
    codimension_one_faces,
    codimension_one_index,
    maximal_ray_sets,
)
    vertices = cofaces_containing(face, codimension_one_faces, codimension_one_index)
    vertices = [
        vertex for vertex in vertices if
        any(cell -> face_contained_in(vertex, maximal_ray_sets[cell]), owners)
    ]
    length(vertices) >= 2 || return false
    adjacency = [Int[] for _ in eachindex(vertices)]
    degrees = zeros(Int, length(vertices))
    for cell in owners
        endpoints = [
            i for i in eachindex(vertices) if
            face_contained_in(vertices[i], maximal_ray_sets[cell])
        ]
        length(endpoints) == 2 || return false
        left, right = endpoints
        left == right && return false
        push!(adjacency[left], right)
        push!(adjacency[right], left)
        degrees[left] += 1
        degrees[right] += 1
    end
    all(==(2), degrees) || return false
    return connected_graph(adjacency)
end

function two_sphere_link(
    face,
    owners,
    vertex_faces,
    edge_faces,
    vertex_index,
    edge_index,
    maximal_ray_sets,
)
    vertices = cofaces_containing(face, vertex_faces, vertex_index)
    edges = cofaces_containing(face, edge_faces, edge_index)
    vertices = [
        vertex for vertex in vertices if
        any(cell -> face_contained_in(vertex, maximal_ray_sets[cell]), owners)
    ]
    edges = [
        edge for edge in edges if
        any(cell -> face_contained_in(edge, maximal_ray_sets[cell]), owners)
    ]
    isempty(vertices) && return false
    isempty(edges) && return false

    edge_endpoints = Vector{Tuple{Int, Int}}(undef, length(edges))
    edge_cells = Vector{Tuple{Int, Int}}(undef, length(edges))
    cell_adjacency = [Int[] for _ in owners]
    owner_position = Dict(cell => i for (i, cell) in enumerate(owners))
    for (edge_number, edge) in enumerate(edges)
        endpoints = [i for i in eachindex(vertices) if face_contained_in(vertices[i], edge)]
        length(endpoints) == 2 || return false
        edge_endpoints[edge_number] = (endpoints[1], endpoints[2])

        incident_cells = [
            cell for cell in owners if face_contained_in(edge, maximal_ray_sets[cell])
        ]
        length(incident_cells) == 2 || return false
        edge_cells[edge_number] = (incident_cells[1], incident_cells[2])
        left = owner_position[incident_cells[1]]
        right = owner_position[incident_cells[2]]
        push!(cell_adjacency[left], right)
        push!(cell_adjacency[right], left)
    end
    connected_graph(cell_adjacency) || return false
    length(vertices) - length(edges) + length(owners) == 2 || return false

    # A closed connected Euler-two pseudomanifold can still be pinched at a
    # vertex. Require every vertex link to be a connected 2-regular graph.
    for vertex in eachindex(vertices)
        incident_cells = [
            cell for cell in owners if
            face_contained_in(vertices[vertex], maximal_ray_sets[cell])
        ]
        cell_position = Dict(cell => i for (i, cell) in enumerate(incident_cells))
        adjacency = [Int[] for _ in incident_cells]
        degrees = zeros(Int, length(incident_cells))
        for edge_number in eachindex(edges)
            vertex in edge_endpoints[edge_number] || continue
            left_cell, right_cell = edge_cells[edge_number]
            haskey(cell_position, left_cell) && haskey(cell_position, right_cell) ||
                return false
            left = cell_position[left_cell]
            right = cell_position[right_cell]
            push!(adjacency[left], right)
            push!(adjacency[right], left)
            degrees[left] += 1
            degrees[right] += 1
        end
        all(==(2), degrees) || return false
        connected_graph(adjacency) || return false
    end
    return true
end

function spherical_link(
    face,
    codimension::Int,
    owners,
    faces_by_dimension,
    face_indices,
    face_dimension::Int,
    maximal_ray_sets,
)
    codimension == 1 && return length(owners) == 2
    codimension == 2 && return circle_link(
        face,
        owners,
        faces_by_dimension[face_dimension + 1],
        face_indices[face_dimension + 1],
        maximal_ray_sets,
    )
    codimension == 3 && return two_sphere_link(
        face,
        owners,
        faces_by_dimension[face_dimension + 1],
        faces_by_dimension[face_dimension + 2],
        face_indices[face_dimension + 1],
        face_indices[face_dimension + 2],
        maximal_ray_sets,
    )
    return false
end

function cached_block_geometry(
    component,
    ray_sets,
    rays,
    quotient,
    quotient_dimension::Int,
)
    length(component) > 1 || return (
        nothing,
        Dict{Int, QQFieldElem}(),
        Dict{Int, Vector{Vector{Int}}}(),
    )
    union_rays = sort!(unique(vcat((ray_sets[cell] for cell in component)...)))
    projected = rays[union_rays, :] * quotient
    rank(matrix(QQ, projected)) == quotient_dimension || return (
        nothing,
        Dict{Int, QQFieldElem}(),
        Dict{Int, Vector{Vector{Int}}}(),
    )
    pivots = independent_coordinate_columns(projected, quotient_dimension)
    coordinates = projected[:, pivots]
    grading = positive_grading(coordinates)
    isnothing(grading) && return (
        nothing,
        Dict{Int, QQFieldElem}(),
        Dict{Int, Vector{Vector{Int}}}(),
    )
    heights = coordinates * reshape(grading, :, 1)
    all(>(0), heights) || return (
        nothing,
        Dict{Int, QQFieldElem}(),
        Dict{Int, Vector{Vector{Int}}}(),
    )
    drop = findfirst(!iszero, grading)
    kept = setdiff(1:quotient_dimension, [drop])
    affine = Matrix{QQFieldElem}(undef, nrows(rays), quotient_dimension - 1)
    fill!(affine, QQ(0))
    for (local_ray, input_ray) in enumerate(union_rays)
        affine[input_ray, :] = coordinates[local_ray, kept] ./ heights[local_ray, 1]
    end
    geometry = SpanBlockGeometry(affine)
    volumes = Dict{Int, QQFieldElem}()
    cell_facets = Dict{Int, Vector{Vector{Int}}}()
    for cell in component
        hull = exact_affine_hull(affine[ray_sets[cell], :])
        volumes[cell] = hull.volume
        cell_facets[cell] = [
            sort!([ray_sets[cell][local_ray] for local_ray in facet])
            for facet in hull.facets
        ]
    end
    return geometry, volumes, cell_facets
end

"""
Precompute the invariants that do not change under any accepted merge.

Every convex full-dimensional union is connected through codimension-one
faces, so it lies in one equal-span facet component.  That component ID is a
constant-time dimension obstruction.  When the whole component is pointed,
one rational cross-section serves every later convexity test in that block.
"""
function initialize_iterative_cache(
    input_ray_sets,
    input_rays,
    lineality,
    quotient_dimension::Int,
)
    initial_facets = facet_adjacencies(
        input_ray_sets,
        input_rays,
        lineality,
        quotient_dimension,
    )
    components = equal_span_components(
        input_ray_sets,
        initial_facets,
        input_rays,
        lineality,
        quotient_dimension,
    )
    block_of_cell = zeros(Int, length(input_ray_sets))
    for (block, component) in enumerate(components), cell in component
        block_of_cell[cell] = block
    end

    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    geometries = Vector{Union{Nothing, SpanBlockGeometry}}(undef, length(components))
    fill!(geometries, nothing)
    initial_volumes = Vector{Union{Nothing, QQFieldElem}}(undef, length(input_ray_sets))
    fill!(initial_volumes, nothing)
    initial_facets_by_cell = [Vector{Vector{Int}}() for _ in input_ray_sets]
    for (block, component) in enumerate(components)
        geometry, volumes, cell_facets = cached_block_geometry(
            component,
            input_ray_sets,
            input_rays,
            quotient,
            quotient_dimension,
        )
        geometries[block] = geometry
        for (cell, cell_volume) in volumes
            initial_volumes[cell] = cell_volume
            initial_facets_by_cell[cell] = cell_facets[cell]
        end
    end
    cells = [
        IterativeCell(
            [i],
            copy(input_ray_sets[i]),
            initial_volumes[i],
            initial_volumes[i],
            block_of_cell[i],
            initial_facets_by_cell[i],
        ) for i in eachindex(input_ray_sets)
    ]
    return cells, geometries, initial_facets, components, initial_facets_by_cell
end

function affine_hyperplane_values(coordinates, facet_positions)
    dimension = ncols(coordinates)
    anchor = first(facet_positions)
    if dimension == 1
        return [coordinates[i, 1] - coordinates[anchor, 1] for i in 1:nrows(coordinates)],
            QQFieldElem[1]
    end
    if dimension == 2
        second = facet_positions[2]
        direction = [
            coordinates[second, 1] - coordinates[anchor, 1],
            coordinates[second, 2] - coordinates[anchor, 2],
        ]
        normal = QQFieldElem[direction[2], -direction[1]]
        values = [
            sum(normal[column] * (coordinates[i, column] - coordinates[anchor, column])
                for column in 1:2; init=QQ(0))
            for i in 1:nrows(coordinates)
        ]
        return values, normal
    end
    if dimension > 3
        differences = matrix(QQ, [
            coordinates[position, column] - coordinates[anchor, column]
            for position in facet_positions[2:end], column in 1:dimension
        ])
        normal_space = kernel(differences; side=:right)
        ncols(normal_space) == 1 || error("a cached boundary facet is not a hyperplane")
        normal = vec(Matrix{QQFieldElem}(normal_space))
        values = [
            sum(normal[column] * (coordinates[i, column] - coordinates[anchor, column])
                for column in 1:dimension; init=QQ(0))
            for i in 1:nrows(coordinates)
        ]
        return values, normal
    end
    normal = nothing
    for left in 2:(length(facet_positions) - 1), right in (left + 1):length(facet_positions)
        b, c = facet_positions[left], facet_positions[right]
        u = [coordinates[b, column] - coordinates[anchor, column] for column in 1:3]
        v = [coordinates[c, column] - coordinates[anchor, column] for column in 1:3]
        candidate = QQFieldElem[
            u[2] * v[3] - u[3] * v[2],
            u[3] * v[1] - u[1] * v[3],
            u[1] * v[2] - u[2] * v[1],
        ]
        if any(!iszero, candidate)
            normal = candidate
            break
        end
    end
    isnothing(normal) && error("a cached boundary facet is not a hyperplane")
    values = [
        sum(normal[column] * (coordinates[i, column] - coordinates[anchor, column])
            for column in 1:3; init=QQ(0))
        for i in 1:nrows(coordinates)
    ]
    return values, normal
end

function supporting_boundary_certificate(boundary_facets, volume_sum, geometry, union_rays)
    isempty(boundary_facets) && return nothing

    coordinates = geometry.affine_coordinates[union_rays, :]
    local_index = Dict(ray => i for (i, ray) in enumerate(union_rays))
    supporting_planes = Dict{Tuple{Vararg{Int}}, Vector{QQFieldElem}}()
    for facet in boundary_facets
        positions = [local_index[ray] for ray in facet]
        values, normal = affine_hyperplane_values(coordinates, positions)
        any(>(0), values) && any(<(0), values) && return nothing
        zero_positions = Tuple(findall(iszero, values))
        supporting_planes[zero_positions] = normal
    end

    merged_facets = Vector{Vector{Int}}()
    extreme_rays = Set{Int}()
    dimension = ncols(coordinates)
    for (zero_positions_tuple, normal) in supporting_planes
        zero_positions = collect(zero_positions_tuple)
        facet_extremes = if dimension == 1
            [first(zero_positions)]
        elseif dimension == 2
            direction_column::Int = findfirst(column -> begin
                values = coordinates[zero_positions, column]
                maximum(values) != minimum(values)
            end, 1:2)
            values = coordinates[zero_positions, direction_column]
            [zero_positions[argmin(values)], zero_positions[argmax(values)]]
        elseif dimension == 3
            drop::Int = findfirst(!iszero, normal)
            kept = setdiff(1:3, [drop])
            zero_positions[convex_polygon_order(coordinates[zero_positions, kept])]
        else
            drop = findfirst(!iszero, normal)
            kept = setdiff(1:dimension, [drop])
            zero_positions[exact_affine_hull(coordinates[zero_positions, kept]).vertices]
        end
        facet_rays = sort!(unique(union_rays[facet_extremes]))
        push!(merged_facets, facet_rays)
        union!(extreme_rays, facet_rays)
    end
    isempty(extreme_rays) && return nothing
    return (
        convex=true,
        volume_sum=volume_sum,
        hull_volume=volume_sum,
        union_rays=union_rays,
        extreme_rays=sort!(collect(extreme_rays)),
        boundary_facets=sort!(unique!(merged_facets); by=Tuple),
    )
end

"""Certify a spherical star convex by checking every boundary facet supports it."""
function boundary_convexity_certificate(owners, cells, geometry, union_rays)
    facet_counts = Dict{Tuple{Vararg{Int}}, Int}()
    for cell in owners, facet in cells[cell].boundary_facets
        key = Tuple(facet)
        facet_counts[key] = get(facet_counts, key, 0) + 1
    end
    boundary_facets = [collect(facet) for (facet, count) in facet_counts if count == 1]
    volume_sum = sum((cells[cell].hull_volume for cell in owners); init=QQ(0))
    return supporting_boundary_certificate(
        boundary_facets,
        volume_sum,
        geometry,
        union_rays,
    )
end

"""Certify a complete span block, using a hull volume only as a fallback."""
function span_block_certificate(
    owners,
    cells,
    geometry,
    original_ray_sets,
    original_boundary_facets,
)
    isnothing(geometry) && return nothing
    any(cell -> isnothing(cells[cell].hull_volume), owners) && return nothing
    original_cells = sort!(unique(vcat(
        (cells[cell].original_cells for cell in owners)...,
    )))
    union_rays = sort!(unique(vcat(
        (original_ray_sets[cell] for cell in original_cells)...,
    )))
    facet_counts = Dict{Tuple{Vararg{Int}}, Int}()
    for cell in original_cells, facet in original_boundary_facets[cell]
        key = Tuple(facet)
        facet_counts[key] = get(facet_counts, key, 0) + 1
    end
    boundary_facets = [collect(facet) for (facet, count) in facet_counts if count == 1]
    volume_sum = sum((cells[cell].hull_volume for cell in owners); init=QQ(0))
    boundary_certificate = supporting_boundary_certificate(
        boundary_facets,
        volume_sum,
        geometry,
        union_rays,
    )
    isnothing(boundary_certificate) || return boundary_certificate

    println("      boundary certificate inconclusive; using volume fallback")
    hull = exact_affine_hull(geometry.affine_coordinates[union_rays, :])
    volume_sum == hull.volume || return nothing
    return (
        convex=true,
        volume_sum=volume_sum,
        hull_volume=hull.volume,
        union_rays=union_rays,
        extreme_rays=sort!(union_rays[hull.vertices]),
        boundary_facets=sort!(unique!([
            sort!(union_rays[facet]) for facet in hull.facets
        ]); by=Tuple),
    )
end

function candidate_certificate(
    owners,
    cells,
    geometries,
    input_rays,
    lineality,
    quotient_dimension,
    certificate_cache;
    require_volume::Bool=false,
)
    block = cells[first(owners)].span_block
    union_rays = sort!(unique(vcat(
        (cells[cell].original_extreme_rays for cell in owners)...,
    )))
    original_cells = sort!(unique(vcat(
        (cells[cell].original_cells for cell in owners)...,
    )))
    key = (block, Tuple(original_cells), require_volume)
    haskey(certificate_cache, key) && return certificate_cache[key]
    geometry = geometries[block]
    certificate = if !isnothing(geometry) && require_volume
        hull = exact_affine_hull(geometry.affine_coordinates[union_rays, :])
        volume_sum = sum((cells[cell].hull_volume for cell in owners); init=QQ(0))
        volume_sum == hull.volume ? (
            convex=true,
            volume_sum=volume_sum,
            hull_volume=hull.volume,
            union_rays=union_rays,
            extreme_rays=sort!(union_rays[hull.vertices]),
            boundary_facets=sort!(unique!([
                sort!(union_rays[facet]) for facet in hull.facets
            ]); by=Tuple),
        ) : nothing
    elseif !isnothing(geometry)
        boundary_convexity_certificate(owners, cells, geometry, union_rays)
    else
        local_ray_sets = [cells[cell].original_extreme_rays for cell in owners]
        try
            fallback = convex_component_certificate(
                collect(eachindex(local_ray_sets)),
                local_ray_sets,
                input_rays,
                lineality,
                quotient_dimension,
            )
            fallback.convex ? fallback : nothing
        catch exception
            message = sprint(showerror, exception)
            if occursin("not pointed", message) ||
                    occursin("positive cross-section grading", message)
                nothing
            else
                rethrow()
            end
        end
    end
    certificate_cache[key] = certificate
    return certificate
end

function apply_merge_batch(
    selected,
    current_cells,
    input_rays,
    lineality,
    input_signature_to_ray,
    input_quotient_basis,
)
    consumed = Set{Int}()
    next_cells = IterativeCell[]
    for candidate in selected
        union!(consumed, candidate.owners)
        extremes = if hasproperty(candidate.certificate, :extreme_rays)
            candidate.certificate.extreme_rays
        else
            cone = positive_hull(
                input_rays[candidate.certificate.union_rays, :],
                lineality;
                non_redundant=false,
            )
            extreme_original_rays(
                cone,
                collect(1:nrows(input_rays)),
                input_signature_to_ray,
                input_quotient_basis,
            )
        end
        originals = sort!(unique(vcat(
            (current_cells[cell].original_cells for cell in candidate.owners)...,
        )))
        push!(next_cells, IterativeCell(
            originals,
            extremes,
            candidate.certificate.volume_sum,
            candidate.certificate.hull_volume,
            current_cells[first(candidate.owners)].span_block,
            hasproperty(candidate.certificate, :boundary_facets) ?
                candidate.certificate.boundary_facets : Vector{Vector{Int}}(),
        ))
    end
    for cell in eachindex(current_cells)
        cell in consumed || push!(next_cells, current_cells[cell])
    end
    sort!(next_cells; by=cell -> first(cell.original_cells))
    return next_cells
end

function one_combinatorial_round(
    cells,
    codimension::Int,
    fan_dimension::Int,
    quotient_dimension::Int,
    input_rays,
    lineality,
    input_signature_to_ray,
    input_quotient_basis,
    geometries,
    certificate_cache;
    cached_codimension_one_owners=nothing,
)
    face_dimension = fan_dimension - codimension
    maximal_ray_sets = [cell.original_extreme_rays for cell in cells]
    owners_by_ray = ray_to_cells(maximal_ray_sets, nrows(input_rays))
    faces_by_dimension = Dict{Int, Vector{Vector{Int}}}()
    codimension_one_owners = Dict{Tuple{Vararg{Int}}, Vector{Int}}()
    if codimension == 1
        codimension_one_owners = isnothing(cached_codimension_one_owners) ?
            facet_adjacencies(
                maximal_ray_sets,
                input_rays,
                lineality,
                quotient_dimension,
            ) : cached_codimension_one_owners
        faces_by_dimension[face_dimension] = [
            collect(face) for face in sort!(collect(keys(codimension_one_owners)))
        ]
    else
        codimension_one_owners = isnothing(cached_codimension_one_owners) ?
            facet_adjacencies(
                maximal_ray_sets,
                input_rays,
                lineality,
                quotient_dimension,
            ) : cached_codimension_one_owners
        faces_by_dimension[fan_dimension - 1] = [
            collect(face) for face in sort!(collect(keys(codimension_one_owners)))
        ]
        faces_by_dimension[fan_dimension - 2] = codimension_two_faces(
            maximal_ray_sets,
            codimension_one_owners,
            input_rays,
            lineality,
            quotient_dimension,
        )
        for level in 3:codimension
            faces_by_dimension[fan_dimension - level] = if level == quotient_dimension - 1
                [[ray] for ray in eachindex(owners_by_ray) if !isempty(owners_by_ray[ray])]
            else
                lower_shared_faces(
                    maximal_ray_sets, codimension_one_owners, input_rays,
                    lineality, quotient_dimension, level,
                )
            end
        end
    end
    face_indices = Dict(
        dimension => ray_to_cells(faces, nrows(input_rays))
        for (dimension, faces) in faces_by_dimension
    )

    candidate_faces = 0
    dimension_candidates = 0
    spherical_links = 0
    convex_candidates = 0
    spherical_candidates = NamedTuple[]
    seen_owner_sets = Set{Tuple{Vararg{Int}}}()
    for face in faces_by_dimension[face_dimension]
        global_owners = codimension == 1 ?
            codimension_one_owners[Tuple(face)] :
            owners_of_face(face, owners_by_ray, length(maximal_ray_sets))
        length(global_owners) > 1 || continue
        owners_by_block = Dict{Int, Vector{Int}}()
        for cell in global_owners
            push!(get!(owners_by_block, cells[cell].span_block, Int[]), cell)
        end
        for block in sort!(collect(keys(owners_by_block)))
            owners = owners_by_block[block]
            length(owners) > 1 || continue
            owner_key = Tuple(owners)
            owner_key in seen_owner_sets && continue
            push!(seen_owner_sets, owner_key)
            candidate_faces += 1
            dimension_candidates += 1
            if codimension <= 3
                spherical_link(
                    face,
                    codimension,
                    owners,
                    faces_by_dimension,
                    face_indices,
                    face_dimension,
                    maximal_ray_sets,
                ) || continue
                spherical_links += 1
            end
            push!(spherical_candidates, (face=copy(face), owners=owners))
        end
    end

    # Prefer the stars eliminating the most cells.  This deterministic greedy
    # packing reduces the number of later face-complex scans without spending
    # convexity tests on overlapping candidates in the same batch.
    sort!(spherical_candidates; by=candidate -> (
        -length(candidate.owners),
        Tuple(candidate.face),
    ))
    selected = NamedTuple[]
    used_cells = Set{Int}()
    for candidate in spherical_candidates
        owners = candidate.owners
        any(cell -> cell in used_cells, owners) && continue
        certificate = candidate_certificate(
            owners,
            cells,
            geometries,
            input_rays,
            lineality,
            quotient_dimension,
            certificate_cache;
            require_volume=codimension > 3,
        )
        certificate === nothing && continue
        convex_candidates += 1
        union!(used_cells, owners)
        push!(selected, (
            face=candidate.face,
            owners=owners,
            certificate=certificate,
        ))
    end

    if isempty(selected)
        round = MergeRound(
            codimension,
            length(cells),
            candidate_faces,
            dimension_candidates,
            spherical_links,
            convex_candidates,
            0,
            length(cells),
        )
        return cells, round, codimension_one_owners
    end

    next_cells = apply_merge_batch(
        selected,
        cells,
        input_rays,
        lineality,
        input_signature_to_ray,
        input_quotient_basis,
    )
    round = MergeRound(
        codimension,
        length(cells),
        candidate_faces,
        dimension_candidates,
        spherical_links,
        convex_candidates,
        length(selected),
        length(next_cells),
    )
    return next_cells, round, codimension_one_owners
end

"""
Merge each remaining invariant span block whose complete union is convex.

Spherical links remove interior subdivisions.  This closure also removes
subdivisions of a block boundary, whose relative links are balls rather than
spheres and therefore cannot be reached by the preceding scans alone.  An
exact supporting-boundary certificate is tried first; an exact rational
cross-section volume comparison is retained as a fallback.
"""
function one_span_block_round(
    cells,
    input_rays,
    lineality,
    input_signature_to_ray,
    input_quotient_basis,
    geometries,
    original_ray_sets,
    original_boundary_facets,
)
    owners_by_block = Dict{Int, Vector{Int}}()
    for cell in eachindex(cells)
        push!(get!(owners_by_block, cells[cell].span_block, Int[]), cell)
    end
    candidates = [
        owners_by_block[block] for block in sort!(collect(keys(owners_by_block)))
        if length(owners_by_block[block]) > 1
    ]
    selected = NamedTuple[]
    for (candidate_number, owners) in enumerate(candidates)
        println(
            "    certifying span block $candidate_number/$(length(candidates)) " *
            "with $(length(owners)) current cells",
        )
        block = cells[first(owners)].span_block
        certificate = span_block_certificate(
            owners,
            cells,
            geometries[block],
            original_ray_sets,
            original_boundary_facets,
        )
        certificate === nothing || push!(
            selected,
            (face=Int[], owners=owners, certificate=certificate),
        )
        candidate_number % 5 == 0 && GC.gc()
    end
    if isempty(selected)
        round = MergeRound(
            0,
            length(cells),
            length(candidates),
            length(candidates),
            0,
            0,
            0,
            length(cells),
        )
        return cells, round
    end
    next_cells = apply_merge_batch(
        selected,
        cells,
        input_rays,
        lineality,
        input_signature_to_ray,
        input_quotient_basis,
    )
    round = MergeRound(
        0,
        length(cells),
        length(candidates),
        length(candidates),
        0,
        length(selected),
        length(selected),
        length(next_cells),
    )
    return next_cells, round
end

"The open equal-span walls of `fan`, as pairs of maximal cones, with its rays and lineality."
function _open_equal_span_walls(fan)
    data = rays_modulo_lineality(fan)
    rays = Matrix{QQFieldElem}(matrix(QQ, data.rays_modulo_lineality))
    lineality = Matrix{QQFieldElem}(matrix(QQ, data.lineality_basis))
    projected = rays * Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    incidence = maximal_cones(IncidenceMatrix, fan)
    ray_sets = [sort!(collect(row(incidence, i))) for i in 1:nrows(incidence)]
    quotient_dimensions = [rank(matrix(QQ, projected[ray_set, :])) for ray_set in ray_sets]
    owners = ray_to_cells(ray_sets, size(rays, 1))
    walls = Tuple{Int, Int}[]
    for quotient_dimension in unique(quotient_dimensions)
        # Rays modulo lineality cannot merge without enlarging the lineality.
        quotient_dimension < 2 && continue
        members = findall(==(quotient_dimension), quotient_dimensions)
        member_sets = ray_sets[members]
        for (a, b) in facet_candidate_pairs(member_sets, quotient_dimension)
            common = intersect(member_sets[a], member_sets[b])
            rank(matrix(QQ, projected[common, :])) == quotient_dimension - 1 || continue
            rank(matrix(QQ, projected[union(member_sets[a], member_sets[b]), :])) ==
                quotient_dimension || continue
            i, j = members[a], members[b]
            any(k -> k != i && k != j && issubset(common, ray_sets[k]), owners[first(common)]) &&
                continue
            push!(walls, (i, j))
        end
    end
    return (walls=walls, rays=rays, lineality=lineality, ray_sets=ray_sets)
end

"""
    open_equal_span_walls(fan) -> Int

Count walls that might still be erased: two maximal cones of equal dimension
meet in a common facet, span the same linear space, and no further maximal cone
contains that facet. Any convex union of maximal cones has such a wall in its
relative interior, so zero certifies that no coarsening of `fan` exists. A
positive count is inconclusive: convexity and faces attached elsewhere can still
block a merge. `open_wall_certificate` decides the walls one by one.
"""
open_equal_span_walls(fan) = length(_open_equal_span_walls(fan).walls)

"All connected sets of at least two cones of a component of the open-wall graph, or nothing past `limit`."
function _connected_cone_sets(neighbours, component, limit)
    level = unique!([sort!([i, j]) for i in component for j in neighbours[i] if i < j])
    sets = copy(level)
    while !isempty(level)
        next = Set{Vector{Int}}()
        for cells in level, cell in cells, neighbour in neighbours[cell]
            neighbour in cells && continue
            push!(next, sort!([cells; neighbour]))
            length(sets) + length(next) > limit && return nothing
        end
        level = collect(next)
        append!(sets, level)
    end
    return sets
end

"""
Is the union of `cells`, connected through shared facets and of equal span,
convex? Exactly when every facet of one of them that no other one shares has all
of them on its side: such facets make up the boundary of the union, so the union
is then the intersection of their half-spaces.
"""
function _cone_union_is_convex(cells, ray_sets, facets_of)
    used = unique(vcat(ray_sets[cells]...))
    for cell in cells, facet in facets_of(cell)
        any(other -> other != cell && issubset(facet.rays, ray_sets[other]), cells) && continue
        all(facet.values[ray] <= 0 for ray in used) || return false
    end
    return true
end

"""
    open_wall_certificate(fan; max_unions=100_000)

Decide whether some maximal cones of `fan` have a convex union, by checking every
candidate its open equal-span walls leave. Cones merged into one convex cone are
connected through the walls inside it, and each of these walls is open: another
maximal cone containing it would meet the interior of the merged cone. So only
the connected sets of at least two cones in the graph whose edges are the open
walls can merge. When none of them has a convex union, `certified` is true: no
coarsening of `fan` exists.

`certified` stays false, which is inconclusive, when a convex union is found
(`convex_union` lists its cones; attached lower-dimensional cones may still block
the merge) or when a component has more than `max_unions` connected sets.
"""
function open_wall_certificate(fan; max_unions::Int=100_000)
    data = _open_equal_span_walls(fan)
    neighbours = Dict{Int, Vector{Int}}()
    for (i, j) in data.walls
        push!(get!(neighbours, i, Int[]), j)
        push!(get!(neighbours, j, Int[]), i)
    end
    facet_cache = Dict{Int, Any}()
    facets_of(cell) = get!(facet_cache, cell) do
        cone = positive_hull(matrix(QQ, data.rays[data.ray_sets[cell], :]), matrix(QQ, data.lineality))
        map(facets(cone)) do facet
            # Facets of a cone are half-spaces a.x <= 0.
            values = data.rays * collect(QQFieldElem, normal_vector(facet))
            (values=values, rays=filter(ray -> iszero(values[ray]), data.ray_sets[cell]))
        end
    end
    components = Vector{Int}[]
    seen = Set{Int}()
    for start in sort!(collect(keys(neighbours)))
        start in seen && continue
        push!(components, [start]); push!(seen, start)
        queue = [start]
        while !isempty(queue)
            for neighbour in neighbours[popfirst!(queue)]
                neighbour in seen && continue
                push!(seen, neighbour); push!(components[end], neighbour); push!(queue, neighbour)
            end
        end
    end
    checked = 0
    result(convex_union, certified) = (
        open_walls=length(data.walls), components=length(components),
        largest_component=maximum(length.(components); init=0),
        checked_unions=checked, convex_union=convex_union, certified=certified,
    )
    for component in components
        sets = _connected_cone_sets(neighbours, component, max_unions)
        isnothing(sets) && return result(nothing, false)
        for cells in sets
            checked += 1
            _cone_union_is_convex(cells, data.ray_sets, facets_of) && return result(cells, false)
        end
    end
    return result(nothing, true)
end

function write_link_report(
    path,
    label,
    input_path,
    output_path,
    provenance_path,
    result,
)
    mkpath(dirname(path))
    open(path, "w") do io
        heading = result.strategy == :equal_span ? "Equal-span coarsening" : "Iterative local-link coarsening"
        println(io, "# $heading of `$label`")
        println(io)
        println(io, "Input: `$(relpath(input_path, dirname(path)))`")
        println(io)
        println(io, "- input maximal cones: `$(result.input_cones)`")
        println(io, "- input rays modulo lineality: `$(result.input_rays)`")
        println(io, "- merged maximal cones: `$(result.output_cones)`")
        println(io, "- merged rays modulo lineality: `$(result.output_rays)`")
        println(io, "- ambient / fan / lineality dimensions: `$(result.ambient_dimension) / $(result.fan_dimension) / $(result.lineality_dimension)`")
        println(io, "- open equal-span walls: `$(result.open_equal_span_walls)`")
        println(io, "- unions across open walls checked not convex: `$(result.checked_open_wall_unions)`")
        println(io, "- coarsest certified: `$(result.coarsest_certified)`")
        println(io)
        println(io, "## Coarseness")
        println(io)
        if result.coarsest_certified && iszero(result.open_equal_span_walls)
            println(io, "No two maximal cones of equal dimension share a facet that lies in no other")
            println(io, "maximal cone while both span the same linear space. Every convex union of")
            println(io, "maximal cones would contain such a facet, so no coarsening of the output exists.")
        elseif result.coarsest_certified
            println(io, "Cones merged into a convex cone are connected through open equal-span walls,")
            println(io, "facets in no other maximal cone shared by two cones of equal span. Every")
            println(io, "connected set of at least two cones across these walls was checked, and none")
            println(io, "has a convex union, so no coarsening of the output exists.")
        else
            println(io, "Some maximal cones of equal dimension still share a facet in no other maximal")
            println(io, "cone while spanning the same linear space. The merges tried above found no")
            println(io, "convex union there, but this certificate cannot exclude one: coarseness is")
            println(io, "not certified.")
        end
        println(io)
        println(io, "## Algorithm")
        println(io)
        if result.strategy == :equal_span
            println(io, "Cells in the same quotient-linear span are grouped, and a group is merged")
            println(io, "only when exact supporting-boundary tests certify that its complete union")
            println(io, "is convex. Rational cross-section volume is retained as a fallback certificate.")
        else
            println(io, "After each batch the shared-facet complex is rebuilt. Lower faces are derived")
            println(io, "from ray pairs verified by shared-facet intersections (codimension two) and rays")
            println(io, "(when their quotient dimension is one); in larger quotients exact intersections")
            println(io, "of shared facets inside individual cells index the lower faces. A candidate first")
            println(io, "passes an invariant equal-span-block dimension obstruction. Through codimension")
            println(io, "three, its link must be a combinatorial sphere (two")
            println(io, "points, a polygon, or a closed combinatorial 2-sphere). Disjoint candidates")
            println(io, "eliminating the most cells are tried first. A spherical star is merged only when every")
            println(io, "boundary facet supports the whole star; for a spherical star this exact")
            println(io, "rational sign test certifies that its union is convex. Quotient coordinates,")
            println(io, "cross-sections, cell boundary facets, facet incidences, and completed convexity")
            println(io, "tests are cached. After the spherical scans stabilize, a convex span-block")
            println(io, "closure removes boundary subdivisions whose relative links are balls. It uses")
            println(io, "exact supporting-boundary tests, with rational cross-section volumes only as")
            println(io, "a fallback. Any unresolved block restarts the scans at codimension one.")
            if result.fan_dimension - result.lineality_dimension > 4
                println(io, "In codimension four and higher no sphere-recognition claim is made: a star")
                println(io, "is accepted only when the exact rational cross-section volumes of its")
                println(io, "constituent cells sum to its convex-hull volume. Since these cells form a fan,")
                println(io, "their interiors are disjoint, so this equality certifies complete convex coverage.")
            end
        end
        println(io)
        println(io, "## Rounds")
        println(io)
        println(io, "| round | codim | cells before | faces | dimension | spheres | convex | merged stars | cells after |")
        println(io, "|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for (number, round) in enumerate(result.rounds)
            codimension = round.codimension == 0 ? "block" : string(round.codimension)
            println(io, "| $number | $codimension | $(round.cells_before) | $(round.candidate_faces) | $(round.dimension_candidates) | $(round.spherical_links) | $(round.convex_candidates) | $(round.merged_stars) | $(round.cells_after) |")
        end
        println(io)
        println(io, "Output: `$(relpath(output_path, dirname(path)))`")
        println(io, "Provenance: `$(relpath(provenance_path, dirname(path)))`")
    end
end

"""
    coarsen_fan(fan; kwargs...) -> LinkCoarseningResult

Iteratively merge cells using spherical links through codimension three and
exact volume coverage in higher codimensions. Every accepted union is certified
convex over the rationals. Scans restart at codimension one after each merge batch.
`coarsest_certified` is true when no set of output cones connected through open
equal-span walls has a convex union (see `open_wall_certificate`); false means
inconclusive, not that a coarser fan exists. `certify=false` only counts the open
walls, and certifies when there are none.
"""
function coarsen_fan(
    input_fan;
    algorithm::Symbol=:support_link,
    label::AbstractString="tropical variety",
    input_path::Union{Nothing, AbstractString}=nothing,
    output_path::Union{Nothing, AbstractString}=nothing,
    provenance_path::Union{Nothing, AbstractString}=nothing,
    report_path::Union{Nothing, AbstractString}=nothing,
    certify::Bool=true,
)
    algorithm in (:support_link, :equal_span) ||
        throw(ArgumentError("coarsening algorithm must be support_link or equal_span"))
    input_data = rays_modulo_lineality(input_fan)
    input_rays = Matrix{QQFieldElem}(matrix(QQ, input_data.rays_modulo_lineality))
    lineality = Matrix{QQFieldElem}(matrix(QQ, input_data.lineality_basis))
    input_incidence = maximal_cones(IncidenceMatrix, input_fan)
    input_ray_sets = [
        sort!(collect(row(input_incidence, i))) for i in 1:nrows(input_incidence)
    ]
    ambient_dimension = ncols(input_rays)
    lineality_dimension = rank(matrix(QQ, lineality))
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    isempty(input_ray_sets) && error("input fan has no maximal cones")
    quotient_dimension = rank(matrix(QQ, input_rays[first(input_ray_sets), :] * quotient))
    fan_dimension = lineality_dimension + quotient_dimension
    quotient_dimension >= 2 || error("coarsening expects quotient dimension at least two")

    signature_to_ray, quotient_basis = signature_map(input_rays, lineality)
    println("Precomputing invariant span blocks and cross-sections")
    cells, geometries, initial_facets, span_blocks, original_boundary_facets =
        initialize_iterative_cache(
        input_ray_sets,
        input_rays,
        lineality,
        quotient_dimension,
        )
    cached_codimension_one_owners = initial_facets
    certificate_cache = Dict{Any, Any}()
    println(
        "  span blocks=$(length(span_blocks)), cached pointed blocks=" *
        "$(count(x -> !isnothing(x), geometries))",
    )
    rounds = MergeRound[]
    while true
        if algorithm == :support_link
            codimension = 1
            while codimension <= quotient_dimension - 1
                println("Scanning codimension $codimension with $(length(cells)) current cells")
                next_cells, round, current_facets = one_combinatorial_round(
                    cells,
                    codimension,
                    fan_dimension,
                    quotient_dimension,
                    input_rays,
                    lineality,
                    signature_to_ray,
                    quotient_basis,
                    geometries,
                    certificate_cache;
                    cached_codimension_one_owners=cached_codimension_one_owners,
                )
                cached_codimension_one_owners = round.merged_stars > 0 ? nothing : current_facets
                push!(rounds, round)
                println(
                    "  candidates=$(round.candidate_faces), dimension=$(round.dimension_candidates), " *
                    "spherical=$(round.spherical_links), convex=$(round.convex_candidates), " *
                    "merged=$(round.merged_stars), cells=$(round.cells_after)",
                )
                cells = next_cells
                codimension = round.merged_stars > 0 ? 1 : codimension + 1
            end
        end

        println("Scanning convex span-block closure with $(length(cells)) current cells")
        next_cells, closure = one_span_block_round(
            cells,
            input_rays,
            lineality,
            signature_to_ray,
            quotient_basis,
            geometries,
            input_ray_sets,
            original_boundary_facets,
        )
        push!(rounds, closure)
        println(
            "  candidates=$(closure.candidate_faces), convex=$(closure.convex_candidates), " *
            "merged=$(closure.merged_stars), cells=$(closure.cells_after)",
        )
        cells = next_cells
        closure.merged_stars == 0 && break
        algorithm == :equal_span && continue
        cached_codimension_one_owners = nothing
        length(unique(cell.span_block for cell in cells)) == length(cells) && break
    end

    fan = validated_output_fan(cells, input_rays, lineality)
    order, _ = output_order(fan, cells, signature_to_ray, quotient_basis)
    cells = cells[order]
    output_data = rays_modulo_lineality(fan)
    output_lineality = Matrix{QQFieldElem}(matrix(QQ, output_data.lineality_basis))
    ncols(output_lineality) == ambient_dimension || error("ambient dimension changed")
    rank(matrix(QQ, output_lineality)) == lineality_dimension ||
        error("lineality dimension changed")
    certificate = certify ? open_wall_certificate(fan) :
        (open_walls=open_equal_span_walls(fan), checked_unions=0)
    certified = certify ? certificate.certified : iszero(certificate.open_walls)
    result = LinkCoarseningResult(
        fan,
        algorithm,
        length(input_ray_sets),
        nrows(input_rays),
        n_maximal_cones(fan),
        length(output_data.rays_modulo_lineality),
        ambient_dimension,
        fan_dimension,
        lineality_dimension,
        rounds,
        [copy(cell.original_cells) for cell in cells],
        certificate.open_walls,
        certificate.checked_unions,
        certified,
    )

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fan)
    end
    if !isnothing(provenance_path)
        write_provenance(provenance_path, collect(eachindex(cells)), cells)
    end
    if !isnothing(report_path)
        isnothing(input_path) && error("report output requires input_path")
        isnothing(output_path) && error("report output requires output_path")
        isnothing(provenance_path) && error("report output requires provenance_path")
        write_link_report(
            report_path,
            label,
            input_path,
            output_path,
            provenance_path,
            result,
        )
    end
    return result
end

"""
    merge_cells(fan; kwargs...) -> PolyhedralFan

Iteratively merge complete stars with exact rational convexity certificates,
using combinatorial spheres through codimension three and exact volume coverage
in higher codimensions, then close any remaining convex equal-span blocks.
Use `coarsen_fan` when the round-by-round
certificate summary is also needed.
"""
merge_cells(input_fan; kwargs...) = coarsen_fan(input_fan; kwargs...).fan

function coarsen_fan_file(
    input_path::AbstractString,
    output_path::AbstractString;
    algorithm::Symbol=:support_link,
    label::AbstractString="tropical variety",
    provenance_path::AbstractString=splitext(output_path)[1] * "_provenance.csv",
    report_path::AbstractString=splitext(output_path)[1] * "_report.md",
)
    println("Loading $input_path")
    input_fan = load(input_path)
    return coarsen_fan(
        input_fan;
        algorithm=algorithm,
        label=label,
        input_path=input_path,
        output_path=output_path,
        provenance_path=provenance_path,
        report_path=report_path,
    )
end

merge_cells_file(input_path::AbstractString, output_path::AbstractString; kwargs...) =
    coarsen_fan_file(input_path, output_path; kwargs...).fan

end # module

using .CoarsestFan: open_equal_span_walls, open_wall_certificate

"Strategy marker for exact coarsening algorithms."
abstract type CoarseningStrategy end

"Iterative support-link merging with exact rational convexity certificates."
struct SupportLinkCoarsening <: CoarseningStrategy end

"Equal-span closure, using the same certified coarsener's block phase."
struct EqualSpanCoarsening <: CoarseningStrategy end

function coarsen_fan(
    fan;
    strategy::CoarseningStrategy=SupportLinkCoarsening(),
    kwargs...,
)
    strategy isa Union{SupportLinkCoarsening, EqualSpanCoarsening} ||
        throw(ArgumentError("unsupported coarsening strategy"))
    algorithm = strategy isa SupportLinkCoarsening ? :support_link : :equal_span
    return CoarsestFan.coarsen_fan(fan; algorithm=algorithm, kwargs...)
end

function coarsen_fan_file(
    input_path::AbstractString,
    output_path::AbstractString;
    strategy::CoarseningStrategy=SupportLinkCoarsening(),
    kwargs...,
)
    strategy isa Union{SupportLinkCoarsening, EqualSpanCoarsening} ||
        throw(ArgumentError("unsupported coarsening strategy"))
    algorithm = strategy isa SupportLinkCoarsening ? :support_link : :equal_span
    return CoarsestFan.coarsen_fan_file(
        input_path,
        output_path;
        algorithm=algorithm,
        kwargs...,
    )
end

function _prevariety_cell(original_cells, extreme_rays, input_ray_sets, dimension)
    originals = sort!(unique!(collect(original_cells)))
    return (
        original_cells=originals,
        extreme_rays=sort!(unique!(collect(extreme_rays))),
        covered_rays=sort!(unique(vcat((input_ray_sets[i] for i in originals)...))),
        dimension=dimension,
    )
end

"The original rays on each supporting facet of an exact merged cone."
function _prevariety_facet_ray_sets(cell, rays, quotient, lineality_dimension)
    pointed_dimension = cell.dimension - lineality_dimension
    pointed_dimension == 0 && return Vector{Set{Int}}()
    projected = rays[cell.extreme_rays, :] * quotient
    pivots = CoarsestFan.independent_coordinate_columns(projected, pointed_dimension)
    cone = positive_hull(matrix(QQ, projected[:, pivots]))
    coordinates = rays[cell.covered_rays, :] * quotient[:, pivots]
    return [
        Set(cell.covered_rays[i] for i in eachindex(cell.covered_rays) if
            iszero(sum(coordinates[i, j] * normal[j] for j in 1:pointed_dimension; init=QQ(0))))
        for normal in (collect(normal_vector(facet)) for facet in facets(Halfspace, cone))
    ]
end

"Test whether an intersection generated by original fan rays is a face."
function _prevariety_intersection_is_face(cell, common_rays, facet_ray_sets)
    common = Set(common_rays)
    smallest_face_extremes = Set(cell.extreme_rays)
    for facet in facet_ray_sets
        issubset(common, facet) && intersect!(smallest_face_extremes, facet)
    end
    return issubset(smallest_face_extremes, common)
end

"""
Restore original cells wherever independent convex unions cease to meet in faces.

For two unions of cells of the original fan, their intersection is generated
by the original rays covered by both unions. The smallest supporting face
containing these rays must have no further extreme rays. This exact test also
protects boundaries between strata of different dimensions. Offending merges
are restored simultaneously, so whole-block coarsening remains equivariant.
"""
function _protect_prevariety_faces(cells, originals, rays, quotient, lineality_dimension)
    facet_cache = Dict{Tuple{Vararg{Int}}, Vector{Set{Int}}}()
    restored = 0
    while any(cell -> length(cell.original_cells) > 1, cells)
        owners_by_ray = CoarsestFan.ray_to_cells(
            [cell.covered_rays for cell in cells], size(rays, 1),
        )
        candidates = Set{Tuple{Int, Int}}()
        for owners in owners_by_ray
            for a in 1:(length(owners) - 1), b in (a + 1):length(owners)
                i, j = owners[a], owners[b]
                length(cells[i].original_cells) == length(cells[j].original_cells) == 1 && continue
                push!(candidates, minmax(i, j))
            end
        end
        bad_merges = Set{Int}()
        for (i, j) in candidates
            left, right = cells[i], cells[j]
            common = intersect(left.covered_rays, right.covered_rays)
            left_facets = get!(facet_cache, Tuple(left.original_cells)) do
                _prevariety_facet_ray_sets(left, rays, quotient, lineality_dimension)
            end
            right_facets = get!(facet_cache, Tuple(right.original_cells)) do
                _prevariety_facet_ray_sets(right, rays, quotient, lineality_dimension)
            end
            if !_prevariety_intersection_is_face(left, common, left_facets) ||
                    !_prevariety_intersection_is_face(right, common, right_facets)
                length(left.original_cells) > 1 && push!(bad_merges, i)
                length(right.original_cells) > 1 && push!(bad_merges, j)
            end
        end
        isempty(bad_merges) && break
        restored += length(bad_merges)
        next_cells = eltype(cells)[]
        for (i, cell) in enumerate(cells)
            if i in bad_merges
                append!(next_cells, originals[cell.original_cells])
            else
                push!(next_cells, cell)
            end
        end
        cells = sort!(next_cells; by=cell -> first(cell.original_cells))
    end
    return cells, restored
end

function _remove_contained_prevariety_cells(cells)
    retained = trues(length(cells))
    removed_originals = Vector{Vector{Int}}()
    owners = Dict{Int, Vector{Int}}()
    for (i, cell) in enumerate(cells), ray in cell.covered_rays
        push!(get!(owners, ray, Int[]), i)
    end
    for i in sort!(collect(eachindex(cells)); by=i -> (cells[i].dimension, i))
        retained[i] || continue
        cell = cells[i]
        candidates = isempty(cell.extreme_rays) ? collect(eachindex(cells)) :
            get(owners, first(cell.extreme_rays), Int[])
        for j in candidates
            i == j && continue
            retained[j] || continue
            cells[j].dimension < cell.dimension && continue
            # Equal cones have the same extreme rays; retain the earlier one.
            cells[j].dimension == cell.dimension && j > i && continue
            all(ray -> ray in cells[j].covered_rays, cell.extreme_rays) || continue
            retained[i] = false
            push!(removed_originals, copy(cell.original_cells))
            break
        end
    end
    return cells[retained], removed_originals
end

"The permutations of the rays induced by coordinate permutations `p`, acting by `v -> v[p]`."
function _prevariety_ray_actions(symmetry, rays, lineality, signature_to_ray, quotient_basis, cone_index)
    n = size(rays, 2)
    lineality_rank = rank(matrix(QQ, lineality))
    return map(symmetry) do generator
        p = Int.(collect(generator))
        sort(p) == 1:n || throw(ArgumentError("$p is not a permutation of the $n coordinates"))
        rank(matrix(QQ, vcat(lineality, lineality[:, p]))) == lineality_rank ||
            throw(ArgumentError("$p does not preserve the lineality space"))
        action = map(1:size(rays, 1)) do i
            image = CoarsestFan.canonical_ray_signature(CoarsestFan.matrix_row(rays, i)[p], quotient_basis)
            get(() -> throw(ArgumentError("$p does not send ray $i to a ray of the fan")),
                signature_to_ray, image)
        end
        all(haskey(cone_index, sort(action[ray_set])) for ray_set in keys(cone_index)) ||
            throw(ArgumentError("$p does not permute the maximal cones of the fan"))
        action
    end
end

"""
The other span partitions in the orbit of partition `start`, each with a ray
permutation sending `start` onto it.
"""
function _span_group_orbit(start, ordered_groups, group_of_cone, cone_index, input_ray_sets, ray_actions)
    isempty(ray_actions) && return Pair{Int, Vector{Int}}[]
    orbit = Dict(start => collect(eachindex(first(ray_actions))))
    queue = [start]
    while !isempty(queue)
        transport = orbit[popfirst!(queue)]
        for action in ray_actions
            image = action[transport]
            members = sort!([cone_index[sort(image[input_ray_sets[i]])] for i in ordered_groups[start]])
            target = group_of_cone[first(members)]
            haskey(orbit, target) && continue
            members == ordered_groups[target] ||
                throw(ArgumentError("the symmetry does not permute the span partitions; are the labels invariant?"))
            orbit[target] = image
            push!(queue, target)
        end
    end
    delete!(orbit, start)
    return collect(orbit)
end

"""
    coarsen_prevariety_fan(fan; strategy=EqualSpanCoarsening(), symmetry=nothing, verbose=false)

Coarsen a possibly nonpure rational prevariety without changing its support or
common lineality. Maximal cones are partitioned by their actual quotient span
and dimension, and the existing exact pure-fan coarsener is applied separately
to each partition. The default merges complete convex equal-span components;
`SupportLinkCoarsening()` additionally tries local link merges.

Recombination verifies that intersections remain faces and restores any
incompatible merge. Contained lower-dimensional cones are then discarded.
The result reports the actual reductions and original-cell membership.
`coarsest_certified` is true when no set of output cones connected through open
equal-span walls has a convex union (see `CoarsestFan.open_wall_certificate`), so no
maximal cones can be merged; false only means the certificate is inconclusive, not
that a coarser fan exists. Whole-block closure and
simultaneous compatibility rollback respect any coordinate symmetry of the
input fan. A caller using local link merges should verify its desired group.
Maximal cones with different `labels` (for instance weights) are never merged.

`symmetry` takes coordinate permutations `p`, acting by `v -> v[p]`, that
generate a group preserving the fan, its lineality and the labels. A
permutation sends a span partition with all its structure onto another, so only
one partition per orbit is coarsened, and its cells are carried to the others
by the induced permutation of the rays. Each generator is checked to permute
the rays and maximal cones; the face checks and the certificate above still
run on the whole result.
"""
function coarsen_prevariety_fan(
    input_fan;
    strategy::CoarseningStrategy=EqualSpanCoarsening(),
    labels::Union{Nothing, AbstractVector}=nothing,
    symmetry::Union{Nothing, AbstractVector}=nothing,
    verbose::Bool=false,
)
    strategy isa Union{SupportLinkCoarsening, EqualSpanCoarsening} ||
        throw(ArgumentError("unsupported coarsening strategy"))
    input_data = rays_modulo_lineality(input_fan)
    rays = Matrix{QQFieldElem}(matrix(QQ, input_data.rays_modulo_lineality))
    lineality = Matrix{QQFieldElem}(matrix(QQ, input_data.lineality_basis))
    incidence = maximal_cones(IncidenceMatrix, input_fan)
    input_ray_sets = [sort!(collect(row(incidence, i))) for i in 1:nrows(incidence)]
    isempty(input_ray_sets) && throw(ArgumentError("the prevariety fan has no maximal cones"))
    labels === nothing || length(labels) == length(input_ray_sets) ||
        throw(ArgumentError("one label per maximal cone"))
    lineality_dimension = rank(matrix(QQ, lineality))
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    projected = rays * quotient
    groups = Dict{Any, Vector{Int}}()
    dimensions = Int[]
    for (i, ray_set) in enumerate(input_ray_sets)
        quotient_dimension, reduced = rref(matrix(QQ, projected[ray_set, :]))
        push!(dimensions, lineality_dimension + quotient_dimension)
        key = (labels === nothing ? nothing : labels[i], quotient_dimension,
               Tuple(vec(Matrix{QQFieldElem}(reduced)[1:quotient_dimension, :])))
        push!(get!(groups, key, Int[]), i)
    end
    originals = [
        _prevariety_cell([i], input_ray_sets[i], input_ray_sets, dimensions[i])
        for i in eachindex(input_ray_sets)
    ]
    cells = eltype(originals)[]
    ordered_groups = sort!(collect(values(groups)); by=first)
    stratum_rounds = Dict(dimension => CoarsestFan.MergeRound[] for dimension in unique(dimensions))
    signature_to_ray, quotient_basis = CoarsestFan.signature_map(rays, lineality)
    cone_index = Dict(ray_set => i for (i, ray_set) in enumerate(input_ray_sets))
    group_of_cone = zeros(Int, length(input_ray_sets))
    for (group_number, members) in enumerate(ordered_groups)
        group_of_cone[members] .= group_number
    end
    ray_actions = isnothing(symmetry) ? Vector{Int}[] :
        _prevariety_ray_actions(symmetry, rays, lineality, signature_to_ray, quotient_basis, cone_index)
    transported = falses(length(ordered_groups))
    for (group_number, members) in enumerate(ordered_groups)
        transported[group_number] && continue
        dimension = dimensions[first(members)]
        quotient_dimension = dimension - lineality_dimension
        if length(members) == 1 || quotient_dimension < 2
            append!(cells, originals[members])
            continue
        end
        verbose && println("Coarsening prevariety span $group_number/$(length(ordered_groups)): " *
                           "dimension=$dimension, cells=$(length(members))")
        used_rays = sort!(unique(vcat(input_ray_sets[members]...)))
        local_index = Dict(ray => i for (i, ray) in enumerate(used_rays))
        group_fan = polyhedral_fan(
            IncidenceMatrix([[local_index[ray] for ray in input_ray_sets[i]] for i in members]),
            matrix(QQ, rays[used_rays, :]), matrix(QQ, lineality);
            non_redundant=true,
        )
        # The whole output is certified once below.
        result = if verbose
            coarsen_fan(group_fan; strategy=strategy, certify=false)
        else
            redirect_stdout(devnull) do
                coarsen_fan(group_fan; strategy=strategy, certify=false)
            end
        end
        append!(stratum_rounds[dimension], result.rounds)
        result_rays = matrix(QQ, rays_modulo_lineality(result.fan).rays_modulo_lineality)
        result_to_original = [
            signature_to_ray[CoarsestFan.canonical_ray_signature(
                CoarsestFan.matrix_row(result_rays, i), quotient_basis,
            )] for i in 1:nrows(result_rays)
        ]
        result_incidence = maximal_cones(IncidenceMatrix, result.fan)
        group_cells = [
            _prevariety_cell(
                members[result.original_cells[i]],
                result_to_original[collect(row(result_incidence, i))],
                input_ray_sets, dimension,
            ) for i in 1:nrows(result_incidence)
        ]
        append!(cells, group_cells)
        orbit = _span_group_orbit(group_number, ordered_groups, group_of_cone, cone_index,
                                  input_ray_sets, ray_actions)
        verbose && !isempty(orbit) && println("  carried to $(length(orbit)) partitions of its orbit")
        for (image, action) in orbit
            transported[image] = true
            for cell in group_cells
                push!(cells, _prevariety_cell(
                    [cone_index[sort(action[input_ray_sets[i]])] for i in cell.original_cells],
                    action[cell.extreme_rays], input_ray_sets, dimension,
                ))
            end
        end
    end
    cells = sort!(cells; by=cell -> first(cell.original_cells))
    cells, rolled_back_merges = _protect_prevariety_faces(
        cells, originals, rays, quotient, lineality_dimension,
    )
    cells, removed_original_cells = _remove_contained_prevariety_cells(cells)
    used_rays = sort!(unique(vcat((cell.extreme_rays for cell in cells)...)))
    compact_index = Dict(ray => i for (i, ray) in enumerate(used_rays))
    fan = polyhedral_fan(
        IncidenceMatrix([[compact_index[ray] for ray in cell.extreme_rays] for cell in cells]),
        matrix(QQ, rays[used_rays, :]), matrix(QQ, lineality);
        non_redundant=true,
    )
    n_maximal_cones(fan) == length(cells) || error("mixed coarsening changed cell count on construction")
    output_rays = matrix(QQ, rays_modulo_lineality(fan).rays_modulo_lineality)
    output_to_original = [
        signature_to_ray[CoarsestFan.canonical_ray_signature(
            CoarsestFan.matrix_row(output_rays, i), quotient_basis,
        )] for i in 1:nrows(output_rays)
    ]
    cells_by_rays = Dict(Tuple(cell.extreme_rays) => cell for cell in cells)
    output_incidence = maximal_cones(IncidenceMatrix, fan)
    cells = [cells_by_rays[Tuple(sort!(output_to_original[collect(row(output_incidence, i))]))]
             for i in 1:nrows(output_incidence)]
    output_lineality = matrix(QQ, rays_modulo_lineality(fan).lineality_basis)
    output_lineality_rank, output_lineality_reduced = rref(output_lineality)
    input_lineality_rank, input_lineality_reduced = rref(matrix(QQ, lineality))
    output_lineality_rank == input_lineality_rank == lineality_dimension &&
        Matrix(output_lineality_reduced)[1:lineality_dimension, :] ==
        Matrix(input_lineality_reduced)[1:lineality_dimension, :] ||
        error("mixed coarsening changed common lineality")
    strata = [
        (
            dimension=dimension,
            input_cones=count(==(dimension), dimensions),
            output_cones=count(cell -> cell.dimension == dimension, cells),
            rounds=stratum_rounds[dimension],
        ) for dimension in sort!(unique(dimensions))
    ]
    certificate = CoarsestFan.open_wall_certificate(fan)
    return (
        fan=fan,
        strategy=strategy isa EqualSpanCoarsening ? :equal_span : :support_link,
        input_cones=length(input_ray_sets),
        input_rays=size(rays, 1),
        output_cones=length(cells),
        output_rays=nrows(output_rays),
        ambient_dimension=size(rays, 2),
        fan_dimension=maximum(dimensions),
        lineality_dimension=lineality_dimension,
        strata=strata,
        original_cells=[copy(cell.original_cells) for cell in cells],
        removed_redundant_cones=length(removed_original_cells),
        removed_original_cells=removed_original_cells,
        rolled_back_merges=rolled_back_merges,
        transported_span_partitions=count(transported),
        support_preserved=true,
        open_equal_span_walls=certificate.open_walls,
        checked_open_wall_unions=certificate.checked_unions,
        coarsest_certified=certificate.certified,
    )
end

"""
    prevariety_dimension_skeleton(fan, dimension_bound)

Construct the complete subfan of faces of dimension at most `dimension_bound`.
Every lower-dimensional maximal cell is retained. Higher-dimensional cells
are replaced by all faces at the bound, using exact supporting-facet zero
sets recursively. The result uses only original rays and preserves the full
common lineality. Deduplication and removal of contained faces change no
skeleton support. This operation does not preserve overdimensional interiors;
their algebraic exclusion must be established by the caller.
"""
function prevariety_dimension_skeleton(input_fan, dimension_bound::Integer)
    input_data = rays_modulo_lineality(input_fan)
    rays = Matrix{QQFieldElem}(matrix(QQ, input_data.rays_modulo_lineality))
    lineality = Matrix{QQFieldElem}(matrix(QQ, input_data.lineality_basis))
    lineality_dimension = rank(matrix(QQ, lineality))
    dimension_bound >= lineality_dimension ||
        throw(ArgumentError("the dimension bound is smaller than common lineality"))
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    projected = rays * quotient
    incidence = maximal_cones(IncidenceMatrix, input_fan)
    input_cells = [sort!(collect(row(incidence, i))) for i in 1:nrows(incidence)]
    isempty(input_cells) && throw(ArgumentError("the input fan has no maximal cones"))
    dimensions = [lineality_dimension + rank(matrix(QQ, projected[cell, :])) for cell in input_cells]
    face_dimensions = Dict{Tuple,Int}()
    face_parents = Dict{Tuple,Vector{Int}}()
    generated = Dict{Tuple,Bool}()
    facet_cache = Dict{Tuple,Vector{Vector{Int}}}()
    enumerated_faces = 0
    function add_face!(face, dimension, parent, is_generated)
        get(face_dimensions, face, dimension) == dimension || error("a face changed dimension")
        face_dimensions[face] = dimension
        push!(get!(face_parents, face, Int[]), parent)
        generated[face] = get(generated, face, false) || is_generated
    end
    for (parent, (cell, dimension)) in enumerate(zip(input_cells, dimensions))
        if dimension <= dimension_bound
            add_face!(Tuple(cell), dimension, parent, false)
            continue
        end
        pending = Tuple{Tuple,Int}[(Tuple(cell), dimension)]
        visited = Set{Tuple}()
        while !isempty(pending)
            face, face_dimension = pop!(pending)
            face in visited && continue
            push!(visited, face)
            if face_dimension == dimension_bound
                add_face!(face, face_dimension, parent, true)
                enumerated_faces += 1
                continue
            end
            facet_sets = get!(facet_cache, face) do
                record = (dimension=face_dimension, extreme_rays=collect(face), covered_rays=collect(face))
                [sort!(collect(facet)) for facet in _prevariety_facet_ray_sets(
                    record, rays, quotient, lineality_dimension,
                )]
            end
            for facet in facet_sets
                facet_dimension = lineality_dimension + rank(matrix(QQ, projected[facet, :]))
                facet_dimension == face_dimension - 1 ||
                    error("supporting facet zero set has incorrect codimension")
                facet_dimension >= dimension_bound && push!(pending, (Tuple(facet), facet_dimension))
            end
        end
    end
    # Distinct faces of the same dimension cannot properly contain each other.
    # Search only already retained higher-dimensional cells through the rarest
    # incident ray; this avoids quadratic work for a pure skeleton.
    ordered = sort!(collect(keys(face_dimensions)); by=face -> (-face_dimensions[face], face))
    retained = Tuple[]
    owners = Dict{Int,Vector{Int}}()
    removed_contained_faces = 0
    largest_dimension = maximum(values(face_dimensions))
    for face in ordered
        candidates = if face_dimensions[face] == largest_dimension
            Int[]
        elseif isempty(face)
            collect(eachindex(retained))
        else
            ray = argmin(index -> length(get(owners, index, Int[])), face)
            get(owners, ray, Int[])
        end
        if any(index -> face_dimensions[retained[index]] > face_dimensions[face] &&
                        all(ray in retained[index] for ray in face), candidates)
            removed_contained_faces += 1
            continue
        end
        push!(retained, face)
        for ray in face
            push!(get!(owners, ray, Int[]), length(retained))
        end
    end
    used_rays = sort!(unique(Int[ray for face in retained for ray in face]))
    compact_index = Dict(ray => index for (index, ray) in enumerate(used_rays))
    fan = polyhedral_fan(
        IncidenceMatrix([[compact_index[ray] for ray in face] for face in retained]),
        matrix(QQ, rays[used_rays, :]), matrix(QQ, lineality);
        non_redundant=true,
    )
    n_maximal_cones(fan) == length(retained) || error("skeleton construction changed its cell count")
    signature_to_ray, quotient_basis = CoarsestFan.signature_map(rays, lineality)
    output_data = rays_modulo_lineality(fan)
    output_rays = matrix(QQ, output_data.rays_modulo_lineality)
    output_to_original = [signature_to_ray[CoarsestFan.canonical_ray_signature(
        CoarsestFan.matrix_row(output_rays, i), quotient_basis,
    )] for i in 1:nrows(output_rays)]
    output_incidence = maximal_cones(IncidenceMatrix, fan)
    provenance = [begin
        face = Tuple(sort!(output_to_original[collect(row(output_incidence, i))]))
        (output_cell=i, original_cells=sort!(unique(face_parents[face])),
         original_ray_indices=collect(face), dimension=face_dimensions[face],
         generated_from_overdimensional_cone=generated[face])
    end for i in 1:nrows(output_incidence)]
    output_lineality_rank, output_lineality_reduced = rref(matrix(QQ, output_data.lineality_basis))
    input_lineality_rank, input_lineality_reduced = rref(matrix(QQ, lineality))
    output_lineality_rank == input_lineality_rank == lineality_dimension &&
        Matrix(output_lineality_reduced)[1:lineality_dimension, :] ==
        Matrix(input_lineality_reduced)[1:lineality_dimension, :] ||
        error("skeleton changed common lineality")
    original_low_faces = Set(Tuple(input_cells[i]) for i in eachindex(input_cells) if dimensions[i] <= dimension_bound)
    counts = (
        input_cones=length(input_cells), input_rays=size(rays, 1),
        retained_original_cones=count(<=(dimension_bound), dimensions),
        removed_overdimensional_cones=count(>(dimension_bound), dimensions),
        enumerated_faces=enumerated_faces,
        added_faces=count(face -> !(face in original_low_faces), retained),
        removed_contained_faces=removed_contained_faces,
        output_cones=nrows(output_incidence), output_rays=nrows(output_rays),
    )
    return (; fan, counts, provenance, dimension_bound=Int(dimension_bound),
             lineality_dimension, skeleton_certified=true, support_included=true,
             input_support_preserved=all(<=(dimension_bound), dimensions))
end
