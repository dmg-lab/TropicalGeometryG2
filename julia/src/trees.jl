"An exact unrooted phylogenetic tree with rational edge lengths."
struct TreeMetric
    ntaxa::Int
    adjacency::Dict{Int, Vector{Tuple{Int, QQFieldElem}}}
end

"Gauge data for Pluecker tree metrics; the all-ones direction is kept separate."
struct TreeGauge
    ntaxa::Int
    orthogonal_leaf_directions::Matrix{QQFieldElem}

    function TreeGauge(ntaxa::Integer, directions::AbstractMatrix)
        ntaxa >= 3 || throw(ArgumentError("a tree gauge needs at least three taxa"))
        size(directions, 2) == ntaxa ||
            throw(ArgumentError("gauge directions must have one column per taxon"))
        rows = QQ.(Matrix(directions))
        all(iszero(sum(rows[i, :]; init=QQ(0))) for i in axes(rows, 1)) ||
            throw(ArgumentError("orthogonal gauge directions must preserve the all-ones direction"))
        new(Int(ntaxa), rows)
    end
end

TreeGauge(ntaxa::Integer) = TreeGauge(ntaxa, zeros(Int, 0, ntaxa))

pluecker_pairs(ntaxa::Integer) =
    [(i, j) for i in 1:(ntaxa - 1) for j in (i + 1):ntaxa]

"Map leaf-star lengths to min-convention Pluecker coordinates."
function star_pluecker(lengths::AbstractVector)
    return QQ.([lengths[i] + lengths[j] for (i, j) in pluecker_pairs(length(lengths))])
end

"Clear denominators and divide by the gcd without changing the ray sign."
function primitive_integer(vector::AbstractVector)
    q = QQ.(collect(vector))
    isempty(q) && throw(ArgumentError("cannot normalize an empty vector"))
    all(iszero, q) && throw(ArgumentError("the zero vector does not generate a ray"))
    denominator_lcm = foldl(lcm, denominator.(q); init=ZZ(1))
    integral = denominator_lcm .* q
    divisor = foldl(
        gcd,
        (abs(numerator(value)) for value in integral if !iszero(value));
        init=ZZ(0),
    )
    return integral ./ divisor
end

function number_of_taxa(ncoordinates::Integer)
    discriminant = 1 + 8 * ncoordinates
    root = isqrt(discriminant)
    ntaxa = (1 + root) ÷ 2
    root^2 == discriminant && binomial(ntaxa, 2) == ncoordinates ||
        throw(ArgumentError("$ncoordinates is not binomial(n, 2)"))
    return ntaxa
end

function pair_dictionary(values::AbstractVector)
    ntaxa = number_of_taxa(length(values))
    distances = Dict{Tuple{Int, Int}, QQFieldElem}()
    index = 0
    for i in 1:(ntaxa - 1), j in (i + 1):ntaxa
        index += 1
        distances[i, j] = distances[j, i] = QQ(values[index])
    end
    return distances, ntaxa
end

"Test the min-convention three-term tropical Pluecker relations."
function satisfies_tropical_pluecker(weight::AbstractVector)
    values, ntaxa = pair_dictionary(weight)
    for i in 1:(ntaxa - 3), j in (i + 1):(ntaxa - 2),
            k in (j + 1):(ntaxa - 1), l in (k + 1):ntaxa
        terms = (
            values[i, j] + values[k, l],
            values[i, k] + values[j, l],
            values[i, l] + values[j, k],
        )
        minimum_term = min(terms...)
        count(==(minimum_term), terms) >= 2 || return false
    end
    return true
end

function neighbor_joining(
    input_distances::Dict{Tuple{Int, Int}, QQFieldElem},
    ntaxa::Int,
)
    distances = copy(input_distances)
    edges = Tuple{Int, Int, QQFieldElem}[]
    active = collect(1:ntaxa)
    next_label = ntaxa + 1
    while length(active) > 2
        count_active = length(active)
        row_sums = Dict(
            i => sum((distances[i, k] for k in active if k != i); init=QQ(0))
            for i in active
        )
        best_i, best_j = active[1], active[2]
        best_value = (count_active - 2) * distances[best_i, best_j] -
                     row_sums[best_i] - row_sums[best_j]
        for a in 1:(count_active - 1), b in (a + 1):count_active
            i, j = active[a], active[b]
            value = (count_active - 2) * distances[i, j] - row_sums[i] - row_sums[j]
            if value < best_value
                best_i, best_j, best_value = i, j, value
            end
        end
        vertex = next_label
        next_label += 1
        i, j = best_i, best_j
        length_i = distances[i, j] / 2 +
                   (row_sums[i] - row_sums[j]) / (2 * (count_active - 2))
        length_j = distances[i, j] - length_i
        push!(edges, (i, vertex, length_i), (j, vertex, length_j))
        for k in active
            k in (i, j) && continue
            distances[vertex, k] = distances[k, vertex] =
                (distances[i, k] + distances[j, k] - distances[i, j]) / 2
        end
        active = push!([k for k in active if k != i && k != j], vertex)
    end
    i, j = active
    push!(edges, (i, j, distances[i, j]))
    adjacency = Dict(
        vertex => Tuple{Int, QQFieldElem}[] for vertex in 1:(next_label - 1)
    )
    for (left, right, edge_length) in edges
        push!(adjacency[left], (right, edge_length))
        push!(adjacency[right], (left, edge_length))
    end
    return TreeMetric(ntaxa, adjacency)
end

function contract_zero_internal_edges(tree::TreeMetric)
    nvertices = maximum(keys(tree.adjacency))
    parent = collect(1:nvertices)
    function representative(vertex)
        while parent[vertex] != vertex
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        end
        return vertex
    end
    function unite!(left, right)
        left, right = representative(left), representative(right)
        left == right || (parent[max(left, right)] = min(left, right))
    end
    for left in 1:nvertices, (right, edge_length) in tree.adjacency[left]
        left < right || continue
        left > tree.ntaxa && right > tree.ntaxa && iszero(edge_length) &&
            unite!(left, right)
    end
    roots = sort!(unique(representative(v) for v in (tree.ntaxa + 1):nvertices))
    labels = Dict(root => tree.ntaxa + index for (index, root) in enumerate(roots))
    relabel(vertex) = vertex <= tree.ntaxa ? vertex : labels[representative(vertex)]
    adjacency = Dict(
        vertex => Tuple{Int, QQFieldElem}[] for vertex in 1:(tree.ntaxa + length(roots))
    )
    for left in 1:nvertices, (right, edge_length) in tree.adjacency[left]
        left < right || continue
        new_left, new_right = relabel(left), relabel(right)
        new_left == new_right && continue
        push!(adjacency[new_left], (new_right, edge_length))
        push!(adjacency[new_right], (new_left, edge_length))
    end
    return TreeMetric(tree.ntaxa, adjacency)
end

"Reconstruct the exact tree whose dissimilarities are `-weight`."
function tree_from_min_pluecker(weight::AbstractVector)
    satisfies_tropical_pluecker(weight) ||
        throw(ArgumentError("the vector is not in Trop(Gr(2,n))"))
    distances, ntaxa = pair_dictionary(-QQ.(collect(weight)))
    return contract_zero_internal_edges(neighbor_joining(distances, ntaxa))
end

function tree_distance(tree::TreeMetric, source::Int, target::Int)
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
    error("tree is disconnected")
end

pendant_lengths(tree::TreeMetric) =
    [only(tree.adjacency[taxon])[2] for taxon in 1:tree.ntaxa]

tree_edge_lengths(tree::TreeMetric) = [
    edge_length
    for left in sort!(collect(keys(tree.adjacency)))
    for (right, edge_length) in tree.adjacency[left]
    if left < right
]

function realizes(tree::TreeMetric, distances::AbstractVector)
    expected, ntaxa = pair_dictionary(distances)
    ntaxa == tree.ntaxa || return false
    return all(
        tree_distance(tree, i, j) == expected[i, j]
        for i in 1:(ntaxa - 1) for j in (i + 1):ntaxa
    )
end

function leaves_on_side(tree::TreeMetric, start::Int, blocked::Int)
    leaves = Int[]
    stack = [(start, blocked)]
    while !isempty(stack)
        vertex, parent = pop!(stack)
        vertex <= tree.ntaxa && push!(leaves, vertex)
        for (neighbor, _) in tree.adjacency[vertex]
            neighbor == parent || push!(stack, (neighbor, vertex))
        end
    end
    return sort!(leaves)
end

"Positive internal edges as canonical unordered taxon splits."
function internal_splits(tree::TreeMetric)
    result = NamedTuple[]
    taxa = collect(1:tree.ntaxa)
    for left in sort!(collect(keys(tree.adjacency)))
        for (right, edge_length) in tree.adjacency[left]
            left < right && left > tree.ntaxa && right > tree.ntaxa && edge_length > 0 ||
                continue
            side = leaves_on_side(tree, left, right)
            complement = setdiff(taxa, side)
            if (length(side), Tuple(side)) > (length(complement), Tuple(complement))
                side, complement = complement, side
            end
            push!(result, (side=side, complement=complement, length=edge_length))
        end
    end
    sort!(result; by=split -> (length(split.side), Tuple(split.side)))
    return result
end

inner_product(left, right) =
    sum((QQ(left[i]) * QQ(right[i]) for i in eachindex(left)); init=QQ(0))

function independent_rows(rows::Matrix{QQFieldElem})
    isempty(rows) && return rows
    rank_value, reduced = rref(matrix(QQ, rows))
    return Matrix{QQFieldElem}(reduced)[1:rank_value, :]
end

function orthogonal_projection(vector, directions)
    basis = independent_rows(QQ.(Matrix(directions)))
    isempty(basis) && return (coefficients=QQFieldElem[], shift=zeros(QQFieldElem, length(vector)))
    gram = basis * transpose(basis)
    products = [inner_product(vector, basis[row, :]) for row in axes(basis, 1)]
    coefficients = inv(gram) * products
    shift = vec(transpose(basis) * coefficients)
    return (coefficients=coefficients, shift=shift)
end

"Normalize a quotient ray while reserving the all-ones shift for nonnegativity."
function normalize_ray(weight::AbstractVector, gauge::TreeGauge)
    length(weight) == binomial(gauge.ntaxa, 2) ||
        throw(ArgumentError("weight has the wrong number of Pluecker coordinates"))
    original = QQ.(collect(weight))
    original_tree = tree_from_min_pluecker(original)
    outer = pendant_lengths(original_tree)
    projection = orthogonal_projection(outer, gauge.orthogonal_leaf_directions)

    orthogonal_weight = original + star_pluecker(projection.shift)
    orthogonal_tree = tree_from_min_pluecker(orthogonal_weight)
    orthogonal_outer = pendant_lengths(orthogonal_tree)
    @assert all(
        iszero(inner_product(orthogonal_outer, direction))
        for direction in eachrow(independent_rows(gauge.orthogonal_leaf_directions))
    )

    integer_scale = foldl(lcm, denominator.(tree_edge_lengths(orthogonal_tree)); init=ZZ(1))
    integral_weight = integer_scale .* orthogonal_weight
    diagonal_leaf_shift = minimum(pendant_lengths(tree_from_min_pluecker(integral_weight)))
    shifted_weight = integral_weight + star_pluecker(fill(diagonal_leaf_shift, gauge.ntaxa))
    shifted_tree = tree_from_min_pluecker(shifted_weight)
    primitive_divisor = foldl(
        gcd,
        (abs(numerator(value)) for value in tree_edge_lengths(shifted_tree));
        init=ZZ(0),
    )
    primitive_divisor > 0 || error("normalization produced the zero tree")
    normalized = shifted_weight ./ primitive_divisor
    tree = tree_from_min_pluecker(normalized)
    @assert minimum(pendant_lengths(tree)) == 0
    @assert all(value >= 0 && denominator(value) == 1 for value in tree_edge_lengths(tree))
    @assert realizes(tree, -normalized)
    return (
        original=original,
        projection_coefficients=projection.coefficients,
        lineality_leaf_shift=projection.shift,
        integer_scale=integer_scale,
        diagonal_leaf_shift=diagonal_leaf_shift,
        primitive_divisor=primitive_divisor,
        normalized=normalized,
        balanced_generator=normalized,
        pendant_lengths=pendant_lengths(tree),
        tree=tree,
    )
end

function row_combination_solution(rows, target)
    equations = matrix(QQ, transpose(QQ.(Matrix(rows))))
    right_hand_side = matrix(QQ, reshape(QQ.(collect(target)), :, 1))
    solvable, solution = can_solve_with_solution(equations, right_hand_side; side=:right)
    solvable || error("requested shift is not in the supplied lineality space")
    return vec(Matrix{QQFieldElem}(solution))
end

"""
Normalize a vector containing a Pluecker block. The first projection uses the
full ambient lineality; the tree gauge is then lifted back through that same
lineality. An optional anchor coordinate fixes the residual projective shift.
"""
function normalize_ray(
    vector::AbstractVector,
    gauge::TreeGauge,
    lineality::AbstractMatrix;
    tree_coordinates::AbstractVector{<:Integer},
    anchor_coordinate::Union{Nothing, Integer}=nothing,
)
    rows = QQ.(Matrix(lineality))
    input = QQ.(collect(vector))
    length(input) == size(rows, 2) || throw(ArgumentError("ambient dimensions differ"))
    length(tree_coordinates) == binomial(gauge.ntaxa, 2) ||
        throw(ArgumentError("tree coordinate block has the wrong size"))
    gram = rows * transpose(rows)
    orthogonal = input - vec(transpose(rows) * (inv(gram) * (rows * input)))
    @assert all(iszero, rows * orthogonal)
    tree_part = orthogonal[tree_coordinates]
    tree_lineality = rows[:, tree_coordinates]
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, tree_lineality); side=:right))
    vertical = all(iszero, vec(transpose(tree_part) * quotient))

    tree_record = vertical ? nothing : normalize_ray(tree_part, gauge)
    scale = vertical ? QQ(1) :
            QQ(tree_record.integer_scale) / QQ(tree_record.primitive_divisor)
    target_tree = vertical ? zeros(QQFieldElem, length(tree_coordinates)) :
                  QQ.(tree_record.normalized)
    coefficients = row_combination_solution(tree_lineality, target_tree - scale .* tree_part)
    normalized = scale .* orthogonal + vec(transpose(rows) * coefficients)
    normalized[tree_coordinates] == target_tree || error("tree-block lift failed")

    if !isnothing(anchor_coordinate)
        anchor = Int(anchor_coordinate)
        constraint_rows = vcat(
            transpose(tree_lineality),
            reshape(rows[:, anchor], 1, :),
        )
        target = vcat(zeros(QQFieldElem, length(tree_coordinates)), -normalized[anchor])
        solvable, correction = can_solve_with_solution(
            matrix(QQ, constraint_rows),
            matrix(QQ, reshape(target, :, 1));
            side=:right,
        )
        solvable || error("cannot impose the requested anchor through residual lineality")
        normalized .+= vec(transpose(rows) * vec(Matrix{QQFieldElem}(correction)))
        iszero(normalized[anchor]) || error("anchor normalization failed")
        normalized[tree_coordinates] == target_tree || error("anchor changed the tree block")
    end
    return (
        original=input,
        orthogonal=orthogonal,
        normalized=normalized,
        balanced_generator=normalized,
        tree=vertical ? nothing : tree_record.tree,
        tree_record=tree_record,
        vertical=vertical,
        pluecker_scale=scale,
    )
end

function gauge_lineality(gauge::TreeGauge)
    directions = independent_rows(gauge.orthogonal_leaf_directions)
    rows = Vector{Vector{QQFieldElem}}([ones(QQFieldElem, binomial(gauge.ntaxa, 2))])
    append!(rows, [star_pluecker(collect(row)) for row in eachrow(directions)])
    return reduce(vcat, permutedims.(rows))
end

"Normalize every ray and rebuild a fan with the same maximal-cell incidence."
function normalize_fan(fan, gauge::TreeGauge)
    ray_data = collect(rays_modulo_lineality(fan).rays_modulo_lineality)
    records = [normalize_ray(primitive_integer(collect(ray)), gauge) for ray in ray_data]
    rays = reduce(vcat, permutedims.([record.normalized for record in records]))
    incidence = maximal_cones(IncidenceMatrix, fan)
    normalized_fan = polyhedral_fan(
        incidence,
        matrix(QQ, rays),
        matrix(QQ, gauge_lineality(gauge));
        non_redundant=true,
    )
    return (fan=normalized_fan, rays=matrix(QQ, rays), records=records)
end

"""
    normalize_fan(fan)

Without a tree gauge, represent each ray by the primitive integer vector of its
orthogonal projection onto the complement of the lineality space. A coordinate
permutation preserving the lineality space then permutes these rays exactly.
"""
function normalize_fan(fan)
    data = rays_modulo_lineality(fan)
    n = ambient_dim(fan)
    lineality = matrix(QQ, data.lineality_basis)
    projection = identity_matrix(QQ, n)
    nrows(lineality) == 0 ||
        (projection -= transpose(lineality) * inv(lineality * transpose(lineality)) * lineality)
    rows = [primitive_integer(vec(collect(matrix(QQ, 1, n, collect(ray)) * projection))) for ray in data.rays_modulo_lineality]
    rays = isempty(rows) ? zero_matrix(QQ, 0, n) : rows_matrix(rows; ring=QQ)
    normalized_fan = polyhedral_fan(maximal_cones(IncidenceMatrix, fan), rays, lineality; non_redundant=true)
    return (fan=normalized_fan, rays=rays)
end

function normalize_fan(
    fan,
    gauge::TreeGauge,
    tree_coordinates::AbstractVector{<:Integer};
    anchor_coordinate::Union{Nothing, Integer}=nothing,
)
    fan_data = rays_modulo_lineality(fan)
    raw_rays = Matrix{QQFieldElem}(matrix(QQ, fan_data.rays_modulo_lineality))
    lineality = Matrix{QQFieldElem}(matrix(QQ, fan_data.lineality_basis))
    records = [
        normalize_ray(
            collect(row),
            gauge,
            lineality;
            tree_coordinates=tree_coordinates,
            anchor_coordinate=anchor_coordinate,
        )
        for row in eachrow(raw_rays)
    ]
    rays = reduce(vcat, permutedims.([record.normalized for record in records]))
    normalized_fan = polyhedral_fan(
        maximal_cones(IncidenceMatrix, fan),
        matrix(QQ, rays),
        matrix(QQ, lineality);
        non_redundant=true,
    )
    return (fan=normalized_fan, rays=matrix(QQ, rays), records=records)
end

"Return normalized ray trees, and optionally canonical maximal-cell trees."
function trees_from_fan(fan, gauge::TreeGauge; maximal_cells::Bool=false)
    normalized = normalize_fan(fan, gauge)
    ray_trees = getproperty.(normalized.records, :tree)
    maximal_cells || return (rays=ray_trees, normalized=normalized)
    weights = canonical_cell_weights(normalized.fan, normalized.rays)
    cell_trees = [tree_from_min_pluecker(weight) for weight in weights]
    return (rays=ray_trees, cells=cell_trees, weights=weights, normalized=normalized)
end

"Exact equal-positive sums of normalized extreme rays, one per maximal cell."
function canonical_cell_weights(fan, normalized_rays)
    incidence = maximal_cones(IncidenceMatrix, fan)
    rays = Matrix{QQFieldElem}(normalized_rays)
    return [
        vec(sum(rays[collect(row(incidence, cell)), :]; dims=1))
        for cell in 1:nrows(incidence)
    ]
end

"""
The normalized tree of each canonical cell weight, in the gauge of the rays.
The gauge is linear on a cell, so a ray sum can only fail to have a pendant
edge of length zero and coprime edge lengths.
"""
normalized_cell_trees(weights, gauge::TreeGauge) =
    [normalize_ray(weight, gauge).normalized for weight in weights]

function normalized_cell_trees(
    weights,
    gauge::TreeGauge,
    lineality::AbstractMatrix,
    tree_coordinates::AbstractVector{<:Integer};
    anchor_coordinate::Union{Nothing, Integer}=nothing,
)
    return [
        normalize_ray(weight, gauge, lineality;
            tree_coordinates=tree_coordinates, anchor_coordinate=anchor_coordinate).normalized
        for weight in weights
    ]
end
