Base.@kwdef struct DegenerationConfig
    expected_orbits::Union{Nothing, Int} = nothing
    compute_bodies::Bool = true
    # Coordinate blocks of a multigrading; `nothing` for a single grading.
    body_gradings::Union{Nothing, Vector{Vector{Int}}} = nothing
    compute_saturation::Bool = true
    compute_degree::Bool = true
    certified_tropical::Bool = false
    coordinate_images::Function = pluecker_coordinate_images
    weight_coordinate_indices::Function = coordinate_permutation
    chevalley_test::Function = body -> missing
end

"Images of Pluecker variables under the exterior-square action of a permutation."
function pluecker_coordinate_images(permutation, variables; leaf_signs=ones(Int, length(permutation)))
    ntaxa = length(permutation)
    pairs = pluecker_pairs(ntaxa)
    length(variables) == length(pairs) ||
        throw(ArgumentError("the ring is not a full Pluecker coordinate ring"))
    index = Dict(pair => position for (position, pair) in enumerate(pairs))
    return [
        begin
            left, right = permutation[i], permutation[j]
            orientation = left < right ? 1 : -1
            orientation * leaf_signs[i] * leaf_signs[j] *
                variables[index[(min(left, right), max(left, right))]]
        end
        for (i, j) in pairs
    ]
end

function verify_ideal_preservation(ideal_value, permutations, coordinate_images)
    ring = base_ring(ideal_value)
    variables = gens(ring)
    generators = gens(ideal_value)
    normalize_generator(polynomial) = iszero(polynomial) ? polynomial :
        inv(leading_coefficient(polynomial)) * polynomial
    normalized_generators = Set(normalize_generator.(generators))
    for permutation in permutations
        images = coordinate_images(permutation, variables)
        length(images) == length(variables) ||
            throw(ArgumentError("coordinate action returned the wrong number of images"))
        action = hom(ring, ring, images)
        image_generators = [action(generator) for generator in generators]
        Set(normalize_generator.(image_generators)) == normalized_generators && continue
        image_ideal = ideal(ring, image_generators)
        image_ideal == ideal_value || return false
    end
    return true
end

function quotient_positive_equivalent(left, right, lineality)
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    projected_left = vec(transpose(QQ.(left)) * quotient)
    projected_right = vec(transpose(QQ.(right)) * quotient)
    pivot = findfirst(!iszero, projected_right)
    pivot === nothing && return all(iszero, projected_left)
    scale = projected_left[pivot] / projected_right[pivot]
    return scale > 0 && projected_left == scale .* projected_right
end

"Certify that canonical cell weights are equivariant modulo lineality and positive scale."
function verify_canonical_weight_action(
    weights, symmetry, lineality;
    coordinate_indices=coordinate_permutation,
)
    # In the compact model the quotient is simply w - w[end]. Normalize once
    # with a positive gcd. Permuting and subtracting another coordinate keeps
    # that gcd unchanged, so every later comparison uses exact integers.
    if size(lineality, 1) == 1 && size(lineality, 2) > 0 &&
       !iszero(lineality[1, 1]) && all(==(lineality[1, 1]), lineality)
        canonical = [_primitive_integer_row(QQ.(weight) .- QQ(last(weight))) for weight in weights]
        limit = BigInt(typemax(Int)) ÷ 2
        keys = all(all(value -> abs(value) <= limit, key) for key in canonical) ?
            [Int.(key) for key in canonical] : canonical
        return _verify_all_ones_weight_action(keys, symmetry, size(lineality, 2), coordinate_indices)
    end
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    projected = [vec(transpose(QQ.(weight)) * quotient) for weight in weights]
    for (permutation, cell_action) in zip(symmetry.permutations, symmetry.cell_actions)
        indices = coordinate_indices(permutation)
        for cell in eachindex(weights)
            left = vec(transpose(QQ.(weights[cell][indices])) * quotient)
            right = projected[cell_action[cell]]
            pivot = findfirst(!iszero, right)
            if pivot === nothing
                all(iszero, left) || return false
            else
                scale = left[pivot] / right[pivot]
                scale > 0 && left == scale .* right || return false
            end
        end
    end
    return true
end

function _verify_all_ones_weight_action(keys::Vector{Vector{T}}, symmetry, ambient, coordinate_indices) where T
    all(length(key) == ambient for key in keys) || throw(ArgumentError("canonical weight has the wrong ambient dimension"))
    length(symmetry.permutations) == length(symmetry.cell_actions) ||
        throw(ArgumentError("missing canonical weight actions"))
    for (permutation, cell_action) in zip(symmetry.permutations, symmetry.cell_actions)
        indices = coordinate_indices(permutation)
        sort(collect(indices)) == collect(1:ambient) || throw(ArgumentError("invalid weight coordinate permutation"))
        length(cell_action) == length(keys) || throw(ArgumentError("weight action omits cells"))
        for cell in eachindex(keys)
            target_cell = cell_action[cell]
            1 <= target_cell <= length(keys) || throw(ArgumentError("weight action references a missing cell"))
            source, target = keys[cell], keys[target_cell]
            shift = source[indices[end]]
            for coordinate in 1:ambient
                source[indices[coordinate]] - shift == target[coordinate] || return false
            end
        end
    end
    return true
end

function integer_matrix(rows)
    rational = QQ.(Matrix(rows))
    all(denominator(value) == 1 for value in rational) ||
        throw(ArgumentError("value matrix must be integral"))
    return matrix(ZZ, numerator.(rational))
end

function distinct_integer_rows(matrix_value)
    seen = Set{Tuple}()
    rows = Vector{Vector{Int}}()
    for row_index in 1:nrows(matrix_value)
        values = [Int(matrix_value[row_index, column]) for column in 1:ncols(matrix_value)]
        key = Tuple(values)
        key in seen || (push!(seen, key); push!(rows, values))
    end
    return rows
end

"""
Values of the coordinates on a cell: row `i` is the image of coordinate `i` in
the saturated lattice of the span of the lineality and the rays of the cell.
Rational rows are scaled to integers, which keeps that span.
"""
function cell_values(lineality, normalized_rays, ray_indices)
    spanning_rows = vcat(QQ.(Matrix(lineality)), QQ.(Matrix(normalized_rays)[ray_indices, :]))
    for row in eachrow(spanning_rows)
        row .*= foldl(lcm, denominator.(row); init=ZZ(1))
    end
    weights = integer_matrix(spanning_rows)
    relation_columns = kernel(weights; side=:right)
    relation_rows = transpose(relation_columns)
    return kernel(relation_rows; side=:right)
end

"""
The polytope of the rows of `values` in the affine lattice generated by their
differences, with its exact invariants.
"""
function affine_lattice_body(values; expected_dimension=nothing)
    differences = matrix(
        ZZ,
        nrows(values),
        ncols(values),
        [values[i, j] - values[1, j] for i in 1:nrows(values) for j in 1:ncols(values)],
    )
    hermite = hnf(differences)
    basis_indices = [
        i for i in 1:nrows(hermite)
        if any(!iszero(hermite[i, j]) for j in 1:ncols(hermite))
    ]
    isnothing(expected_dimension) || length(basis_indices) == expected_dimension ||
        error("degree-one values do not span the expected affine hyperplane")
    affine_basis = hermite[basis_indices, :]
    points = zero_matrix(ZZ, nrows(values), length(basis_indices))
    for i in 1:nrows(values)
        solvable, coordinates = can_solve_with_solution(
            affine_basis,
            differences[i:i, :];
            side=:left,
        )
        solvable || error("value is outside the intrinsic affine lattice")
        points[i, :] = coordinates
    end
    polytope = convex_hull(points)
    return (
        points=points,
        unique_points=distinct_integer_rows(points),
        polytope=polytope,
        f_vector=Int.(f_vector(polytope)),
        normalized_volume=Int(normalized_volume(polytope)),
        lattice_points=length(lattice_points(polytope)),
    )
end

"""
    intrinsic_cell_body(lineality, normalized_rays, ray_indices; gradings=nothing)

Build the saturated intrinsic Newton--Okounkov point configuration of a cell,
from the values of the coordinates, and its convex hull. `gradings`
partitions the coordinates into the blocks of a multigrading, all coordinates
by default. The body is that of the sum of the gradings: the Minkowski sum of
the `slices`, the bodies of each grading, given only for several gradings.
"""
function intrinsic_cell_body(lineality, normalized_rays, ray_indices; gradings=nothing)
    values = cell_values(lineality, normalized_rays, ray_indices)
    blocks = isnothing(gradings) ? [collect(1:nrows(values))] : collect.(gradings)
    sort!(reduce(vcat, blocks)) == collect(1:nrows(values)) ||
        throw(ArgumentError("the gradings must partition the coordinates"))
    block_rows(block) = [[values[i, j] for j in 1:ncols(values)] for i in block]
    # The degree-one values of the sum of the gradings are the sums of one
    # degree-one value of each grading.
    sums = vec([reduce(+, combination) for combination in Iterators.product(block_rows.(blocks)...)])
    body = affine_lattice_body(rows_matrix(sums); expected_dimension=ncols(values) - length(blocks))
    slices = length(blocks) == 1 ? NamedTuple[] :
        [affine_lattice_body(rows_matrix(block_rows(block))) for block in blocks]
    return merge(body, (slices=slices,))
end

function affine_basis_indices(points)
    count_points, dimension = size(points)
    basis = [1]
    for index in 2:count_points
        trial = vcat(basis, [index])
        differences = matrix(
            QQ,
            length(trial) - 1,
            dimension,
            [points[trial[i], j] - points[trial[1], j]
             for i in 2:length(trial) for j in 1:dimension],
        )
        rank(differences) == length(trial) - 1 && push!(basis, index)
        length(basis) == dimension + 1 && return basis
    end
    error("point configuration is not full-dimensional")
end

function affine_map_is_unimodular(source, target, mapping)
    dimension = ncols(source)
    basis = affine_basis_indices(source)
    source_edges = matrix(
        QQ,
        dimension,
        dimension,
        [source[basis[i], j] - source[basis[1], j]
         for i in 2:length(basis) for j in 1:dimension],
    )
    target_edges = matrix(
        QQ,
        dimension,
        dimension,
        [target[mapping[basis[i]], j] - target[mapping[basis[1]], j]
         for i in 2:length(basis) for j in 1:dimension],
    )
    linear = inv(source_edges) * target_edges
    all(denominator(linear[i, j]) == 1 for i in 1:dimension for j in 1:dimension) ||
        return false
    abs(det(linear)) == 1 || return false
    translation = [
        target[mapping[basis[1]], j] -
        sum(source[basis[1], k] * linear[k, j] for k in 1:dimension)
        for j in 1:dimension
    ]
    target_set = Set(Tuple(target[i, j] for j in 1:dimension) for i in 1:nrows(target))
    mapped_set = Set(
        Tuple(sum(source[i, k] * linear[k, j] for k in 1:dimension) + translation[j]
              for j in 1:dimension)
        for i in 1:nrows(source)
    )
    return mapped_set == target_set
end

"Vertex sets of the facets, in the vertex order of `vertices(polytope)`."
facet_vertex_sets(polytope) =
    Set(sort!(collect(row(vertex_indices(facets(polytope)), i))) for i in 1:n_facets(polytope))

"""
Exact affine-unimodular equivalence, reduced to combinatorial vertex maps.

A unimodular map induces a face-lattice isomorphism, and every such isomorphism
is one fixed isomorphism followed by a combinatorial automorphism of the
target, so testing those compositions is complete. A vertex-edge graph
isomorphism cannot serve as the fixed one: for non-simple polytopes it need
not preserve faces.
"""
function affine_unimodular_equivalent(source_polytope, target_polytope)
    dim(source_polytope) == dim(target_polytope) || return false
    f_vector(source_polytope) == f_vector(target_polytope) || return false
    source = matrix(QQ, vertices(source_polytope))
    target = matrix(QQ, vertices(target_polytope))
    size(source) == size(target) || return false
    permutations = Polymake.polytope.find_facet_vertex_permutations(
        Oscar.pm_object(source_polytope), Oscar.pm_object(target_polytope),
    )
    isnothing(permutations) && return false
    vertex_permutation = Int.(collect(last(permutations))) .+ 1
    # Polymake does not document the direction; keep the one that maps every
    # source facet onto a target facet.
    target_facets = facet_vertex_sets(target_polytope)
    maps_facets(seed) = all(
        facet -> sort!(seed[facet]) in target_facets,
        facet_vertex_sets(source_polytope),
    )
    seed = findfirst(maps_facets, (vertex_permutation, invperm(vertex_permutation)))
    isnothing(seed) && error("polymake returned a vertex map that does not preserve facets")
    initial = (vertex_permutation, invperm(vertex_permutation))[seed]
    target_group = automorphism_group(target_polytope; action=:on_vertices)
    for automorphism in target_group
        action = [Int(automorphism(index)) for index in 1:nrows(target)]
        affine_map_is_unimodular(source, target, action[initial]) && return true
    end
    return false
end

"Index of a binomial exponent-difference lattice in its saturation, by exact Smith form."
function binomial_lattice_index(reduced::AbstractVector)
    isempty(reduced) && return BigInt(1)
    nvariables = ngens(parent(first(reduced)))
    differences = Vector{Vector{Int}}()
    for polynomial in reduced
        exponent_vectors = collect(exponents(polynomial))
        isempty(exponent_vectors) && continue
        length(exponent_vectors) == 2 ||
            throw(ArgumentError("lattice index requires a monomial-free binomial basis"))
        push!(differences, exponent_vectors[1] - exponent_vectors[2])
    end
    isempty(differences) && return BigInt(1)
    lattice = matrix(ZZ, length(differences), nvariables, reduce(vcat, differences))
    smith = snf(lattice)
    return prod((abs(BigInt(smith[index, index])) for index in 1:min(size(smith)...)
                 if !iszero(smith[index, index])); init=BigInt(1))
end

"Compute algebraic invariants at one exact canonical cell-interior weight."
function canonical_initial_ideal(
    ideal_value, weight;
    compute_saturation::Bool=true,
    compute_degree::Bool=true,
    certified_tropical::Bool=false,
)
    initial_ideal = initial(
        ideal_value, tropical_semiring_map(coefficient_ring(ideal_value), min), QQ.(weight),
    )
    reduced = gens(groebner_basis(initial_ideal; complete_reduction=true))
    term_counts = length.(terms.(reduced))
    binomial = isempty(term_counts) || maximum(term_counts) <= 2
    ring = base_ring(initial_ideal)
    compute_saturation || certified_tropical ||
        throw(ArgumentError("skipping saturation requires a certified tropical weight"))
    saturated = compute_saturation ?
        saturation(initial_ideal, ideal(ring, prod(gens(ring)))) : nothing
    monomial_free = compute_saturation ? !isone(saturated) : true
    prime = is_prime(initial_ideal)
    lattice_index = binomial && monomial_free ? binomial_lattice_index(reduced) : missing
    # QQ primality can conceal conjugate components after extending scalars.
    # A prime binomial torus ideal stays prime over the algebraic closure
    # exactly when its exponent-difference lattice is saturated.
    geometrically_prime = !prime ? false :
        (binomial && monomial_free ? lattice_index == 1 : missing)
    return (
        ideal=initial_ideal,
        saturation=saturated,
        prime=prime,
        binomial_lattice_index=lattice_index,
        geometrically_prime=geometrically_prime,
        binomial=binomial,
        monomial_free=monomial_free,
        torus_saturated=compute_saturation ? saturated == initial_ideal :
            (prime ? true : missing),
        toric=prime && binomial && monomial_free && geometrically_prime === true,
        dimension=dim(initial_ideal),
        degree=compute_degree ? degree(initial_ideal) : missing,
        groebner_generators=length(reduced),
        maximum_terms=isempty(term_counts) ? 0 : maximum(term_counts),
    )
end

"""
    classify_degenerations(ideal, fan, normalized_rays, symmetry; config)

Classify one exact equal-ray-sum weight for each coarsest-cell orbit. The
routine intentionally makes no assertion that an initial ideal is constant on
the whole coarsest support cell.
"""
function classify_degenerations(
    ideal_value,
    fan,
    normalized_rays,
    symmetry;
    config::DegenerationConfig=DegenerationConfig(),
    on_record::Union{Nothing, Function}=nothing,
    skip_orbits=Set{Int}(),
)
    incidence = maximal_cones(IncidenceMatrix, fan)
    cells = [sort!(collect(row(incidence, index))) for index in 1:nrows(incidence)]
    weights = canonical_cell_weights(fan, normalized_rays)
    lineality = Matrix{QQFieldElem}(
        matrix(QQ, rays_modulo_lineality(fan).lineality_basis),
    )
    isnothing(config.expected_orbits) || length(symmetry.cell_orbits) == config.expected_orbits ||
        error("expected $(config.expected_orbits) coarsest-cell orbits, found $(length(symmetry.cell_orbits))")
    sort!(reduce(vcat, symmetry.cell_orbits)) == collect(eachindex(cells)) ||
        error("symmetry orbits do not cover all coarsest cells")
    verify_canonical_weight_action(
        weights, symmetry, lineality;
        coordinate_indices=config.weight_coordinate_indices,
    ) ||
        error("canonical cell weights are not equivariant modulo lineality")
    verify_ideal_preservation(ideal_value, symmetry.permutations, config.coordinate_images) ||
        error("the supplied taxon action does not preserve the defining ideal")

    representatives = Any[]
    cell_to_orbit = zeros(Int, length(cells))
    progress = TaskProgress("degenerations",
        count(index -> !(index in skip_orbits), eachindex(symmetry.cell_orbits)))
    for (orbit_index, orbit) in enumerate(symmetry.cell_orbits)
        cell_to_orbit[orbit] .= orbit_index
        orbit_index in skip_orbits && continue
        representative = first(orbit)
        algebra = canonical_initial_ideal(
            ideal_value, weights[representative];
            compute_saturation=config.compute_saturation,
            compute_degree=config.compute_degree,
            certified_tropical=config.certified_tropical,
        )
        # Only a prime initial ideal has a Newton-Okounkov body.
        body = config.compute_bodies && algebra.prime ?
            intrinsic_cell_body(lineality, normalized_rays, cells[representative];
                                gradings=config.body_gradings) : nothing
        chevalley = body === nothing ? missing : config.chevalley_test(body)
        polytope_type(P) = RecursivePolytopeTypes.describe_type(
            RecursivePolytopeTypes.determine_polytope(P))
        typed = body !== nothing && algebra.toric
        body_type = typed ? polytope_type(body.polytope) : missing
        slice_types = typed && !isempty(body.slices) ?
            join((polytope_type(slice.polytope) for slice in body.slices), " + ") : missing
        record = (
            orbit=orbit_index,
            representative=representative,
            cells=orbit,
            rays=cells[representative],
            canonical_weight=weights[representative],
            algebra=algebra,
            body=body,
            chevalley_equivalent=chevalley,
            body_type=body_type,
            body_slice_types=slice_types,
        )
        isnothing(on_record) ? push!(representatives, record) : on_record(record)
        advance!(progress)
    end
    return (
        policy=:canonical_coarsest_cell_interior,
        representatives=representatives,
        cell_to_orbit=cell_to_orbit,
        canonical_weights=weights,
    )
end
