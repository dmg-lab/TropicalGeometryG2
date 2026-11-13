#!/usr/bin/env julia

using LinearAlgebra
using Oscar
using TropicalGeometryG2

const DATA_DIRECTORY = data_path("g2", "grassmannian")
const REPORT_PATH = joinpath(@__DIR__, "REPORT.md")
const PROJECTION_PATH = joinpath(@__DIR__, "projected_rays.csv")
const HEXAGON = [1, 5, 3, 4, 2, 6]

compose(left, right) = left[right]

function permutation_power(permutation, exponent::Integer)
    exponent >= 0 || throw(ArgumentError("the exponent must be nonnegative"))
    result = collect(eachindex(permutation))
    for _ in 1:exponent
        result = compose(permutation, result)
    end
    return result
end

function cycle_type(permutation)
    seen = falses(length(permutation))
    lengths = Int[]
    for seed in eachindex(permutation)
        seen[seed] && continue
        length = 0
        current = seed
        while !seen[current]
            seen[current] = true
            length += 1
            current = permutation[current]
        end
        push!(lengths, length)
    end
    sort!(lengths; rev=true)
    return join(lengths, "+")
end

"Matrix for the row-vector pullback used by `coordinate_permutation`."
function coordinate_action_matrix(permutation)
    indices = coordinate_permutation(permutation)
    action = zeros(Int, length(indices), length(indices))
    for column in eachindex(indices)
        action[indices[column], column] = 1
    end
    return action
end

function fixed_codimension(action::AbstractMatrix)
    identity = Matrix{Int}(I, size(action)...)
    return rank(matrix(QQ, identity - action))
end

function quotient_fixed_codimension(action::AbstractMatrix, lineality::AbstractMatrix)
    identity = Matrix{Int}(I, size(action)...)
    difference = matrix(QQ, identity - action)
    return rank(vcat(matrix(QQ, lineality), difference)) - rank(matrix(QQ, lineality))
end

function fan_span_fixed_codimension(action::AbstractMatrix, lineality::AbstractMatrix, rays)
    identity = Matrix{Int}(I, size(action)...)
    moved_rays = matrix(QQ, rays) * matrix(QQ, identity - action)
    return rank(vcat(matrix(QQ, lineality), moved_rays)) - rank(matrix(QQ, lineality))
end

function row_span_trace(action::AbstractMatrix, spanning_rows)
    basis = TropicalGeometryG2.canonical_row_span(spanning_rows)
    pivots = [findfirst(!iszero, collect(row)) for row in eachrow(basis)]
    transformed = basis * QQ.(action)
    restriction = transformed[:, pivots]
    return sum(restriction[index, index] for index in axes(restriction, 1); init=QQ(0))
end

const CLASS_LABELS = ["1", "r,+r^5", "r^2,+r^4", "r^3", "vertex refl.", "edge refl."]
const CLASS_SIZES = [1, 2, 2, 1, 3, 3]
const IRREDUCIBLE_CHARACTERS = [
    (name="A1", dimension=1, values=[1, 1, 1, 1, 1, 1]),
    (name="A2", dimension=1, values=[1, 1, 1, 1, -1, -1]),
    (name="B1", dimension=1, values=[1, -1, 1, -1, 1, -1]),
    (name="B2", dimension=1, values=[1, -1, 1, -1, -1, 1]),
    (name="E1", dimension=2, values=[2, 1, -1, -2, 0, 0]),
    (name="E2", dimension=2, values=[2, -1, -1, 2, 0, 0]),
]

function decompose_character(character)
    multiplicities = Dict{String, Int}()
    for irreducible in IRREDUCIBLE_CHARACTERS
        numerator = sum(CLASS_SIZES .* character .* irreducible.values)
        @assert iszero(numerator % 12)
        multiplicities[irreducible.name] = numerator ÷ 12
    end
    @assert sum(
        irreducible.dimension * multiplicities[irreducible.name]
        for irreducible in IRREDUCIBLE_CHARACTERS
    ) == first(character)
    return multiplicities
end

function format_decomposition(multiplicities)
    terms = String[]
    for irreducible in IRREDUCIBLE_CHARACTERS
        multiplicity = multiplicities[irreducible.name]
        iszero(multiplicity) && continue
        push!(terms, multiplicity == 1 ? irreducible.name : "$(multiplicity)$(irreducible.name)")
    end
    return join(terms, " + ")
end

function character_table(characters)
    lines = [
        "| representation | " * join(CLASS_LABELS, " | ") * " | decomposition |",
        "|---|" * join(fill("---:", length(CLASS_LABELS)), "|") * "|---|",
    ]
    for representation in characters
        decomposition = format_decomposition(decompose_character(representation.values))
        push!(lines,
            "| $(representation.name) | " * join(representation.values, " | ") * " | `$decomposition` |")
    end
    return join(lines, '\n')
end

function orthogonal_projector(spanning_rows)
    basis = TropicalGeometryG2.canonical_row_span(spanning_rows)
    gram_inverse = Matrix{QQFieldElem}(inv(matrix(QQ, basis * transpose(basis))))
    return transpose(basis) * gram_inverse * basis
end

function projective_key(coordinates; oriented=true)
    pivot = findfirst(!iszero, coordinates)
    pivot === nothing && return nothing
    denominator = oriented ? abs(coordinates[pivot]) : coordinates[pivot]
    return Tuple(coordinates ./ denominator)
end

function row_action(vector, action)
    return [
        sum(vector[index] * action[index, column] for index in eachindex(vector); init=QQ(0))
        for column in axes(action, 2)
    ]
end

function direction_orbits(directions, actions)
    unseen = Set(directions)
    orbits = Vector{Vector{eltype(directions)}}()
    while !isempty(unseen)
        seed = first(unseen)
        orbit = sort!(collect(Set(
            projective_key(row_action(collect(seed), action); oriented=true)
            for action in actions
        )))
        foreach(direction -> delete!(unseen, direction), orbit)
        push!(orbits, orbit)
    end
    sort!(orbits; by=orbit -> (length(orbit), first(orbit)))
    return orbits
end

"The primitive integer vector of the rational direction `direction`, as a tuple of `Int`s."
primitive_integer_direction(direction) = Tuple(Int.(primitive_integer(collect(direction))))


function lattice_primitive_direction(direction, plane_basis)
    primitive = primitive_integer_direction(direction)
    common_denominator = foldl(lcm, Int.(denominator.(plane_basis)); init=1)
    for multiplier in 1:common_denominator
        candidate = multiplier .* collect(primitive)
        ambient_vector = row_action(candidate, plane_basis)
        all(denominator(coordinate) == 1 for coordinate in ambient_vector) &&
            return Tuple(candidate)
    end
    error("failed to find an ambient-integral vector on a rational direction")
end

metric_pairing(left, right, metric) = sum(
    left[i] * metric[i, j] * right[j]
    for i in eachindex(left), j in eachindex(right);
    init=QQ(0),
)

function induced_linear_fan_action(action, rays, cells, lineality)
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    signatures = [quotient_ray_signature(row, quotient) for row in eachrow(rays)]
    signature_to_ray = Dict(signature => index for (index, signature) in enumerate(signatures))
    ray_action = Int[]
    for row in eachrow(rays)
        image = row_action(collect(row), action)
        signature = quotient_ray_signature(image, quotient)
        haskey(signature_to_ray, signature) || return nothing
        push!(ray_action, signature_to_ray[signature])
    end
    cell_index = Dict(Tuple(cell) => index for (index, cell) in enumerate(cells))
    cell_action = TropicalGeometryG2.induced_cell_action(ray_action, cells, cell_index)
    cell_action === nothing && return nothing
    return (rays=ray_action, cells=cell_action)
end

function markdown_table(rows)
    lines = [
        "| kind | exponent | order | taxon cycle type | codim Fix in R^21 | codim Fix on fan span/L | codim Fix on R^21/L |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    append!(lines, [
        "| $(row.kind) | $(row.exponent) | $(row.order) | `$(row.cycle_type)` | $(row.ambient_codimension) | $(row.fan_span_codimension) | $(row.quotient_codimension) |"
        for row in rows
    ])
    return join(lines, '\n')
end

fan = load(joinpath(DATA_DIRECTORY, "03_normalized_coarsest_fan.mrdi"))
rays = load(joinpath(DATA_DIRECTORY, "03_normalized_rays.mrdi"))
stored_permutations = [
    Int.(collect(row))
    for row in eachrow(load(joinpath(DATA_DIRECTORY, "05_symmetry_taxon_permutations.mrdi")))
]
lineality = Matrix{QQFieldElem}(matrix(QQ, rays_modulo_lineality(fan).lineality_basis))
cells = TropicalGeometryG2.fan_cells(fan)

# The six-cycle follows the intrinsic hexagon 1-5-3-4-2-6.  The second
# generator is reflection in the axis through the opposite vertices 1 and 4.
rotation = [5, 6, 4, 2, 3, 1, 7]
reflection = [1, 3, 2, 4, 6, 5, 7]
identity = collect(1:7)

@assert permutation_order(rotation) == 6
@assert permutation_order(reflection) == 2
@assert compose(reflection, compose(rotation, reflection)) == invperm(rotation)

group = generated_subgroup([rotation, reflection])
@assert length(group) == 12
@assert Set(Tuple.(group)) == Set(Tuple.(stored_permutations))

# This is an independent exact check against the normalized fan, not merely a
# check that the generators occur in the saved permutation list.
fan_action = taxon_symmetry(
    rays,
    cells;
    candidates=group,
    lineality=rays_modulo_lineality(fan).lineality_basis,
)
@assert length(fan_action.permutations) == 12
@assert all(preserves_lineality(lineality, coordinate_permutation(element)) for element in group)

abstract_reflections = [
    compose(permutation_power(rotation, exponent), reflection)
    for exponent in 0:5
]
@assert length(Set(Tuple.(abstract_reflections))) == 6
@assert all(permutation_order(element) == 2 for element in abstract_reflections)

rows = NamedTuple[]
for exponent in 0:5
    element = permutation_power(rotation, exponent)
    action = coordinate_action_matrix(element)
    push!(rows, (
        kind=exponent == 0 ? "identity" : "rotation",
        exponent=exponent,
        order=permutation_order(element),
        cycle_type=cycle_type(element),
        ambient_codimension=fixed_codimension(action),
        fan_span_codimension=fan_span_fixed_codimension(action, lineality, rays),
        quotient_codimension=quotient_fixed_codimension(action, lineality),
    ))
end
for exponent in 0:5
    element = abstract_reflections[exponent + 1]
    action = coordinate_action_matrix(element)
    push!(rows, (
        kind="hexagon reflection",
        exponent=exponent,
        order=2,
        cycle_type=cycle_type(element),
        ambient_codimension=fixed_codimension(action),
        fan_span_codimension=fan_span_fixed_codimension(action, lineality, rays),
        quotient_codimension=quotient_fixed_codimension(action, lineality),
    ))
end

nonidentity_codimensions = [
    fixed_codimension(coordinate_action_matrix(element))
    for element in group if element != identity
]
ray_span_modulo_lineality = rank(vcat(matrix(QQ, lineality), matrix(QQ, rays))) - rank(matrix(QQ, lineality))
quotient_dimension = size(rays, 2) - rank(matrix(QQ, lineality))

class_representatives = [
    identity,
    rotation,
    permutation_power(rotation, 2),
    permutation_power(rotation, 3),
    reflection,
    compose(rotation, reflection),
]
class_actions = coordinate_action_matrix.(class_representatives)
fan_and_lineality = vcat(lineality, Matrix{QQFieldElem}(rays))
characters = [
    (
        name="R^21",
        values=Int[sum(action[index, index] for index in axes(action, 1)) for action in class_actions],
    ),
    (
        name="lineality L",
        values=[Int(row_span_trace(action, lineality)) for action in class_actions],
    ),
    (
        name="fan span modulo L",
        values=[
            Int(row_span_trace(action, fan_and_lineality) - row_span_trace(action, lineality))
            for action in class_actions
        ],
    ),
    (
        name="R^21/L",
        values=[
            Int(sum(action[index, index] for index in axes(action, 1)) - row_span_trace(action, lineality))
            for action in class_actions
        ],
    ),
]

@assert minimum(nonidentity_codimensions) > 1
@assert Set(row.ambient_codimension for row in rows if row.kind == "hexagon reflection") == Set([8, 9])
@assert minimum(row.fan_span_codimension for row in rows if row.kind == "hexagon reflection") > 1
@assert format_decomposition(decompose_character(characters[2].values)) == "A1 + E1"

# The E1 central idempotent is especially simple because its character
# vanishes on reflections.  Removing its lineality summand leaves the unique
# faithful two-dimensional D6 summand in the span of the fan modulo L.
rotation_actions = [QQ.(coordinate_action_matrix(permutation_power(rotation, exponent))) for exponent in 0:5]
e1_projector = (
    2rotation_actions[1] + rotation_actions[2] - rotation_actions[3] -
    2rotation_actions[4] - rotation_actions[5] + rotation_actions[6]
) ./ 6
lineality_projector = orthogonal_projector(lineality)
identity_q = Matrix{QQFieldElem}(I, 21, 21)
residual_e1_projector = e1_projector * (identity_q - lineality_projector)
projected_rays = Matrix{QQFieldElem}(rays) * residual_e1_projector
plane_basis = TropicalGeometryG2.canonical_row_span(projected_rays)
plane_projector = orthogonal_projector(projected_rays)
coxeter_taxon_reflections = [reflection, compose(rotation, reflection)]
ambient_reflections = [
    identity_q - plane_projector +
    plane_projector * QQ.(coordinate_action_matrix(element))
    for element in coxeter_taxon_reflections
]
@assert all(action * action == identity_q for action in ambient_reflections)
@assert all(rank(matrix(QQ, identity_q - action)) == 1 for action in ambient_reflections)
ambient_reflection_fan_actions = [
    induced_linear_fan_action(action, rays, cells, lineality)
    for action in ambient_reflections
]
plane_pivots = [findfirst(!iszero, collect(row)) for row in eachrow(plane_basis)]
plane_coordinates = projected_rays[:, plane_pivots]
nonzero_projected_rays = [collect(row) for row in eachrow(plane_coordinates) if any(!iszero, row)]
oriented_directions = Set(projective_key(row; oriented=true) for row in nonzero_projected_rays)
unoriented_directions = Set(projective_key(row; oriented=false) for row in nonzero_projected_rays)
plane_actions = [
    (plane_basis * QQ.(coordinate_action_matrix(element)))[:, plane_pivots]
    for element in group
]
projected_direction_orbits = direction_orbits(oriented_directions, plane_actions)
direction_counts = Dict(direction => 0 for direction in oriented_directions)
for row in nonzero_projected_rays
    direction_counts[projective_key(row; oriented=true)] += 1
end
direction_multiplicities = sort!(collect(values(direction_counts)))
plane_metric = plane_basis * transpose(plane_basis)
root_direction_orbits = [
    sort!([lattice_primitive_direction(direction, plane_basis) for direction in orbit])
    for orbit in projected_direction_orbits
]
orbit_squared_lengths = [
    only(Set(metric_pairing(root, root, plane_metric) for root in orbit))
    for orbit in root_direction_orbits
]
length_order = sortperm(orbit_squared_lengths)
short_roots = root_direction_orbits[length_order[1]]
long_roots = root_direction_orbits[length_order[2]]
short_squared_length = orbit_squared_lengths[length_order[1]]
long_squared_length = orbit_squared_lengths[length_order[2]]
simple_short = (-1, -1)
simple_long = (3, 0)
cartan_matrix = [
    2  2metric_pairing(simple_short, simple_long, plane_metric) / short_squared_length
    2metric_pairing(simple_long, simple_short, plane_metric) / long_squared_length  2
]

@assert rank(matrix(QQ, residual_e1_projector)) == 4
@assert size(plane_basis, 1) == 2
@assert sort(length.(projected_direction_orbits)) == [6, 6]
@assert direction_multiplicities == fill(3, 12)
@assert long_squared_length / short_squared_length == 3
@assert cartan_matrix == [2 -3; -1 2]

plane_denominator = foldl(lcm, Int.(denominator.(plane_basis)); init=1)
admissible_residues = Set(
    (a, b) for a in 0:(plane_denominator - 1), b in 0:(plane_denominator - 1)
    if all(denominator(coordinate) == 1 for coordinate in row_action((a, b), plane_basis))
)
plane_lattice_basis = [1 1; 2 -1]
@assert plane_denominator == 3
@assert admissible_residues == Set([(0, 0), (1, 1), (2, 2)])
@assert abs(round(Int, det(plane_lattice_basis))) == plane_denominator^2 ÷ length(admissible_residues)
@assert all(denominator(coordinate) == 1 for coordinate in plane_lattice_basis * plane_basis)

root_type = Dict(root => "short" for root in short_roots)
merge!(root_type, Dict(root => "long" for root in long_roots))

# Test the affine-root heuristic.  The A1-isotypic part is fixed by D6 and is
# therefore the only place where a null-root coordinate can live.  Normalize
# every nonzero E1 projection to its primitive root and retain its A1
# coordinates at the same scale.
group_actions = [QQ.(coordinate_action_matrix(element)) for element in group]
a1_projector = reduce(+, group_actions) ./ length(group_actions)
residual_a1_projector = a1_projector * (identity_q - lineality_projector)
projected_invariants = Matrix{QQFieldElem}(rays) * residual_a1_projector
invariant_basis = TropicalGeometryG2.canonical_row_span(projected_invariants)
invariant_pivots = [findfirst(!iszero, collect(row)) for row in eachrow(invariant_basis)]
invariant_coordinates = projected_invariants[:, invariant_pivots]
@assert size(invariant_basis, 1) == 2

invariant_denominator = foldl(lcm, Int.(denominator.(invariant_basis)); init=1)
invariant_admissible_residues = Set(
    (a, b)
    for a in 0:(invariant_denominator - 1), b in 0:(invariant_denominator - 1)
    if all(denominator(coordinate) == 1 for coordinate in row_action((a, b), invariant_basis))
)
invariant_lattice_basis = [1 0; 0 2]
@assert invariant_denominator == 2
@assert invariant_admissible_residues == Set([(0, 0), (1, 0)])
@assert all(denominator(coordinate) == 1 for coordinate in invariant_lattice_basis * invariant_basis)

affine_coordinates = Dict{Tuple{Int, Int}, Vector{Tuple{QQFieldElem, QQFieldElem}}}()
normalized_invariant_by_ray = Dict{Int, Tuple{QQFieldElem, QQFieldElem}}()
for index in axes(plane_coordinates, 1)
    finite_part = collect(plane_coordinates[index, :])
    all(iszero, finite_part) && continue
    root = lattice_primitive_direction(projective_key(finite_part; oriented=true), plane_basis)
    pivot = findfirst(!iszero, root)
    scale = finite_part[pivot] / root[pivot]
    invariant_part = Tuple(collect(invariant_coordinates[index, :]) ./ scale)
    normalized_invariant_by_ray[index] = invariant_part
    push!(get!(affine_coordinates, root, Tuple{QQFieldElem, QQFieldElem}[]), invariant_part)
end
@assert all(length(values) == 3 for values in values(affine_coordinates))

normalized_affine_lifts = Set(
    (root, invariant_part)
    for (root, invariant_parts) in affine_coordinates
    for invariant_part in invariant_parts
)
distinct_invariant_parts = collect(Set(last(lift) for lift in normalized_affine_lifts))
@assert length(normalized_affine_lifts) == 24
@assert length(distinct_invariant_parts) == 4
@assert rank(matrix(QQ, reduce(vcat, permutedims.(collect.(distinct_invariant_parts))))) == 2

zero_finite_invariants = [
    Tuple(collect(invariant_coordinates[index, :]))
    for index in axes(plane_coordinates, 1)
    if all(iszero, plane_coordinates[index, :])
]
zero_finite_directions = Set(
    projective_key(collect(coordinates); oriented=false)
    for coordinates in zero_finite_invariants
)
zero_lattice_directions = sort!([
    lattice_primitive_direction(direction, invariant_basis)
    for direction in zero_finite_directions
])
@assert length(zero_finite_invariants) == 15
@assert length(zero_finite_directions) == 3

rank_one_reflections_preserve_fan = .!isnothing.(ambient_reflection_fan_actions)
@assert rank_one_reflections_preserve_fan == [false, false]

# Incidence test for the three possible fixed directions.  A useful affine
# null direction should meet both root lengths uniformly and its star should
# look like a real-root interface rather than a zero-root boundary stratum.
root_by_ray = Dict(
    index => lattice_primitive_direction(
        projective_key(collect(plane_coordinates[index, :]); oriented=true),
        plane_basis,
    )
    for index in axes(plane_coordinates, 1)
    if any(!iszero, plane_coordinates[index, :])
)
zero_direction_indices = Dict{Tuple{Int, Int}, Vector{Int}}()
for index in axes(plane_coordinates, 1)
    all(iszero, plane_coordinates[index, :]) || continue
    direction = primitive_integer_direction(
        projective_key(collect(invariant_coordinates[index, :]); oriented=false),
    )
    push!(get!(zero_direction_indices, direction, Int[]), index)
end

ray_cell_degrees = [count(cell -> index in cell, cells) for index in axes(rays, 1)]
function restricted_orbit_sizes(indices, actions)
    unseen = Set(indices)
    sizes = Int[]
    while !isempty(unseen)
        seed = first(unseen)
        orbit = Set(action[seed] for action in actions)
        @assert orbit ⊆ Set(indices)
        foreach(index -> delete!(unseen, index), orbit)
        push!(sizes, length(orbit))
    end
    return sort!(sizes)
end

function candidate_incidence(direction, indices)
    touched_cells = [cell for cell in cells if !isempty(intersect(cell, indices))]
    neighboring_roots = Set(
        root_by_ray[index]
        for cell in touched_cells for index in cell
        if haskey(root_by_ray, index)
    )
    short_directions = count(root -> root_type[root] == "short", neighboring_roots)
    long_directions = count(root -> root_type[root] == "long", neighboring_roots)
    real_counts = [count(index -> haskey(root_by_ray, index), cell) for cell in touched_cells]
    return (
        direction=direction,
        indices=sort(indices),
        orbit_sizes=restricted_orbit_sizes(indices, fan_action.ray_actions),
        degrees=sort(ray_cell_degrees[indices]),
        cells=length(touched_cells),
        zero_only=count(iszero, real_counts),
        all_two_real=all(==(2), real_counts),
        short_directions=short_directions,
        long_directions=long_directions,
    )
end
candidate_summaries = sort!([
    candidate_incidence(direction, indices)
    for (direction, indices) in zero_direction_indices
]; by=summary -> summary.direction)

balanced_incidence_direction = only([
    summary.direction for summary in candidate_summaries
    if summary.short_directions == 6 && summary.long_directions == 6 && summary.all_two_real
])
@assert balanced_incidence_direction == (1, 8)
fixed_fan_ray = only([
    index for index in axes(rays, 1)
    if all(action[index] == index for action in fan_action.ray_actions)
])
fixed_ray_cells = [cell for cell in cells if fixed_fan_ray in cell]
@assert fixed_fan_ray == 37
@assert 37 in zero_direction_indices[(1, -6)]
@assert length(fixed_ray_cells) == 9
@assert all(cell -> all(index -> !haskey(root_by_ray, index), cell), fixed_ray_cells)

function lift_positions(root_kind)
    return sort!(collect(Set(
        invariant_part
        for (root, invariant_parts) in affine_coordinates
        if root_type[root] == root_kind
        for invariant_part in invariant_parts
    )))
end
add_tuple(left, right) = Tuple(left[index] + right[index] for index in eachindex(left))
subtract_tuple(left, right) = Tuple(left[index] - right[index] for index in eachindex(left))
scale_tuple(scale, coordinates) = Tuple(scale * coordinate for coordinate in coordinates)
collinear(left, right) = left[1] * right[2] == left[2] * right[1]

short_positions = lift_positions("short")
long_positions = lift_positions("long")
@assert length(short_positions) == 2
@assert length(long_positions) == 2
short_lift_step = primitive_integer_direction(projective_key(
    collect(subtract_tuple(short_positions[2], short_positions[1])); oriented=false,
))
long_lift_step = primitive_integer_direction(projective_key(
    collect(subtract_tuple(long_positions[2], long_positions[1])); oriented=false,
))
@assert short_lift_step == balanced_incidence_direction
@assert long_lift_step == (3, 10)
@assert long_lift_step ∉ keys(zero_direction_indices)

# The unique globally fixed ray singles out c as the homogenizing direction.
# Quotienting the two-dimensional A1 part by Qc makes the other two zero-root
# directions equal, and cone incidence singles out their common class as delta.
homogenizing_direction = (1, -6)
delta_direction = balanced_incidence_direction
other_delta_direction = (3, -4)
@assert other_delta_direction == add_tuple(delta_direction, scale_tuple(2, homogenizing_direction))

short_level_zero = only(position for position in short_positions if collinear(position, homogenizing_direction))
short_level_one = only(position for position in short_positions if position != short_level_zero)
long_level_zero = only(position for position in long_positions if collinear(position, homogenizing_direction))
long_level_minus_three = only(position for position in long_positions if position != long_level_zero)
delta_vector = scale_tuple(QQ(5) / QQ(21), delta_direction)

@assert short_level_zero == scale_tuple(QQ(10) / QQ(21), homogenizing_direction)
@assert short_level_one == add_tuple(short_level_zero, delta_vector)
@assert long_level_zero == scale_tuple(-QQ(5) / QQ(14), homogenizing_direction)
@assert long_level_minus_three == add_tuple(
    scale_tuple(-QQ(5) / QQ(7), homogenizing_direction),
    scale_tuple(-QQ(5) / QQ(7), delta_direction),
)

equivalent_mod_homogenizing(left, right) =
    collinear(subtract_tuple(left, right), homogenizing_direction)
affine_level(position) = QQ(3) / QQ(10) * (6position[1] + position[2])
@assert equivalent_mod_homogenizing(short_level_zero, (0, 0))
@assert equivalent_mod_homogenizing(short_level_one, delta_vector)
@assert equivalent_mod_homogenizing(long_level_zero, (0, 0))
@assert equivalent_mod_homogenizing(long_level_minus_three, scale_tuple(-3, delta_vector))
@assert affine_level(homogenizing_direction) == 0
@assert affine_level(delta_vector) == 1
@assert affine_level.([
    short_level_zero, short_level_one, long_level_zero, long_level_minus_three,
]) == [0, 1, 0, -3]

# A twisted affine simple system.  The finite parts give D4^(3), while their
# A1 coordinates satisfy the null-root relation after quotienting by Qc.
twisted_simple_roots = [(-1, 2), (-1, -1), (3, 0)]
twisted_simple_invariants = [short_level_one, short_level_zero, long_level_zero]
twisted_cartan_matrix = [
    2metric_pairing(twisted_simple_roots[i], twisted_simple_roots[j], plane_metric) /
    metric_pairing(twisted_simple_roots[i], twisted_simple_roots[i], plane_metric)
    for i in 1:3, j in 1:3
]
twisted_kac_labels = [1, 2, 1]
twisted_finite_relation = Tuple(
    sum(twisted_kac_labels[index] * twisted_simple_roots[index][coordinate] for index in 1:3)
    for coordinate in 1:2
)
twisted_invariant_relation = Tuple(
    sum(twisted_kac_labels[index] * twisted_simple_invariants[index][coordinate] for index in 1:3)
    for coordinate in 1:2
)
@assert twisted_cartan_matrix == [2 -1 0; -1 2 -3; 0 -1 2]
@assert twisted_finite_relation == (0, 0)
@assert equivalent_mod_homogenizing(twisted_invariant_relation, delta_vector)
@assert subtract_tuple(twisted_invariant_relation, delta_vector) ==
        scale_tuple(QQ(15) / QQ(14), homogenizing_direction)

function rays_at_lift(root_kind, position)
    return sort!([
        index for (index, invariant_part) in normalized_invariant_by_ray
        if root_type[root_by_ray[index]] == root_kind && invariant_part == position
    ])
end
short_level_zero_rays = rays_at_lift("short", short_level_zero)
short_level_one_rays = rays_at_lift("short", short_level_one)
long_level_zero_rays = rays_at_lift("long", long_level_zero)
long_level_minus_three_rays = rays_at_lift("long", long_level_minus_three)
@assert length(short_level_zero_rays) == 6
@assert length(short_level_one_rays) == 12
@assert length(long_level_zero_rays) == 12
@assert length(long_level_minus_three_rays) == 6
for root in short_roots
    @assert count(index -> root_by_ray[index] == root, short_level_zero_rays) == 1
    @assert count(index -> root_by_ray[index] == root, short_level_one_rays) == 2
end
for root in long_roots
    @assert count(index -> root_by_ray[index] == root, long_level_zero_rays) == 2
    @assert count(index -> root_by_ray[index] == root, long_level_minus_three_rays) == 1
end

function rays_for_root_lift(root, position)
    return sort!([
        index for (index, invariant_part) in normalized_invariant_by_ray
        if root_by_ray[index] == root && invariant_part == position
    ])
end
alpha_zero_rays = rays_for_root_lift(twisted_simple_roots[1], short_level_one)
alpha_one_rays = rays_for_root_lift(twisted_simple_roots[2], short_level_zero)
alpha_two_rays = rays_for_root_lift(twisted_simple_roots[3], long_level_zero)
@assert alpha_zero_rays == [25, 46]
@assert alpha_one_rays == [45]
@assert alpha_two_rays == [1, 39]
affine_simple_realizations = collect(Iterators.product(
    alpha_zero_rays, alpha_one_rays, alpha_two_rays,
))
cofacial_affine_simple_realizations = [
    realization for realization in affine_simple_realizations
    if any(cell -> all(index -> index in cell, realization), cells)
]
level_zero_mixed_cells = [
    cell for cell in cells
    if !isempty(intersect(cell, short_level_zero_rays)) &&
       !isempty(intersect(cell, long_level_zero_rays))
]
@assert length(affine_simple_realizations) == 4
@assert isempty(cofacial_affine_simple_realizations)
@assert isempty(level_zero_mixed_cells)

# The affine reflections do generate the expected infinite Coxeter group on
# root-level coordinates (x,y,n), but the affine reflection s0 does not
# preserve the finite 24-root window cut out by the fan.
affine_reflections = [
    [1 0 0; 1 -1 -1; 0 0 1],
    [0 -1 0; -1 0 0; 0 0 1],
    [-1 0 0; -1 1 0; 0 0 1],
]
affine_identity = Matrix{Int}(I, 3, 3)
integer_row_action(vector, action) = Tuple(
    sum(vector[index] * action[index, column] for index in eachindex(vector); init=0)
    for column in axes(action, 2)
)
function matrix_order(action, bound)
    power = copy(affine_identity)
    for order in 1:bound
        power *= action
        power == affine_identity && return order
    end
    return nothing
end
@assert all(action * action == affine_identity for action in affine_reflections)
@assert matrix_order(affine_reflections[1] * affine_reflections[2], 6) == 3
@assert matrix_order(affine_reflections[1] * affine_reflections[3], 6) == 2
@assert matrix_order(affine_reflections[2] * affine_reflections[3], 6) == 6
coxeter_translation = (affine_reflections[1] * affine_reflections[2] * affine_reflections[3])^2
@assert coxeter_translation == [1 0 1; 0 1 0; 0 0 1]
affine_root_window = Set(vcat(
    [(root[1], root[2], level) for root in short_roots for level in (0, 1)],
    [(root[1], root[2], level) for root in long_roots for level in (-3, 0)],
))
s_zero_preserved_roots = count(
    root -> integer_row_action(root, affine_reflections[1]) in affine_root_window,
    affine_root_window,
)
@assert length(affine_root_window) == 24
@assert s_zero_preserved_roots == 12

# The delta rays have exact K_3,3 real-root links.  Each side consists of two
# long roots and one short root, with lambda_1 + lambda_2 = -3 alpha.
delta_ray_indices = sort(zero_direction_indices[delta_direction])
@assert delta_ray_indices == [18, 35, 36]
@assert all(==(9), ray_cell_degrees[delta_ray_indices])
@assert all(cell -> count(index -> index in delta_ray_indices, cell) <= 1, cells)
expected_delta_link_parts = Dict(
    18 => ([6, 7, 17], [13, 40, 41]),
    35 => ([6, 12, 44], [11, 20, 40]),
    36 => ([7, 11, 45], [12, 19, 41]),
)
for delta_ray in delta_ray_indices
    star = [cell for cell in cells if delta_ray in cell]
    @assert sort(length.(star)) == vcat(fill(3, 8), [4])
    @assert all(cell -> count(index -> haskey(root_by_ray, index), cell) == 2, star)
    link_edges = Set(
        Tuple(sort([index for index in cell if haskey(root_by_ray, index)]))
        for cell in star
    )
    left, right = expected_delta_link_parts[delta_ray]
    expected_edges = Set(Tuple(sort([a, b])) for a in left for b in right)
    @assert link_edges == expected_edges
    @assert Set(root_by_ray[index] for index in left) ==
            Set(scale_tuple(-1, root_by_ray[index]) for index in right)
    for part in (left, right)
        short = only(root_by_ray[index] for index in part if root_type[root_by_ray[index]] == "short")
        longs = [root_by_ray[index] for index in part if root_type[root_by_ray[index]] == "long"]
        @assert length(longs) == 2
        @assert add_tuple(add_tuple(longs[1], longs[2]), scale_tuple(3, short)) == (0, 0)
    end
end

delta_four_cells = sort!([
    cell for cell in cells
    if length(cell) == 4 && count(index -> index in delta_ray_indices, cell) == 1
])
@assert delta_four_cells == [[13, 14, 17, 18], [19, 21, 36, 45], [20, 22, 35, 44]]
rational_rays = Matrix{QQFieldElem}(rays)
lineality_rank = rank(matrix(QQ, lineality))
for cell in delta_four_cells
    delta_ray = only(index for index in cell if index in delta_ray_indices)
    real_rays = [index for index in cell if haskey(root_by_ray, index)]
    baseline_ray = only(index for index in cell if !haskey(root_by_ray, index) && index != delta_ray)
    relation = (
        rational_rays[real_rays[1], :] + rational_rays[real_rays[2], :] -
        6rational_rays[baseline_ray, :] - 6rational_rays[delta_ray, :]
    )
    @assert rank(matrix(QQ, vcat(lineality, permutedims(relation)))) == lineality_rank
end

homogenizing_rays = sort!([
    index for (index, position) in normalized_invariant_by_ray
    if collinear(position, homogenizing_direction)
])
@assert homogenizing_rays == sort(vcat(short_level_zero_rays, long_level_zero_rays))

candidate_by_direction = Dict(summary.direction => summary for summary in candidate_summaries)
@assert candidate_by_direction[(1, -6)].indices == [8, 9, 10, 30, 33, 34, 37]
@assert candidate_by_direction[(1, -6)].orbit_sizes == [1, 3, 3]
@assert candidate_by_direction[(1, -6)].degrees == vcat([9], fill(21, 6))
@assert candidate_by_direction[(3, -4)].indices == [14, 15, 21, 22, 29]
@assert candidate_by_direction[(3, -4)].orbit_sizes == [2, 3]
@assert candidate_by_direction[(1, 8)].orbit_sizes == [3]

candidate_table = join(vcat(
    [
        "| fixed direction | ray indices | D6 orbit sizes | degrees | maximal cells | zero-only | short/long directions | always two real roots? |",
        "|---|---|---|---|---:|---:|---:|---|",
    ],
    [
        "| `$(summary.direction)` | `$(summary.indices)` | `$(summary.orbit_sizes)` | `$(summary.degrees)` | $(summary.cells) | $(summary.zero_only) | $(summary.short_directions)/$(summary.long_directions) | $(summary.all_two_real) |"
        for summary in candidate_summaries
    ],
), '\n')

delta_link_table = join(vcat(
    [
        "| delta ray | K3,3 bipartition of its real-root link |",
        "|---:|---|",
    ],
    [
        "| $index | `$(expected_delta_link_parts[index][1])` versus `$(expected_delta_link_parts[index][2])` |"
        for index in delta_ray_indices
    ],
), '\n')

open(PROJECTION_PATH, "w") do io
    println(io, "ray,x,y,direction_x,direction_y,root_type,a1_x,a1_y,normalized_a1_x,normalized_a1_y")
    for (index, coordinates) in enumerate(eachrow(plane_coordinates))
        vector = collect(coordinates)
        invariant_part = collect(invariant_coordinates[index, :])
        if all(iszero, vector)
            println(io, join((
                index, 0, 0, "", "", "zero",
                TropicalGeometryG2.rational_key(invariant_part[1]),
                TropicalGeometryG2.rational_key(invariant_part[2]),
                "", "",
            ), ','))
            continue
        end
        direction = lattice_primitive_direction(projective_key(vector; oriented=true), plane_basis)
        pivot = findfirst(!iszero, direction)
        scale = vector[pivot] / direction[pivot]
        normalized_invariant_part = invariant_part ./ scale
        println(io, join((
            index,
            TropicalGeometryG2.rational_key(vector[1]),
            TropicalGeometryG2.rational_key(vector[2]),
            direction[1],
            direction[2],
            root_type[direction],
            TropicalGeometryG2.rational_key(invariant_part[1]),
            TropicalGeometryG2.rational_key(invariant_part[2]),
            TropicalGeometryG2.rational_key(normalized_invariant_part[1]),
            TropicalGeometryG2.rational_key(normalized_invariant_part[2]),
        ), ','))
    end
end

open(REPORT_PATH, "w") do io
    println(io, """# Dihedral symmetries of the G2/P2 fan

## Result

The coarsest fan has a faithful dihedral taxon symmetry

```text
D6 = <r,s | r^6 = s^2 = 1, srs = r^(-1)>,    |D6| = 12,
```

where `r` rotates the taxon hexagon `$(join(HEXAGON, '-'))` and `s` reflects it
through the opposite vertices 1 and 4.  Both fix taxon 7.  The induced
permutations of the 21 Pluecker coordinates preserve the lineality space, all
51 quotient rays, and all 248 maximal cones.  Thus they are exact linear fan
automorphisms of R^21 and give the desired abstract dihedral action on taxa.

There is an important qualification: the six elements called *reflections* in
the abstract dihedral group are not Euclidean hyperplane reflections in this
natural 21-dimensional representation.  Their fixed spaces have codimension 8
or 9 (and not 1).  Their restrictions to the span of the fan, modulo
lineality, still have fixed-space codimension greater than 1, so changing only
their action on an ambient complement cannot repair this.  In fact, no
nonidentity element of this D6 image fixes a hyperplane.  The group is
therefore an orthogonal permutation representation of an abstract reflection
group, but it is not itself a real reflection group inside GL(21) in the
strict Coxeter sense.

Allowing signs in a monomial lift does not change this obstruction.  Each
2-cycle of coordinates in an involution contributes one `-1` eigenvalue,
independently of the compatible signs.  The two taxon-reflection types induce
8 and 9 such coordinate 2-cycles.

The fan rays span $ray_span_modulo_lineality of the $quotient_dimension
dimensions modulo the 3-dimensional lineality space.

## What remains outside lineality

The representation-theoretic decomposition explains precisely where the
literal G2 reflection representation lives.  Here `E1` is the faithful
2-dimensional representation in which `r` is rotation through 60 degrees and
the abstract reflections fix lines; `E2` is the 2-dimensional representation
that factors through a rotation of order 3.  The four `A/B` representations
are 1-dimensional.

$(character_table(characters))

Consequently the lineality space is exactly `A1 + E1`: a fixed scaling
direction plus the expected rank-2 G2 reflection plane.  Outside lineality the
fan itself spans

```text
$(format_decomposition(decompose_character(characters[3].values))),
```

and the full ambient quotient is

```text
$(format_decomposition(decompose_character(characters[4].values))).
```

The entire complement is therefore not a reflection representation: on the
fan span each abstract reflection has a 5-dimensional `-1` eigenspace.
Nevertheless, its lone `E1` summand is canonical.  It is obtained by applying
the `E1` character projector and taking the orthogonal complement of
lineality.

This plane gives an exact second G2 picture.  In the rational basis selected
by row reduction, its intrinsic lattice is

```text
E intersect Z^21 = {(a,b) in Z^2 : a = b mod 3},
```

with basis `(1,1), (2,-1)`.  The primitive lattice vectors on the projected
directions are

```text
short: $short_roots
long:  $long_roots
```

The inherited Gram matrix is

```text
[$(TropicalGeometryG2.rational_key(plane_metric[1, 1]))  $(TropicalGeometryG2.rational_key(plane_metric[1, 2]))]
[$(TropicalGeometryG2.rational_key(plane_metric[2, 1]))  $(TropicalGeometryG2.rational_key(plane_metric[2, 2]))],
```

so the squared long/short length ratio is exactly 3.  Taking the short root
`(-1,-1)` and long root `(3,0)` gives the G2 Cartan matrix

```text
[ 2 -3]
[-1  2].
```

Of the 51 fan rays, $(length(nonzero_projected_rays)) have nonzero image and
15 vanish.  The nonzero images occupy $(length(oriented_directions)) oriented
root directions, split into the short and long hexagons above, and every root
direction is hit by exactly three fan rays.  The ray-by-ray exact coordinates
are recorded in `projected_rays.csv`.

## Hidden affine G2 structure

The next incidence test gives a clean answer: there is an exact finite window
of the twisted affine root system `D4^(3)`, whose affine Weyl group is of G2
type.  It is a root pattern in a quotient of the fan's canonical isotypic
projection, not an affine Weyl symmetry of the finite fan itself.

The relevant four-dimensional projection is

```text
U = 2A1 + E1
```

where `E1` is the finite G2 root plane and `2A1` supplies a homogenizing
coordinate and an affine level.  After normalizing the `E1` part, the 36
nonzero rays give 24 distinct lifts: two lifts over each of the 12 finite G2
roots.

Cone incidence distinguishes three different roles among the fixed
projective directions:

$candidate_table

Set

```text
c = (1,-6),      d = (1,8),      delta = (5/21)[d] in A1^2 / Qc.
```

The line `Qc` is canonically the homogenizing line: it contains the unique fan
ray fixed by all of D6, ray $fixed_fan_ray, and the nine maximal cones in that
ray's star contain only zero-finite rays.  By contrast, the three rays on
`Qd`, namely `$delta_ray_indices`, have 27 incident maximal cones and see all
12 finite roots.  Every such cone contains exactly two real-root rays.  Their
real-root links are exactly:

$delta_link_table

Each link is `K3,3`.  Each side contains two long roots and one short root,
and the exact relation on each side is

```text
lambda_1 + lambda_2 = -3 alpha_short.
```

The three nonsimplicial cones in these stars are

```text
$(join(string.(delta_four_cells), "\n"))
```

and exact arithmetic gives, modulo lineality, the circuit relation

```text
[r_+] + [r_-] = 6([z_+] + [z_-]),
```

where `r_+` and `r_-` are the antipodal real-root rays and `z_+,z_-` are the
two zero-finite rays in the cone.  This is why the `Qd` incidence is an affine
interface, whereas ray $fixed_fan_ray is a vertical apex.

The third direction is not independent after removing the homogenizing line:

```text
(3,-4) = d + 2c.
```

Thus `(3,-4)` and `(1,8)` define the same one-dimensional quotient direction.
In the chosen A1 coordinates, the intrinsic fixed-plane lattice is

```text
Lambda = {(a,b) in Z^2 : b is even}.
```

The primitive quotient coordinate on this lattice is

```text
q(a,b) = 3a + b/2,        ker(q) = Zc.
```

The normalized level functional

```text
n(a,b) = (3/5)q(a,b) = (3/10)(6a+b)
```

kills `c` and sends `delta` to 1.  The four normalized lift positions are

| finite roots | A1 coordinate | level | rays per root |
|---|---|---:|---:|
| short | `(10/21)c` | 0 | 1 |
| short | `(10/21)c + (5/21)d` | 1 | 2 |
| long | `(-5/14)c` | 0 | 2 |
| long | `(-5/7)c - (5/7)d` | -3 | 1 |

Consequently the 24 distinct lifts are exactly

```text
{alpha, alpha+delta : alpha short}
    union {beta, beta-3delta : beta long}.
```

Short-root levels are spaced by 1 and long-root levels by 3.  This is the
twisted affine system `D4^(3)`; untwisted `G2^(1)` would require a long-root
level `+1` or `-1`, and no such lift occurs.

A concrete affine simple system is

```text
alpha_0 = (-1,2) + delta,
alpha_1 = (-1,-1),
alpha_2 = (3,0).
```

It is represented by rays `$(alpha_zero_rays)`, `$(alpha_one_rays)`, and
`$(alpha_two_rays)`, respectively.  Its Cartan matrix is

```text
[ 2 -1  0]
[-1  2 -3]
[ 0 -1  2],
```

and the Kac-label relation is exactly

```text
alpha_0 + 2 alpha_1 + alpha_2 = delta  modulo Qc.
```

On abstract level coordinates `(x,y,n)`, the three corresponding affine
reflections are

```text
s_0(x,y,n) = (x+y, -y, n-y),
s_1(x,y,n) = (-y, -x, n),
s_2(x,y,n) = (-x-y, y, n).
```

They are involutions with Coxeter orders `m_01=3`, `m_02=2`, `m_12=6`.
Moreover, `(s_0 s_1 s_2)^2` is a nontrivial unipotent translation, so they
generate the expected infinite affine Coxeter group.

The finite fan is only a window into that group.  No realization of the three
simple roots above is cofacial; in fact no level-zero short and level-zero
long lift share a maximal cone.  Also `s_0` preserves only
$s_zero_preserved_roots of the 24 displayed real roots.  Hence the finite fan
does not carry the full affine Weyl action, even though its projected rays and
incidence encode an exact `D4^(3)`/affine-G2 fragment.

## Exact certificate

- Rotation: `r = $rotation` (taxon `i` maps to `r[i]`).
- Reflection: `s = $reflection`.
- Exact relations checked: `order(r)=6`, `order(s)=2`, and `s*r*s=r^(-1)`.
- Generated group: 12 elements, exactly the stored fan-symmetry group.
- Fan action checked from the rational ray and maximal-cone data.

$(markdown_table(rows))

Here `codim Fix in R^21` is `rank(I-g)` on the pair-coordinate permutation
representation.  A genuine linear reflection would have value 1.  The last
two columns perform the analogous exact computations on the 12-dimensional
span of the fan rays modulo lineality and on the whole 18-dimensional ambient
quotient, respectively.

## Reproduce

From julia/:

```sh
julia --project=. studies/g2/fan_symmetries/analyze_dihedral_action.jl
```
""")
end

println("verified D6 fan action (order 12)")
println("abstract reflections have ambient fixed codimensions ",
        sort([row.ambient_codimension for row in rows if row.kind == "hexagon reflection"]))
println("residual E1 plane: ", length(nonzero_projected_rays), " nonzero projected rays, ",
        length(oriented_directions), " oriented directions in orbits ",
        sort(length.(projected_direction_orbits)), ", ", length(unoriented_directions), " lines")
println("G2 squared-length ratio ", long_squared_length / short_squared_length,
        ", Cartan matrix ", cartan_matrix)
println("affine test: ", length(normalized_affine_lifts),
        " normalized real-root lifts; fixed-part rank ",
        rank(matrix(QQ, reduce(vcat, permutedims.(collect.(distinct_invariant_parts))))),
        "; rank-one reflections preserve fan ", rank_one_reflections_preserve_fan)
println("intrinsic zero-finite directions: ", zero_lattice_directions)
println("hidden D4^(3) window: homogenizing line ", homogenizing_direction,
        "; delta direction ", delta_direction,
        "; levels short {0,1}, long {0,-3}; fixed fan ray ", fixed_fan_ray)
println("wrote ", REPORT_PATH)
println("wrote ", PROJECTION_PATH)
