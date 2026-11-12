###############################################################################
# Exact G2 Grassmannian data shared by the atlas generators.
#
# Included by to_tex/render_maximal_cell_tree.jl, compute_label_orbits.jl,
# and the flag atlas. It loads the tree library of the fine Groebner reference
# study, reads the W(G2) orbit CSVs written by compute_label_orbits.jl, builds
# the 248 coarsest support cells (the 245 variety triangles plus the three
# merged squares), and writes the censuses of both fields in atlas rays.
###############################################################################

include(joinpath(
    @__DIR__, "..", "..", "reference", "fine_groebner",
    "build_balanced_triangle_atlas.jl",
))

const LABEL_GROUP_CSV = joinpath(@__DIR__, "label_permutation_group.csv")
const RAY_ORBITS_CSV = joinpath(@__DIR__, "ray_label_orbits.csv")
const TRIANGLE_ORBITS_CSV = joinpath(@__DIR__, "triangle_label_orbits.csv")
const PREVARIETY_FAN_PATH =
    TropicalCensus.data_path("g2", "dressian", "01_prevariety_fan.mrdi")

# These are written in the ray indexing of the tropical variety, as in
# observations.md. The prevariety's raw ray order and scaling are recovered
# and certified by prevariety_only_cone_data. The first member is the
# representative.
const EXPECTED_PREVARIETY_ONLY_CONE_ORBITS = (
    (
        kind=:tetrahedron,
        quotient_dimension=4,
        representative=(8, 11, 12, 30),
        members=(
            (8, 11, 12, 30),
            (7, 10, 33, 41),
            (6, 9, 34, 40),
        ),
    ),
    (
        kind=:square,
        quotient_dimension=3,
        representative=(13, 14, 17, 18),
        members=(
            (13, 14, 17, 18),
            (20, 22, 35, 44),
            (19, 21, 36, 45),
        ),
    ),
)

const SUPPORT_SQUARES = EXPECTED_PREVARIETY_ONLY_CONE_ORBITS[2].members

function row_as_vector(matrix, row::Int)
    return [QQ(matrix[row, column]) for column in 1:ncols(matrix)]
end

###############################################################################
# The censuses of both fields in atlas rays
###############################################################################

"""The stage-3 file `name` of the Grassmannian census over `field`, `:qq` or `:gf2`."""
function grassmannian_census_path(field::Symbol, name::AbstractString)
    field in (:qq, :gf2) || throw(ArgumentError("field must be :qq or :gf2, not :$field"))
    return TropicalCensus.data_path(
        "g2", field == :qq ? "grassmannian" : "grassmannian_gf2", name,
    )
end

"""
The rays of the Grassmannian census over `field` in atlas numbering. Atlas
rays 1:51 are the rows of `balanced_rays`; they are the rays of the QQ census,
in its order. The rays that occur only over GF(2) continue the numbering, by
their balanced generators. Returns the ray matrix and the atlas index of every
census ray.
"""
function census_atlas_rays(field::Symbol, balanced_rays)
    census_rays = Matrix{QQFieldElem}(
        load(grassmannian_census_path(field, "03_normalized_rays.mrdi")),
    )
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, balanced_lineality_basis()); side=:right))
    signature(rays, row) =
        TropicalCensus.quotient_ray_signature(row_as_vector(rays, row), quotient)
    atlas_index = Dict(signature(balanced_rays, row) => row for row in 1:nrows(balanced_rays))
    length(atlas_index) == nrows(balanced_rays) || error("the atlas rays are not distinct")

    ray_index = zeros(Int, nrows(census_rays))
    new_rows = Vector{QQFieldElem}[]
    for ray in 1:nrows(census_rays)
        ray_index[ray] = get(atlas_index, signature(census_rays, ray)) do
            push!(new_rows, QQ.(
                balance_outer_edges(row_as_vector(census_rays, ray)).balanced_generator,
            ))
            nrows(balanced_rays) + length(new_rows)
        end
    end
    rays = vcat(
        Matrix{QQFieldElem}(balanced_rays),
        reduce(vcat, permutedims.(new_rows); init=zeros(QQFieldElem, 0, ncols(balanced_rays))),
    )
    sort(ray_index) == collect(1:nrows(rays)) ||
        error("the census rays over $field do not extend the atlas rays bijectively")
    field == :qq && ray_index != collect(1:nrows(balanced_rays)) &&
        error("the QQ census rays are not the atlas rays in their order")
    return (rays=rays, ray_index=ray_index)
end

"""
The maximal cells of the Grassmannian census over `field`, in census order
and written in atlas rays through the `ray_index` of `census_atlas_rays`.
"""
function census_atlas_cells(field::Symbol, ray_index)
    fan = load(grassmannian_census_path(field, "03_normalized_coarsest_fan.mrdi"))
    lineality = matrix(QQ, balanced_lineality_basis())
    census_lineality = matrix(QQ, rays_modulo_lineality(fan).lineality_basis)
    rank(lineality) == rank(census_lineality) == rank(vcat(lineality, census_lineality)) ||
        error("the census fan over $field has another lineality space than the atlas")
    return [Tuple(sort!(ray_index[cell])) for cell in TropicalCensus.fan_cells(fan)]
end

###############################################################################
# Orbit CSVs
###############################################################################

function load_index_orbits(path::AbstractString, expected_kind::Symbol)
    lines = readlines(path)
    isempty(lines) && error("empty orbit CSV: $path")
    header = split(chomp(first(lines)), ',')
    if expected_kind == :ray
        header[1:3] == ["orbit", "size", "rays"] || error("unexpected ray-orbit CSV")
        expected_total = 51
    elseif expected_kind == :triangle
        header[1:3] == ["orbit", "size", "triangles"] || error("unexpected triangle-orbit CSV")
        expected_total = 251
    else
        error("unknown orbit kind: $expected_kind")
    end

    result = Vector{Vector{Int}}()
    for (offset, line) in enumerate(lines[2:end])
        line_number = offset + 1
        fields = split(chomp(line), ','; keepempty=true)
        length(fields) >= 3 || error("short CSV line $line_number in $path")
        orbit_index = parse(Int, fields[1])
        orbit_index == length(result) + 1 || error("nonconsecutive orbit numbering")
        indices = parse.(Int, split(fields[3], ';'))
        parse(Int, fields[2]) == length(indices) || error("wrong orbit size on line $line_number")
        push!(result, indices)
    end
    flattened = reduce(vcat, result)
    sort(flattened) == collect(1:expected_total) || error("orbit CSV does not partition 1:$expected_total")
    length(unique(flattened)) == expected_total || error("duplicate orbit member")
    return result
end

function certify_declared_representatives(path::AbstractString, orbits)
    lines = readlines(path)
    length(lines) == length(orbits) + 1 || error("unexpected orbit CSV length: $path")
    for (orbit_number, line) in enumerate(lines[2:end])
        fields = split(chomp(line), ','; keepempty=true)
        declared = parse(Int, fields[4])
        declared == first(orbits[orbit_number]) || error(
            "orbit $orbit_number declares representative $declared but atlas chose " *
            "$(first(orbits[orbit_number]))",
        )
    end
    return true
end

function label_ray_actions(path::AbstractString=LABEL_GROUP_CSV)
    lines = readlines(path)
    isempty(lines) && error("empty label-action CSV: $path")
    split(chomp(first(lines)), ',') == [
        "element",
        "label_images",
        "cycle_notation",
        "order",
        "ray_permutation",
        "triangle_permutation",
    ] || error("unexpected label-action CSV header: $path")

    actions = Vector{Vector{Int}}()
    for (offset, line) in enumerate(lines[2:end])
        fields = split(chomp(line), ','; keepempty=true)
        length(fields) == 6 || error("label-action CSV line $(offset + 1) is malformed")
        action = parse.(Int, split(fields[5], ';'))
        sort(action) == collect(1:51) || error("ray action $(offset + 1) is not a permutation")
        push!(actions, action)
    end
    length(actions) == 12 || error("expected the 12 elements of W(G2)")
    return actions
end

###############################################################################
# Rooted tree signatures
###############################################################################

"""
Root a seven-taxon tree at taxon 7. Every positive internal edge then defines
the cluster on the side away from the root. This removes all dependence on
the arbitrary internal-vertex labels.
"""
function rooted_internal_edges(tree::TreeMetric)
    result = NamedTuple{(:cluster, :length), Tuple{Vector{Int}, QQFieldElem}}[]
    for split in internal_splits(tree)
        cluster = 7 in split.side ? split.complement : split.side
        @assert !(7 in cluster)
        push!(result, (cluster=sort(cluster), length=split.length))
    end
    sort!(result; by=edge -> (length(edge.cluster), join(edge.cluster, ",")))
    return result
end

cluster_key(cluster) = isempty(cluster) ? "-" : join(cluster, "")

function rational_key(value)
    value = QQ(value)
    denominator(value) == 1 && return string(numerator(value))
    return "$(numerator(value))/$(denominator(value))"
end

"""Scale-invariant exact signature of a rooted metric tree."""
function rooted_metric_signature(tree::TreeMetric, permutation::Vector{Int}=collect(1:7))
    edges = rooted_internal_edges(tree)
    isempty(edges) && error("a fan ray has no positive internal edge")
    total_length = sum((edge.length for edge in edges); init=QQ(0))
    total_length > 0 || error("nonpositive internal metric")
    entries = [
        (
            cluster=sort([permutation[taxon] for taxon in edge.cluster]),
            length=edge.length / total_length,
        )
        for edge in edges
    ]
    sort!(entries; by=entry -> (length(entry.cluster), join(entry.cluster, ",")))
    return join(
        ("$(cluster_key(entry.cluster)):$(rational_key(entry.length))" for entry in entries),
        ";",
    )
end

###############################################################################
# Coarsest support cells and prevariety-only cones
###############################################################################

"""Merge the two variety triangles contained in each declared square."""
function support_maximal_cells(variety_fan)
    triangles = Tuple.(maximal_triangle_ray_indices(variety_fan))
    triangle_sets = Set.(triangles)
    merges = Any[]
    removed_triangle_indices = Set{Int}()

    for square in SUPPORT_SQUARES
        square_set = Set(square)
        indices = [index for (index, triangle) in enumerate(triangle_sets)
                   if triangle <= square_set]
        length(indices) == 2 || error(
            "square $square contains $(length(indices)) variety triangles, expected two",
        )
        left, right = triangles[indices]
        union(Set(left), Set(right)) == square_set || error(
            "the two triangles do not cover square $square",
        )
        diagonal = Tuple(sort!(collect(intersect(Set(left), Set(right)))))
        length(diagonal) == 2 || error("square $square has no common diagonal edge")
        foreach(index -> push!(removed_triangle_indices, index), indices)
        push!(merges, (
            square=square,
            triangle_indices=Tuple(indices),
            triangles=(left, right),
            diagonal=diagonal,
        ))
    end

    length(removed_triangle_indices) == 6 || error("expected six merged triangles")
    cells = Tuple[triangle for (index, triangle) in enumerate(triangles)
                  if !(index in removed_triangle_indices)]
    append!(cells, SUPPORT_SQUARES)
    length(cells) == 248 || error("support fan must have 248 maximal cells")
    return (
        triangles=triangles,
        cells=cells,
        merges=merges,
        removed_triangle_indices=removed_triangle_indices,
    )
end

function support_cell_orbit_specs(construction)
    old_triangle_orbits = load_index_orbits(TRIANGLE_ORBITS_CSV, :triangle)
    removed = construction.removed_triangle_indices
    specs = Any[]
    for (old_orbit_number, indices) in enumerate(old_triangle_orbits)
        removed_here = [index for index in indices if index in removed]
        if !isempty(removed_here)
            length(removed_here) == length(indices) || error(
                "triangle orbit T$old_orbit_number is only partially merged",
            )
            continue
        end
        members = Tuple(construction.triangles[index] for index in indices)
        push!(specs, (
            kind=:triangle,
            old_orbits=(old_orbit_number,),
            sort_key=first(indices),
            representative=first(members),
            members=members,
        ))
    end

    merged_orbits = sort!(unique(
        orbit_number
        for (orbit_number, indices) in enumerate(old_triangle_orbits)
        if any(index -> index in removed, indices)
    ))
    merged_orbits == [17, 19] || error("expected old triangle orbits T17 and T19 to merge")
    push!(specs, (
        kind=:square,
        old_orbits=Tuple(merged_orbits),
        sort_key=minimum(removed),
        representative=first(SUPPORT_SQUARES),
        members=SUPPORT_SQUARES,
    ))
    sort!(specs; by=spec -> spec.sort_key)

    actions = label_ray_actions()
    for spec in specs
        images = Set(
            Tuple(sort(action[collect(spec.representative)]))
            for action in actions
        )
        images == Set(spec.members) || error(
            "support cell representative $(spec.representative) has the wrong W(G2) orbit",
        )
    end
    length(specs) == 31 || error("expected 31 support maximal-cell orbits")
    count(spec -> spec.kind == :triangle, specs) == 30 || error("expected 30 triangle orbits")
    count(spec -> spec.kind == :square, specs) == 1 || error("expected one square orbit")
    sum(spec -> length(spec.members), specs) == 248 || error("support orbits do not cover 248 cells")
    Set(member for spec in specs for member in spec.members) == Set(construction.cells) ||
        error("support maximal-cell orbits do not partition the support fan")
    return specs
end

"""
Match every raw prevariety ray to the primitive-integral ray normalization used
by the atlases. This simultaneously corrects the prevariety's ray order and
its arbitrary positive row scalings.
"""
function prevariety_to_atlas_ray_indices(prevariety, balanced_ray_matrix)
    atlas_rows = [row_vector(balanced_ray_matrix, index)
                  for index in 1:nrows(balanced_ray_matrix)]
    atlas_index = Dict(Tuple(row) => index for (index, row) in enumerate(atlas_rows))
    length(atlas_index) == 51 || error("atlas ray normalizations are not unique")

    raw_rays = rays_modulo_lineality(prevariety).rays_modulo_lineality
    length(raw_rays) == 51 || error("expected 51 prevariety rays")
    result = Int[]
    for (prevariety_index, raw_ray) in enumerate(raw_rays)
        normalized = balance_outer_edges(collect(raw_ray)).balanced_generator
        key = Tuple(normalized)
        haskey(atlas_index, key) || error(
            "renormalized prevariety ray $prevariety_index has no atlas ray",
        )
        push!(result, atlas_index[key])
    end
    sort(result) == collect(1:51) || error("prevariety-to-atlas ray map is not bijective")
    return result
end

function quotient_ray_dimension(rays, balanced_ray_matrix)
    lineality = matrix(QQ, balanced_lineality_basis())
    ray_rows = matrix(
        QQ,
        vectors_as_matrix([row_vector(balanced_ray_matrix, index) for index in rays]),
    )
    return rank(vcat(lineality, ray_rows)) - rank(lineality)
end

function normalized_prevariety_cone_record(rays, balanced_ray_matrix)
    ray_sum = zeros(QQFieldElem, ncols(balanced_ray_matrix))
    for ray_index in rays
        # These rows, rather than rays_modulo_lineality(prevariety), are the
        # primitive-integral Coxeter-balanced normalizations used in the atlas.
        ray_sum .+= row_vector(balanced_ray_matrix, ray_index)
    end
    record = balance_outer_edges(ray_sum)
    roots = [rational_row(COXETER_STAR_ROOTS, index)
             for index in axes(COXETER_STAR_ROOTS, 1)]
    @assert all(iszero, record.root_coefficients)
    @assert all(iszero, record.lineality_leaf_shift)
    @assert record.integer_scale == 1
    @assert minimum(record.pendant_lengths) == 0
    @assert all(iszero(inner_product(record.pendant_lengths, root)) for root in roots)
    @assert all(edge_length >= 0 && denominator(edge_length) == 1
                for edge_length in tree_edge_lengths(record.tree))
    return merge((rays=collect(rays), ray_sum=ray_sum), record)
end

"""Certify and normalize the six maximal prevariety cones absent from the variety fan."""
function prevariety_only_cone_data(variety_fan, balanced_ray_matrix)
    prevariety = load(PREVARIETY_FAN_PATH)
    ray_map = prevariety_to_atlas_ray_indices(prevariety, balanced_ray_matrix)
    incidence = Matrix{Bool}(maximal_cones(IncidenceMatrix, prevariety))
    prevariety_cones = [
        Tuple(sort(ray_map[findall(identity, collect(row))]))
        for row in eachrow(incidence)
    ]

    variety_triangles = Set(Tuple(triangle) for triangle in maximal_triangle_ray_indices(variety_fan))
    all(cone -> length(cone) == 4 || cone in variety_triangles, prevariety_cones) ||
        error("a maximal prevariety triangle is absent from the variety fan")
    exceptional = Set(cone for cone in prevariety_cones if length(cone) == 4)
    expected = Set(member for orbit in EXPECTED_PREVARIETY_ONLY_CONE_ORBITS
                   for member in orbit.members)
    exceptional == expected || error(
        "prevariety-only maximal cones differ from the six declared cones",
    )

    actions = label_ray_actions()
    records = Dict{Tuple{Vararg{Int}}, Any}()
    for orbit in EXPECTED_PREVARIETY_ONLY_CONE_ORBITS
        images = Set(
            Tuple(sort(action[collect(orbit.representative)]))
            for action in actions
        )
        images == Set(orbit.members) || error(
            "declared $(orbit.kind) members do not form one W(G2) orbit",
        )
        for member in orbit.members
            quotient_ray_dimension(member, balanced_ray_matrix) == orbit.quotient_dimension ||
                error("wrong quotient dimension for exceptional cone $member")
            records[member] = normalized_prevariety_cone_record(
                member,
                balanced_ray_matrix,
            )
        end
    end
    length(records) == 6 || error("expected six normalized exceptional-cone records")
    return (
        prevariety=prevariety,
        ray_map=ray_map,
        orbits=EXPECTED_PREVARIETY_ONLY_CONE_ORBITS,
        records=records,
    )
end
