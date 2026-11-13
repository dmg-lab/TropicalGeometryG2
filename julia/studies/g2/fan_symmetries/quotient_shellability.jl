###############################################################################
# Shellability of the tropical G2-Grassmannian and G2 flag variety modulo D12.
#
# Modulo lineality, Trop(G2/P2) and Trop(G2/B) are polyhedral complexes Delta,
# of dimensions 2 and 3, on which the Weyl group D12 acts. The action is not
# regular: some stabilizers fold their cell, and some cells have two vertices
# in one ray orbit, so the cell orbits are not the cells of a regular CW
# complex. On the barycentric subdivision the action is regular, since an
# element fixing a chain of faces fixes each face of it. The orbits of chains
# are therefore the simplices of a regular CW complex homeomorphic to
# |Delta|/D12, whose closed cells are simplices. The script searches that
# complex exhaustively for a shelling, checks the shelling it finds
# independently, and computes the integral homology of Delta and of the
# quotient.
#
# Run from julia/:
#   julia --project=. studies/g2/fan_symmetries/quotient_shellability.jl \
#       [--target=grassmannian|flag] [--fan=coarsest|raw] [--output=FILE]
# --fan=raw takes the simplicial Groebner structure of the Grassmannian (251
# triangles) instead of the coarsest one (245 triangles, 3 quadrilaterals);
# the raw flag fan has rays of its own, with no stored action. --output=FILE
# writes the shelling as CSV.
###############################################################################

using Oscar
using TropicalGeometryG2

const TARGETS = Dict("grassmannian" => "Trop(G2/P2)", "flag" => "Trop(G2/B)")
const FAN_FILES = Dict(
    "coarsest" => "03_normalized_coarsest_fan.mrdi",
    "raw" => "01_raw_fan.mrdi",
)

"Faces modulo lineality as sorted ray sets, by dimension: rays first, maximal cells last."
function fan_faces(fan)
    lineality = lineality_dim(fan)
    return map((lineality + 1):dim(fan)) do dimension
        incidence = IncidenceMatrix(cones(fan, dimension))
        [sort!(collect(row(incidence, index))) for index in 1:nrows(incidence)]
    end
end

"The stored action permutes the normalized rays, so the fan must list them in that order."
function check_ray_order(fan, normalized_rays)
    data = rays_modulo_lineality(fan)
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, data.lineality_basis); side=:right))
    rays = Matrix{QQFieldElem}(matrix(QQ, data.rays_modulo_lineality))
    size(rays, 1) == size(normalized_rays, 1) ||
        error("the fan has $(size(rays, 1)) rays, not $(size(normalized_rays, 1))")
    for index in axes(rays, 1)
        quotient_ray_signature(rays[index, :], quotient) ==
            quotient_ray_signature(normalized_rays[index, :], quotient) ||
            error("ray $index of the fan is not normalized ray $index")
    end
end

"For each face, the faces of one dimension less that it contains; empty for the rays."
function covered_faces(faces)
    # A cone of a fan containing the rays of a smaller cone contains it as a face.
    return [level == 1 ? [Int[] for _ in faces[1]] :
            [[j for (j, small) in enumerate(faces[level - 1]) if issubset(small, big)] for big in faces[level]]
            for level in eachindex(faces)]
end

"The permutation of the faces of each dimension induced by each ray permutation."
function face_actions(faces, ray_actions)
    return map(eachindex(faces)) do level
        index = Dict{Vector{Int}, Int}(face => position for (position, face) in enumerate(faces[level]))
        actions = [TropicalGeometryG2.induced_cell_action(action, faces[level], index) for action in ray_actions]
        any(isnothing, actions) && error("the group does not permute the faces of dimension $(level - 1)")
        actions
    end
end

"Per dimension: face orbits, those whose stabilizer moves one of their rays, and those with two rays in one orbit."
function orbit_cell_defects(faces, actions, ray_actions)
    ray_orbit = zeros(Int, length(first(ray_actions)))
    for (number, orbit) in enumerate(action_orbits(ray_actions, length(ray_orbit)))
        ray_orbit[orbit] .= number
    end
    return map(eachindex(faces)) do level
        orbits = action_orbits(actions[level], length(faces[level]))
        folds(face) = any(eachindex(ray_actions)) do g
            actions[level][g][face] == face && any(ray -> ray_actions[g][ray] != ray, faces[level][face])
        end
        (dimension=level - 1, orbits=length(orbits),
         folded=count(orbit -> folds(first(orbit)), orbits),
         pinched=count(orbit -> !allunique(ray_orbit[faces[level][first(orbit)]]), orbits))
    end
end

"The maximal chains of faces, entry `i` indexing a face of dimension `i - 1`."
function maximal_chains(faces, below)
    top = length(faces)
    chains = Vector{Int}[]
    chain = zeros(Int, top)
    function descend(level, face)
        chain[level] = face
        level == 1 && return push!(chains, copy(chain))
        foreach(smaller -> descend(level - 1, smaller), below[level][face])
    end
    foreach(face -> descend(top, face), eachindex(faces[top]))
    return chains
end

"The barycentric subdivision, as a simplicial complex on the faces numbered by dimension."
function barycentric_subdivision(faces, below)
    offsets = cumsum([0; length.(faces)[1:(end - 1)]])
    return simplicial_complex([chain .+ offsets for chain in maximal_chains(faces, below)])
end

"""
    barycentric_quotient(faces, below, actions)

The orbits of chains of faces, as a pure complex whose closed cells are
simplices. Facet `j` is the orbit of the maximal chain `facets[j]`, whose entry
`i` is a face of dimension `i - 1`; `face_ids[mask, j]` is the cell of its face
on the chain positions in the bitmask `mask`.
"""
function barycentric_quotient(faces, below, actions)
    top = length(faces)
    full = (1 << top) - 1
    cell_of = Dict{Tuple{Int, Vector{Int}}, Int}()
    function cell(mask, chain)
        positions = [i for i in 1:top if isodd(mask >> (i - 1))]
        # The group preserves dimensions, so the orbit is that of the faces at the same positions.
        key = minimum([actions[i][g][chain[i]] for i in positions] for g in eachindex(first(actions)))
        return get!(() -> length(cell_of) + 1, cell_of, (mask, key))
    end
    facets = Vector{Int}[]
    for chain in maximal_chains(faces, below)
        cell(full, chain) > length(facets) && push!(facets, chain)
    end
    face_ids = [cell(mask, chain) for mask in 1:full, chain in facets]
    # The faces of a chain have distinct dimensions, so a closed cell embeds.
    all(allunique, eachcol(face_ids)) || error("a closed cell of the quotient does not embed")
    cell_dimension = zeros(Int, length(cell_of))
    for ((mask, _), id) in cell_of
        cell_dimension[id] = count_ones(mask) - 1
    end
    return (dimension=top - 1, facets=facets, face_ids=face_ids, cell_dimension=cell_dimension)
end

f_vector_of(complex) = [count(==(k), complex.cell_dimension) for k in 0:complex.dimension]

struct SearchLimit <: Exception end

"""
    find_shelling(complex; max_states=10^6)

Search the facets of `complex` for a shelling: an order in which each facet
after the first meets the earlier ones in a pure subcomplex of codimension one
of its boundary. Returns `(status, order, states)` with `status` one of
`:shellable`, `:not_shellable` (the search was exhaustive) and `:unknown`
(more than `max_states` dead ends).
"""
function find_shelling(complex; max_states=10^6)
    d, face_ids = complex.dimension, complex.face_ids
    full = size(face_ids, 1)
    nfacets = size(face_ids, 2)
    containing = [Int[] for _ in complex.cell_dimension]
    for j in 1:nfacets, mask in 1:(full - 1)
        push!(containing[face_ids[mask, j]], j)
    end
    covering = zeros(Int, length(complex.cell_dimension))
    placed = falses(nfacets)
    order = Int[]
    failed = Set{BitVector}()
    present(mask, j) = mask == 0 ? !isempty(order) : covering[face_ids[mask, j]] > 0
    # Facet j can come next exactly when its faces missing from the earlier
    # facets form an interval [R, j]; then R spans the positions whose
    # opposite ridge is present. Returns the bitmask R, or -1.
    function restriction(j)
        r = 0
        for i in 0:d
            present(full ⊻ (1 << i), j) && (r |= 1 << i)
        end
        for mask in 0:(full - 1)
            present(mask, j) == (mask & r != r) || return -1
        end
        return r
    end
    # Only the facets containing a face that becomes present can change status.
    function place!(j, queue)
        for mask in 1:full
            id = face_ids[mask, j]
            (covering[id] += 1) == 1 && append!(queue, containing[id])
        end
        placed[j] = true
        push!(order, j)
    end
    function unplace!()
        j = pop!(order)
        foreach(mask -> covering[face_ids[mask, j]] -= 1, 1:full)
        placed[j] = false
    end
    function search(queue)
        # A facet meeting the earlier ones in all its ridges, or all but one,
        # adds at most one ridge, whose boundary is already present. That can
        # only complete the boundary of other facets, so if any shelling
        # continues the current one, one continues it with this facet next:
        # place these at once.
        forced = 0
        while !isempty(queue)
            j = pop!(queue)
            placed[j] && continue
            r = restriction(j)
            r >= 0 && count_ones(r) >= d || continue
            place!(j, queue)
            forced += 1
        end
        length(order) == nfacets && return true
        if !(placed in failed)
            candidates = Tuple{Int, Int}[]
            for j in 1:nfacets
                placed[j] && continue
                r = restriction(j)
                r >= 0 && push!(candidates, (-count_ones(r), j))
            end
            # Facets meeting the earlier ones in more ridges add fewer new faces: try them first.
            for (_, j) in sort!(candidates)
                queue = Int[]
                place!(j, queue)
                search(queue) && return true
                unplace!()
            end
            # Whether a set of facets begins a shelling does not depend on its order.
            push!(failed, copy(placed))
            length(failed) > max_states && throw(SearchLimit())
        end
        foreach(_ -> unplace!(), 1:forced)
        return false
    end
    # The recursion can be as deep as there are facets, so it gets a large stack.
    task = Task(() -> search(Int[]), 2^30)
    status = try
        schedule(task)
        fetch(task) ? :shellable : :not_shellable
    catch error
        error isa TaskFailedException && task.exception isa SearchLimit || rethrow()
        :unknown
    end
    return (status=status, order=status == :shellable ? copy(order) : Int[], states=length(failed))
end

"""
    shelling_h_vector(complex, order)

Check that `order` is a shelling from the definition, and return its h-vector:
entry `k + 1` counts the facets meeting the earlier ones in exactly `k` ridges.
"""
function shelling_h_vector(complex, order)
    d, face_ids = complex.dimension, complex.face_ids
    full = size(face_ids, 1)
    sort(order) == 1:size(face_ids, 2) || error("the order does not list every facet once")
    seen = Set{Int}()
    h = zeros(Int, d + 2)
    for (step, j) in enumerate(order)
        shared = [mask for mask in 1:(full - 1) if face_ids[mask, j] in seen]
        ridges = filter(mask -> count_ones(mask) == d, shared)
        step == 1 || !isempty(ridges) || error("facet $step meets the earlier ones in no ridge")
        all(mask -> any(ridge -> mask & ridge == mask, ridges), shared) ||
            error("facet $step meets the earlier ones in a face outside their common ridges")
        h[length(ridges) + 1] += 1
        union!(seen, face_ids[:, j])
    end
    return h
end

"""
    cell_subdivision(complex)

The barycentric subdivision of `complex`, i.e. the order complex of its cells.
The quotient is not a simplicial complex, as distinct cells can have the same
vertices, but this subdivision is one and is homeomorphic to it.
"""
function cell_subdivision(complex)
    face_ids = complex.face_ids
    # Within a facet, the flags of faces are the cumulative unions of its chain positions.
    flags = [face_ids[cumsum(1 .<< permutation), j]
             for j in axes(face_ids, 2) for permutation in all_permutations(0:complex.dimension)]
    return simplicial_complex(flags)
end

"Reduced integral homology computed by polymake, as the free rank and the torsion coefficients in each dimension."
function reduced_homology(complex::SimplicialComplex)
    return map(0:dim(complex)) do k
        divisors = elementary_divisors(homology(complex, k))
        (dimension=k, rank=count(iszero, divisors) - (k == 0),
         torsion=filter(divisor -> !iszero(divisor) && !isone(divisor), divisors))
    end
end

function describe_homology(groups)
    parts = String[]
    for group in groups
        summands = vcat(iszero(group.rank) ? String[] : ["Z^$(group.rank)"], ["Z/$t" for t in group.torsion])
        isempty(summands) || push!(parts, "H_$(group.dimension) = " * join(summands, " + "))
    end
    return isempty(parts) ? "zero" : join(parts, ", ")
end

function write_shelling(path, faces, complex, order)
    full = size(complex.face_ids, 1)
    open(path, "w") do io
        println(io, "step,shared_ridges," * join(("face_dimension_$k" for k in 0:complex.dimension), ','))
        seen = Set{Int}()
        for (step, j) in enumerate(order)
            shared = count(i -> complex.face_ids[full ⊻ (1 << i), j] in seen, 0:complex.dimension)
            chain = complex.facets[j]
            println(io, "$step,$(step == 1 ? 0 : shared)," *
                join(("\"" * join(faces[level][chain[level]], ';') * "\"" for level in eachindex(chain)), ','))
            union!(seen, complex.face_ids[:, j])
        end
    end
end

function main(arguments=ARGS)
    options = Dict("target" => "grassmannian", "fan" => "coarsest", "output" => "")
    for argument in arguments
        parsed = match(r"^--(target|fan|output)=(.+)$", argument)
        isnothing(parsed) && error("unknown argument $argument; the options are " *
                                   "--target=grassmannian|flag, --fan=coarsest|raw and --output=FILE")
        options[parsed[1]] = parsed[2]
    end
    haskey(TARGETS, options["target"]) || error("--target must be grassmannian or flag")
    haskey(FAN_FILES, options["fan"]) || error("--fan must be coarsest or raw")
    options["target"] == "flag" && options["fan"] == "raw" &&
        error("the raw flag fan has rays of its own, with no stored D12 action")
    directory = data_path("g2", options["target"])
    fan = load(joinpath(directory, FAN_FILES[options["fan"]]))
    is_pure(fan) || error("the fan is not pure")
    normalized_rays = Matrix{QQFieldElem}(load(joinpath(directory, "03_normalized_rays.mrdi")))
    check_ray_order(fan, normalized_rays)
    symmetry = TropicalGeometryG2.load_symmetry_artifacts(joinpath(directory, "05_symmetry"))
    length(symmetry.ray_actions) == 12 || error("expected the order-12 group D12")

    faces = fan_faces(fan)
    below = covered_faces(faces)
    actions = face_actions(faces, symmetry.ray_actions)
    println("$(TARGETS[options["target"]]) modulo lineality, $(options["fan"]) structure, f-vector ",
            Tuple(length.(faces)))
    defects = orbit_cell_defects(faces, actions, symmetry.ray_actions)
    println("D12 orbits of faces: ", Tuple(defect.orbits for defect in defects))
    for defect in defects[2:end]
        println("  dimension $(defect.dimension): $(defect.folded) orbits folded by their stabilizer, ",
                "$(defect.pinched) with two rays in one ray orbit")
    end
    any(defect -> defect.folded + defect.pinched > 0, defects) &&
        println("The face orbits are not the cells of a regular CW complex; we shell the barycentric quotient.")

    subdivision = barycentric_subdivision(faces, below)
    quotient = barycentric_quotient(faces, below, actions)
    println("Barycentric subdivision sd(Delta): f-vector ", Tuple(f_vector(subdivision)))
    println("  reduced integral homology: ", describe_homology(reduced_homology(subdivision)))
    println("Quotient sd(Delta)/D12: f-vector ", Tuple(f_vector_of(quotient)))
    homology = reduced_homology(cell_subdivision(quotient))
    println("  reduced integral homology: ", describe_homology(homology))

    result = @timed find_shelling(quotient)
    search = result.value
    println("Shelling search: $(search.status), $(search.states) dead-end states, ",
            round(result.time; digits=1), " s")
    search.status == :shellable || return search
    h = shelling_h_vector(quotient, search.order)
    println("  verified shelling, h-vector ", Tuple(h))
    # A shellable complex is a wedge of spheres, one per facet meeting the earlier ones in its whole boundary.
    last(h) == homology[end].rank && all(group -> iszero(group.rank) && isempty(group.torsion), homology[1:(end - 1)]) ||
        error("the shelling contradicts the homology")
    println("  so |Delta|/D12 is homotopy equivalent to a wedge of $(last(h)) spheres of dimension $(quotient.dimension)")
    if !isempty(options["output"])
        write_shelling(options["output"], faces, quotient, search.order)
        println("  shelling written to ", options["output"])
    end
    return search
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
