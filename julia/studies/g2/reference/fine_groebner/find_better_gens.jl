using Oscar
import TropicalCensus

const DEFAULT_FAN_PATH =
    TropicalCensus.data_path("g2", "grassmannian", "01_raw_fan.mrdi")
const PLUECKER_PAIRS = [(i, j) for i in 1:7 for j in (i + 1):7]

# A row s represents a star tree with leaf-edge lengths s_1, ..., s_7.  Its
# image in Pluecker coordinates is (s_i + s_j)_{1 <= i < j <= 7}.
const LINEALITY_SPACE_STAR = [
     1  1  1  1  1  1  1
     1 -1  0 -1  1  0  0
    -1  0  1  1  0 -1  0
     0  1 -1  0 -1  1  0
]
const COXETER_STAR_ROOTS = LINEALITY_SPACE_STAR[2:4, :]

"""
    TreeMetric

An exact, unrooted phylogenetic tree.  Taxa are vertices `1:ntaxa`; larger
vertex labels are internal nodes.  Edge lengths are rational and pendant-edge
lengths are allowed to be negative, as is natural before quotienting tropical
Grassmannian tree space by its full star-tree lineality space.
"""
struct TreeMetric
    ntaxa::Int
    adjacency::Dict{Int, Vector{Tuple{Int, QQFieldElem}}}
end

"""Map leaf-edge lengths to the corresponding 21 Pluecker coordinates."""
function star_pluecker(s::AbstractVector)
    length(s) == 7 || throw(ArgumentError("a seven-taxon star must have 7 edge lengths"))
    return QQ.([s[i] + s[j] for (i, j) in PLUECKER_PAIRS])
end

"""The four supplied star-tree rows, expressed in Pluecker coordinates."""
function star_lineality_matrix()
    M = Matrix{QQFieldElem}(undef, size(LINEALITY_SPACE_STAR, 1), length(PLUECKER_PAIRS))
    for i in axes(LINEALITY_SPACE_STAR, 1)
        M[i, :] = star_pluecker(view(LINEALITY_SPACE_STAR, i, :))
    end
    return M
end

"""
    primitive_integer(v)

Clear denominators and divide by the gcd of the resulting entries.  No sign
change is made, since a ray has a distinguished positive direction.
"""
function primitive_integer(v::AbstractVector)
    q = QQ.(collect(v))
    isempty(q) && throw(ArgumentError("cannot normalize an empty vector"))
    all(iszero, q) && throw(ArgumentError("the zero vector does not generate a ray"))

    denominator_lcm = foldl(lcm, denominator.(q); init=ZZ(1))
    integral = denominator_lcm .* q
    entry_gcd = foldl(gcd, (abs(numerator(x)) for x in integral if !iszero(x)); init=ZZ(0))
    return integral ./ entry_gcd
end

function number_of_taxa(ncoordinates::Integer)
    discriminant = 1 + 8 * ncoordinates
    root = isqrt(discriminant)
    n = (1 + root) ÷ 2
    root^2 == discriminant && binomial(n, 2) == ncoordinates ||
        throw(ArgumentError("$ncoordinates is not binomial(n, 2) for any integer n"))
    return n
end

function pair_dictionary(values::AbstractVector)
    n = number_of_taxa(length(values))
    D = Dict{Tuple{Int, Int}, QQFieldElem}()
    k = 0
    for i in 1:(n - 1), j in (i + 1):n
        k += 1
        D[i, j] = D[j, i] = QQ(values[k])
    end
    return D, n
end

"""
    satisfies_tropical_pluecker(w) -> Bool

Check the three-term tropical Pluecker relations for `Gr(2,n)` in the min
convention.  Equivalently, `-w` obeys the four-point condition for additive
tree dissimilarities.
"""
function satisfies_tropical_pluecker(w::AbstractVector)
    W, n = pair_dictionary(w)
    for i in 1:(n - 3), j in (i + 1):(n - 2), k in (j + 1):(n - 1), l in (k + 1):n
        terms = (W[i, j] + W[k, l], W[i, k] + W[j, l], W[i, l] + W[j, k])
        minimum_term = min(terms...)
        count(==(minimum_term), terms) >= 2 || return false
    end
    return true
end

"""
    neighbor_joining(distances, ntaxa) -> TreeMetric

Exact neighbor joining for an additive dissimilarity.  For a tropical
Pluecker vector this reconstructs the tree without numerical tolerances.
Unresolved vertices may appear as zero-length internal edges in the returned
binary refinement.
"""
function neighbor_joining(
    distances::Dict{Tuple{Int, Int}, QQFieldElem},
    ntaxa::Int,
)
    D = copy(distances)
    edges = Tuple{Int, Int, QQFieldElem}[]
    active = collect(1:ntaxa)
    next_label = ntaxa + 1

    while length(active) > 2
        m = length(active)
        row_sums = Dict(
            i => sum((D[i, k] for k in active if k != i); init=QQ(0))
            for i in active
        )

        best_i, best_j = active[1], active[2]
        best_q = (m - 2) * D[best_i, best_j] - row_sums[best_i] - row_sums[best_j]
        for a in 1:(m - 1), b in (a + 1):m
            i, j = active[a], active[b]
            q = (m - 2) * D[i, j] - row_sums[i] - row_sums[j]
            if q < best_q
                best_i, best_j, best_q = i, j, q
            end
        end

        u = next_label
        next_label += 1
        i, j = best_i, best_j
        length_iu = D[i, j] / 2 + (row_sums[i] - row_sums[j]) / (2 * (m - 2))
        length_ju = D[i, j] - length_iu
        push!(edges, (i, u, length_iu), (j, u, length_ju))

        for k in active
            (k == i || k == j) && continue
            D[u, k] = D[k, u] = (D[i, k] + D[j, k] - D[i, j]) / 2
        end
        active = push!([k for k in active if k != i && k != j], u)
    end

    i, j = active
    push!(edges, (i, j, D[i, j]))

    adjacency = Dict(
        v => Tuple{Int, QQFieldElem}[]
        for v in 1:(next_label - 1)
    )
    for (u, v, length) in edges
        push!(adjacency[u], (v, length))
        push!(adjacency[v], (u, length))
    end
    return TreeMetric(ntaxa, adjacency)
end

"""
    contract_zero_internal_edges(tree) -> TreeMetric

Remove the zero-length edges between internal vertices that neighbor joining
uses to give a binary refinement of a genuinely unresolved tree.  Pendant
edges, including zero-length pendant edges, are retained so that every taxon
remains a leaf.
"""
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

    function merge!(u, v)
        root_u, root_v = representative(u), representative(v)
        root_u == root_v && return
        parent[max(root_u, root_v)] = min(root_u, root_v)
    end

    for u in 1:nvertices
        for (v, edge_length) in tree.adjacency[u]
            u < v || continue
            u > tree.ntaxa && v > tree.ntaxa && iszero(edge_length) && merge!(u, v)
        end
    end

    internal_roots = sort!(unique(representative(v) for v in (tree.ntaxa + 1):nvertices))
    internal_label = Dict(root => tree.ntaxa + i for (i, root) in enumerate(internal_roots))
    new_label(v) = v <= tree.ntaxa ? v : internal_label[representative(v)]
    new_nvertices = tree.ntaxa + length(internal_roots)
    adjacency = Dict(v => Tuple{Int, QQFieldElem}[] for v in 1:new_nvertices)

    for u in 1:nvertices
        for (v, edge_length) in tree.adjacency[u]
            u < v || continue
            new_u, new_v = new_label(u), new_label(v)
            new_u == new_v && continue
            push!(adjacency[new_u], (new_v, edge_length))
            push!(adjacency[new_v], (new_u, edge_length))
        end
    end
    return TreeMetric(tree.ntaxa, adjacency)
end

"""
    tree_from_min_pluecker(w) -> TreeMetric

Reconstruct the phylogenetic tree encoded by a min-convention tropical
Pluecker vector `w`.  Tree dissimilarities are `d_ij = -w_ij`.
"""
function tree_from_min_pluecker(w::AbstractVector)
    satisfies_tropical_pluecker(w) ||
        throw(ArgumentError("the vector does not satisfy the tropical Pluecker relations"))
    distances, ntaxa = pair_dictionary(-QQ.(collect(w)))
    return contract_zero_internal_edges(neighbor_joining(distances, ntaxa))
end

"""Return the distance between two taxa in a reconstructed tree."""
function tree_distance(tree::TreeMetric, source::Int, target::Int)
    1 <= source <= tree.ntaxa && 1 <= target <= tree.ntaxa ||
        throw(ArgumentError("taxa must lie in 1:$(tree.ntaxa)"))
    source == target && return QQ(0)

    stack = [(source, 0, QQ(0))]
    while !isempty(stack)
        vertex, parent, distance = pop!(stack)
        for (neighbor, length) in tree.adjacency[vertex]
            neighbor == parent && continue
            neighbor == target && return distance + length
            push!(stack, (neighbor, vertex, distance + length))
        end
    end
    error("the reconstructed graph is disconnected")
end

"""Return the seven pendant-edge lengths, ordered by taxon."""
function pendant_lengths(tree::TreeMetric)
    result = Vector{QQFieldElem}(undef, tree.ntaxa)
    for taxon in 1:tree.ntaxa
        length(tree.adjacency[taxon]) == 1 ||
            error("taxon $taxon is not a leaf of the reconstructed tree")
        result[taxon] = only(tree.adjacency[taxon])[2]
    end
    return result
end

"""Return every edge length once."""
function tree_edge_lengths(tree::TreeMetric)
    return [
        edge_length
        for u in sort!(collect(keys(tree.adjacency)))
        for (v, edge_length) in tree.adjacency[u]
        if u < v
    ]
end

"""Verify that the reconstructed tree has precisely the requested distances."""
function realizes(tree::TreeMetric, distances::AbstractVector)
    D, n = pair_dictionary(distances)
    n == tree.ntaxa || return false
    return all(tree_distance(tree, i, j) == D[i, j] for i in 1:(n - 1) for j in (i + 1):n)
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

"""
    internal_splits(tree)

Return the positive-length internal edges as canonical leaf splits.  These
splits and their lengths are unchanged by every star-tree lineality shift;
only pendant-edge lengths change.
"""
function internal_splits(tree::TreeMetric)
    result = NamedTuple{(:side, :complement, :length), Tuple{Vector{Int}, Vector{Int}, QQFieldElem}}[]
    all_taxa = collect(1:tree.ntaxa)
    for u in sort!(collect(keys(tree.adjacency)))
        for (v, edge_length) in tree.adjacency[u]
            u < v || continue
            u > tree.ntaxa && v > tree.ntaxa && edge_length > 0 || continue
            side = leaves_on_side(tree, u, v)
            complement = setdiff(all_taxa, side)
            if length(side) > length(complement) ||
                    (length(side) == length(complement) && join(side, ",") > join(complement, ","))
                side, complement = complement, side
            end
            push!(result, (side=side, complement=complement, length=edge_length))
        end
    end
    sort!(result; by=split -> (length(split.side), join(split.side, ",")))
    return result
end

"""Convert the exact reconstruction to Oscar's `PhylogeneticTree` type."""
function oscar_phylogenetic_tree(tree::TreeMetric)
    vertices = sort!(collect(keys(tree.adjacency)))
    vertices == collect(1:length(vertices)) || error("tree vertices must be consecutively labelled")
    root_index = findfirst(>(tree.ntaxa), vertices)
    root_index === nothing && error("the tree has no internal vertex at which to root")
    root = vertices[root_index]

    graph = Graph{Directed}(length(vertices))
    edge_lengths = Dict{NTuple{2, Int}, QQFieldElem}()
    visited = falses(length(vertices))
    visited[root] = true
    queue = [root]
    while !isempty(queue)
        vertex = popfirst!(queue)
        for (neighbor, length) in tree.adjacency[vertex]
            visited[neighbor] && continue
            visited[neighbor] = true
            add_edge!(graph, vertex, neighbor)
            edge_lengths[vertex, neighbor] = length
            push!(queue, neighbor)
        end
    end
    label!(graph, edge_lengths, nothing; name=:distance)
    label!(graph, nothing, Dict(i => string(i) for i in 1:tree.ntaxa); name=:leaves)
    return Oscar.phylogenetic_tree(graph)
end

function rational_row(M, i)
    return QQ.(collect(view(M, i, :)))
end

function inner_product(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || throw(ArgumentError("vectors have different lengths"))
    return sum((QQ(x[i]) * QQ(y[i]) for i in eachindex(x)); init=QQ(0))
end

"""
    coxeter_balance(outer_lengths)

Orthogonally project the pendant-edge vector onto the span of the last three
star-tree lineality generators.  If the roots are `r_1,r_2,r_3`, their Gram
matrix has diagonal 4 and off-diagonal -2, so they form a tight frame with
frame constant 6.  Hence the projection and its canonical coefficients are

    shift = sum_i <outer_lengths,r_i>/6 * r_i.

The three coefficients sum to zero.  Subtracting `shift` gives the unique
representative perpendicular to the Coxeter-root plane and therefore the one
with minimum squared pendant-edge length (equivalently minimum variance,
since every root has coordinate sum zero).
"""
function coxeter_balance(outer_lengths::AbstractVector)
    length(outer_lengths) == 7 || throw(ArgumentError("expected seven pendant-edge lengths"))
    outer = QQ.(collect(outer_lengths))
    roots = [rational_row(COXETER_STAR_ROOTS, i) for i in axes(COXETER_STAR_ROOTS, 1)]
    coefficients = [inner_product(outer, root) / 6 for root in roots]
    @assert sum(coefficients; init=QQ(0)) == 0

    shift = zeros(QQFieldElem, 7)
    for (coefficient, root) in zip(coefficients, roots)
        shift .+= coefficient .* root
    end
    balanced = outer - shift
    @assert all(iszero(inner_product(balanced, root)) for root in roots)
    @assert sum(balanced; init=QQ(0)) == sum(outer; init=QQ(0))
    return (coefficients=coefficients, shift=shift, balanced=balanced)
end

"""
    balance_outer_edges(generator)

First use the last three (Coxeter-root) star generators to choose the
least-squares-balanced representative of a min-convention Pluecker ray. Then
clear the smallest common denominator of all tree edges, use the diagonal
star direction to translate the shortest pendant edge to zero, and divide by
the gcd of all shifted tree-edge lengths. The final tree is primitive integral,
with nonnegative pendant lengths and positive internal lengths.
"""
function balance_outer_edges(generator::AbstractVector)
    original_generator = QQ.(collect(generator))
    original_tree = tree_from_min_pluecker(original_generator)
    original_outer = pendant_lengths(original_tree)
    balance = coxeter_balance(original_outer)

    # Adding star_pluecker(s) to min-convention coordinates subtracts s_i
    # from pendant edge i of the associated distance tree.
    coxeter_balanced_generator = original_generator + star_pluecker(balance.shift)
    coxeter_balanced_tree = tree_from_min_pluecker(coxeter_balanced_generator)
    coxeter_balanced_outer = pendant_lengths(coxeter_balanced_tree)

    @assert coxeter_balanced_outer == balance.balanced
    @assert internal_splits(coxeter_balanced_tree) == internal_splits(original_tree)
    @assert realizes(coxeter_balanced_tree, -coxeter_balanced_generator)
    @assert inner_product(coxeter_balanced_outer, coxeter_balanced_outer) <=
        inner_product(original_outer, original_outer)

    # Scaling preserves the ray.  Use the smallest scale that makes every
    # edge integral, not merely every Pluecker coordinate integral.
    integer_scale = foldl(
        lcm,
        denominator.(tree_edge_lengths(coxeter_balanced_tree));
        init=ZZ(1),
    )
    integral_generator = integer_scale .* coxeter_balanced_generator
    integral_tree = tree_from_min_pluecker(integral_generator)
    integral_outer = pendant_lengths(integral_tree)

    # If a is added to all seven leaf-star lengths, then star_pluecker(a*1)
    # is the diagonal Pluecker shift 2a*(1,...,1), and the associated pendant
    # lengths change by -a. Taking a=min(l) makes the new minimum exactly zero.
    # This is still a translation by lineality, not a scaling operation.
    diagonal_leaf_shift = minimum(integral_outer)
    diagonal_shift = fill(diagonal_leaf_shift, 7)
    shifted_generator = integral_generator + star_pluecker(diagonal_shift)
    shifted_tree = tree_from_min_pluecker(shifted_generator)
    shifted_outer = pendant_lengths(shifted_tree)

    @assert minimum(shifted_outer) == 0
    @assert all(edge_length >= 0 && denominator(edge_length) == 1
                for edge_length in tree_edge_lengths(shifted_tree))

    # Only after all lineality translations have been made do we choose the
    # primitive lattice generator. Division by this gcd is a positive scaling,
    # so it preserves the ray in the quotient fan.
    primitive_divisor = foldl(
        gcd,
        abs.(numerator.(tree_edge_lengths(shifted_tree)));
        init=ZZ(0),
    )
    @assert primitive_divisor > 0
    balanced_generator = shifted_generator ./ primitive_divisor
    balanced_tree = tree_from_min_pluecker(balanced_generator)
    balanced_outer = pendant_lengths(balanced_tree)

    roots = [rational_row(COXETER_STAR_ROOTS, i) for i in axes(COXETER_STAR_ROOTS, 1)]
    @assert minimum(balanced_outer) == 0
    @assert all(edge_length >= 0 && denominator(edge_length) == 1
                for edge_length in tree_edge_lengths(balanced_tree))
    @assert all(split.length > 0 for split in internal_splits(balanced_tree))
    @assert all(denominator(entry) == 1 for entry in balanced_generator)
    @assert all(iszero(inner_product(balanced_outer, root)) for root in roots)
    @assert realizes(balanced_tree, -balanced_generator)
    @assert foldl(
        gcd,
        abs.(numerator.(tree_edge_lengths(balanced_tree)));
        init=ZZ(0),
    ) == 1

    # Explicit quotient-space certificate:
    #   g*w_new = m*w_old + star_pluecker(m*root_shift + diagonal_shift).
    # The difference between the two scaled generators is therefore in the
    # stated star-tree lineality space.
    @assert primitive_divisor .* balanced_generator ==
        integer_scale .* original_generator +
        star_pluecker(integer_scale .* balance.shift + diagonal_shift)

    original_splits = internal_splits(original_tree)
    final_splits = internal_splits(balanced_tree)
    @assert [(s.side, s.complement) for s in final_splits] ==
        [(s.side, s.complement) for s in original_splits]
    quotient_scale = QQ(integer_scale) / QQ(primitive_divisor)
    @assert [s.length for s in final_splits] ==
        quotient_scale .* [s.length for s in original_splits]

    return (
        generator=original_generator,
        root_coefficients=balance.coefficients,
        lineality_leaf_shift=balance.shift,
        coxeter_balanced_generator=coxeter_balanced_generator,
        coxeter_balanced_pendant_lengths=coxeter_balanced_outer,
        integer_scale=integer_scale,
        diagonal_leaf_shift=diagonal_leaf_shift,
        diagonal_pluecker_shift=2 * diagonal_leaf_shift,
        primitive_divisor=primitive_divisor,
        balanced_generator=balanced_generator,
        original_pendant_lengths=original_outer,
        pendant_lengths=balanced_outer,
        tree=balanced_tree,
        oscar_tree=oscar_phylogenetic_tree(balanced_tree),
    )
end

"""
    ray_tree_data(fan)

Primitive-normalize every ray representative of `fan` and reconstruct its
exact seven-taxon tree.  Each result has fields `generator`, `tree`, and
`oscar_tree`.
"""
function ray_tree_data(fan)
    rays = rays_modulo_lineality(fan).rays_modulo_lineality
    return [
        begin
            generator = primitive_integer(collect(ray))
            tree = tree_from_min_pluecker(generator)
            @assert realizes(tree, -generator)
            @assert all(
                length >= 0
                for (u, neighbors) in tree.adjacency
                for (v, length) in neighbors
                if u > tree.ntaxa && v > tree.ntaxa
            ) "ray $i reconstructed with a negative internal edge"
            (generator=generator, tree=tree, oscar_tree=oscar_phylogenetic_tree(tree))
        end
        for (i, ray) in enumerate(rays)
    ]
end

"""Return just the 51 Oscar tree objects, one for each ray of the fan."""
phylogenetic_trees(fan) = [record.oscar_tree for record in ray_tree_data(fan)]

"""Balance all fan rays using the Coxeter-root lineality directions."""
function balanced_ray_tree_data(fan)
    rays = rays_modulo_lineality(fan).rays_modulo_lineality
    return [balance_outer_edges(primitive_integer(collect(ray))) for ray in rays]
end

"""Return the 51 balanced Oscar tree objects, one for each ray of the fan."""
balanced_phylogenetic_trees(fan) =
    [record.oscar_tree for record in balanced_ray_tree_data(fan)]

function matrix_from_rows(rows)
    M = Matrix{QQFieldElem}(undef, length(rows), length(first(rows)))
    for (i, row) in enumerate(rows)
        M[i, :] = QQ.(collect(row))
    end
    return M
end

"""Check that the stated star trees span exactly the fan's lineality space."""
function verify_star_lineality(fan)
    fan_rows = rays_modulo_lineality(fan).lineality_basis
    A = matrix(QQ, star_lineality_matrix())
    B = matrix(QQ, matrix_from_rows(fan_rows))
    expected_rank = rank(A)
    return expected_rank == 3 && rank(B) == expected_rank && rank(vcat(A, B)) == expected_rank
end

function format_split(split)
    left = join(split.side, ",")
    right = join(split.complement, ",")
    return "{$left}|{$right}:$(split.length)"
end

format_vector(v) = "[" * join(string.(v), ", ") * "]"

function main(io::IO=stdout; fan_path::AbstractString=DEFAULT_FAN_PATH)
    fan = load(fan_path)
    verify_star_lineality(fan) || error("the supplied star trees do not span the fan lineality space")
    records = balanced_ray_tree_data(fan)

    println(io, "# $(length(records)) primitive integral Coxeter-balanced ray trees from $(basename(fan_path))")
    println(io, "# convention: d_ij = -w_ij; input ray representatives are primitive integral vectors")
    for (i, record) in enumerate(records)
        splits = internal_splits(record.tree)
        split_string = isempty(splits) ? "star" : join(format_split.(splits), "; ")
        println(io, "ray $i")
        println(io, "  root coefficients: ", format_vector(record.root_coefficients))
        println(io, "  pendant edges before: ", format_vector(record.original_pendant_lengths))
        println(io, "  pendant edges Coxeter-balanced: ",
                format_vector(record.coxeter_balanced_pendant_lengths))
        println(io, "  integer scale: ", record.integer_scale)
        println(io, "  diagonal leaf shift: ", record.diagonal_leaf_shift,
                " (Pluecker shift ", record.diagonal_pluecker_shift, ")")
        println(io, "  primitive divisor: ", record.primitive_divisor)
        println(io, "  primitive nonnegative pendant edges: ", format_vector(record.pendant_lengths))
        println(io, "  internal splits: $split_string")
        println(io, "  balanced generator: ", format_vector(record.balanced_generator))
        println(io, "  Newick: ", newick(record.oscar_tree))
    end
    return records
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
