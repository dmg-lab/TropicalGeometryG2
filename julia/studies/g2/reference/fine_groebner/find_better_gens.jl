using Oscar
import TropicalGeometryG2
# The exact trees of the package; G2Symmetry holds the G2 roots and gauge.
using TropicalGeometryG2: TreeMetric, star_pluecker, primitive_integer, tree_from_min_pluecker,
    tree_distance, pendant_lengths, tree_edge_lengths, realizes, leaves_on_side,
    internal_splits, inner_product, normalize_ray
isdefined(@__MODULE__, :G2Symmetry) ||
    include(joinpath(@__DIR__, "..", "..", "ideals", "G2Symmetry.jl"))
using .G2Symmetry: G2_ROOT_DIRECTIONS, G2_GAUGE

const DEFAULT_FAN_PATH =
    TropicalGeometryG2.data_path("g2", "grassmannian", "01_raw_fan.mrdi")
const PLUECKER_PAIRS = TropicalGeometryG2.pluecker_pairs(7)

# A row s represents a star tree with leaf-edge lengths s_1, ..., s_7.  Its
# image in Pluecker coordinates is (s_i + s_j)_{1 <= i < j <= 7}.  The last
# three rows are the Coxeter roots G2_ROOT_DIRECTIONS.
const LINEALITY_SPACE_STAR = vcat(ones(Int, 1, 7), G2_ROOT_DIRECTIONS)
const COXETER_STAR_ROOTS = LINEALITY_SPACE_STAR[2:4, :]

"""The four supplied star-tree rows, expressed in Pluecker coordinates."""
function star_lineality_matrix()
    M = Matrix{QQFieldElem}(undef, size(LINEALITY_SPACE_STAR, 1), length(PLUECKER_PAIRS))
    for i in axes(LINEALITY_SPACE_STAR, 1)
        M[i, :] = star_pluecker(view(LINEALITY_SPACE_STAR, i, :))
    end
    return M
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

`normalize_ray(generator, G2_GAUGE)`: first use the last three (Coxeter-root)
star generators to choose the least-squares-balanced representative of a
min-convention Pluecker ray. Then clear the smallest common denominator of all
tree edges, use the diagonal star direction to translate the shortest pendant
edge to zero, and divide by the gcd of all shifted tree-edge lengths. The final
tree is primitive integral, with nonnegative pendant lengths and positive
internal lengths. The record adds the coefficients of `coxeter_balance`, the
Coxeter-balanced generator and the Oscar tree.
"""
function balance_outer_edges(generator::AbstractVector)
    record = normalize_ray(generator, G2_GAUGE)
    original_outer = pendant_lengths(tree_from_min_pluecker(record.original))
    balance = coxeter_balance(original_outer)
    @assert balance.shift == record.lineality_leaf_shift
    return (
        generator=record.original,
        root_coefficients=balance.coefficients,
        lineality_leaf_shift=balance.shift,
        coxeter_balanced_generator=record.original + star_pluecker(balance.shift),
        coxeter_balanced_pendant_lengths=balance.balanced,
        integer_scale=record.integer_scale,
        diagonal_leaf_shift=record.diagonal_leaf_shift,
        diagonal_pluecker_shift=2 * record.diagonal_leaf_shift,
        primitive_divisor=record.primitive_divisor,
        balanced_generator=record.balanced_generator,
        original_pendant_lengths=original_outer,
        pendant_lengths=record.pendant_lengths,
        tree=record.tree,
        oscar_tree=oscar_phylogenetic_tree(record.tree),
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
