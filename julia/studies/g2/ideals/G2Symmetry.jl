module G2Symmetry

using TropicalGeometryG2

export G2_FORM_SUPPORT,
       G2_GROUP_CANDIDATES,
       G2_PERMUTATIONS,
       g2_leaf_signs,
       g2_coordinate_images,
       g2_flag_coordinate_images,
       G2_ROOT_DIRECTIONS,
       G2_GAUGE

# The support of the three-form Phi = e123 + e456 + e147 + e257 + e367.
const G2_FORM_SUPPORT = Set([(1, 2, 3), (4, 5, 6), (1, 4, 7), (2, 5, 7), (3, 6, 7)])

# The taxon permutations fixing coordinate seven. The symmetry stage keeps
# those that act on the fan.
const G2_GROUP_CANDIDATES = [vcat(permutation, [7]) for permutation in all_permutations(collect(1:6))]

# The order-12 split G2 taxon group fixes coordinate seven and preserves the
# five triples of the three-form. Use the full group, not only generators.
const G2_PERMUTATIONS = [vcat(permutation, [7])
    for permutation in all_permutations(1:6)
    if Set(Tuple(sort(vcat(permutation, [7])[collect(triple)]))
           for triple in G2_FORM_SUPPORT) == G2_FORM_SUPPORT]
length(G2_PERMUTATIONS) == 12 || error("expected the order-12 split G2 group")
TropicalGeometryG2.validate_tropical_permutation_group(G2_PERMUTATIONS)

function _permutation_orientation(values)
    inversions = count(i < j && values[i] > values[j]
                       for i in eachindex(values) for j in eachindex(values))
    return isodd(inversions) ? -1 : 1
end

"""
Diagonal signs for a taxon lift preserving the split three-form up to sign:
signs s such that e_i -> s_i e_pi(i) maps Phi to a multiple of itself.
"""
function g2_leaf_signs(permutation)
    length(permutation) == 7 && sort(collect(permutation)) == collect(1:7) ||
        throw(ArgumentError("the split G2 action requires a permutation of 1:7"))
    for mask in 0:(2^7 - 1)
        signs = [iszero(mask & (1 << (index - 1))) ? 1 : -1 for index in 1:7]
        coefficients = Dict{NTuple{3, Int}, Int}()
        for triple in G2_FORM_SUPPORT
            image = [permutation[index] for index in triple]
            key = Tuple(sort(image))
            coefficients[key] = _permutation_orientation(image) *
                                prod(signs[index] for index in triple)
        end
        Set(keys(coefficients)) == G2_FORM_SUPPORT || continue
        length(unique(values(coefficients))) == 1 && return signs
    end
    error("taxon permutation has no signed lift preserving the split G2 three-form")
end

"Signed exterior-square images of the 21 Pluecker coordinates."
g2_coordinate_images(permutation, variables) = TropicalGeometryG2.pluecker_coordinate_images(
    permutation,
    variables;
    leaf_signs=g2_leaf_signs(permutation),
)

"Signed point and exterior-square images in the full 28-coordinate ring."
function g2_flag_coordinate_images(permutation, variables)
    length(variables) == 28 || throw(ArgumentError("expected 28 flag coordinates"))
    signs = g2_leaf_signs(permutation)
    return vcat(
        [signs[index] * variables[permutation[index]] for index in 1:7],
        TropicalGeometryG2.pluecker_coordinate_images(
            permutation, variables[8:28]; leaf_signs=signs,
        ),
    )
end

# One-parameter subgroups of the maximal torus of G2 on e1, ..., e7: exponents
# (a, b, c, -a, -b, -c, 0) with a + b + c = 0. Stage 3 normalizes a ray so
# that its pendant lengths are orthogonal to them.
const G2_ROOT_DIRECTIONS = [
     1  -1   0  -1   1   0  0
    -1   0   1   1   0  -1  0
     0   1  -1   0  -1   1  0
]
const G2_GAUGE = TreeGauge(7, G2_ROOT_DIRECTIONS)

end # module
