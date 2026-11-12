module G2GrassmannianSlongIdeal

using Oscar

include(joinpath(@__DIR__, "G2GrassmannianIdeal.jl"))
using .G2GrassmannianIdeal

export build_g2_grassmannian_slong_ideal, SLONG_POSITIONS, slong_substitution,
    SLONG_WEYL_GENERATORS, slong_coordinate_images, slong_signed_permutation, g2_leaf_signs,
    build_g2_grassmannian_short_long_ideal, SHORT_LONG_WEYL_GENERATORS, short_long_signed_permutation

const PAIRS = [(i, j) for i in 1:6 for j in (i + 1):7]
const PAIR_INDEX = Dict(pair => k for (k, pair) in enumerate(PAIRS))

# The support of the three-form Phi = e123 + e456 + e147 + e257 + e367.
const FORM_SUPPORT = Set([(1, 2, 3), (4, 5, 6), (1, 4, 7), (2, 5, 7), (3, 6, 7)])

orientation(values) = isodd(count(values[i] > values[j] for i in eachindex(values) for j in eachindex(values) if i < j)) ? -1 : 1

"Signs s such that e_i -> s_i e_pi(i) maps Phi to a multiple of itself, as in run_census.jl."
function g2_leaf_signs(permutation)
    for mask in 0:(2^7 - 1)
        signs = [iszero(mask & (1 << (index - 1))) ? 1 : -1 for index in 1:7]
        coefficients = Dict{NTuple{3, Int}, Int}()
        for triple in FORM_SUPPORT
            image = [permutation[index] for index in triple]
            coefficients[Tuple(sort(image))] = orientation(image) * prod(signs[index] for index in triple)
        end
        Set(keys(coefficients)) == FORM_SUPPORT || continue
        length(unique(values(coefficients))) == 1 && return signs
    end
    error("taxon permutation has no signed lift preserving the G2 three-form")
end

"""
Positions of the zero-weight coordinates u = p25 - p36, v = p36 - p14 and
w = p14 - p25 of the long tripod S_long. Each takes the slot of the Pluecker
coordinate on the short root line that it is orthogonal to: u that of p14, v
that of p25, w that of p36, so the Weyl group permutes the 21 positions as it
permutes the Pluecker coordinates.
"""
const SLONG_POSITIONS = (
    findfirst(==((1, 4)), PAIRS), findfirst(==((2, 5)), PAIRS), findfirst(==((3, 6)), PAIRS),
)

"""
The images of the 21 Pluecker coordinates in the S_long ring: the identity on
the 18 others, and on the plane p14 + p25 + p36 = 0, where (u, v, w) determines
(p14, p25, p36) since char 0 is not 3, p14 = (w - v)/3, p25 = (u - w)/3,
p36 = (v - u)/3.
"""
function slong_substitution(variables)
    u, v, w = (variables[k] for k in SLONG_POSITIONS)
    third = inv(coefficient_ring(parent(u))(3))
    images = collect(variables)
    images[SLONG_POSITIONS[1]] = third * (w - v)
    images[SLONG_POSITIONS[2]] = third * (u - w)
    images[SLONG_POSITIONS[3]] = third * (v - u)
    return images
end

"""
The signed lift of the taxon permutation `permutation` (fixing 7) on the S_long
coordinates, gfan-style: variable k goes to signs[k] times variable p[k]. On
the 18 Pluecker coordinates it is the lift e_i -> s_i e_pi(i) of
`g2_leaf_signs`. The lift permutes p14, p25, p36 up to one common sign c,
since it maps p14 + p25 + p36 to a multiple of itself, so it maps the
difference p_a - p_b to c (p_pi(a) - p_pi(b)), which is plus or minus one of
u, v, w.
"""
function slong_signed_permutation(permutation)
    target, signs = pluecker_signed_permutation(permutation)
    images = [difference_image(k, target, signs) for k in SLONG_POSITIONS]
    for (k, (position, sign)) in zip(SLONG_POSITIONS, images)
        target[k], signs[k] = position, sign
    end
    return target => signs
end

"The signed lift on the 21 Pluecker coordinates, gfan-style."
function pluecker_signed_permutation(permutation)
    signs_leaf = g2_leaf_signs(permutation)
    target = Vector{Int}(undef, length(PAIRS))
    signs = Vector{Int}(undef, length(PAIRS))
    for (k, (i, j)) in enumerate(PAIRS)
        a, b = permutation[i], permutation[j]
        target[k] = PAIR_INDEX[(min(a, b), max(a, b))]
        signs[k] = (a < b ? 1 : -1) * signs_leaf[i] * signs_leaf[j]
    end
    return target, signs
end

"""
The image (slot, sign) of the zero-weight form in the slot `k` of p14, p25 or
p36 (u, v or w) under the Pluecker lift (`target`, `signs`).
"""
function difference_image(k, target, signs)
    p14, p25, p36 = SLONG_POSITIONS
    # The difference of Pluecker coordinates each of u, v, w is, by its slot.
    difference = Dict(p14 => (p25, p36), p25 => (p36, p14), p36 => (p14, p25))
    slot = Dict(pair => position for (position, pair) in difference)
    a, b = difference[k]
    signs[a] == signs[b] || error("the lift does not permute p14, p25, p36 with a common sign")
    image = (target[a], target[b])
    return haskey(slot, image) ? (slot[image], signs[a]) : (slot[reverse(image)], -signs[a])
end

"""
The signed lift on the 24 coordinates of S_short and S_long together: the 21
Pluecker coordinates, then u, v, w as coordinates 22, 23, 24.
"""
function short_long_signed_permutation(permutation)
    target, signs = pluecker_signed_permutation(permutation)
    extra = [difference_image(k, target, signs) for k in SLONG_POSITIONS]
    shift = Dict(k => length(PAIRS) + index for (index, k) in enumerate(SLONG_POSITIONS))
    return vcat(target, [shift[position] for (position, _) in extra]) =>
        vcat(signs, [sign for (_, sign) in extra])
end

"""
The rotation rho = (1 5 3 4 2 6) and the reflection tau = (2 3)(5 6) of
Section 8, which generate the Weyl group, as signed variable permutations of
the S_long ring (`symmetry_action=:variables` of `run_tropicalization`).
"""
const WEYL_GENERATOR_TAXA = [[5, 6, 4, 2, 3, 1, 7], [1, 3, 2, 4, 6, 5, 7]]
const SLONG_WEYL_GENERATORS = slong_signed_permutation.(WEYL_GENERATOR_TAXA)
const SHORT_LONG_WEYL_GENERATORS = short_long_signed_permutation.(WEYL_GENERATOR_TAXA)

"""
The ideal of G2/P2 in the 24 coordinates of S_short and S_long together: the
Pluecker ideal in p12, ..., p67 and u - (p25 - p36), v - (p36 - p14),
w - (p14 - p25) in u, v, w. Both tropicalizations are coordinate projections
of its tropicalization.
"""
function build_g2_grassmannian_short_long_ideal(; coefficient_field=QQ)
    pluecker = build_g2_grassmannian_ideal(; coefficient_field=coefficient_field)
    names = vcat(["p$i$j" for (i, j) in PAIRS], ["u", "v", "w"])
    ring, variables = polynomial_ring(coefficient_field, names)
    embed = hom(pluecker.ring, ring, variables[1:length(PAIRS)])
    p14, p25, p36 = (variables[k] for k in SLONG_POSITIONS)
    u, v, w = variables[(length(PAIRS) + 1):end]
    generators = vcat(embed.(pluecker.generators), [u - (p25 - p36), v - (p36 - p14), w - (p14 - p25)])
    return (
        ring=ring,
        ideal=ideal(ring, generators),
        generators=generators,
    )
end

"The images of the variables under the signed lift of a taxon permutation, for stage 6."
function slong_coordinate_images(permutation, variables)
    target, signs = slong_signed_permutation(collect(permutation))
    return [signs[k] * variables[target[k]] for k in eachindex(variables)]
end

"""
The ideal of G2/P2 in the coordinates of S_long: the 18 Pluecker coordinates
p_ij with ij not in {14, 25, 36}, and u, v, w in their positions. It is the
image of the Pluecker ideal under `slong_substitution`, plus u + v + w; the
substitution kills the linear form p14 + p25 + p36.
"""
function build_g2_grassmannian_slong_ideal(; coefficient_field=QQ)
    characteristic(coefficient_field) == 3 &&
        throw(ArgumentError("S_long does not span the zero weight space in characteristic 3"))
    pluecker = build_g2_grassmannian_ideal(; coefficient_field=coefficient_field)
    names = ["p$i$j" for (i, j) in PAIRS]
    for (name, position) in zip(("u", "v", "w"), SLONG_POSITIONS)
        names[position] = name
    end
    ring, variables = polynomial_ring(coefficient_field, names)
    phi = hom(pluecker.ring, ring, slong_substitution(variables))
    u, v, w = (variables[k] for k in SLONG_POSITIONS)
    generators = vcat([u + v + w], filter(!iszero, phi.(pluecker.generators)))
    return (
        ring=ring,
        ideal=ideal(ring, generators),
        generators=generators,
    )
end

end
