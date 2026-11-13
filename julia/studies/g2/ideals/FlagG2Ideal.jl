module FlagG2Ideal

using Oscar
using TropicalGeometryG2
include(joinpath(@__DIR__, "G2Symmetry.jl"))
using .G2Symmetry: G2_PERMUTATIONS, g2_flag_coordinate_images

export build_flag_g2_ideal,
       build_full_flag_g2_ideal,
       G2_FLAG_PERMUTATIONS,
       flag_coordinate_images,
       original_variable_names,
       reduced_variable_names

const PAIRS = [
    (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7),
    (2, 3), (2, 4), (2, 5), (2, 6), (2, 7),
    (3, 4), (3, 5), (3, 6), (3, 7),
    (4, 5), (4, 6), (4, 7),
    (5, 6), (5, 7),
    (6, 7),
]

# The signed G2 action of G2Symmetry.jl, under the names the census problems use.
const G2_FLAG_PERMUTATIONS = G2_PERMUTATIONS
const flag_coordinate_images = g2_flag_coordinate_images

# Six of the seven G2 linear equations are binomials.  Eliminating their
# left-hand variables is an isomorphism on the dense torus.  This is the main
# inexpensive reduction: 28 ambient coordinates become 22.
const ELIMINATED_P = Set([(2, 3), (1, 3), (1, 2), (5, 6), (4, 6), (4, 5)])
const KEPT_PAIRS = [ij for ij in PAIRS if !(ij in ELIMINATED_P)]

original_variable_names() = vcat(
    ["x$i" for i in 1:7],
    ["p$i$j" for (i, j) in PAIRS],
)

reduced_variable_names() = vcat(
    ["x$i" for i in 1:7],
    ["p$i$j" for (i, j) in KEPT_PAIRS],
)

"""Return a scalar-normalized polynomial, so that `f` and `-f` deduplicate."""
function _normalize_polynomial(f)
    iszero(f) && return f
    c = leading_coefficient(f)
    return inv(c) * f
end

function _deduplicate_polynomials(polynomials)
    result = eltype(polynomials)[]
    seen = Set{String}()
    for f in polynomials
        iszero(f) && continue
        f = _normalize_polynomial(f)
        key = TropicalGeometryG2.gfan_polynomial(f)
        key in seen && continue
        push!(seen, key)
        push!(result, f)
    end
    return result
end

"""
    build_flag_g2_ideal()

Construct the homogeneous flag ideal in the reduced 22-variable torus chart.
The returned named tuple also contains the coordinate map needed to lift
tropical weights back to the original 28 coordinates.
"""
function build_flag_g2_ideal()
    names = reduced_variable_names()
    R, variables = polynomial_ring(QQ, names)

    x = collect(variables[1:7])
    kept_p = Dict(KEPT_PAIRS[i] => variables[7 + i] for i in eachindex(KEPT_PAIRS))

    # Substitution dictated by
    #   p23=-p47, p13=p57, p12=-p67,
    #   p56=p17, p46=-p27, p45=p37.
    p = copy(kept_p)
    p[(2, 3)] = -kept_p[(4, 7)]
    p[(1, 3)] = kept_p[(5, 7)]
    p[(1, 2)] = -kept_p[(6, 7)]
    p[(5, 6)] = kept_p[(1, 7)]
    p[(4, 6)] = -kept_p[(2, 7)]
    p[(4, 5)] = kept_p[(3, 7)]

    pluecker = eltype(x)[]
    for i in 1:4, j in (i + 1):5, k in (j + 1):6, l in (k + 1):7
        push!(
            pluecker,
            p[(i, j)] * p[(k, l)] -
            p[(i, k)] * p[(j, l)] +
            p[(i, l)] * p[(j, k)],
        )
    end

    # The first six contracted 3-form equations vanish after substitution.
    g2_linear = p[(1, 4)] + p[(2, 5)] + p[(3, 6)]

    incidence = eltype(x)[]
    for i in 1:5, j in (i + 1):6, k in (j + 1):7
        push!(
            incidence,
            x[i] * p[(j, k)] - x[j] * p[(i, k)] + x[k] * p[(i, j)],
        )
    end

    quadric = x[1] * x[4] + x[2] * x[5] + x[3] * x[6] - x[7]^2
    generators = _deduplicate_polynomials(vcat(pluecker, [g2_linear], incidence, [quadric]))

    # One entry for each original coordinate: the index of its reduced
    # coordinate.  Coefficient signs disappear under valuation.
    reduced_index = Dict(name => i for (i, name) in enumerate(names))
    lift_indices = Int[]
    append!(lift_indices, 1:7)
    replacement = Dict(
        (2, 3) => (4, 7),
        (1, 3) => (5, 7),
        (1, 2) => (6, 7),
        (5, 6) => (1, 7),
        (4, 6) => (2, 7),
        (4, 5) => (3, 7),
    )
    for ij in PAIRS
        kept_ij = get(replacement, ij, ij)
        push!(lift_indices, reduced_index["p$(kept_ij[1])$(kept_ij[2])"])
    end

    return (
        ring=R,
        ideal=ideal(R, generators),
        generators=generators,
        original_variable_names=original_variable_names(),
        lift_indices=lift_indices,
    )
end

"""
    build_full_flag_g2_ideal()

Build the 78 defining equations in all 28 flag coordinates: seven split G2
linear contractions, 35 Pluecker relations, 35 point-plane incidences, and
the invariant point quadric. No coordinate is eliminated, so the signed G2
group acts directly on the original coordinate ring during tropicalization.
The coefficients are reduced in `coefficient_field`.
"""
function build_full_flag_g2_ideal(; coefficient_field=QQ)
    names = original_variable_names()
    ring, variables = polynomial_ring(coefficient_field, names)
    x = collect(variables[1:7])
    p = Dict(pair => variable for (pair, variable) in zip(PAIRS, variables[8:28]))
    linears = [
        p[(2, 3)] + p[(4, 7)],
        -p[(1, 3)] + p[(5, 7)],
        p[(1, 2)] + p[(6, 7)],
        p[(5, 6)] - p[(1, 7)],
        -p[(4, 6)] - p[(2, 7)],
        p[(4, 5)] - p[(3, 7)],
        p[(1, 4)] + p[(2, 5)] + p[(3, 6)],
    ]
    pluecker = eltype(variables)[]
    for i in 1:4, j in (i + 1):5, k in (j + 1):6, l in (k + 1):7
        push!(pluecker, p[(i, j)] * p[(k, l)] -
                       p[(i, k)] * p[(j, l)] +
                       p[(i, l)] * p[(j, k)])
    end
    incidence = eltype(variables)[]
    for i in 1:5, j in (i + 1):6, k in (j + 1):7
        push!(incidence, x[i] * p[(j, k)] - x[j] * p[(i, k)] + x[k] * p[(i, j)])
    end
    quadric = x[1] * x[4] + x[2] * x[5] + x[3] * x[6] - x[7]^2
    generators = vcat(linears, pluecker, incidence, [quadric])
    return (
        ring=ring,
        ideal=ideal(ring, generators),
        generators=generators,
    )
end

include(joinpath(@__DIR__, "FlagPurityVerifier.jl"))

end # module
