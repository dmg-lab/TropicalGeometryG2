using Oscar
using Random

"""
    tropical_hyperplane(a, convention=min; fan_structure=:coarse)

The tropicalization of the linear hyperplane `H = {x : a ⋅ x = 0}` in `QQ^n`,
with `a` its rational normal vector. `H` realizes the dual of the rank-one
matroid of the 1×n matrix `a`, so its tropicalization is the Bergman fan of
that dual: the points where the minimum (or maximum) over the support of `a`
is attained at least twice, free in the coordinates where `a` vanishes.
"""
function tropical_hyperplane(a::AbstractVector, convention::Union{typeof(min),typeof(max)}=min;
                             fan_structure::Symbol=:coarse)
    @req count(!iszero, a) >= 2 "the normal vector needs a support of size at least 2, otherwise H misses the torus"
    normal = matroid_from_matrix_columns(matrix(QQ, 1, length(a), a))
    return bergman_fan(dual_matroid(normal), convention; fan_structure)
end

"""
    random_shifted_tropical_hyperplane(n, convention=min; rng, height=10, shift=10, fan_structure=:coarse,
                                       coordinates=1:n, support_sizes=2:length(coordinates), shift_basis=nothing)

A random tropical hyperplane in `QQ^n`, shifted by a random integer vector.
The normal vector has a random support among `coordinates`, of a size drawn
from `support_sizes`, with nonzero entries `p//q` where `0 < |p|, q <= height`.
The shift has entries in `-shift:shift`, or, with an integral matrix
`shift_basis` of `n` columns, is the combination of its rows with random
coefficients in `-shift:shift`. Returns the shifted polyhedral complex
together with the normal vector, the shift and the coefficients (`nothing`
without a basis).
"""
function random_shifted_tropical_hyperplane(n::Int, convention::Union{typeof(min),typeof(max)}=min;
                                            rng::AbstractRNG=Random.default_rng(), height::Int=10,
                                            shift::Int=10, fan_structure::Symbol=:coarse,
                                            coordinates::AbstractVector{Int}=1:n,
                                            support_sizes::AbstractUnitRange{Int}=2:length(coordinates),
                                            shift_basis::Union{Nothing,AbstractMatrix{Int}}=nothing)
    @req n >= 2 "a tropical hyperplane needs ambient dimension at least 2"
    @req length(coordinates) >= 2 && allunique(coordinates) && all(in(1:n), coordinates) "the support needs at least 2 coordinates of 1:n"
    @req first(support_sizes) >= 2 && last(support_sizes) <= length(coordinates) "support sizes must lie in 2:$(length(coordinates))"
    support = coordinates[randperm(rng, length(coordinates))[1:rand(rng, support_sizes)]]
    normal = zeros(QQFieldElem, n)
    for i in support
        normal[i] = rand(rng, [-height:-1; 1:height]) // rand(rng, 1:height)
    end
    if isnothing(shift_basis)
        coefficients, translation = nothing, rand(rng, -shift:shift, n)
    else
        @req size(shift_basis, 2) == n "the shift basis needs $n columns"
        coefficients = rand(rng, -shift:shift, size(shift_basis, 1))
        translation = vec(permutedims(coefficients) * shift_basis)
    end
    fan = tropical_hyperplane(normal, convention; fan_structure)
    return (complex=polyhedral_complex(fan) + QQ.(translation), normal=normal, shift=translation, coefficients=coefficients)
end
