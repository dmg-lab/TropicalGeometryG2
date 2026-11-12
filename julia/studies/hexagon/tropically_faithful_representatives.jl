using Oscar

# Finding a tropically faithful representative of a homogeneous polynomial f modulo a homogeneous ideal I, at a point
# w of trop(V(I)). Min convention throughout.
#
# A representative f~ ∈ f + I is tropically faithful at w if in_w(f~) vanishes on no component of V(in_w I) in the torus
# (Proposition prop:tropically_faithful). Then trop(f~)(w) = val(f(x)) for a generic lift x of w.
#
# The minimal primes in step 2 are taken over K instead of the algebraic closure, for K perfect. As f has coefficients
# in K, the components of a K-prime are Galois conjugate: in_w(f~) lies in the K-prime if and only if it vanishes
# on one, equivalently all, of its geometric components. Both tests of the algorithm therefore have the same answer.
#
# For the monomial x^m of step 3 we take a divisor of x_1^k ⋯ x_n^k that is as small as possible, since the polynomials then
# stay much smaller: for the ideal of Remark 9.2 and 28 variables, x_1 ⋯ x_28 gives normal forms with 60000 terms. A single
# variable is preferred (`loop_multiplier`): the choice can change the number of multiplications drastically.

"An integer vector on the ray of a rational vector."
function integer_weight(w)
    v = QQ.(w)
    d = reduce(lcm, denominator.(v); init=ZZ(1))
    return [Int(numerator(c * d)) for c in v]
end

"Largest weight passed to Singular: weight orders with entries of about 60000 and more crashed it (segmentation fault)."
const SINGULAR_MAX_WEIGHT = 30000

"""
A monomial order on `S` refining the weight `w` (min convention) on homogeneous polynomials: leading terms have minimal
`w`-weight, ties are broken by degrevlex. The weights are shifted by a constant to be positive, which does not matter
for homogeneous polynomials.
"""
function min_weight_ordering(S, w)
    wi = integer_weight(w)
    shifted = (maximum(wi) + 1) .- wi
    maximum(shifted) <= SINGULAR_MAX_WEIGHT ||
        throw(ArgumentError("the weight is too large for Singular; pass an ordering of a smaller weight in the same Gröbner cone"))
    return weight_ordering(shifted, degrevlex(S))
end

"""
Run `f()` with Julia's garbage collector switched off, then collect. Singular computations for many different weight
orders crashed (segmentation fault while reading a freshly computed Gröbner basis), which looks like objects being
finalised while in use.
"""
function without_gc(f)
    previous = GC.enable(false)
    try
        return f()
    finally
        GC.enable(previous)
        GC.gc(false)
    end
end

"The data of steps 1-3 of Algorithm 2. It only depends on the Gröbner cone of `weight`, not on `f`."
struct FaithfulSetup
    I                      # the ideal
    weight                 # the weight it was computed for
    ordering               # the monomial order ≺ refining the weight
    groebner_basis         # G, the reduced Gröbner basis of I ∩ K[x_1, ..., x_n] for ≺
    initial_ideal          # in_w(<G>)
    saturation             # J = in_w(<G>) : <x_1 ⋯ x_n>^∞
    primes                 # the minimal primes P_i of J
    k::Int                 # the smallest k with <x_1 ⋯ x_n>^k J ⊂ in_w(<G>)
    multiplier::Vector{Int}  # the exponent vector of the monomial m of the loop, see above
    tropical_map           # the trivial valuation on K, min convention
    normal_forms::Dict     # the normal forms computed so far
end

"""
The exponent vector of a divisor m of x_1^k ⋯ x_n^k with m J ⊂ `inG`: a power of a single variable if one works, smallest
exponent first, else variables are dropped from x_1^k ⋯ x_n^k one by one if possible. Products of variables found this way
made Algorithm 2 need unboundedly many multiplications near some walls of the 44 cones of half_hexagon_cover_fix/, where
a single variable needs one.
"""
function loop_multiplier(inG, J, k)
    S = base_ring(inG)
    e = fill(k, ngens(S))
    k == 0 && return e
    works(v) = (m = prod(gens(S) .^ v); all(g -> m * g in inG, gens(J)))
    for d in 1:k, i in 1:ngens(S)
        single = [j == i ? d : 0 for j in 1:ngens(S)]
        works(single) && return single
    end
    for i in eachindex(e)
        trial = [j == i ? 0 : e[j] for j in eachindex(e)]
        m = prod(gens(S) .^ trial)
        all(g -> m * g in inG, gens(J)) && (e = trial)
    end
    return e
end

"""
    faithful_setup(I, w; ordering=min_weight_ordering(base_ring(I), w), saturate=true)

Steps 1-3 of Algorithm 2 for a homogeneous ideal `I` of a polynomial ring, standing for the ideal it generates in the
Laurent polynomial ring, and a point `w` of trop(V(I)). With `saturate`, `I` is replaced by I ∩ K[x_1, ..., x_n], its
saturation by the product of the variables; pass `saturate=false` if `I` is known to be saturated.
"""
function faithful_setup(I, w; ordering=min_weight_ordering(base_ring(I), w), saturate::Bool=true)
    S = base_ring(I)
    all(is_homogeneous, gens(I)) || throw(ArgumentError("the ideal must be homogeneous"))
    torus = ideal(S, [prod(gens(S))])
    Isat = saturate ? saturation(I, torus) : I
    nu = tropical_semiring_map(base_ring(S))
    G = without_gc(() -> collect(groebner_basis(Isat; ordering, complete_reduction=true)))
    inG = ideal(S, [initial(g, nu, QQ.(w)) for g in G])
    J, k = saturation_with_index(inG, torus)
    is_one(J) && throw(ArgumentError("w is not a point of the tropical variety of I"))
    return FaithfulSetup(Isat, QQ.(w), ordering, G, inG, J, minimal_primes(J), k, loop_multiplier(inG, J, k), nu, Dict())
end

"The normal form of `f` for the reduced Gröbner basis of the setup."
reduce_modulo(T::FaithfulSetup, f) = without_gc(() -> normal_form(f, T.I; ordering=T.ordering))

"""
The value trop(g)(w) of a polynomial g, for the trivial valuation. Oscar's `tropical_polynomial` adds its terms one at a
time, which takes minutes for the polynomials occurring after a multiplication by x_1^k ⋯ x_n^k.
"""
tropical_value(g, w) = minimum(sum(e .* w) for e in AbstractAlgebra.exponent_vectors(g))

"""
    tropically_faithful_representative(f, T::FaithfulSetup; w=T.weight, max_steps=20)
    tropically_faithful_representative(I, w, f; max_steps=20, kwargs...)

The loop of Algorithm 2 for a homogeneous polynomial `f`, at the weight `w`, which must lie in the relative interior of
the Gröbner cone of `T.weight`. Returns a named tuple with

- `status`: `:faithful` if a tropically faithful representative was found, `:no_representative` if none exists and the
  representative is the certificate of Lemma lem:no_faithful_representative, or `:in_ideal` if `f` lies in `I`, which
  the algorithm excludes,
- `representative`, `shift`: the polynomial f~ and the exponent vector u of the representative x^{-u} f~ of `f`,
- `value`: the value of the representative at `w`, i.e. trop(f~)(w) - <u, w> (`PosInf` for `:in_ideal`). For `:faithful`
  this is val(f(x)) of a generic lift x of `w`, in all cases it is a lower bound of val(f(x)),
- `steps`: the number of iterations of the loop that multiplied by a monomial.
"""
function tropically_faithful_representative(f, T::FaithfulSetup; w=T.weight, max_steps::Integer=20)
    is_homogeneous(f) || throw(ArgumentError("f must be homogeneous"))
    S = parent(f)
    w = QQ.(w)
    u = zeros(Int, ngens(S))
    g = get!(() -> reduce_modulo(T, f), T.normal_forms, f)
    for step in 0:max_steps
        if iszero(g)
            return (status=:in_ideal, representative=g, shift=u, value=PosInf(), steps=step)
        end
        h = initial(g, T.tropical_map, w)
        value = tropical_value(g, w) - sum(u .* w)
        any(P -> h in P, T.primes) || return (status=:faithful, representative=g, shift=u, value=value, steps=step)
        h in T.saturation || return (status=:no_representative, representative=g, shift=u, value=value, steps=step)
        g = reduce_modulo(T, prod(gens(S) .^ T.multiplier) * g)
        u = u .+ T.multiplier
    end
    error("no decision after $max_steps steps of Algorithm 2")
end

function tropically_faithful_representative(I, w, f; max_steps::Integer=20, kwargs...)
    return tropically_faithful_representative(f, faithful_setup(I, w; kwargs...); max_steps)
end
