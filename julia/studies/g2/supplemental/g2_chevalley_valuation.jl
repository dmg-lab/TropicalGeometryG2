###############################################################################
#
#   The Chevalley valuation of X_{G_2} in Oscar's split G_2, shared by the
#   scripts of this directory.
#
#   Oscar's lie_algebra(QQ, :G, 2) is the split form over Q, and the weight
#   basis b_1..b_7 of V_7 = simple_module(L, [1,0]) is a split basis, but not
#   ../main.tex's: its invariant 3-form is e147+e156+e237+e246+e345, so the
#   distinguished index is 4 where ../main.tex has 7.  SIGMA is the fixed
#   dictionary between the two; it is derived and checked in
#   g2_oscar_notation.jl.
#
###############################################################################

using TropicalGeometryG2
using Oscar

L = lie_algebra(QQ, :G, 2); V7 = simple_module(L, [1, 0]); W = exterior_power(V7, 2)[1]
pos, neg, hs = chevalley_basis(L)
prs = [(i, j) for i in 1:7 for j in (i+1):7]    # lex order = Oscar's basis of /\^2 V_7

"The basis vector of the module `V` of weight `weight` killed by e_1 and e_2."
function highest_vector(V, weight)
    weights = basis_weights(L, V)
    return only([t for t in 1:dim(V) if weights[t] == weight &&
                 iszero(pos[1] * basis(V, t)) && iszero(pos[2] * basis(V, t))])
end

"""
The torus expansion of the highest vector of weight `weight` of the module `V`
for the paper word `word`, and the Chevalley valuation of each basis vector.
"""
function chevalley_valuations(V, weight, word)
    expansion = torus_expansion(L, V, highest_vector(V, weight), word)
    return expansion, [minimal_exponent(f, expansion.ordering) for f in expansion.vector]
end

const SIGMA = [1, 5, 6, 7, 3, 2, 4]              # Oscar's labels -> ../main.tex's
"The Pluecker coordinate `p` of Oscar's labels in the labels `sigma`."
relabel(p, sigma=SIGMA) = Tuple(sort([sigma[p[1]], sigma[p[2]]]))

# --- rank-one reductions of the valuation ---------------------------------
# degrevlex compares exponent vectors in ring variable order (a_l, ..., a_1).
ringorder(e) = reverse(e)
"degrevlex on Z^l: -1 if e < f, 0 if e = f, +1 if e > f."
function dcmp(e, f)
    u, w = ringorder(e), ringorder(f)
    su, sw = sum(u), sum(w)
    su != sw && return su < sw ? -1 : 1
    d = u .- w
    i = findlast(!iszero, d)
    i === nothing && return 0
    return d[i] < 0 ? 1 : -1
end
"""
lambda: Z^l -> Z, strictly order-preserving for degrevlex on a bounded set when
`B` and `K` are large enough; it turns nu into a rank-one valuation.
"""
lam(e, B, K) = K * sum(e) - sum(B^(j-1) * ringorder(e)[j] for j in eachindex(e))

"The degree of the monomial `m` in the grading of its variables by the vectors `A`."
function Adeg(m, A)
    e = first(AbstractAlgebra.exponent_vectors(m))
    d = zeros(Int, length(first(A)))
    for t in eachindex(A); d .+= e[t] .* A[t]; end
    return d
end
