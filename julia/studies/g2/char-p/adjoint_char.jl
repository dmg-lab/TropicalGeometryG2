# Characteristic checks for the definition of G2 in the introduction of main.tex
# (the SL_7-stabilizer of Phi is smooth in every characteristic) and for Remark
# rem:adjoint-char: the 14-dimensional G2-module E = ker(wedge_2 K^7 -> (K^7)^*, w -> Phi(w, .)) is simple,
# and isomorphic to Lie(G2), exactly when char K != 3.
#
# Method.  Lie(G2) is the stabilizer of the trilinear form Phi in gl_7(QQ).  Its twelve root
# vectors X_a are taken primitive in ZZ^{7x7}; their divided powers X_a^k/k! are integral, so
# the root subgroups x_a(r) = sum_k r^k X_a^k/k! are defined over ZZ.  The coefficients of r^k
# (k >= 1) of x_a(r) acting on a lattice M give the hyperalgebra operators on M.  Over GF(p):
# - dim of the U^+-invariants (common kernel of the operators for positive roots) is 1
#   iff every nonzero submodule contains the highest weight vector;
# - spinning that vector under all operators gives the submodule it generates.
# Both numbers are invariant under field extension, so GF(p) speaks for all K of char p.
# Lattices: E over ZZ (saturated, as the map is onto over ZZ) and the Chevalley lattice of
# Lie(G2) spanned by the X_a and [X_a, X_{-a}], with the adjoint action.
#
# Run (Oscar master, as in ../README.md; about 20 s):
#   env JULIA_DEPOT_PATH=$HOME/local/oscar-dev-depot: julia-1.13 --project=$HOME/oscar-system/Oscar.jl \
#     julia/studies/g2/char-p/adjoint_char.jl
# Expected (logs/adjoint_char.log): the SL_7-stabilizer of Phi has a 14-dimensional Lie algebra
# in every characteristic (elementary divisors all 1; so it is smooth); p = 2 and p = 5 give (1, 14) for both modules;
# p = 3 gives 2-dimensional invariants.
using Oscar

const TRIPLES = [(1, 2, 3), (4, 5, 6), (1, 4, 7), (2, 5, 7), (3, 6, 7)]  # Phi = e123+e456+e147+e257+e367
const PAIRS = [(i, j) for i in 1:6 for j in (i + 1):7]

function phi(i, j, k)
  t = [i, j, k]
  allunique(t) || return 0
  s = sort(t)
  Tuple(s) in TRIPLES || return 0
  return sign(perm(sortperm(t)))   # sign of the permutation sorting t
end

# weight of the matrix unit E_ij is eps_i - eps_j, eps_{i+3} = -eps_i, eps_7 = 0;
# encoded in ZZ^3 modulo (1,1,1)
const EPS = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0], [0, -1, 0], [0, 0, -1], [0, 0, 0]]
normalize_weight(w) = Tuple(w .- minimum(w))
weight(i, j) = normalize_weight(EPS[i] .- EPS[j])
negate(w) = normalize_weight(-collect(w))
height(w) = 7 * w[1] + 3 * w[2] - 10 * w[3]   # generic functional on ZZ^3/(1,1,1)

idx(a, i) = (i - 1) * 7 + a   # column-major position of entry A[a,i]

"The 35x49 integer system A.Phi = 0 for A in gl_7, columns indexed by matrix entries (a,i)."
function stabilizer_system()
  S = zero_matrix(ZZ, 35, 49)
  for (r, (i, j, k)) in enumerate([(i, j, k) for i in 1:5 for j in (i + 1):6 for k in (j + 1):7])
    for a in 1:7
      S[r, idx(a, i)] += phi(a, j, k)
      S[r, idx(a, j)] += phi(i, a, k)
      S[r, idx(a, k)] += phi(i, j, a)
    end
  end
  return S
end

"Lie(G2) over QQ as the kernel of the stabilizer system."
function lie_algebra_basis()
  K = kernel(change_base_ring(QQ, stabilizer_system()); side = :right)
  return [matrix(QQ, [K[idx(a, i), c] for a in 1:7, i in 1:7]) for c in 1:ncols(K)]
end

# Smoothness of the stabilizer scheme of Phi.  Its Lie algebra over a field F of char p is the
# kernel of the stabilizer system over F (plus trace 0 inside sl_7).  Dimension 14 = dim G2 for
# the SL_7-stabilizer means that the scheme is smooth, hence equal to the reduced group G2.
# In GL_7 the identity lies in the Lie algebra in char 3 (I.Phi = 3 Phi): the factor mu_3.
function check_stabilizer_lie_algebra()
  S = stabilizer_system()
  St = vcat(S, matrix(ZZ, 1, 49, [a == i ? 1 : 0 for i in 1:7 for a in 1:7]))
  for p in (0, 2, 3, 5, 7)
    F = p == 0 ? QQ : GF(p)
    dSL = 49 - rank(change_base_ring(F, St))
    dGL = 49 - rank(change_base_ring(F, S))
    println("p = $p: dim Lie(Stab_SL7(Phi)) = $dSL, dim Lie(Stab_GL7(Phi)) = $dGL")
    @assert dSL == 14 && dGL == (p == 3 ? 15 : 14)
  end
  # all primes at once: the nonzero elementary divisors over ZZ
  eSL = filter(!iszero, diagonal(snf(St)))
  eGL = filter(!iszero, diagonal(snf(S)))
  println("elementary divisors: SL_7 system ", unique(eSL), ", GL_7 system ", unique(eGL))
  @assert length(eSL) == 35 && all(isone, eSL)   # rank 35 mod every p: dimension 14 for all p
  @assert length(eGL) == 35 && unique(eGL) == [1, 3]
end

"Primitive integral root vectors, one per nonzero weight."
function root_vectors()
  B = lie_algebra_basis()
  @assert length(B) == 14
  roots = Dict{NTuple{3,Int},ZZMatrix}()
  for w in unique(weight(i, j) for i in 1:7 for j in 1:7)
    w == (0, 0, 0) && continue
    comps = [matrix(QQ, 7, 7, [weight(i, j) == w ? X[i, j] : zero(QQ) for i in 1:7, j in 1:7]) for X in B]
    comps = filter(!iszero, comps)
    isempty(comps) && continue   # e.g. 2 eps_1 = weight(1, 4) is not a root
    @assert rank(matrix(QQ, permutedims(reduce(hcat, [vec(Matrix(C)) for C in comps])))) == 1
    C = comps[1]
    C = C * lcm([denominator(c) for c in C])
    C = C / gcd([numerator(c) for c in C])
    roots[w] = map_entries(c -> ZZ(numerator(c)), C)
  end
  @assert length(roots) == 12
  return roots
end

"x_a(r) = sum_k r^k X^k/k! over ZZ[r] (asserting integrality) and its inverse x_a(-r)."
function root_element(X::ZZMatrix, R, r, s)
  M = identity_matrix(R, 7)
  P = X
  k = 1
  while !iszero(P)
    D = divexact(P, factorial(ZZ(k)))
    @assert D * factorial(ZZ(k)) == P
    M += (s * r)^k * change_base_ring(R, D)
    P *= X
    k += 1
  end
  return M
end

"Action of a 7x7 matrix on wedge_2 in the basis e_i ^ e_j, i < j."
function wedge2(M)
  return matrix(base_ring(M), 21, 21,
    [M[a, i] * M[b, j] - M[b, i] * M[a, j] for (a, b) in PAIRS, (i, j) in PAIRS])
end

"Integer matrices of the coefficients of r^k, k >= 1, of a matrix over ZZ[r]."
function coefficient_matrices(W)
  d = maximum(degree(x) for x in W if !iszero(x))
  return [map_entries(x -> coeff(x, k), W) for k in 1:d]
end

"Hyperalgebra operators of every root on the lattice spanned by the columns of L (basis of M)."
function operators(roots, L, act)
  R, r = polynomial_ring(ZZ, "r")
  ops = Dict{NTuple{3,Int},Vector{ZZMatrix}}()
  for (w, X) in roots
    W = act(root_element(X, R, r, 1), root_element(X, R, r, -1), change_base_ring(R, L))
    ops[w] = [solve(L, Wk; side = :right) for Wk in coefficient_matrices(W)]   # integral: L is saturated
  end
  return ops
end

function test_module(name, ops, p)
  F = GF(p)
  pos = [w for w in keys(ops) if height(w) > 0]
  A = reduce(vcat, [change_base_ring(F, X) for w in pos for X in ops[w]])
  inv = kernel(A; side = :right)
  n = ncols(inv)
  gen = 0
  if n == 1
    allops = [change_base_ring(F, X) for w in keys(ops) for X in ops[w]]
    S = inv
    while true
      S2 = reduce(hcat, [S; [O * S for O in allops]])
      rank(S2) == rank(S) && break
      S = S2
    end
    gen = rank(S)
  end
  println("p = $p, $name: dim U^+-invariants = $n", n == 1 ? ", highest weight vector generates dim $gen" : "")
  return (n, gen)
end

check_stabilizer_lie_algebra()

roots = root_vectors()

# E over ZZ: kernel of the surjection wedge_2 ZZ^7 -> ZZ^7
Map = matrix(ZZ, 7, 21, [phi(i, j, k) for k in 1:7, (i, j) in PAIRS])
EZ = kernel(Map; side = :right)
@assert ncols(EZ) == 14
opsE = operators(roots, EZ, (g, gi, L) -> wedge2(g) * L)

# Chevalley lattice of Lie(G2) in ZZ^49 (column-major vectorization)
vecm(X) = matrix(base_ring(X), 49, 1, vec(Matrix(X)))
unvec(v) = matrix(base_ring(v), [v[(j - 1) * 7 + i, 1] for i in 1:7, j in 1:7])
gens_lie = [collect(values(roots)); [roots[w] * roots[negate(w)] - roots[negate(w)] * roots[w] for w in keys(roots)]]
H = hnf(reduce(vcat, [transpose(vecm(X)) for X in gens_lie]))
LZ = transpose(H[1:rank(H), :])
@assert ncols(LZ) == 14
adjoint_act(g, gi, L) = reduce(hcat, [vecm(g * unvec(L[:, c:c]) * gi) for c in 1:ncols(L)])
opsL = operators(roots, LZ, adjoint_act)

for p in (2, 3, 5)
  e = test_module("E", opsE, p)
  l = test_module("Lie(G2)", opsL, p)
  if p != 3
    @assert e == (1, 14) && l == (1, 14)
  else
    @assert e[1] > 1 && l[1] > 1
  end
end
println("E and Lie(G2) are simple with highest weight theta, hence isomorphic, for p = 2, 5; not for p = 3.")
