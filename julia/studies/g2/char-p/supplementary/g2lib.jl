# Shared code for the characteristic-p study of the G2-Grassmannian.
# MIN convention throughout (as Oscar and the stored fans in julia/data/g2).
using Oscar
const PAIRS = [(i, j) for i in 1:6 for j in (i + 1):7]          # p12, ..., p67
const IDX = Dict(p => k for (k, p) in enumerate(PAIRS))
const DATA = joinpath(@__DIR__, "../../../../data/g2/grassmannian")

# same generators as julia/studies/g2/ideals/G2GrassmannianIdeal.jl, over any field
function g2_ideal(K)
  R, x = polynomial_ring(K, ["p$i$j" for (i, j) in PAIRS])
  P(i, j) = x[IDX[(i, j)]]
  lin = [P(2,3) + P(4,7), -P(1,3) + P(5,7), P(1,2) + P(6,7), P(5,6) - P(1,7),
         -P(4,6) - P(2,7), P(4,5) - P(3,7), P(1,4) + P(2,5) + P(3,6)]
  quad = elem_type(R)[]
  for (i, j, k, l) in [(i, j, k, l) for i in 1:7 for j in i+1:7 for k in j+1:7 for l in k+1:7]
    push!(quad, P(i,j)*P(k,l) - P(i,k)*P(j,l) + P(i,l)*P(j,k))
  end
  return R, ideal(R, vcat(lin, quad))
end

function initial_form(f, w)
  t = [(c, e, sum(e .* w)) for (c, e) in zip(coefficients(f), exponents(f))]
  m = minimum(x[3] for x in t); R = parent(f)
  sum(c * prod(gens(R) .^ e) for (c, e, v) in t if v == m)
end

# J : (prod x)^infty == 1 ?  via Bayer's trick (homogeneous J)
function bayer_monomial_free(J)
  R = base_ring(J); x = gens(R); n = length(x)
  for i in 1:n
    G = groebner_basis(J; ordering = degrevlex([x[[j for j in 1:n if j != i]]; x[i]]))
    H = [divexact(g, x[i]^minimum(m[i] for m in exponents(g))) for g in G]
    any(is_constant, H) && return false
    J = ideal(R, H)
  end
  return !is_one(J)
end

# min convention, integral weight w
function in_trop(I, w::Vector{Int})
  R = base_ring(I); c = maximum(w) + 1
  G = groebner_basis(I; ordering = weight_ordering(c .- w, degrevlex(R)))
  inG = [initial_form(g, w) for g in G]
  any(g -> length(g) == 1, inG) && return false
  return bayer_monomial_free(ideal(R, inG))
end

field(name) = Dict("QQ" => QQ, "GF2" => GF(2), "GF3" => GF(3), "GF5" => GF(5), "GF7" => GF(7),
                   "GF11" => GF(11), "GF13" => GF(13))[name]

# Small interior point of a cone with lineality: sum of the primitive integer ray
# vectors (modulo lineality).  Keep weights small: huge weights (~1e5) crashed Singular
# (packed weighted degrees), see PLAN.md progress log.
function interior_point(C)
  v = zeros(ZZRingElem, 21)
  for r in rays_modulo_lineality(C).rays_modulo_lineality
    q = Vector{QQFieldElem}(r)
    d = lcm([denominator(a) for a in q]); z = [numerator(a * d) for a in q]
    g = gcd(z); v .+= [a ÷ g for a in z]
  end
  return [Int(a) for a in v]
end
