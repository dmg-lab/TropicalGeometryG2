using Oscar
const PAIRS = [[i,j] for i in 1:7 for j in i+1:7]
const IDX = Dict(p => k for (k, p) in enumerate(PAIRS))
function g2_ideal(K)
  R0, _ = polynomial_ring(K, ["p$(i)$(j)" for (i,j) in PAIRS])
  I0 = grassmann_pluecker_ideal(R0, 2, 7); R = base_ring(I0); p = gens(R)
  P(i,j) = p[IDX[[i,j]]]
  lin = [P(2,3)+P(4,7), P(1,3)-P(5,7), P(1,2)+P(6,7), P(1,7)-P(5,6), P(2,7)+P(4,6), P(3,7)-P(4,5), P(1,4)+P(2,5)+P(3,6)]
  return R, I0 + ideal(R, lin), P
end
initial_form(f, w) = begin
  t = [(c, e, sum(e .* w)) for (c, e) in zip(coefficients(f), exponents(f))]
  m = minimum(x[3] for x in t); R = parent(f)
  sum(c * prod(gens(R) .^ e) for (c, e, v) in t if v == m)
end
function bayer_free(J)
  R = base_ring(J); x = gens(R); n = length(x)
  for i in 1:n
    G = groebner_basis(J; ordering = degrevlex([x[[j for j in 1:n if j != i]]; x[i]]))
    H = [divexact(g, x[i]^minimum(m[i] for m in exponents(g))) for g in G]
    any(is_constant, H) && return false
    J = ideal(R, H)
  end
  !is_one(J)
end
# max convention W  ->  min convention w = -W
function in_trop(I, W)
  R = base_ring(I); w = -W; c = maximum(w) + 1
  G = groebner_basis(I; ordering = weight_ordering(c .- w, degrevlex(R)))
  inG = [initial_form(g, w) for g in G]
  any(g -> length(g) == 1, inG) && return false
  bayer_free(ideal(R, inG))
end
treemetric(a,b,c,d,x) = [x[1]+x[2]+a+b+c, x[1]+x[3]+a+b+d, x[1]+x[4], x[1]+x[5]+a+b+d, x[1]+x[6]+a+b+c, x[1]+x[7]+a,
  x[2]+x[3]+c+d, x[2]+x[4]+a+b+c, x[2]+x[5]+c+d, x[2]+x[6], x[2]+x[7]+b+c, x[3]+x[4]+a+b+d, x[3]+x[5], x[3]+x[6]+c+d,
  x[3]+x[7]+b+d, x[4]+x[5]+a+b+d, x[4]+x[6]+a+b+c, x[4]+x[7]+a, x[5]+x[6]+c+d, x[5]+x[7]+b+d, x[6]+x[7]+b+c]
T = treemetric(1,1,2,2,[1,0,0,1,0,0,2])

println("--- the ideal in various characteristics")
for (nm, K) in [("QQ", QQ), ("GF2", GF(2)), ("GF3", GF(3)), ("GF5", GF(5))]
  R, I, P = g2_ideal(K)
  println(nm, ": dim ", dim(I), ", degree ", degree(quo(R, I)[1]))
end

println("--- Bernd's witness")
R, I, P = g2_ideal(QQ)
wit = 2*P(1,7)*P(2,7)*P(3,7) + P(3,6)*P(4,7)*P(5,6) - P(1,4)*P(4,5)*P(6,7) - P(2,6)*P(4,5)*P(5,7) - P(3,5)*P(4,6)*P(6,7)
println("witness in I over QQ: ", wit in I)
R2, I2, P2 = g2_ideal(GF(2))
wit2 = P2(3,6)*P2(4,7)*P2(5,6) + P2(1,4)*P2(4,5)*P2(6,7) + P2(2,6)*P2(4,5)*P2(5,7) + P2(3,5)*P2(4,6)*P2(6,7)
println("its reduction mod 2 (four terms) in I over GF2: ", wit2 in I2)
termW(W) = [W[IDX[[1,7]]]+W[IDX[[2,7]]]+W[IDX[[3,7]]], W[IDX[[3,6]]]+W[IDX[[4,7]]]+W[IDX[[5,6]]],
            W[IDX[[1,4]]]+W[IDX[[4,5]]]+W[IDX[[6,7]]], W[IDX[[2,6]]]+W[IDX[[4,5]]]+W[IDX[[5,7]]], W[IDX[[3,5]]]+W[IDX[[4,6]]]+W[IDX[[6,7]]]]
println("term weights at T: ", termW(T))

println("--- T in trop?")
for (nm, K) in [("QQ", QQ), ("GF2", GF(2)), ("GF3", GF(3)), ("GF5", GF(5))]
  R, I, P = g2_ideal(K)
  println(nm, ": ", in_trop(I, T))
end
