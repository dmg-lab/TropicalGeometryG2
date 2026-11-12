# Direct construction of the tropical G2 flag variety from the Grassmannian (cross-check of
# the census, see thm:flag in main.tex).  MIN convention, as the stored fans.
#   trop(F) = {(a,b) : b in trop(G), a in L_b},
#   L_b = {a : min(a_i+b_jk, a_j+b_ik, a_k+b_ij) attained at least twice for all i<j<k}.
# For each maximal cell tau of the coarsest Grassmannian fan, with interior point b0:
#  * trace the tree L_b0 along the tropical segments a(t) = min(t + v_i, v_j) between leaves
#    i and j (v_i = (b_i1,...,b_i7), v_i[i] = +infinity); one sample point per edge
#    (midpoints between consecutive breakpoints, plus both unbounded ends);
#  * the tie pattern of a sample (for each triple, the set of minimizing terms) defines the
#    cell  {(a,b) : b in tau, ties of the pattern are equal and <= the other terms};
#  * check dim 8, that (a0,b0) is in its relative interior, and find the census cell of
#    julia/data/g2/flag/02_coarsest_fan.mrdi equal to it (mutual containment).
# Rays of each cell (projected orthogonally to the flag lineality, primitive integral) go
# to the csv for direct_fvector.jl.
# usage: julia direct.jl CHUNK NCHUNKS   -> direct-CHUNK.csv
using Oscar
chunk, nchunks = parse(Int, ARGS[1]), parse(Int, ARGS[2])
D = joinpath(@__DIR__, "../../../../../data/g2")
GG = load(joinpath(D, "grassmannian/02_coarsest_fan.mrdi"))
FF = load(joinpath(D, "flag/02_coarsest_fan.mrdi"))
const PAIRS = [(i, j) for i in 1:6 for j in (i + 1):7]
const IDX = Dict(p => k for (k, p) in enumerate(PAIRS))
const n = 28
bidx(i, j) = IDX[(min(i, j), max(i, j))]
const TRIPLES = [(i, j, k) for i in 1:7 for j in (i + 1):7 for k in (j + 1):7]
# term weight vectors (on R^28) of the trinomial of a triple
function termvecs(t)
  i, j, k = t
  T = [zeros(Int, n) for _ in 1:3]
  T[1][i] = 1; T[1][7 + bidx(j, k)] = 1
  T[2][j] = 1; T[2][7 + bidx(i, k)] = 1
  T[3][k] = 1; T[3][7 + bidx(i, j)] = 1
  T
end
const TV = [termvecs(t) for t in TRIPLES]
function termvals(a, b, t)
  i, j, k = t
  [a[i] + b[bidx(j, k)], a[j] + b[bidx(i, k)], a[k] + b[bidx(i, j)]]
end
function pattern(a, b)
  P = Vector{Vector{Int}}()
  for t in TRIPLES
    v = termvals(a, b, t); m = minimum(v)
    push!(P, findall(==(m), v))
  end
  P
end
function interior_point(C, dimension)
  v = zeros(ZZRingElem, dimension)
  for r in rays_modulo_lineality(C).rays_modulo_lineality
    q = Vector{QQFieldElem}(r); d = lcm([denominator(x) for x in q]); z = [numerator(x * d) for x in q]
    g = gcd(z); v .+= [x ÷ g for x in z]
  end
  [Rational{BigInt}(BigInt(x)) for x in v]
end
# samples on L_b, one per edge (possibly repeated)
function samples(b)
  out = Vector{Vector{Rational{BigInt}}}()
  for i in 1:7, j in 1:7
    i == j && continue
    a_of(t) = [k == i ? b[bidx(i, j)] : k == j ? t + b[bidx(i, j)] : min(t + b[bidx(i, k)], b[bidx(j, k)]) for k in 1:7]
    bp = sort(unique([b[bidx(j, k)] - b[bidx(i, k)] for k in 1:7 if k != i && k != j]))
    ts = [bp[1] - 1; [(bp[m] + bp[m + 1]) / 2 for m in 1:length(bp) - 1]; bp[end] + 1]
    append!(out, [a_of(t) for t in ts])
  end
  out
end
rowmat(vs) = matrix(QQ, reduce(vcat, [permutedims(QQ.(v)) for v in vs]))
function cell_of(P, tau)
  A = Vector{Vector{QQFieldElem}}(); E = Vector{Vector{QQFieldElem}}()
  for h in eachrow(Matrix{QQFieldElem}(linear_inequality_matrix(facets(tau))))
    push!(A, [zeros(QQFieldElem, 7); Vector{QQFieldElem}(h)])
  end
  for e in eachrow(Matrix{QQFieldElem}(linear_equation_matrix(linear_span(tau))))
    push!(E, [zeros(QQFieldElem, 7); Vector{QQFieldElem}(e)])
  end
  for (s, T) in zip(P, TV)
    p = s[1]
    for q in s[2:end]; push!(E, QQ.(T[p] .- T[q])); end
    for r in setdiff(1:3, s); push!(A, QQ.(T[p] .- T[r])); end
  end
  cone_from_inequalities(rowmat(A), rowmat(E))
end
# census cells with H-representation for fast point tests
fcells = collect(maximal_cones(FF))
fH = [(Matrix{QQFieldElem}(linear_inequality_matrix(facets(C))),
       Matrix{QQFieldElem}(linear_equation_matrix(linear_span(C)))) for C in fcells]
contains_pt(H, x) = all(H[1] * x .<= 0) && all(iszero, H[2] * x)
L = matrix(QQ, Matrix{QQFieldElem}(generator_matrix(lineality_space(FF))))
Pr = identity_matrix(QQ, n) - transpose(L) * inv(L * transpose(L)) * L
function prim(v)
  w = Vector{QQFieldElem}((Pr * matrix(QQ, n, 1, v))[:, 1])
  d = lcm([denominator(x) for x in w]); z = [numerator(x * d) for x in w]
  g = gcd(z); [BigInt(x ÷ g) for x in z]
end
gcells = collect(maximal_cones(GG))
out = open(joinpath(@__DIR__, "direct-$chunk.csv"), "w")
println(out, "tau,ncells,cell,dim,relint_ok,census_cell,rays"); flush(out)
for k in chunk:nchunks:length(gcells)
  tau = gcells[k]
  b0 = interior_point(tau, 21)
  pats = Dict{Vector{Vector{Int}},Vector{Rational{BigInt}}}()
  for a in samples(b0)
    P = pattern(a, b0)
    @assert all(length(s) >= 2 for s in P)       # a lies on L_b0
    haskey(pats, P) || (pats[P] = a)
  end
  m = 0
  for (P, a0) in pats
    m += 1
    C = cell_of(P, tau)
    x = QQ.([a0; b0])
    # (a0,b0) in relint C: strict inequalities for the non-minimal terms and for tau's facets
    strict = all(QQ(sum(T[s[1]] .* [a0; b0]) - sum(T[r] .* [a0; b0])) < 0 for (s, T) in zip(P, TV) for r in setdiff(1:3, s))
    cand = [j for (j, H) in enumerate(fH) if contains_pt(H, x)]
    hit = [j for j in cand if issubset(C, fcells[j]) && issubset(fcells[j], C)]
    rs = [prim(Vector{QQFieldElem}(r)) for r in rays_modulo_lineality(C).rays_modulo_lineality]
    println(out, k, ",", length(pats), ",", m, ",", dim(C), ",", strict, ",\"", join(hit, ";"), "\",\"",
            join([join(r, " ") for r in rs], ";"), "\"")
  end
  flush(out)
end
println(out, "# done"); close(out)
