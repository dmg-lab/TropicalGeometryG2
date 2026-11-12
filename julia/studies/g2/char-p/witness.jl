# First inclusion in the proof of Theorem thm:charp (main.tex): the tetrahedra of the
# G2-Dressian are not in trop over any field of characteristic != 2.  Computes, stores and
# verifies:
#   * the generators of I have coefficients +-1 (so their prevariety, the Dressian, does not
#     depend on the field);
#   * the 12 taxon permutations of the census act by signed permutations
#     p_ij -> s_i s_j p_{sigma(i) sigma(j)} that map the linear forms into their ZZ-span and
#     the Pluecker quadrics into their ZZ-span (so they preserve I over ZZ);
#   * the witness f = 2 p17 p27 p37 + p36 p47 p56 - p14 p45 p67 - p26 p45 p57 - p35 p46 p67
#     and its Weyl images lie in I over ZZ (explicit cofactors, identity checked);
#   * for each of the 3 tetrahedra (7-dim cones of 01_prevariety_fan.mrdi) a Weyl image of f
#     whose term +-2 p17 p27 p37 is, at each of the 4 rays, <= every other term (MIN
#     convention) and strictly < at some ray for each other term; all terms of f have the
#     same weight along the lineality space.
# Output: certificates/witness.mrdi, results/witness.txt.
include("common.jl")
mkpath(CERT); mkpath(joinpath(@__DIR__, "results"))
io = open(joinpath(@__DIR__, "results", "witness.txt"), "w")
say(a...) = (println(io, a...); println(a...); flush(io))
R, gens_Z = ideal_ZZ(); x = gens(R)
I = ideal(R, gens_Z)
P(i, j) = i < j ? x[IDX[(i, j)]] : -x[IDX[(j, i)]]
say("generators: ", length(gens_Z), ", all coefficients +-1: ", all(c -> abs(c) == 1, [c for g in gens_Z for c in coefficients(g)]))
lin = filter(g -> total_degree(g) == 1, gens_Z); quad = filter(g -> total_degree(g) == 2, gens_Z)
say("linear forms: ", length(lin), ", quadrics: ", length(quad))

# ---- signed Weyl action --------------------------------------------------------------
# f is in the ZZ-span of fs  <=>  integral solution of the linear system on coefficients
function in_ZZ_span(f, fs)
  mons = unique(vcat([collect(monomials(h)) for h in [fs; f]]...)); mi = Dict(m => k for (k, m) in enumerate(mons))
  A = zero_matrix(ZZ, length(mons), length(fs)); b = zero_matrix(ZZ, length(mons), 1)
  for (j, h) in enumerate(fs), (c, m) in zip(coefficients(h), monomials(h)); A[mi[m], j] = c; end
  for (c, m) in zip(coefficients(f), monomials(f)); b[mi[m], 1] = c; end
  can_solve(A, b; side = :right)
end
act(sigma, s, f) = evaluate(f, [s[i] * s[j] * P(sigma[i], sigma[j]) for (i, j) in PAIRS])
pm = load(joinpath(DATA, "05_symmetry_taxon_permutations.mrdi"))
perms = [[Int(pm[k, j]) for j in 1:7] for k in 1:nrows(pm)]
weyl = Tuple{Vector{Int},Vector{Int}}[]
for sigma in perms
  found = Int[]
  for bits in 0:127
    s = [((bits >> (k - 1)) & 1 == 1) ? -1 : 1 for k in 1:7]
    if all(f -> in_ZZ_span(act(sigma, s, f), lin), lin)
      found = s; break
    end
  end
  @assert !isempty(found)
  @assert all(f -> in_ZZ_span(act(sigma, found, f), quad), quad)
  push!(weyl, (sigma, found))
  say("sigma = ", sigma, "  signs = ", found, ": linear forms -> ZZ-span of linear forms, quadrics -> ZZ-span of quadrics")
end
say("all ", length(weyl), " taxon permutations act by signed permutations preserving I over ZZ")

# ---- witness and its images --------------------------------------------------------------
f = 2*P(1,7)*P(2,7)*P(3,7) + P(3,6)*P(4,7)*P(5,6) - P(1,4)*P(4,5)*P(6,7) - P(2,6)*P(4,5)*P(5,7) - P(3,5)*P(4,6)*P(6,7)
cof = vec(collect(coordinates(f, I)))
say("witness f = ", f, ";  f = sum_j c_j gens_j over ZZ: ", sum(cof[j] * gens_Z[j] for j in eachindex(gens_Z)) == f)

PV = load(joinpath(DATA, "..", "dressian", "01_prevariety_fan.mrdi"))
tets = filter(C -> dim(C) == 7, collect(maximal_cones(PV)))
say("tetrahedra (7-dim cones of the prevariety): ", length(tets))
images = unique([act(s[1], s[2], f) for s in weyl])
lin_PV = [Vector{QQFieldElem}(l) for l in eachrow(Matrix{QQFieldElem}(generator_matrix(lineality_space(PV))))]
say("all terms of f have equal weight along the lineality space: ",
    all(v -> length(unique([sum(QQ(e[i]) * v[i] for i in 1:21) for e in exponents(f)])) == 1, lin_PV))
chosen = elem_type(R)[]; chosen_cof = Vector{elem_type(R)}[]
for (t, C) in enumerate(tets)
  rs = [Vector{QQFieldElem}(r) for r in rays_modulo_lineality(C).rays_modulo_lineality]
  good = nothing
  for g in images
    es = collect(exponents(g)); cs = collect(coefficients(g))
    two = findall(c -> abs(c) == 2, cs); length(two) == 1 || continue
    i2 = two[1]
    vals = [[sum(QQ(e[k]) * r[k] for k in 1:21) for e in es] for r in rs]
    if all(v -> all(v[i2] <= v[j] for j in eachindex(v)), vals) && all(j -> j == i2 || any(v -> v[i2] < v[j], vals), eachindex(es))
      good = g; break
    end
  end
  @assert good !== nothing
  c = vec(collect(coordinates(good, I)))
  ok = sum(c[j] * gens_Z[j] for j in eachindex(gens_Z)) == good
  push!(chosen, good); push!(chosen_cof, c)
  say("tetrahedron ", t, " (", length(rs), " rays): ", good, "  [2-term strictly minimal on the open cell; in I over ZZ with cofactors: ", ok, "]")
end
save(joinpath(CERT, "witness.mrdi"), Dict{String,Any}("f" => f, "f_cofactors" => cof,
     "weyl_sigma" => [w[1] for w in weyl], "weyl_signs" => [w[2] for w in weyl],
     "tetrahedron_witnesses" => chosen, "tetrahedron_witness_cofactors" => chosen_cof))
say("saved certificates/witness.mrdi")
close(io)
