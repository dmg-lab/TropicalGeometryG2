# First inclusion in the proof of Theorem thm:flagchar2 (main.tex), odd p:
# trop_p(F) is contained in trop(F).  Computes, stores and checks:
#   * the 78 generators of J have coefficients +-1, so trop_p(F) lies in their prevariety,
#     the flag Dressian D of Theorem thm:flag_dressian, for every p;
#   * the 12 taxon permutations of the census act by signed permutations
#     x_i -> s_i x_sigma(i), p_ij -> s_i s_j p_sigma(i)sigma(j) that map every generator into
#     the ZZ-span of the generators of its bidegree and of (variable) * (linear form), hence
#     preserve J over ZZ (also used by compute.jl / verify.jl to certify one cone per orbit);
#   * the three cubics of eq:threecubics lie in J over ZZ (explicit cofactors, identity
#     checked), hence so do their 12 images under the signed action; coefficients 2 and +-1;
#   * the cells of D (coarsest fan, julia/data/g2/flag_dressian/02_coarsest_fan.mrdi): every
#     maximal cone of trop(F) is a cell of D, every cell of dimension <= 7 lies in trop(F), and
#     for each of the remaining cells E (the 15 + 33 cells of dimensions 8 and 9 outside
#     trop(F)) one of the 12 cubics has a term that is strictly minimal on relint(E): its weight
#     is <= that of every other term at each ray of E, < at some ray for each other term, and
#     all terms have equal weight on the lineality space.
# Then for w in D outside trop(F), the cubic mod p (p odd) is in J over F_p and in_w of it is a
# monomial, so w is not in trop_p(F).
# Output: certificates/basis.mrdi, results/basis.txt.     usage: julia --project=julia basis.jl
include("common_flag.jl")
mkpath(FCERT); mkpath(FRES)
io = open(joinpath(FRES, "basis.txt"), "w")
say(a...) = (println(io, a...); println(a...); flush(io))

R, gens_Z = flag_ideal_ZZ(); x = gens(R)
J = ideal(R, gens_Z)
say("generators: ", length(gens_Z), ", all coefficients +-1: ", all(c -> abs(c) == 1, [c for g in gens_Z for c in coefficients(g)]))
bideg(f) = (e = first(exponents(f)); (sum(e[1:7]), sum(e[8:end])))
@assert all(g -> all(e -> (sum(e[1:7]), sum(e[8:end])) == bideg(g), exponents(g)), gens_Z)
kinds = Dict{Tuple{Int,Int},Int}(); for g in gens_Z; kinds[bideg(g)] = get(kinds, bideg(g), 0) + 1; end
say("bidegrees (x, p) of the generators: ", sort(collect(kinds)))

# ---- signed Weyl action ------------------------------------------------------------------
function in_ZZ_span(f, fs)
  mons = unique(vcat([collect(monomials(h)) for h in [fs; f]]...)); mi = Dict(m => k for (k, m) in enumerate(mons))
  A = zero_matrix(ZZ, length(mons), length(fs)); b = zero_matrix(ZZ, length(mons), 1)
  for (j, h) in enumerate(fs), (c, m) in zip(coefficients(h), monomials(h)); A[mi[m], j] = c; end
  for (c, m) in zip(coefficients(f), monomials(f)); b[mi[m], 1] = c; end
  can_solve(A, b; side = :right)
end
lin = filter(g -> bideg(g) == (0, 1), gens_Z)
spanning = Dict(d => [filter(g -> bideg(g) == d, gens_Z);
                      d == (0, 2) ? [v * l for v in x[8:end] for l in lin] : elem_type(R)[];
                      d == (1, 1) ? [v * l for v in x[1:7] for l in lin] : elem_type(R)[]] for d in keys(kinds))
pm = load(joinpath(FDATA, "05_symmetry_taxon_permutations.mrdi"))
perms = [[Int(pm[k, j]) for j in 1:7] for k in 1:nrows(pm)]
weyl = Tuple{Vector{Int},Vector{Int}}[]
for sigma in perms
  found = Int[]
  for bits in 0:127
    s = [((bits >> (k - 1)) & 1 == 1) ? -1 : 1 for k in 1:7]
    img = signed_action(R, sigma, s)
    if all(g -> in_ZZ_span(evaluate(g, img), spanning[bideg(g)]), gens_Z)
      found = s; break
    end
  end
  @assert !isempty(found) "no signs for $sigma"
  push!(weyl, (sigma, found))
  say("sigma = ", sigma, ", signs = ", found, ": every generator maps into J over ZZ")
end
say("all ", length(weyl), " taxon permutations act by signed permutations preserving J over ZZ")

# ---- the cubics of eq:threecubics and their orbit ------------------------------------------
P(i, j) = x[7 + IDX[(i, j)]]
cubics = [2*P(4,5)*P(4,6)*P(5,6) + P(2,3)*P(3,6)*P(5,6) - P(1,4)*P(1,2)*P(4,5) + P(2,6)*P(1,3)*P(4,5) - P(3,5)*P(1,2)*P(4,6),
          2*P(2,3)*P(1,2)*P(1,3) - P(4,6)*P(1,4)*P(1,3) + P(2,5)*P(5,6)*P(2,3) + P(1,6)*P(4,5)*P(2,3) - P(3,4)*P(5,6)*P(1,2),
          2*P(4,6)*P(5,6)*P(4,5) + P(1,2)*P(2,5)*P(4,5) - P(3,6)*P(1,3)*P(4,6) + P(1,5)*P(2,3)*P(4,6) - P(2,4)*P(1,3)*P(5,6)]
cubic_cof = Vector{elem_type(R)}[]
for (c, f) in enumerate(cubics)
  cof = vec(collect(coordinates(f, J)))
  ok = sum(cof[j] * gens_Z[j] for j in eachindex(gens_Z)) == f
  @assert ok
  push!(cubic_cof, cof)
  say("cubic ", c, ": coefficients ", sort(Int.(collect(coefficients(f)))), ", in J over ZZ with integer cofactors: ", ok)
end
normalize_sign(f) = leading_coefficient(f) < 0 ? -f : f
orbit = unique([normalize_sign(evaluate(f, signed_action(R, s...))) for f in cubics for s in weyl])
say("images of the cubics under the signed action: ", length(orbit), " cubics up to sign, coefficients ",
    sort(unique([Int(c) for f in orbit for c in coefficients(f)])))

# ---- cells of the flag Dressian outside trop(F) --------------------------------------------
F = flag_fan()
D = load(joinpath(FDATA, "..", "flag_dressian", "02_coarsest_fan.mrdi"))
linmat(X) = matrix(QQ, Matrix{QQFieldElem}(generator_matrix(lineality_space(X))))
LF, LD = linmat(F), linmat(D)
say("D: ", n_maximal_cones(D), " maximal cones; lineality of D = lineality of trop(F): ",
    rank(vcat(LF, LD)) == rank(LF) == rank(LD))
lin_hom = all(f -> all(v -> length(unique([sum(QQ(e[i]) * v[i] for i in 1:NV) for e in exponents(f)])) == 1,
                       [LD[i, :] for i in 1:nrows(LD)]), orbit)
say("all terms of each cubic have equal weight on the lineality space: ", lin_hom)
canon = canonical_ray_map(F)
cset(C) = Set(canon(r) for r in primitive_rays(C))
trop_sets = [cset(C) for C in maximal_cones(F)]
maxD = collect(maximal_cones(D))
dsets = [cset(C) for C in maxD]
say("maximal cones of D by dimension: ", sort(collect(Dict(k => count(C -> dim(C) == k, maxD) for k in Set(dim.(maxD))))))
in_trop(s) = any(t -> issubset(s, t), trop_sets)
# each maximal cone of trop(F) is a cell of D: a maximal 8-cone of D, or a facet of a 9-cone of D
cells8 = Set(s for (C, s) in zip(maxD, dsets) if dim(C) == 8)
for (C, s) in zip(maxD, dsets)
  dim(C) == 9 || continue
  @assert length(s) == 5                       # simplicial: the facets are the 4-subsets of the rays
  for r in s; push!(cells8, setdiff(s, Set([r]))); end
end
say("cells of D of dimension 8: ", length(cells8), "; every maximal cone of trop(F) is one of them: ", issubset(Set(trop_sets), cells8))
# the cells of D outside trop(F) are faces of maximal cones of D outside trop(F)
bad_max = [k for k in eachindex(maxD) if !in_trop(dsets[k])]
say("maximal cones of D outside trop(F): ", [d => count(k -> dim(maxD[k]) == d, bad_max) for d in (8, 9)])
outside = Set{Tuple{Int,Set{Vector{Int}}}}()
for k in bad_max, d in 5:dim(maxD[k])
  for E in (d == dim(maxD[k]) ? [maxD[k]] : collect(faces(maxD[k], d)))
    s = cset(E); in_trop(s) || push!(outside, (d, s))
  end
end
outside = sort(collect(outside); by = o -> (o[1], sort(collect(o[2]))))
say("cells of D outside trop(F), by dimension: ", [d => count(o -> o[1] == d, outside) for d in 5:9])

# strictly minimal term of f on relint(E), E spanned by the rays rs (and the lineality space)
function dominant_term(f, rs)
  es = collect(exponents(f)); cs = collect(coefficients(f))
  vals = [[sum(e .* r) for e in es] for r in rs]
  for i in eachindex(es)
    abs(cs[i]) in (1, 2) || continue
    all(v -> all(v[i] <= v[j] for j in eachindex(es)), vals) || continue
    all(j -> j == i || any(v -> v[i] < v[j], vals), eachindex(es)) && return i
  end
  nothing
end
# candidates: the three cubics of eq:threecubics first, then their other images
cands = [cubics; [f for f in orbit if !(f in normalize_sign.(cubics))]]
witnesses = Dict{String,Any}[]
nfail = 0
for (k, s) in outside
  rs = collect(s)
  hit = nothing
  for (c, f) in enumerate(cands)
    i = dominant_term(f, rs); i === nothing && continue
    hit = (c, i); break
  end
  if hit === nothing
    global nfail += 1; say("  NO WITNESS for the cell of dimension ", k, " with rays ", rs); continue
  end
  f = cands[hit[1]]; t = collect(exponents(f))[hit[2]]
  push!(witnesses, Dict{String,Any}("dimension" => k, "rays" => rs, "cubic" => hit[1], "term" => t,
                                   "coefficient" => Int(collect(coefficients(f))[hit[2]])))
end
say("cells of D outside trop(F) with a strictly minimal cubic term: ", length(witnesses), " of ", length(outside),
    "; coefficients of the minimal terms: ", sort(unique([w["coefficient"] for w in witnesses])),
    "; witnessed by the three cubics of eq:threecubics themselves: ", count(w -> w["cubic"] <= 3, witnesses))
for d in (8, 9), c in 1:3
  n = count(w -> w["dimension"] == d && w["cubic"] == c, witnesses)
  n > 0 && say("  cubic ", c, " of eq:threecubics: strictly minimal 2-term on ", n, " cells of dimension ", d)
end
save(joinpath(FCERT, "basis.mrdi"), Dict{String,Any}("weyl_sigma" => [w[1] for w in weyl], "weyl_signs" => [w[2] for w in weyl],
     "cubics" => cubics, "cubic_cofactors" => cubic_cof, "candidates" => cands, "witnesses" => witnesses))
say(nfail == 0 && lin_hom ? "RESULT: trop_p(F) is contained in trop(F) for every odd p" : "RESULT: incomplete")
close(io)
