# Shared code for the characteristic-p certificates of the tropical G2 flag variety
# (Theorem thm:flagchar2 in main.tex, odd p).  See README.md in this directory.
# Builds on ../common.jl (weight orders, initial forms, the checks over ZZ); MIN convention.
include(joinpath(@__DIR__, "..", "common.jl"))

const FDATA = joinpath(@__DIR__, "..", "..", "..", "..", "data", "g2", "flag")
const FCERT = joinpath(@__DIR__, "certificates")
const FRES = joinpath(@__DIR__, "results")
const NAMES = [["x$i" for i in 1:7]; ["p$i$j" for (i, j) in PAIRS]]
const NV = length(NAMES)

"The flag ideal J over ZZ: the 78 generators of julia/data/g2/flag/00_ideal.mrdi (coefficients +-1)."
function flag_ideal_ZZ()
  JQ = load(joinpath(FDATA, "00_ideal.mrdi"))
  @assert string.(symbols(base_ring(JQ))) == NAMES
  R, _ = polynomial_ring(ZZ, NAMES)
  gens_Z = [map_coefficients(c -> (@assert isone(denominator(c)); ZZ(numerator(c))), g; parent = R) for g in gens(JQ)]
  return R, gens_Z
end

"trop(F) over QQ: the census fan (01_raw_fan.mrdi equals 02_coarsest_fan.mrdi, 2714 maximal cones)."
flag_fan() = load(joinpath(FDATA, "01_raw_fan.mrdi"))

"x_i -> s_i x_sigma(i), p_ij -> s_i s_j p_sigma(i)sigma(j) (with p_ji = -p_ij)."
function signed_action(R, sigma, s)
  x = gens(R)
  P(i, j) = i < j ? x[7 + IDX[(i, j)]] : -x[7 + IDX[(j, i)]]
  [[s[i] * x[sigma[i]] for i in 1:7]; [s[i] * s[j] * P(sigma[i], sigma[j]) for (i, j) in PAIRS]]
end

"The coordinate permutation of the signed action on R^28 (signs do not change valuations)."
coordinate_permutation(sigma) = [[sigma[i] for i in 1:7]; [7 + IDX[(min(sigma[i], sigma[j]), max(sigma[i], sigma[j]))] for (i, j) in PAIRS]]

"""
    ray_violation(G, w, rs)

The ray conditions giving in_v(g) = in_w(g) on the relative interior of the cone spanned by
the rays rs (and a lineality space for which every g is homogeneous), w = sum(rs), MIN
convention: for every g, the terms of in_w(g) have equal weight at every ray, and every other
term has weight >= theirs at every ray.  Returns nothing if they hold, else the exponent
difference a of a violating pair of terms; then a takes both signs on rs.
"""
function ray_violation(G, w, rs)
  @assert w == sum(rs)
  for g in G
    es = collect(exponents(g)); ws = [sum(e .* w) for e in es]; m = minimum(ws)
    i0 = findfirst(==(m), ws)
    for i in eachindex(es)
      i == i0 && continue
      a = es[i] .- es[i0]
      vals = [sum(a .* r) for r in rs]
      (ws[i] == m ? all(iszero, vals) : all(>=(0), vals)) || return a
    end
  end
  nothing
end

"""
    slice_volume_function(C, canon)

For subcones P of the cone C (given by their primitive rays modulo lineality): the volume of
{v in P/lineality : c(v) <= 1}, c = -(sum of the facet inequalities of C), which is positive
on C outside its lineality space; volumes in the coordinates of a basis of span(C)/lineality.
"""
function slice_volume_function(C, canon)
  A = Matrix{QQFieldElem}(linear_inequality_matrix(facets(C)))
  c = -vec(sum(A; dims = 1))
  B = matrix(QQ, reduce(hcat, [QQ.(canon(r)) for r in primitive_rays(C)]))
  cols = Int[]; for j in 1:ncols(B); rank(B[:, [cols; j]]) == length(cols) + 1 && push!(cols, j); end
  B = B[:, cols]
  coords(v) = (s = can_solve_with_solution(B, matrix(QQ, length(v), 1, v); side = :right); @assert s[1];
               [s[2][i, 1] for i in 1:ncols(B)])
  return function (rs)
    pts = [[QQ(0) for _ in 1:ncols(B)]]
    for r in rs
      v = QQ.(canon(r)); s = sum(c .* v); @assert s > 0
      push!(pts, coords(v .// s))
    end
    volume(convex_hull(matrix(QQ, reduce(vcat, [permutedims(p) for p in pts]))))
  end
end

flag_cone_file(k) =joinpath(FCERT, "cones", "cone_" * lpad(k, 4, '0') * ".mrdi")

"Orthogonal projection to the complement of the lineality space of F, as a primitive integer vector."
function canonical_ray_map(F)
  L = matrix(QQ, Matrix{QQFieldElem}(generator_matrix(lineality_space(F))))
  PROJ = identity_matrix(QQ, ambient_dim(F)) - transpose(L) * inv(L * transpose(L)) * L
  return function (r)
    v = PROJ * matrix(QQ, length(r), 1, QQ.(r)); q = [v[i, 1] for i in 1:length(r)]
    d = lcm([denominator(a) for a in q]); z = [numerator(a * d) for a in q]; g = gcd(z)
    [Int(a ÷ g) for a in z]
  end
end

"""
    cone_orbits(F, perms)

Orbits of the maximal cones of F under coordinate permutations (which must preserve the
lineality space): cone k maps to the cone whose canonical rays are the permuted canonical
rays of k.  Returns (orbit id of each cone, image table img[g][k], representative = smallest
index of each orbit).
"""
function cone_orbits(F, perms)
  canon = canonical_ray_map(F)
  cones = collect(maximal_cones(F))
  sets = [Set(canon(r) for r in primitive_rays(C)) for C in cones]
  index = Dict(s => k for (k, s) in enumerate(sets))
  @assert length(index) == length(cones)
  img = [[index[Set(r[invperm(g)] for r in s)] for s in sets] for g in perms]
  orbit = zeros(Int, length(cones)); reps = Int[]
  for k in eachindex(cones)
    orbit[k] == 0 || continue
    push!(reps, k)
    for t in img; orbit[t[k]] = length(reps); end
  end
  return orbit, img, reps
end
