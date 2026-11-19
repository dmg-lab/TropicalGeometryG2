# Second inclusion in the proof of Theorem thm:flagchar2 (main.tex): every maximal cone of
# trop(F) lies in trop_p(F) for odd p.  One certificate per D12-orbit of maximal cones of
# julia/data/g2/flag/01_raw_fan.mrdi (the signed action preserves J over ZZ, basis.jl).
#   certificates/base.mrdi: "w0", "G0" (reduced Groebner basis of J over ZZ for
#     weight_order(R, w0), w0 = the interior point of cone 1) and "cofactors"
#     (G0[k] = sum_j cofactors[k][j] * gens_j over ZZ); every other Groebner basis below is shown
#     to lie in J by division by G0 (verify.jl).
#   certificates/cones/cone_NNNN.mrdi (NNNN = index of the orbit representative C, the smallest
#     index in its orbit): "pieces", a list of cones P covering C, each a Dict with
#     "rays"       primitive rays of P modulo lineality (the lineality space is that of C)
#     "w"          interior point of P (MIN convention): the sum of these rays
#     "G"          reduced Groebner basis of J over ZZ, for weight_order(R, w)
#     "pm1_point"  a common zero in {+-1}^28 of all in_w(g) (exhaustive search), or []
#     and, without a +-1 point, the slice certificate as in ../compute.jl:
#     "slice_S", "slice_gens" (in_w(G), z_s - 1 for s in S, y prod_{t not in S} z_t - 1), "slice_G".
# The pieces: P = C if the ray conditions (ray_violation in common_flag.jl) hold at the sum of
# the rays of C; otherwise the Groebner cones meeting the interior of C (see certify).
# Existing certificates are kept unless --force.
# usage: julia --project=julia compute.jl [--force] [CONE ...]      (default: all representatives)
include("common_flag.jl")
force = "--force" in ARGS
only = [parse(Int, a) for a in ARGS if !startswith(a, "--")]
mkpath(joinpath(FCERT, "cones"))
R, gens_Z = flag_ideal_ZZ()
J = ideal(R, gens_Z)
F = flag_fan()
cones = collect(maximal_cones(F))
pm = load(joinpath(FDATA, "05_symmetry_taxon_permutations.mrdi"))
perms = [coordinate_permutation([Int(pm[k, j]) for j in 1:7]) for k in 1:nrows(pm)]
orbit, img, reps = cone_orbits(F, perms)
println("maximal cones: ", length(cones), ", orbits: ", length(reps)); flush(stdout)
todo = isempty(only) ? reps : only

"Coordinates on which the span of the rays and lineality generators projects isomorphically."
function slice_coordinates(rs, lin)
  B = matrix(QQ, reduce(vcat, [permutedims(QQ.(v)) for v in [rs; lin]]))
  r = rank(B); S = Int[]
  for j in 1:NV
    rank(B[:, [S; j]]) == length(S) + 1 && push!(S, j)
    length(S) == r && break
  end
  S
end

"""
A common zero in {+-1}^28 of the forms, with x1 = p12 = 1 (the forms are bihomogeneous), by
depth-first search over the signs; [] if there is none (or after 10^7 nodes)."""
function pm1_point(forms)
  byvar = [Vector{Vector{Tuple{Int,Int}}}() for _ in 0:NV]     # forms by their last variable
  for f in forms
    terms = [(Int(c), sum((isodd(e[j]) ? 1 << (j - 1) : 0) for j in 1:NV; init = 0)) for (c, e) in zip(coefficients(f), exponents(f))]
    m = reduce(|, [t[2] for t in terms]; init = 0)
    push!(byvar[(m == 0 ? 0 : 64 - leading_zeros(m)) + 1], terms)
  end
  ok(bits, j) = all(t -> sum(isodd(count_ones(m & bits)) ? -c : c for (c, m) in t) == 0, byvar[j + 1])
  ok(0, 0) || return Int[]
  nodes = Ref(0)
  function dfs(j, bits)
    j > NV && return bits
    for neg in (j == 1 || j == 8 ? (false,) : (false, true))
      (nodes[] += 1) > 10^7 && return nothing
      b = neg ? bits | (1 << (j - 1)) : bits
      if ok(b, j)
        r = dfs(j + 1, b); r === nothing || return r
      end
    end
    nothing
  end
  r = dfs(1, 0)
  r === nothing ? Int[] : [((r >> (j - 1)) & 1 == 1) ? -1 : 1 for j in 1:NV]
end

canon = canonical_ray_map(F)
hrep(C) = (Matrix{QQFieldElem}(linear_inequality_matrix(facets(C))), Matrix{QQFieldElem}(linear_equation_matrix(linear_span(C))))
inside(H, v) = all(<=(0), H[1] * v) && (isempty(H[2]) || all(iszero, H[2] * v))
cone_of(rs, lin) = positive_hull(matrix(QQ, reduce(vcat, [permutedims(QQ.(r)) for r in rs])),
                                 matrix(QQ, reduce(vcat, [permutedims(v) for v in lin])))
"Reduced Groebner basis of J over ZZ for the weight order of w."
gb_at(w) = collect(groebner_basis(J; ordering = weight_order(R, w), complete_reduction = true))

"Torus point for the piece with rays rs, G a Groebner basis for w = sum(rs)."
function piece_certificate(rs, lin, G)
  w = sum(rs)
  @assert ray_violation(G, w, rs) === nothing
  ini = [initial_form(g, w) for g in G]
  S = slice_coordinates(rs, lin)
  pt = any(f -> length(f) == 1, ini) ? Int[] : pm1_point(ini)
  d = Dict{String,Any}("rays" => rs, "w" => w, "G" => G, "pm1_point" => pt)
  if isempty(pt)
    T, z = polynomial_ring(ZZ, [NAMES; "y"])
    emb = hom(R, T, z[1:NV])
    rest = [j for j in 1:NV if !(j in S)]
    sg = [emb.(ini); [z[j] - 1 for j in S]; z[NV + 1] * prod(z[rest]) - 1]
    d["slice_S"] = S; d["slice_gens"] = sg
    d["slice_G"] = collect(groebner_basis(ideal(T, sg); ordering = degrevlex(T)))
  end
  d
end

"The cone of the v in C with in_v(g) = in_q(g) for all g in G (closure; MIN convention)."
function groebner_piece(C, q, G)
  ineq = Set{Vector{Int}}(); eq = Set{Vector{Int}}()
  for g in G
    es = collect(exponents(g)); ws = [sum(e .* q) for e in es]; m = minimum(ws); i0 = findfirst(==(m), ws)
    for i in eachindex(es)
      i == i0 && continue
      push!(ws[i] == m ? eq : ineq, es[i0] .- es[i])            # (e_i0 - e_i) . v <= 0, or = 0
    end
  end
  M(S) = isempty(S) ? zero_matrix(QQ, 0, NV) : matrix(QQ, reduce(vcat, [permutedims(QQ.(a)) for a in S]))
  intersect(C, cone_from_inequalities(M(ineq), M(eq)))
end

"""
Certificate pieces for the cone C spanned by rs and lin: C itself if the ray conditions hold
for the reduced Groebner basis at the sum of the rays; otherwise the Groebner cones meeting
the interior of C, found by stepping across the facets of the pieces found so far, until
their volumes add up to the volume of C.
"""
function certify(rs, lin)
  w = sum(rs); G = gb_at(w)
  ray_violation(G, w, rs) === nothing && return [piece_certificate(rs, lin, G)]
  C = cone_of(rs, lin); HC = hrep(C)
  vol = slice_volume_function(C, canon); total = vol(rs)
  primes = [2, 3, 5, 7, 11, 13, 17, 19]
  cands = [sum(primes[mod1(i + s, length(primes))] .* rs[i] for i in eachindex(rs)) for s in 0:3]
  pieces = Tuple{Vector{Vector{Int}},Any,Any}[]          # (rays, cone, G)
  covered = QQ(0); t = QQ(1, 4)
  while covered < total
    if isempty(cands)
      t //= 2; @assert t > QQ(1, 4096) "no progress"
      for (prs, P, _) in pieces, Fc in faces(P, dim(P) - 1)
        frs = primitive_rays(Fc)
        x = (1 + t) .* (QQ.(sum(frs)) .// length(frs)) .- t .* (QQ.(sum(prs)) .// length(prs))
        d = lcm([denominator(a) for a in x]); push!(cands, [Int(numerator(a * d)) for a in x])
      end
    end
    q = popfirst!(cands)
    all(<(0), HC[1] * QQ.(q)) || continue                   # q in the interior of C
    any(p -> inside(hrep(p[2]), QQ.(q)), pieces) && continue
    Gq = gb_at(q); P = groebner_piece(C, q, Gq)
    dim(P) == dim(C) || continue
    prs = primitive_rays(P)
    push!(pieces, (prs, P, Gq)); covered += vol(prs)
    println("    piece ", length(pieces), ": ", length(prs), " rays, covered ", covered // total); flush(stdout)
  end
  @assert covered == total
  [piece_certificate(prs, lin, Gq) for (prs, _, Gq) in pieces]
end

base_file = joinpath(FCERT, "base.mrdi")
if !isfile(base_file) || force
  w0 = interior_point(cones[1]); G0 = gb_at(w0)
  save(base_file, Dict{String,Any}("w0" => w0, "G0" => G0, "cofactors" => [collect(coordinates(g, J)) for g in G0]))
  println("base: |G0| = ", length(G0)); flush(stdout)
end
for k in todo
  file = flag_cone_file(k)
  isfile(file) && !force && haskey(load(file), "pieces") && !any(p -> haskey(p, "cofactors"), load(file)["pieces"]) && continue
  t0 = time()
  pieces = certify(primitive_rays(cones[k]), lineality_generators(cones[k]))
  save(file, Dict{String,Any}("pieces" => pieces))
  o(p) = weight_order(R, Vector{Int}(p["w"]))
  lcs = sort(unique([abs(Int(leading_coefficient(g; ordering = o(p)))) for p in pieces for g in p["G"]]))
  println("cone ", k, " (orbit ", orbit[k], "): ", length(pieces), " piece(s), |G| = ", [length(p["G"]) for p in pieces],
          ", lcs ", lcs, ", torus ", [isempty(p["pm1_point"]) ? "slice" : "pm1" for p in pieces],
          ", ", round(time() - t0; digits = 1), " s"); flush(stdout)
  GC.gc()
end
println("# done"); flush(stdout)
