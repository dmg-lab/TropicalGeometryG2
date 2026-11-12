# Characteristic 2 in the proof of Theorem thm:char2 (main.tex).  Computes and stores:
#   * the maximal cones of the GF(2) raw fan (julia/data/g2/grassmannian_gf2/01_raw_fan.mrdi)
#     that coincide with cones of the QQ raw fan, and the QQ cones that are missing over GF(2)
#     (expected: the 12 boundary triangles of the three Dressian tetrahedra);
#   * the three new rays (type epsilon), one in the relative interior of each tetrahedron,
#     and the 18 new cones, each spanned by an epsilon ray and an edge of its tetrahedron;
#   * the three cubics of eq:threecubics: membership in I over ZZ, and their term weights at
#     the four vertices and at epsilon of their tetrahedron (MIN convention: the initial form
#     consists of the terms of minimal weight);
#   * for each of the 18 new cones C: a Groebner basis G of I over GF(2) for the weight order
#     of an interior point w (refined by degrevlex), the ray conditions that give
#     in_v(g) = in_w(g) on relint C, and a torus point: the all-ones vector if it is a common
#     zero of the in_w(g), else a slice certificate as in compute.jl (Groebner basis over
#     GF(2) of <in_w(G), p_s - 1 (s in S), y prod_{t not in S} p_t - 1> without constant).
#     Stored in certificates/gf2/cone_NNN.mrdi (NNN = index in the GF(2) raw fan).
# Output: results/char2.txt, results/char2.csv.    usage: julia --project=julia char2.jl
include("common.jl")
const DATA2 = joinpath(DATA, "..", "grassmannian_gf2")
mkpath(joinpath(CERT, "gf2")); mkpath(joinpath(@__DIR__, "results"))
io = open(joinpath(@__DIR__, "results", "char2.txt"), "w")
say(a...) = (println(io, a...); println(a...); flush(io))

R, gens_Z = ideal_ZZ(); x = gens(R)
I = ideal(R, gens_Z)
P(i, j) = x[IDX[(i, j)]]
F0 = raw_fan(); F2 = load(joinpath(DATA2, "01_raw_fan.mrdi"))
PV = load(joinpath(DATA, "..", "dressian", "01_prevariety_fan.mrdi"))
cones0 = collect(maximal_cones(F0)); cones2 = collect(maximal_cones(F2))
tets = filter(C -> dim(C) == 7, collect(maximal_cones(PV)))
say("QQ raw fan: ", length(cones0), " maximal cones;  GF(2) raw fan: ", length(cones2),
    " maximal cones;  tetrahedra of the Dressian: ", length(tets))

# ---- lineality and canonical rays modulo lineality ------------------------------------
linmat(F) = matrix(QQ, Matrix{QQFieldElem}(generator_matrix(lineality_space(F))))
L = linmat(F0)
samespace(A, B) = rank(vcat(A, B)) == rank(A) == rank(B)
say("lineality spaces of the QQ fan, the GF(2) fan and the Dressian agree: ",
    samespace(L, linmat(F2)) && samespace(L, linmat(PV)))
# orthogonal projection to the complement of the lineality space, scaled to a primitive vector
const PROJ = identity_matrix(QQ, 21) - transpose(L) * inv(L * transpose(L)) * L
function canon(r)
  v = PROJ * matrix(QQ, 21, 1, QQ.(r)); q = [v[i, 1] for i in 1:21]
  d = lcm([denominator(a) for a in q]); z = [numerator(a * d) for a in q]; g = gcd(z)
  [Int(a ÷ g) for a in z]
end
crays(C) = Set(canon(r) for r in primitive_rays(C))
rays0 = union([crays(C) for C in cones0]...); rays2 = union([crays(C) for C in cones2]...)
say("rays modulo lineality: QQ ", length(rays0), ", GF(2) ", length(rays2),
    ", QQ rays contained in GF(2) rays: ", issubset(rays0, rays2))
eps = collect(setdiff(rays2, rays0))
say("new rays over GF(2): ", length(eps))

# ---- matching of maximal cones --------------------------------------------------------
sets0 = [crays(C) for C in cones0]; sets2 = [crays(C) for C in cones2]
match2 = [findfirst(==(s), sets0) for s in sets2]           # GF(2) cone -> QQ cone or nothing
common = [k for k in eachindex(cones2) if match2[k] !== nothing]
new2 = [k for k in eachindex(cones2) if match2[k] === nothing]
missing0 = setdiff(eachindex(cones0), [match2[k] for k in common])
say("GF(2) cones equal to a QQ cone: ", length(common), ";  new GF(2) cones: ", length(new2),
    ";  QQ cones missing over GF(2): ", sort(missing0))

# ---- tetrahedra, epsilon rays, facets, new cones -------------------------------------
# coefficients of r in the basis (rays of T, lineality generators): unique since dim T = 7
function tet_coords(T, r)
  V = [canon(v) for v in primitive_rays(T)]
  A = matrix(QQ, reduce(hcat, [[QQ(a) for a in v] for v in V]))
  A = hcat(A, transpose(L))
  fl, sol = can_solve_with_solution(A, matrix(QQ, 21, 1, QQ.(r)); side = :right)
  fl ? [sol[i, 1] for i in 1:length(V)] : nothing
end
tet_eps = Int[]; tet_verts = Vector{Vector{Int}}[]
for (t, T) in enumerate(tets)
  V = [canon(v) for v in primitive_rays(T)]; push!(tet_verts, V)
  found = [e for e in eachindex(eps) if (c = tet_coords(T, eps[e])) !== nothing && all(>(0), c)]
  @assert length(found) == 1
  push!(tet_eps, found[1])
  say("tetrahedron ", t, ": ", length(V), " vertices; epsilon = ", eps[found[1]],
      " = (vertices) * ", tet_coords(T, eps[found[1]]), " modulo lineality")
  facets = [k for k in missing0 if issubset(sets0[k], Set(V)) && length(sets0[k]) == 3]
  say("   its four facets are the QQ cones ", facets, " (missing over GF(2))")
  newc = [k for k in new2 if eps[found[1]] in sets2[k]]
  @assert all(k -> length(sets2[k]) == 3 && length(intersect(sets2[k], Set(V))) == 2, newc)
  say("   new GF(2) cones with apex epsilon: ", newc, ", each spanned by epsilon and an edge: ",
      length(newc) == 6 && length(unique([intersect(sets2[k], Set(V)) for k in newc])) == 6)
end
say("every new ray and every new cone belongs to a tetrahedron: ",
    length(unique(tet_eps)) == length(eps) == 3 && all(k -> any(e -> eps[e] in sets2[k], tet_eps), new2))

# ---- the cubics of eq:threecubics ------------------------------------------------------
cubics = [2*P(4,5)*P(4,6)*P(5,6) + P(2,3)*P(3,6)*P(5,6) - P(1,4)*P(1,2)*P(4,5) + P(2,6)*P(1,3)*P(4,5) - P(3,5)*P(1,2)*P(4,6),
          2*P(2,3)*P(1,2)*P(1,3) - P(4,6)*P(1,4)*P(1,3) + P(2,5)*P(5,6)*P(2,3) + P(1,6)*P(4,5)*P(2,3) - P(3,4)*P(5,6)*P(1,2),
          2*P(4,6)*P(5,6)*P(4,5) + P(1,2)*P(2,5)*P(4,5) - P(3,6)*P(1,3)*P(4,6) + P(1,5)*P(2,3)*P(4,6) - P(2,4)*P(1,3)*P(5,6)]
wt(e, r) = sum(e .* r)
for (c, f) in enumerate(cubics)
  cof = vec(collect(coordinates(f, I)))
  say("cubic ", c, ": ", f, "\n   in I over ZZ with integer cofactors: ", sum(cof[j] * gens_Z[j] for j in eachindex(gens_Z)) == f)
  es = collect(exponents(f)); cs = collect(coefficients(f)); i2 = findfirst(c -> abs(c) == 2, cs)
  for (t, T) in enumerate(tets)
    V = tet_verts[t]
    vals = [[wt(e, v) for e in es] for v in V]
    dom = all(v -> all(v[i2] <= v[j] for j in eachindex(es)), vals) && all(j -> j == i2 || any(v -> v[i2] < v[j], vals), eachindex(es))
    dom || continue
    say("   dominant 2-term on the open tetrahedron ", t, " (weight <= at every vertex, < at some vertex for each other term)")
    ve = [wt(e, eps[tet_eps[t]]) for e in es]
    others = [j for j in eachindex(es) if j != i2]
    say("   at epsilon: weights ", ve, "; the four other terms tie: ", length(unique(ve[others])) == 1,
        ", strictly above the 2-term: ", all(j -> ve[j] > ve[i2], others))
    big = Int[]
    for (i, v) in enumerate(vals)
      m = minimum(v); hi = [j for j in others if v[j] > m]; lo = [j for j in others if v[j] == m]
      say("   at vertex ", i, ": weights ", v, "; 2-term minimal: ", v[i2] == m, "; one other term strictly larger, three tie with the 2-term: ",
          length(hi) == 1 && length(lo) == 3)
      append!(big, hi)
    end
    say("   the strictly larger terms at the four vertices are the four terms other than the 2-term: ", sort(big) == sort(others))
  end
end

# ---- Groebner certificates over GF(2) for the new cones --------------------------------
K = GF(2); RK, xK = polynomial_ring(K, string.(symbols(R)))
gens_K = [map_coefficients(c -> K(c), g; parent = RK) for g in gens_Z]
IK = ideal(RK, gens_K)
try
  I2 = load(joinpath(DATA2, "00_ideal.mrdi"))
  say("the census ideal over GF(2) equals <the 42 generators mod 2>: ", ideal(RK, [evaluate(g, xK) for g in gens(I2)]) == IK)
catch err
  say("could not compare with the census ideal over GF(2): ", sprint(showerror, err))
end
csv = open(joinpath(@__DIR__, "results", "char2.csv"), "w")
println(csv, "cone,tetrahedron,gb_size,lin_hom,H,I,torus,torus_kind")
lin = [Vector{QQFieldElem}(l) for l in eachrow(Matrix{QQFieldElem}(L))]
for k in new2
  C = cones2[k]; w = interior_point(C); o = weight_order(RK, w)
  G = collect(groebner_basis(IK; ordering = o))
  rs = primitive_rays(C)
  lin_hom = all(g -> all(v -> length(unique([sum(QQ(e[i]) * v[i] for i in 1:21) for e in exponents(g)])) <= 1, lin), G)
  Hok = true; Iok = true
  for g in G
    es = collect(exponents(g)); ws = [sum(e .* w) for e in es]; m = minimum(ws)
    ini = findall(==(m), ws); non = findall(>(m), ws)
    for r in rs
      a = sum(es[ini[1]] .* r)
      Hok &= all(i -> sum(es[i] .* r) == a, ini)
      Iok &= all(i -> sum(es[i] .* r) >= a, non)
    end
  end
  J = [initial_form(g, w) for g in G]
  d = Dict{String,Any}("w" => w, "G" => G)
  if !any(f -> length(f) == 1, J) && all(f -> iszero(evaluate(f, fill(K(1), 21))), J)
    torus = true; kind = "ones"; d["ones_point"] = true
  else
    kind = "slice"
    B = matrix(QQ, reduce(vcat, [permutedims(QQ.(v)) for v in [rs; lin]]))
    r = rank(B); S = Int[]
    for j in 1:21
      rank(B[:, [S; j]]) == length(S) + 1 && push!(S, j)
      length(S) == r && break
    end
    T, z = polynomial_ring(K, [["p$i$j" for (i, j) in PAIRS]; "y"])
    rest = [j for j in 1:21 if !(j in S)]
    sg = [[evaluate(f, z[1:21]) for f in J]; [z[j] - 1 for j in S]; z[22] * prod(z[rest]) - 1]
    H = collect(groebner_basis(ideal(T, sg); ordering = degrevlex(T)))
    torus = !any(is_constant, H)
    d["slice_S"] = S; d["slice_gens"] = sg; d["slice_G"] = H
  end
  t = findfirst(e -> eps[e] in sets2[k], tet_eps)
  save(joinpath(CERT, "gf2", "cone_" * lpad(k, 3, '0') * ".mrdi"), d)
  println(csv, k, ",", t, ",", length(G), ",", lin_hom, ",", Hok, ",", Iok, ",", torus, ",", kind); flush(csv)
  say("new cone ", k, " (tetrahedron ", t, "): |G| = ", length(G), ", lin_hom = ", lin_hom, ", H = ", Hok, ", I = ", Iok, ", torus point: ", kind, " ", torus)
end
println(csv, "# done"); close(csv)
say("done"); close(io)
