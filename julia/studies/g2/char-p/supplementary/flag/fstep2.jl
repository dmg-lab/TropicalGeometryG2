# Flag variety, good reduction: all checks of step2/step2b/step2c/step2d
# in one pass, for every maximal cone C of the raw flag fan (28 coordinates x1..x7, p12..p67):
#  gb_size, lcs : strong Groebner basis G over ZZ for the weight order of w = sum of the
#                 primitive rays of C (mod lineality); leading coefficients
#  lin_homog    : every g in G is homogeneous for the lineality grading of the fan
#  H, I         : in_w(g) homogeneous along the rays of C; at every ray the initial terms of
#                 g have weight <= all other terms (MIN convention, as in the stored fan)
#  torus        : Z' = in_w(G) + (z_S - 1) + (y prod_{j not in S} z_j - 1) over ZZ (S = coordinates
#                 with an invertible minor of the span of C): leading coefficients of its
#                 strong GB, and whether it contains a nonzero integer constant
# usage: julia fstep2.jl CHUNK NCHUNKS   -> fstep2-CHUNK.csv
using Oscar
chunk, nchunks = parse(Int, ARGS[1]), parse(Int, ARGS[2])
FD = joinpath(@__DIR__, "../../../../../data/g2/flag")
F = load(joinpath(FD, "01_raw_fan.mrdi"))
cones = collect(maximal_cones(F))
stored = load(joinpath(FD, "00_ideal.mrdi"))
names = string.(symbols(base_ring(stored)))
n = length(names)
R, x = polynomial_ring(ZZ, names)
I = ideal(R, [map_coefficients(c -> ZZ(numerator(c)), g; parent = R) for g in gens(stored)])
lin = [Vector{QQFieldElem}(l) for l in eachrow(Matrix{QQFieldElem}(generator_matrix(lineality_space(F))))]
wt(e, v) = sum(QQ(e[i]) * v[i] for i in 1:n)
lin_homog(f) = all(v -> length(unique([wt(e, v) for e in exponents(f)])) <= 1, lin)
@assert all(lin_homog, gens(I))

function primitive_rays(C)
  out = Vector{Int}[]
  for r in rays_modulo_lineality(C).rays_modulo_lineality
    q = Vector{QQFieldElem}(r); d = lcm([denominator(a) for a in q]); z = [numerator(a * d) for a in q]
    g = gcd(z); push!(out, [Int(a ÷ g) for a in z])
  end
  out
end

out = open(joinpath(@__DIR__, "fstep2-$chunk.csv"), "w")
println(out, "cone,gb_size,lcs,lin_homog,H,I,torus_lcs,torus_const"); flush(out)
for k in chunk:nchunks:length(cones)
  C = cones[k]; rs = primitive_rays(C); w = sum(rs); c = maximum(w) + 1
  o = weight_ordering(c .- w, degrevlex(R))
  G = collect(groebner_basis(I; ordering = o))
  lcs = sort(unique([Int(leading_coefficient(g; ordering = o)) for g in G]))
  H = true; In = true; J = elem_type(R)[]
  for g in G
    es = collect(exponents(g)); cs = collect(coefficients(g)); ws = [sum(e .* w) for e in es]; m = minimum(ws)
    ini = [i for i in eachindex(es) if ws[i] == m]; non = [i for i in eachindex(es) if ws[i] > m]
    for r in rs
      a = sum(es[ini[1]] .* r)
      H &= all(i -> sum(es[i] .* r) == a, ini)
      In &= all(i -> sum(es[i] .* r) >= a, non)
    end
    push!(J, sum(cs[i] * prod(x .^ es[i]) for i in ini))
  end
  lh = all(lin_homog, G)
  # torus point certificate
  B = matrix(QQ, reduce(vcat, [[permutedims(QQ.(r)) for r in rs]; [permutedims(v) for v in lin]]))
  r = rank(B); S = Int[]
  for j in 1:n
    rank(B[:, [S; j]]) == length(S) + 1 && push!(S, j)
    length(S) == r && break
  end
  T, z = polynomial_ring(ZZ, [names; "y"])
  emb = hom(R, T, z[1:n])
  rest = [j for j in 1:n if !(j in S)]
  Z = ideal(T, [emb.(J); [z[j] - 1 for j in S]; z[n + 1] * prod(z[rest]) - 1])
  Hb = groebner_basis(Z; ordering = degrevlex(T))
  tl = sort(unique([abs(Int(leading_coefficient(h; ordering = degrevlex(T)))) for h in Hb]))
  tc = join([string(h) for h in Hb if is_constant(h)], ";")
  println(out, k, ",", length(G), ",\"", join(lcs, ";"), "\",", lh, ",", H, ",", In, ",\"", join(tl, ";"), "\",\"", tc, "\"")
  flush(out)
end
println(out, "# done"); close(out)
