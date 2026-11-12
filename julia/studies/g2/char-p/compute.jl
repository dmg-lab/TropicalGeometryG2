# Compute and store the Groebner-basis certificates for the second inclusion in the proof
# of Theorem thm:charp (main.tex), one file per maximal cone of the raw Grassmannian fan:
#   certificates/cones/cone_NNN.mrdi  (a Dict{String,Any}) with
#     "w"            interior point (MIN convention): sum of the primitive rays
#     "G"            strong Groebner basis of I over ZZ for weight_order(R, w)
#     "cofactors"    cofactors[k][j] with G[k] = sum_j cofactors[k][j] * gens_j  (G in I)
#     "pm1_point"    a common zero in {+-1}^21 of all in_w(g), or [] if there is none
#   and, if there is no +-1 point (14 cones), the slice certificate of a torus point:
#     "slice_S"          six coordinates on which the span of the cone projects isomorphically
#     "slice_gens"       in_w(G), p_s - 1 (s in S), y * prod_{t not in S} p_t - 1   (ring ZZ[p, y])
#     "slice_G"          strong Groebner basis of the slice ideal (degrevlex)
#   (no cofactors for slice_G: the argument only needs slice_gens in <slice_G>, checked by
#   division in verify.jl, and that slice_G has standard representations and no constant)
# usage: julia compute.jl CHUNK NCHUNKS
include("common.jl")
chunk, nchunks = parse(Int, ARGS[1]), parse(Int, ARGS[2])
mkpath(joinpath(CERT, "cones"))
R, gens_Z = ideal_ZZ()
I = ideal(R, gens_Z)
F = raw_fan()
cones = collect(maximal_cones(F))

# +-1 torus point: p12 = 1 (global sign), search the other 20 signs
function pm1_point(forms)
  data = [[(Int(c), sum((isodd(e[k]) ? 1 << (k - 2) : 0) for k in 2:21; init = 0)) for (c, e) in zip(coefficients(f), exponents(f))] for f in forms]
  for bits in 0:(1 << 20) - 1
    ok = true
    for f in data
      s = 0
      for (c, m) in f
        s += isodd(count_ones(m & bits)) ? -c : c
      end
      if s != 0; ok = false; break; end
    end
    ok && return [1; [((bits >> (k - 2)) & 1 == 1) ? -1 : 1 for k in 2:21]]
  end
  return Int[]
end

for k in chunk:nchunks:length(cones)
  C = cones[k]; w = interior_point(C); o = weight_order(R, w)
  G = collect(groebner_basis(I; ordering = o))
  cof = [collect(coordinates(g, I)) for g in G]
  J = [initial_form(g, w) for g in G]
  pt = any(f -> length(f) == 1, J) ? Int[] : pm1_point(J)
  d = Dict{String,Any}("w" => w, "G" => G, "cofactors" => cof, "pm1_point" => pt)
  if isempty(pt)
    B = matrix(QQ, reduce(vcat, [permutedims(QQ.(v)) for v in [primitive_rays(C); lineality_generators(C)]]))
    r = rank(B); S = Int[]
    for j in 1:21
      rank(B[:, [S; j]]) == length(S) + 1 && push!(S, j)
      length(S) == r && break
    end
    T, z = polynomial_ring(ZZ, [["p$i$j" for (i, j) in PAIRS]; "y"])
    emb = hom(R, T, z[1:21])
    rest = [j for j in 1:21 if !(j in S)]
    sg = [emb.(J); [z[j] - 1 for j in S]; z[22] * prod(z[rest]) - 1]
    SI = ideal(T, sg)
    H = collect(groebner_basis(SI; ordering = degrevlex(T)))
    d["slice_S"] = S; d["slice_gens"] = sg; d["slice_G"] = H
  end
  save(cone_file(k), d)
  println("cone ", k, ": |G| = ", length(G), isempty(pt) ? ", slice certificate" : ", +-1 point"); flush(stdout)
end
