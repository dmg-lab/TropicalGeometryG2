# Does the projection (a,b) -> b map each maximal cell of the coarsest flag fan ONTO a
# maximal cell of the coarsest Grassmannian fan?  (claim in the proof of thm:flag)
using Oscar
D = joinpath(@__DIR__, "../../../../../data/g2")
FF = load(joinpath(D, "flag/02_coarsest_fan.mrdi"))
GG = load(joinpath(D, "grassmannian/02_coarsest_fan.mrdi"))
println("flag: ambient ", ambient_dim(FF), ", dim ", dim(FF), ", maximal ", n_maximal_cones(FF),
        "; grassmannian: ambient ", ambient_dim(GG), ", dim ", dim(GG), ", maximal ", n_maximal_cones(GG))
LG = lineality_space(GG)
rowmat(vs) = matrix(QQ, reduce(vcat, [permutedims(Vector{QQFieldElem}(v)) for v in vs]))
gcones = collect(maximal_cones(GG))
proj(v) = Vector{QQFieldElem}(v)[8:28]
hits = Dict{Int,Int}(); onto = 0; inside = 0; none = 0
for C in maximal_cones(FF)
  R = [proj(r) for r in rays_modulo_lineality(C).rays_modulo_lineality]
  L = [proj(l) for l in eachrow(Matrix{QQFieldElem}(generator_matrix(lineality_space(C))))]
  L = filter(!iszero, L); R = filter(!iszero, R)
  P = positive_hull(rowmat(R), rowmat(L))
  ks = [k for (k, S) in enumerate(gcones) if issubset(P, S)]
  if length(ks) == 1
    global inside += 1
    S = gcones[ks[1]]
    issubset(S, P) && (global onto += 1)
    hits[ks[1]] = get(hits, ks[1], 0) + 1
  else
    global none += 1
  end
end
println("flag maximal cells whose image lies in a single maximal Grassmannian cell: ", inside,
        " (onto it: ", onto, "); in no single cell: ", none,
        "; Grassmannian cells hit: ", length(hits), " of ", length(gcones))
