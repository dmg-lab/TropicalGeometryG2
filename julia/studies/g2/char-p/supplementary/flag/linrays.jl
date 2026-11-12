# The flag rays whose image lies in the lineality space of trop(G): print them, and
# the same for a few rays over rays of trop(G), with a modulo (1,..,1).
using Oscar
D = joinpath(@__DIR__, "../../../../../data/g2")
FF = load(joinpath(D, "flag/02_coarsest_fan.mrdi")); GG = load(joinpath(D, "grassmannian/02_coarsest_fan.mrdi"))
LF = Matrix{QQFieldElem}(generator_matrix(lineality_space(FF)))
println("flag lineality generators:"); foreach(r -> println("  ", r), eachrow(LF))
Lmat = matrix(QQ, Matrix{QQFieldElem}(generator_matrix(lineality_space(GG))))
inlin(v) = rank(vcat(Lmat, matrix(QQ, permutedims(v)))) == rank(Lmat)
for r in rays_modulo_lineality(FF).rays_modulo_lineality
  v = Vector{QQFieldElem}(r)
  inlin(v[8:28]) || continue
  # write b = lineality combination; subtract the flag lineality part that carries it
  M = matrix(QQ, LF)
  c = solve(M[:, 8:28], v[8:28]; side = :left)   # c * M[:,8:28] = b
  w = v .- [sum(c[i] * LF[i, j] for i in 1:size(LF, 1)) for j in 1:28]
  println("ray: a = ", v[1:7], "  reduced by flag lineality: ", w[1:7], "  b-part after: ", all(iszero, w[8:28]))
end
