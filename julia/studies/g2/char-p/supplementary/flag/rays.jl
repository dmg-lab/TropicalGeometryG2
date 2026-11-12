# Rays of the coarsest flag fan: image under (a,b) -> b, modulo the lineality of trop(G).
using Oscar
D = joinpath(@__DIR__, "../../../../../data/g2")
FF = load(joinpath(D, "flag/02_coarsest_fan.mrdi")); GG = load(joinpath(D, "grassmannian/02_coarsest_fan.mrdi"))
Lmat = matrix(QQ, Matrix{QQFieldElem}(generator_matrix(lineality_space(GG))))
grays = [Vector{QQFieldElem}(r) for r in rays_modulo_lineality(GG).rays_modulo_lineality]
inlin(v) = rank(vcat(Lmat, matrix(QQ, permutedims(v)))) == rank(Lmat)
function which_ray(v)
  for (k, r) in enumerate(grays)
    M = vcat(Lmat, matrix(QQ, permutedims(r)))
    if rank(vcat(M, matrix(QQ, permutedims(v)))) == rank(M)
      # v = lambda r + lineality with lambda > 0 ?
      c = solve(transpose(M), v; side = :right)
      c[end] > 0 && return k
    end
  end
  return 0
end
nl = 0; nr = 0; other = 0; hit = Set{Int}()
for r in rays_modulo_lineality(FF).rays_modulo_lineality
  v = Vector{QQFieldElem}(r)[8:28]
  if inlin(v); global nl += 1
  else k = which_ray(v); k > 0 ? (global nr += 1; push!(hit, k)) : (global other += 1) end
end
println("flag rays: ", nl + nr + other, "; image in lineality of trop(G): ", nl, "; on a ray of trop(G): ", nr,
        " (", length(hit), " of ", length(grays), " rays hit); other: ", other)
