# Cell types of the census coarsest flag fan: number of rays (mod lineality) per maximal cell.
using Oscar
D = joinpath(@__DIR__, "../../../../../data/g2")
FF = load(joinpath(D, "flag/02_coarsest_fan.mrdi"))
cnt = Dict{Int,Int}()
for C in maximal_cones(FF)
  k = length(rays_modulo_lineality(C).rays_modulo_lineality); cnt[k] = get(cnt, k, 0) + 1
end
println("census 02_coarsest_fan.mrdi, maximal cells by number of rays: ", sort(collect(cnt)))
println("f-vector: ", f_vector(FF))
