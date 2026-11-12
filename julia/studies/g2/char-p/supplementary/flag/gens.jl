using Oscar
I = load(joinpath(@__DIR__, "../../../../../data/g2/flag/00_ideal.mrdi"))
R = base_ring(I); xs = gens(R)[1:7]; ps = gens(R)[8:28]
kind(g) = begin
  dx = maximum(sum(e[1:7]) for e in exponents(g)); dp = maximum(sum(e[8:28]) for e in exponents(g))
  total_degree(g) == 1 ? "linear in p" : dx == 2 ? "quadric in x" : dx == 1 ? "bilinear (x,p)" : "quadric in p"
end
c = Dict{String,Int}(); for g in gens(I); c[kind(g)] = get(c, kind(g), 0) + 1; end
println(c)
