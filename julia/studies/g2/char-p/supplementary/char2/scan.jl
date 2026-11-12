# Remark rem:char2 of main.tex: test the canonical interior weight of every coarsest cell
# and an interior point of each Dressian tetrahedron over one field.
# usage: julia scan.jl FIELD      -> writes scan-FIELD.csv (one line per point, flushed)
include(joinpath(@__DIR__, "..", "g2lib.jl"))
fname = ARGS[1]
K = field(fname)
W = load(joinpath(DATA, "04_canonical_cell_weights.mrdi"))
pts = [("cell$k", [Int(W[k, j]) for j in 1:21]) for k in 1:size(W, 1)]
F = load(joinpath(DATA, "..", "dressian", "01_prevariety_fan.mrdi"))
for (t, C) in enumerate(filter(C -> dim(C) == 7, maximal_cones(F)))
  v = sum(Vector{QQFieldElem}(r) for r in rays_modulo_lineality(C).rays_modulo_lineality)
  d = lcm([denominator(q) for q in v])
  push!(pts, ("tetra$t", [Int(numerator(q * d)) for q in v]))
end
R, I = g2_ideal(K)
out = open(joinpath(@__DIR__, "scan-$fname.csv"), "w")
println(out, "point,field,in_trop,seconds"); flush(out)
for (name, w) in pts
  t0 = time()
  r = in_trop(I, w)
  println(out, name, ",", fname, ",", r, ",", round(time() - t0; digits = 2)); flush(out)
end
println(out, "# done"); close(out)
