# Check of the fiber description of the tropical flag variety (MIN convention):
#   trop(F) = {(a,b) : b in trop(G), a in L_b},
#   L_b = {a : for all i<j<k, min(a_i+b_jk, a_j+b_ik, a_k+b_ij) is attained at least twice}
# (the incidence trinomials are the circuits of the line).  For b = interior point of each
# raw Grassmannian cone, sample a = min(lambda + v_i, mu + v_j) with v_i = (b_i1,..,b_i7),
# v_i[i] = +infinity (a point of L_b), and a perturbed a' off L_b; Bayer test of (a,b) and
# (a',b) for the flag ideal over the given fields.  Expect: (a,b) in trop, (a',b) not.
# usage: julia fiber.jl CHUNK NCHUNKS FIELD[,FIELD...]   -> fiber-CHUNK.csv
include(joinpath(@__DIR__, "..", "g2lib.jl"))
using Random
chunk, nchunks = parse(Int, ARGS[1]), parse(Int, ARGS[2])
fields = split(ARGS[3], ",")
FD = joinpath(@__DIR__, "../../../../../data/g2/flag")
stored = load(joinpath(FD, "00_ideal.mrdi"))
names = string.(symbols(base_ring(stored)))
@assert names == [["x$i" for i in 1:7]; ["p$i$j" for (i, j) in PAIRS]]
ideals = Dict(f => begin
  R, _ = polynomial_ring(field(f), names)
  ideal(R, [map_coefficients(c -> field(f)(numerator(c)), g; parent = R) for g in gens(stored)])
end for f in fields)
G = load(joinpath(DATA, "01_raw_fan.mrdi"))
cones = collect(maximal_cones(G))
bb(b, i, j) = b[IDX[(min(i, j), max(i, j))]]
on_line(a, b) = all(begin
    t = sort([a[i] + bb(b, j, k), a[j] + bb(b, i, k), a[k] + bb(b, i, j)]); t[1] == t[2]
  end for i in 1:7 for j in i+1:7 for k in j+1:7)
out = open(joinpath(@__DIR__, "fiber-$chunk.csv"), "w"); println(out, "cone,i,j,a,on_line,field,in_trop,perturbed_on_line,perturbed_in_trop"); flush(out)
rng = MersenneTwister(chunk)
for k in chunk:nchunks:length(cones)
  b = interior_point(cones[k])
  i, j = randperm(rng, 7)[1:2]
  lam, mu = rand(rng, -20:20), rand(rng, -20:20)
  a = [min(l == i ? typemax(Int) ÷ 4 : lam + bb(b, i, l), l == j ? typemax(Int) ÷ 4 : mu + bb(b, j, l)) for l in 1:7]
  ol = on_line(a, b)
  a2 = copy(a); a2[rand(rng, 1:7)] += rand(rng, [-1, 1]) * rand(rng, 1:5)
  ol2 = on_line(a2, b)
  for f in fields
    t1 = in_trop(ideals[f], [a; b]); t2 = in_trop(ideals[f], [a2; b])
    println(out, k, ",", i, ",", j, ",\"", join(a, ";"), "\",", ol, ",", f, ",", t1, ",", ol2, ",", t2); flush(out)
  end
end
println(out, "# done"); close(out)
