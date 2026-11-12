include(joinpath(@__DIR__, "g2lib.jl"))
# variables v = (a,b,c,d,x1..x7)
E = QQ[-1 0 1 1 0 1 1 -1 0 0 -1;   # p23-p47
        1 0 0 0 1 0 1 0 -1 0 -1;   # p13-p57
        1 0 0 0 1 1 0 0 0 -1 -1;   # p12-p67
        1 0 -1 -1 1 0 0 0 -1 -1 1; # p17-p56
       -1 0 0 0 0 1 0 -1 0 -1 1;   # p27-p46
       -1 0 0 0 0 0 1 -1 -1 0 1;   # p37-p45
        0 0 0 0 0 1 -1 0 1 -1 0]   # tie x2+x5 = x3+x6
println("rank of the 6 equations: ", rank(E[1:6, :]), ", with the tie: ", rank(E))
N = kernel(E; side = :right)   # 11 x k
println("solution space dim: ", ncols(N))
pts = []
while length(pts) < 60
  v = N * matrix(QQ, ncols(N), 1, rand(-6:6, ncols(N)))
  v = [v[i, 1] for i in 1:11]
  all(>(0), v[1:4]) || continue
  x = v[5:11]
  v[6] + v[9] + v[3] + v[4] >= v[5] + v[8] || continue     # x2+x5+c+d >= x1+x4
  den = lcm([denominator(q) for q in v]); v = [Int(numerator(q * den)) for q in v]
  push!(pts, v)
end
pts = vcat([[1,1,2,2,1,0,0,1,0,0,2]], pts)
fields = [("QQ", g2_ideal(QQ)[2]), ("GF2", g2_ideal(GF(2))[2]), ("GF3", g2_ideal(GF(3))[2])]
bad = Dict("H0" => 0, "H2" => 0, "H3" => 0)
for v in pts
  W = treemetric(v[1], v[2], v[3], v[4], v[5:11])
  tw = termW(W)
  r = Dict(nm => in_trop(I, W) for (nm, I) in fields)
  h0 = count(==(maximum(tw)), tw) >= 2            # all five terms
  h2 = count(==(maximum(tw[2:5])), tw[2:5]) >= 2  # four terms, char 2
  bad["H0"] += r["QQ"] != h0; bad["H3"] += r["GF3"] != h0; bad["H2"] += r["GF2"] != h2
  println(v, "  terms ", tw, "  QQ ", r["QQ"], "  GF3 ", r["GF3"], "  GF2 ", r["GF2"])
end
println("mismatches: ", bad)
