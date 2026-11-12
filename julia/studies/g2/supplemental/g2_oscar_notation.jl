###############################################################################
#
#   OSCAR'S CONVENTIONS FOR SPLIT G_2, AND THE TRANSLATION TO ../main.tex.
#
#   Every other script in this directory reports its results in the labels of
#   ../main.tex.  They can do that because of one permutation, fixed once here.
#   This script is the place where that permutation is derived and checked.
#
#   The point is that OSCAR's weight basis is a split basis of V_7 but not
#   ../main.tex's.  Oscar's lie_algebra(QQ, :G, 2) is the split Chevalley form,
#   simple_module(L, [1,0]) is V_7, and the basis b_1..b_7 that comes out is
#   ordered by weight.  In it the invariant 3-form and quadratic form (each
#   unique up to scale, and solved for below rather than assumed) are
#
#       Phi_osc = e147 + e156 + e237 + e246 + e345,
#       q_osc   = x1x7 - x2x6 - x3x5 + x4^2,
#
#   so the distinguished index -- the one carrying the square of q, equivalently
#   the one lying in three terms of Phi -- is 4, where ../main.tex has 7 with
#
#       Phi = e123 + e456 + e147 + e257 + e367,
#       q   = x1x4 + x2x5 + x3x6 - x7^2.
#
#   Exactly 12 permutations, each lifting to 8 signed maps b_i -> c_i e_sigma(i),
#   carry Phi_osc to Phi.  We fix the representative
#
#       sigma = [1,5,6,7,3,2,4],   c = (-1,+1,+1,+1,+1,+1,+1),
#
#   i.e. b_1 -> -e_1, b_2 -> e_5, b_3 -> e_6, b_4 -> e_7, b_5 -> e_3,
#        b_6 -> e_2, b_7 -> e_4.
#
#   Run with:  julia g2_oscar_notation.jl
#
###############################################################################

using TropicalCensus
using Oscar
const CP = TropicalCensus.ChevalleyPolytopes

L  = lie_algebra(QQ, :G, 2)
V7 = simple_module(L, [1, 0])
W  = exterior_power(V7, 2)[1]
pos, neg, hs = chevalley_basis(L)
prs = [(i, j) for i in 1:7 for j in (i+1):7]
idx = Dict(prs[t] => t for t in 1:21)

wt(V, j) = [Int(numerator(QQ(collect(coefficients(hs[i] * basis(V, j)))[j]))) for i in 1:2]
println("="^75)
println("Oscar's weight basis of V_7")
println("="^75)
println("  weights of b_1..b_7: ", [wt(V7, j) for j in 1:7])

# --- the two invariants of Oscar's basis, solved for ------------------------
A7 = [begin
        M = zero_matrix(QQ, 7, 7)
        for j in 1:7
            c = collect(coefficients(X * basis(V7, j)))
            for i in 1:7; M[i, j] = c[i]; end
        end
        M
      end for X in vcat(pos, neg, hs)]

triples = [(i, j, k) for i in 1:7 for j in (i+1):7 for k in (j+1):7]
tix = Dict(t => n for (n, t) in enumerate(triples))
function alt(F, i, j, k)
    (i == j || j == k || i == k) && return zero(QQ)
    v = [i, j, k]; s = 1
    for _ in 1:3, b in 1:2
        v[b] > v[b+1] && ((v[b], v[b+1]) = (v[b+1], v[b]); s = -s)
    end
    return s * F[tix[(v[1], v[2], v[3])]]
end

rows = Vector{QQFieldElem}[]
for M in A7, (i, j, k) in triples
    row = zeros(QQ, 35)
    for n in 1:35
        F = zeros(QQ, 35); F[n] = one(QQ)
        row[n] = sum(M[m,i]*alt(F,m,j,k) + M[m,j]*alt(F,i,m,k) + M[m,k]*alt(F,i,j,m) for m in 1:7)
    end
    push!(rows, row)
end
KPhi = kernel(matrix(QQ, length(rows), 35, reduce(vcat, rows)); side = :right)
@assert ncols(KPhi) == 1
Phi_osc = [KPhi[n,1] for n in 1:35]
Phi_osc = Phi_osc .* lcm([denominator(x) for x in Phi_osc])
Phi_osc = Phi_osc .// gcd([numerator(x) for x in Phi_osc])
println("\n  invariant 3-form  Phi_osc = ",
        join(["$(c < 0 ? "-" : "+") e$(t[1])$(t[2])$(t[3])"
              for (t, c) in zip(triples, Phi_osc) if !iszero(c)], " "))

sympairs = [(i, j) for i in 1:7 for j in i:7]
six = Dict(p => n for (n, p) in enumerate(sympairs))
bil(F, i, j) = F[six[(min(i,j), max(i,j))]]
rows2 = Vector{QQFieldElem}[]
for M in A7, i in 1:7, j in i:7
    row = zeros(QQ, 28)
    for n in 1:28
        F = zeros(QQ, 28); F[n] = one(QQ)
        row[n] = sum(M[m,i]*bil(F,m,j) + M[m,j]*bil(F,i,m) for m in 1:7)
    end
    push!(rows2, row)
end
KB = kernel(matrix(QQ, length(rows2), 28, reduce(vcat, rows2)); side = :right)
@assert ncols(KB) == 1
Bq = [KB[n,1] for n in 1:28]
Bq = Bq .* lcm([denominator(x) for x in Bq]); Bq = Bq .// gcd([numerator(x) for x in Bq])
println("  invariant form    q_osc   = ",
        join([(p[1] == p[2] ? "$(c//2) x$(p[1])^2" : "$(c) x$(p[1])x$(p[2])")
              for (p, c) in zip(sympairs, Bq) if !iszero(c)], " + "),
        "   (bilinear entries; halve the diagonal)")

# --- the seven linear relations, in Oscar's labels --------------------------
Losc = [Dict((j,k) => alt(Phi_osc, i, j, k) for (j,k) in prs if !iszero(alt(Phi_osc, i, j, k)))
        for i in 1:7]
println("\n  the seven linear forms in Oscar's labels (sum_{j<k} Phi_ijk p_jk):")
for i in 1:7
    println("    i=$i:  ", join([(c > 0 ? "+ " : "- ") * "p$(p[1])$(p[2])"
                                 for (p, c) in sort(collect(Losc[i]))], " "))
end
println("  distinguished index (in three terms of Phi_osc): ",
        only([k for k in 1:7 if count(t -> k in t && !iszero(Phi_osc[tix[t]]), triples) == 3]))

# --- the signed permutations to ../main.tex ---------------------------------
Phi_main = Dict((1,2,3) => 1, (4,5,6) => 1, (1,4,7) => 1, (2,5,7) => 1, (3,6,7) => 1)
allperms(x) = length(x) <= 1 ? [x] :
    [vcat(x[i], p) for i in 1:length(x) for p in allperms(vcat(x[1:i-1], x[i+1:end]))]
function signed_matches()
    out = Tuple{Vector{Int},Vector{Int}}[]
    for s in allperms(collect(1:7))
        Set([Tuple(sort([s[t[1]], s[t[2]], s[t[3]]])) for t in triples if !iszero(Phi_osc[tix[t]])]) ==
            Set(keys(Phi_main)) || continue
        for bits in 0:127, u in (1, -1)
            c = [iszero(bits & (1 << (i-1))) ? 1 : -1 for i in 1:7]
            ok = true
            for t in triples
                iszero(Phi_osc[tix[t]]) && continue
                v = [s[t[1]], s[t[2]], s[t[3]]]; sg = 1
                for _ in 1:3, b in 1:2
                    v[b] > v[b+1] && ((v[b], v[b+1]) = (v[b+1], v[b]); sg = -sg)
                end
                u * Phi_osc[tix[t]] * sg * c[t[1]] * c[t[2]] * c[t[3]] == Phi_main[(v[1],v[2],v[3])] ||
                    (ok = false; break)
            end
            ok && push!(out, (s, u == 1 ? c : -c))
        end
    end
    return out
end
matches = signed_matches()
sigmas = sort(unique(first.(matches)))
println("\n  signed permutations b_i -> c_i e_{sigma(i)} carrying Phi_osc to ../main.tex's Phi: ",
        length(matches), ", on ", length(sigmas), " permutations")
println("  all send Oscar's index 4 to ../main.tex's index 7: ", all(s[4] == 7 for s in sigmas))
const SIGMA = sigmas[1]
signs = first([c for (s, c) in matches if s == SIGMA])
println("  the representative fixed for all scripts: sigma = ", SIGMA, ",  c = ", signs)
println("    ", join(["b$i -> $(signs[i] < 0 ? "-" : "")e$(SIGMA[i])" for i in 1:7], ",  "))

# --- what it does to the Pluecker labels ------------------------------------
println("\n  induced relabelling of the 21 Pluecker coordinates (Oscar -> ../main.tex):")
for row in 1:7
    println("    ", join([begin
                             q = sort([SIGMA[p[1]], SIGMA[p[2]]])
                             "p$(p[1])$(p[2]) -> p$(q[1])$(q[2])"
                          end for p in prs[(3*(row-1)+1):min(3*row, 21)]], "   "))
end

# --- the valuations, in both labellings -------------------------------------
hw = only([t for t in 1:21 if wt(W, t) == [0, 1] &&
           iszero(pos[1] * basis(W, t)) && iszero(pos[2] * basis(W, t))])
r = [2, 1, 2, 1, 2]; l = 5
Fcol = CP._lowering_columns(L, W, 2)
Ra, gs = polynomial_ring(QQ, ["a$k" for k in l:-1:1])
a = Any[gs[l+1-k] for k in 1:l]; o = degrevlex(Ra)
v = Any[zero(Ra) for _ in 1:21]; v[hw] = one(Ra)
for k in 1:l
    global v = CP._exp_f(Fcol[r[k]], a[k], v, Ra)
end
nu = [CP._min_exponent(v[t], o) for t in 1:21]
println("\n  highest weight vector of the g_2-summand of /\\^2 V_7 in Oscar's labels: p",
        prs[hw])
println("\n  the Chevalley valuations, in both labellings:")
println("     Oscar        ../main.tex")
for t in 1:21
    p = prs[t]; q = Tuple(sort([SIGMA[p[1]], SIGMA[p[2]]]))
    println("     p$(p[1])$(p[2])  $(join(nu[t]))     p$(q[1])$(q[2])  $(join(nu[t]))")
end
println("\n  (the vector attached to a coordinate does not change; only its name does)")
