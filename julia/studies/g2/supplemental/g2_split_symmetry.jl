###############################################################################
#
#   The split G_2 model in the labels of ../main.tex, verified.
#
#   ../main.tex fixes the split trilinear form and quadratic form
#
#     Phi = e123 + e456 + e147 + e257 + e367,
#     q   = x1x4 + x2x5 + x3x6 - x7^2,
#
#   on indices 1..7, with the Pluecker coordinates of Gr(2,7) written p_ij.
#   Section 6 of an earlier note, tropical-mukai-g10.tex, used to work in a different
#   basis -- the Zorn vector-matrix basis, indices 0..6 with h' at 3, and a
#   factor sqrt(-1) coming from the rescaling h' = sqrt(-1) h -- and that basis
#   is NOT related to ../main.tex's by a relabelling alone.  This script
#   establishes the dictionary and re-verifies, in ../main.tex's labels, every
#   claim that Section 6 makes.
#
#   Everything is exact: rationals, and Gaussian rationals where the old Zorn
#   normalisation needs them.  No Oscar, no external packages.
#
#   Run from julia/ with:  julia --project=. studies/g2/supplemental/g2_split_symmetry.jl
#
###############################################################################

const Q = Rational{BigInt}
const C = Complex{Rational{BigInt}}

# ---------------------------------------------------------------- 3-forms ---
# sortsign, alt, permsign, subsets, top, contract and allperms
include(joinpath(@__DIR__, "exterior_forms.jl"))

"apply a relabelling to a 3-form"
function relabel(phi::Dict{NTuple{3,Int},T}, sig::Dict{Int,Int}) where {T}
    out = Dict{NTuple{3,Int},T}()
    for (t, c) in phi
        tt, s = sortsign((sig[t[1]], sig[t[2]], sig[t[3]]))
        out[tt] = c * s
    end
    return out
end

# --------------------------------------------------- exact linear algebra ---
"reduced row echelon form, in place; returns the pivot columns"
function rref!(A::Matrix{T}) where {T}
    m, n = size(A); piv = Int[]; r = 1
    for c in 1:n
        r > m && break
        p = findfirst(k -> !iszero(A[k, c]), r:m)
        p === nothing && continue
        p += r - 1
        A[r, :], A[p, :] = A[p, :], A[r, :]
        A[r, :] ./= A[r, c]
        for k in 1:m
            k != r && !iszero(A[k, c]) && (A[k, :] .-= A[k, c] .* A[r, :])
        end
        push!(piv, c); r += 1
    end
    return piv
end

"a basis of the right kernel of A"
function nullbasis(A::Matrix{T}) where {T}
    B = copy(A); piv = rref!(B); n = size(A, 2)
    free = setdiff(1:n, piv); out = Vector{T}[]
    for f in free
        v = zeros(T, n); v[f] = one(T)
        for (r, c) in enumerate(piv); v[c] = -B[r, f]; end
        push!(out, v)
    end
    return out
end

# ------------------------------------------- the objects of ../main.tex -----
idx = collect(1:7)
triples = [(i, j, k) for i in 1:7 for j in (i+1):7 for k in (j+1):7]   # 35
prs = [(i, j) for i in 1:7 for j in (i+1):7]                            # 21

Phi = Dict{NTuple{3,Int},Q}((1,2,3) => 1, (4,5,6) => 1,
                            (1,4,7) => 1, (2,5,7) => 1, (3,6,7) => 1)

println("="^75)
println("The split G_2 model in ../main.tex's labels")
println("="^75)
println("Phi = ", join(["e$(t[1])$(t[2])$(t[3])" for t in triples if haskey(Phi, t)], " + "))

# --- (1) the annihilator of Phi in gl_7 is 14-dimensional -------------------
rows = Vector{Q}[]
for (i, j, k) in triples
    row = zeros(Q, 49)
    for m in 1:7, n in 1:7
        # unknown M[m,n]; contributes M[m,i]*Phi(m,j,k) etc.
        v = zero(Q)
        n == i && (v += alt(Phi, m, j, k))
        n == j && (v += alt(Phi, i, m, k))
        n == k && (v += alt(Phi, i, j, m))
        row[(m-1)*7 + n] = v
    end
    push!(rows, row)
end
Ann = nullbasis(permutedims(reduce(hcat, rows)))
println("\n(1) dim of the annihilator of Phi in gl_7 = ", length(Ann), "   (expect 14 = dim g_2)")
@assert length(Ann) == 14

# --- (2) the invariant quadratic form ---------------------------------------
sym = [(i, j) for i in 1:7 for j in i:7]                                 # 28
six = Dict(p => n for (n, p) in enumerate(sym))
bidx(i, j) = six[(min(i, j), max(i, j))]
rows2 = Vector{Q}[]
for v in Ann
    M = reshape(v, 7, 7) |> permutedims        # M[m,n] = v[(m-1)*7+n]
    for i in 1:7, j in i:7
        row = zeros(Q, 28)
        for m in 1:7
            row[bidx(m, j)] += M[m, i]
            row[bidx(i, m)] += M[m, j]
        end
        push!(rows2, row)
    end
end
Bs = nullbasis(permutedims(reduce(hcat, rows2)))
println("(2) dim of the space of invariant symmetric forms = ", length(Bs), "   (expect 1)")
@assert length(Bs) == 1
b = Bs[1]; b = b .// b[bidx(1, 4)]             # normalise so that the x1x4 entry is 1
qterms = [(sym[n], b[n]) for n in 1:28 if !iszero(b[n])]
println("    q = ", join([(t[1] == t[2] ? "$(c//2) x$(t[1])^2" : "$(c) x$(t[1])x$(t[2])")
                          for (t, c) in qterms], " + "), "   (bilinear entries; halve the diagonal)")
expected = Dict((1,4) => Q(1), (2,5) => Q(1), (3,6) => Q(1), (7,7) => Q(-2))
@assert Dict(qterms) == expected "the invariant form is not x1x4+x2x5+x3x6-x7^2"
println("    matches q = x1x4 + x2x5 + x3x6 - x7^2 of ../main.tex: true")

# --- (3) the seven linear forms ---------------------------------------------
L = [Dict((j, k) => alt(Phi, i, j, k) for (j, k) in prs if !iszero(alt(Phi, i, j, k)))
     for i in 1:7]
println("\n(3) the seven linear forms L_i(p) = sum_{j<k} Phi_{ijk} p_jk:")
for i in 1:7
    println("    L_$i = ", join([(c > 0 ? "+" : "-") * " p$(p[1])$(p[2])"
                                 for (p, c) in sort(collect(L[i]))], " "))
end
mainforms = [Set([(2,3),(4,7)]), Set([(1,3),(5,7)]), Set([(1,2),(6,7)]),
             Set([(1,7),(5,6)]), Set([(2,7),(4,6)]), Set([(3,7),(4,5)]),
             Set([(1,4),(2,5),(3,6)])]
@assert Set(Set(keys(l)) for l in L) == Set(mainforms) "not eq:7linear of ../main.tex"
println("    these are exactly eq:7linear of ../main.tex: true")

# --- (4) the split matroid, the tripod, and the hexagon ---------------------
supports = [Set(keys(l)) for l in L]
binoms = [s for s in supports if length(s) == 2]
trinom = only([s for s in supports if length(s) == 3])
hexset = setdiff(Set(prs), union(supports...))
println("\n(4) matroid: ", length(binoms), " binomials (6 x U_{1,2}), one trinomial (U_{2,3}) on ",
        sort(collect(trinom)), ", and ", length(hexset), " coloops")
@assert length(binoms) == 6 && length(hexset) == 6
@assert 6 * 1 + 6 * 1 + 2 == 14
hex = sort(collect(hexset))
println("    coloops = ", hex)
adj = Dict(i => Int[] for i in 1:7)
for (a, c) in hex; push!(adj[a], c); push!(adj[c], a); end
@assert all(length(adj[i]) == 2 for i in 1:7 if !isempty(adj[i]))
@assert isempty(adj[7]) "index 7 should lie on no hexagon edge"
function cycle(adj, start)
    cyc = [start, adj[start][1]]
    while true
        nxt = only([x for x in adj[cyc[end]] if x != cyc[end-1]])
        nxt == start && return cyc
        push!(cyc, nxt)
    end
end
hexcyc = cycle(adj, hex[1][1])
println("    they form the hexagon ", join(hexcyc, "-"), "-", hexcyc[1],
        "   (a single 6-cycle: ", length(hexcyc) == 6, ")")

# --- (5) the two symmetry groups --------------------------------------------
"permutations with sigma.phi = lambda*phi for a global scalar lambda"
function stab_scalar(phi::Dict{NTuple{3,Int},T}, idx) where {T}
    out = Tuple{Vector{Int},T}[]
    for p in allperms(idx)
        sig = Dict(idx[i] => p[i] for i in 1:length(idx))
        img = relabel(phi, sig)
        Set(keys(img)) == Set(keys(phi)) || continue
        t0 = first(keys(phi)); lam = img[t0] / phi[t0]
        all(lam * phi[t] == img[t] for t in keys(phi)) && push!(out, (p, lam))
    end
    return out
end

"permutations preserving phi after some diagonal sign change c_i = +-1"
function stab_signed(phi::Dict{NTuple{3,Int},T}, idx) where {T}
    n = length(idx); out = Vector{Int}[]
    for p in allperms(idx)
        sig = Dict(idx[i] => p[i] for i in 1:n)
        img = relabel(phi, sig)
        Set(keys(img)) == Set(keys(phi)) || continue
        ok = false
        for bits in 0:(2^n - 1), u in (1, -1)
            c = Dict(idx[i] => (iszero(bits & (1 << (i-1))) ? 1 : -1) for i in 1:n)
            good = true
            for (t, v) in phi
                v * u * c[t[1]] * c[t[2]] * c[t[3]] == img[t] || (good = false; break)
            end
            good && (ok = true; break)
        end
        ok && push!(out, p)
    end
    return out
end

perm_order(p, idx) = begin
    m = Dict(idx[i] => p[i] for i in 1:length(idx)); cur = copy(m); k = 1
    while any(cur[x] != x for x in idx)
        cur = Dict(x => m[cur[x]] for x in idx); k += 1
    end
    k
end
compose(a, b, idx) = begin      # a after b
    ma = Dict(idx[i] => a[i] for i in 1:length(idx))
    mb = Dict(idx[i] => b[i] for i in 1:length(idx))
    [ma[mb[idx[i]]] for i in 1:length(idx)]
end

A = stab_scalar(Phi, idx); Ap = first.(A)
B = stab_signed(Phi, idx)
println("\n(5) permutations of the seven indices preserving Phi")
println("    up to a global scalar   |A| = ", length(A))
println("    up to a diagonal sign   |B| = ", length(B), "     index [B:A] = ", length(B) ÷ length(A))
ab = all(compose(a, b, idx) == compose(b, a, idx) for a in Ap, b in Ap)
ords = sort([perm_order(p, idx) for p in Ap])
println("    element orders of A = ", ords, ",  abelian = ", ab,
        "  ->  A is ", 6 in ords ? "cyclic Z/6" : "S_3")
@assert ords == [1,2,2,2,3,3] && !ab "A should be S_3 in ../main.tex's labels"
for (p, lam) in sort(A, by = x -> perm_order(x[1], idx))
    println("      ", p, "   order ", perm_order(p, idx), "   lambda = ", lam)
end

# how A sits inside the hexagon's dihedral group
pos = Dict(hexcyc[k] => k for k in 1:6)
function hexaction(p, idx)
    m = Dict(idx[i] => p[i] for i in 1:length(idx))
    all(haskey(pos, m[x]) for x in hexcyc) || return nothing
    sh = [pos[m[hexcyc[k]]] for k in 1:6]
    for r in 0:5
        all(sh[k] == mod1(k + r, 6) for k in 1:6) && return ("rotation by $r", r)
        all(sh[k] == mod1(r - k, 6) for k in 1:6) && return ("reflection", r)
    end
    return ("neither", -1)
end
println("    action on the hexagon ", join(hexcyc, "-"), ":")
for p in Ap
    println("      ", p, "  ->  ", hexaction(p, idx)[1])
end
@assert all(hexaction(p, idx) !== nothing for p in B) "B should preserve the hexagon"
rotB = count(p -> hexaction(p, idx)[1][1:3] == "rot", B)
println("    B preserves the hexagon: true;  ", rotB, " rotations + ", length(B) - rotB,
        " reflections = the dihedral group of order ", length(B))

# --- (6) the dictionary from the old Zorn labels ----------------------------
println("\n(6) dictionary from Section 6's old Zorn basis (indices 0..6, h' at 3)")
one_ = C(1); i_ = C(0, 1)
phi_zorn = Dict{NTuple{3,Int},C}((0,1,2) => one_, (0,3,4) => -i_, (1,3,5) => -i_,
                                 (2,3,6) => -i_, (4,5,6) => -one_)
phi_zorn_real = Dict{NTuple{3,Int},C}((0,1,2) => one_, (0,3,4) => -one_, (1,3,5) => -one_,
                                      (2,3,6) => -one_, (4,5,6) => -one_)   # the h'=h variant
sig = Dict(0 => 1, 1 => 2, 2 => 3, 3 => 7, 4 => 4, 5 => 5, 6 => 6)
println("    relabelling 0,1,2,3,4,5,6 -> ", [sig[i] for i in 0:6],
        "   (h' at 3 becomes index 7)")
for (nm, ph) in (("as printed (h' = sqrt(-1) h)", phi_zorn), ("h' = h variant", phi_zorn_real))
    r = relabel(ph, sig)
    println("    ", nm, ":  ",
            join(["$(c) e$(t[1])$(t[2])$(t[3])" for (t, c) in sort(collect(r), by = x -> x[1])], "  "))
    # solve for a diagonal d_1..d_7 with d_i d_j d_k * r_{ijk} = Phi_{ijk}
    sols = C[]
    found = nothing
    for a in (C(1), C(-1), i_, -i_), bb in (C(1), C(-1), i_, -i_), t in (C(1), C(-1), i_, -i_)
        d = Dict(1 => a, 2 => a, 3 => a, 4 => bb, 5 => bb, 6 => bb, 7 => t)
        all(d[t1[1]] * d[t1[2]] * d[t1[3]] * c == C(Phi[t1]) for (t1, c) in r) &&
            (found = d; break)
    end
    if found === nothing
        println("        no diagonal of that shape; searching all sign patterns...")
        for bits in 0:127
            d = Dict(k => C(iszero(bits & (1 << (k-1))) ? 1 : -1) for k in 1:7)
            all(d[t1[1]] * d[t1[2]] * d[t1[3]] * c == C(Phi[t1]) for (t1, c) in r) &&
                (found = d; break)
        end
    end
    println("        rescaling to ../main.tex's Phi: ",
            found === nothing ? "none found" : join(["x$k -> $(found[k])" for k in 1:7], ", "))
    @assert found !== nothing
    println("        is it a relabelling alone (all entries 1)? ",
            all(found[k] == one_ for k in 1:7))
    A0 = stab_scalar(ph, collect(0:6))
    o0 = sort([perm_order(p, collect(0:6)) for p in first.(A0)])
    println("        its own up-to-scalar stabiliser: order ", length(A0),
            ", element orders ", o0, "  ->  ", 6 in o0 ? "cyclic Z/6" : "S_3")
end
println("\n    Conclusion: the old Zorn form is NOT a relabelling of ../main.tex's Phi -- a")
println("    diagonal rescaling is needed as well -- which is why its up-to-scalar stabiliser")
println("    is the cyclic Z/6 while ../main.tex's is the S_3.  In ../main.tex's labels the")
println("    group is S_3, and the note now uses those labels throughout.")

# --- (7) the two ways of recovering q from Phi -------------------------------
println("\n(7) recovering q from Phi")
g = Dict((a, b) => top(contract(Phi, a), contract(Phi, b), Phi) // 6 for a in 1:7 for b in a:7)
nz = sort([(k, v) for (k, v) in g if !iszero(v)], by = x -> x[1])
println("    eq:metric-from-phi, g(e_a,e_b) = (1/6)(i_a Phi)^(i_b Phi)^Phi:")
println("      ", join(["g(e$(k[1]),e$(k[2])) = $v" for (k, v) in nz], ",  "))
# the polarisation of q = x1x4+x2x5+x3x6-x7^2:  g(x,x) = q(x)
qbil = Dict((1,4) => Q(1)//2, (2,5) => Q(1)//2, (3,6) => Q(1)//2, (7,7) => Q(-1))
c = nz[1][2] // qbil[nz[1][1]]
println("    proportional to the polarisation of q = x1x4+x2x5+x3x6-x7^2, factor ", c, ": ",
        all(v == c * get(qbil, k, Q(0)) for (k, v) in g))
@assert all(g[k] == c * get(qbil, k, Q(0)) for k in keys(g))
@assert c == 1 "with ../main.tex's normalisation of Phi the constant should be exactly 1"

# ../main.tex's own formula q(x) = (x _| Phi) ^ (x _| Phi) ^ Phi, as a quadratic form
coef = Dict{Tuple{Int,Int},Q}()
for a in 1:7, b in 1:7
    v = top(contract(Phi, a), contract(Phi, b), Phi)
    iszero(v) || (coef[(min(a,b), max(a,b))] = get(coef, (min(a,b), max(a,b)), Q(0)) + v)
end
println("    ../main.tex's formula (x _| Phi)^(x _| Phi)^Phi gives the quadratic form")
println("      ", join(["$(v) x$(k[1])x$(k[2])" for (k, v) in sort(collect(coef), by = x -> x[1])], " + "))
qquad = Dict((1,4) => Q(1), (2,5) => Q(1), (3,6) => Q(1), (7,7) => Q(-1))
c2 = coef[(1,4)] // qquad[(1,4)]
println("    = ", c2, " * (x1x4 + x2x5 + x3x6 - x7^2): ",
        all(v == c2 * qquad[k] for (k, v) in coef))
@assert all(coef[k] == c2 * qquad[k] for k in keys(coef))

# --- (8) the Zorn structure constants, in ../main.tex's labels ---------------
# Basis of Im(O_split): u_1,u_2,u_3, w_1,w_2,w_3, h, with (writing h = h'/sqrt(-1),
# i.e. undoing the rescaling of Section 6.3)
#     h x u_i = u_i,   h x w_i = -w_i,
#     u_i x u_j = eps_ijk w_k,   w_i x w_j = -eps_ijk u_k,   u_i x w_j = (1/2) delta_ij h,
#     B(u_i,w_j) = delta_ij/2,  B(h,h) = 1,  all other pairings zero.
# The claim to check: in the basis (u_1,u_2,u_3,-w_1,-w_2,-w_3,-h), labelled 1..7,
# twice the 3-form B(x x y, z) is exactly ../main.tex's Phi.
U = [1, 2, 3]; W = [4, 5, 6]; H = 7                     # positions in the Zorn list
zb = [:u1, :u2, :u3, :w1, :w2, :w3, :h]
eps(i, j, k) = permsign([i, j, k])
function cross(a::Int, b::Int)                           # returns Dict position => coeff
    out = Dict{Int,Q}()
    add!(p, v) = (out[p] = get(out, p, Q(0)) + v)
    if a == H && b in U; add!(b, Q(1))
    elseif a == H && b in W; add!(b, Q(-1))
    elseif b == H; for (p, v) in cross(b, a); add!(p, -v); end
    elseif a in U && b in U
        i, j = a, b; for k in 1:3; e = eps(i, j, k); e != 0 && add!(W[k], Q(e)); end
    elseif a in W && b in W
        i, j = a - 3, b - 3; for k in 1:3; e = eps(i, j, k); e != 0 && add!(U[k], Q(-e)); end
    elseif a in U && b in W
        a == b - 3 && add!(H, Q(1)//2)
    elseif a in W && b in U
        for (p, v) in cross(b, a); add!(p, -v); end
    end
    return out
end
Bz(a::Int, b::Int) = (a in U && b == a + 3) || (b in U && a == b + 3) ? Q(1)//2 :
                     (a == H && b == H ? Q(1) : Q(0))
phi_z(a, b, c) = sum(Q(v) * Bz(p, c) for (p, v) in cross(a, b); init = Q(0))

# the new basis v_1..v_7 = u_1,u_2,u_3,-w_1,-w_2,-w_3,-h
newb = [(1, Q(1)), (2, Q(1)), (3, Q(1)), (4, Q(-1)), (5, Q(-1)), (6, Q(-1)), (7, Q(-1))]
phi_new = Dict{NTuple{3,Int},Q}()
for (i, j, k) in triples
    (pi_, ci) = newb[i]; (pj, cj) = newb[j]; (pk, ck) = newb[k]
    v = 2 * ci * cj * ck * phi_z(pi_, pj, pk)
    iszero(v) || (phi_new[(i, j, k)] = v)
end
println("\n(8) the Zorn construction of Section 6.3, in the basis")
println("      (v_1,...,v_7) = (u_1, u_2, u_3, -w_1, -w_2, -w_3, -h)   [note: h, not h' = sqrt(-1) h]")
println("    twice the 3-form B(x x y, z) is  ",
        join(["$(c > 0 ? "+" : "-") e$(t[1])$(t[2])$(t[3])" for (t, c) in sort(collect(phi_new), by = x -> x[1])], " "))
println("    equals ../main.tex's Phi exactly: ", phi_new == Phi)
@assert phi_new == Phi
qz = Dict{Tuple{Int,Int},Q}()
for i in 1:7, j in i:7
    (pi_, ci) = newb[i]; (pj, cj) = newb[j]
    v = ci * cj * Bz(pi_, pj) * (i == j ? 1 : 2)
    iszero(v) || (qz[(i, j)] = v)
end
println("    the associated quadratic form is  ",
        join([(k[1] == k[2] ? "$(v) x$(k[1])^2" : "$(v) x$(k[1])x$(k[2])")
              for (k, v) in sort(collect(qz), by = x -> x[1])], " + "))
println("    = -(x1x4 + x2x5 + x3x6 - x7^2), i.e. q of ../main.tex up to the global sign: ",
        qz == Dict((1,4) => Q(-1), (2,5) => Q(-1), (3,6) => Q(-1), (7,7) => Q(1)))
@assert qz == Dict((1,4) => Q(-1), (2,5) => Q(-1), (3,6) => Q(-1), (7,7) => Q(1))
println("\n    So no sqrt(-1) is needed anywhere: over Q, the Zorn basis rescaled this way")
println("    gives ../main.tex's Phi on the nose, and its quadratic form is -q.  Since")
println("    {q = 0} = {-q = 0} and O(cB) = O(B), the sign is immaterial.")

# --- (9) is A the Weyl group of the GL_3 Levi factor? ------------------------
# The Levi fixing h is GL_3, acting on the span of the u_i (labels 1,2,3) and on the
# span of the w_i (labels 4,5,6) contragradiently.  Its Weyl group therefore acts by
# SIMULTANEOUS permutations of {1,2,3} and {4,5,6}, fixing 7.
levi = Vector{Int}[]
for p in allperms([1, 2, 3])
    push!(levi, vcat(p, [q + 3 for q in p], 7))
end
println("\n(9) the Weyl group of the GL_3 Levi factor, as permutations of the labels")
println("    |W(Levi)| = ", length(levi), ":  ", levi)
Aset = Set(Ap); Lset = Set(levi); Bset = Set(B)
println("    W(Levi) is contained in B (the order-12 group): ", issubset(Lset, Bset))
println("    W(Levi) equals A (the up-to-scalar stabiliser): ", Aset == Lset)
println("    A cap W(Levi) = ", length(intersect(Aset, Lset)), " elements: ",
        sort(collect(intersect(Aset, Lset))))
println("    so A and W(Levi) are two DIFFERENT order-6 subgroups sharing the 3-cycles.")
# which reflections: through opposite vertices (two fixed labels) or opposite edges (none)
function reflkind(p)
    m = Dict(idx[i] => p[i] for i in 1:7)
    fixed = [x for x in hexcyc if m[x] == x]
    return isempty(fixed) ? "through opposite edges (no fixed vertex)" :
           "through the opposite vertices $(fixed)"
end
println("    reflections in A:")
for p in Ap
    perm_order(p, idx) == 2 && println("      ", p, "  ", reflkind(p))
end
println("    reflections in W(Levi):")
for p in levi
    perm_order(p, idx) == 2 && println("      ", p, "  ", reflkind(p))
end
# both are normal of index 2 in B, hence not conjugate
function isnormal(H, G)
    Hs = Set(H)
    inv_(g) = [findfirst(==(k), g) for k in 1:7]
    all(compose(compose(g, h, idx), inv_(g), idx) in Hs for g in G for h in H)
end
println("    A normal in B: ", isnormal(Ap, B), ";  W(Levi) normal in B: ", isnormal(levi, B))
println("    (both of index 2, so the two S_3's are not conjugate in B)")
