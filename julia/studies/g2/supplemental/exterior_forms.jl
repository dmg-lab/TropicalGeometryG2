###############################################################################
#
#   Exact alternating forms on the index set {1,...,7}, as dictionaries from
#   sorted index tuples to coefficients, and permutations of the indices.
#   Pure Julia, no Oscar; shared by the scripts of this directory.
#
###############################################################################

"sort a triple, returning the sorted triple and the sign of the permutation"
function sortsign(t::NTuple{3,Int})
    v = [t...]; s = 1
    for _ in 1:3, b in 1:2
        v[b] > v[b+1] && ((v[b], v[b+1]) = (v[b+1], v[b]); s = -s)
    end
    return (v[1], v[2], v[3]), s
end

"value of an alternating 3-form on an arbitrary (unsorted, possibly repeating) triple"
function alt(phi::Dict{NTuple{3,Int},T}, i, j, k) where {T}
    (i == j || j == k || i == k) && return zero(T)
    t, s = sortsign((i, j, k))
    return s * get(phi, t, zero(T))
end

"sign of the permutation given by a list of distinct indices, or 0 if repeated"
function permsign(v::Vector{Int})
    length(unique(v)) == length(v) || return 0
    w = copy(v); s = 1
    for i in 1:length(w), j in 1:(length(w)-1)
        w[j] > w[j+1] && ((w[j], w[j+1]) = (w[j+1], w[j]); s = -s)
    end
    return s
end

subsets(n, k) = k == 0 ? [Int[]] :
    [vcat(i, t) for i in 1:n for t in subsets(n, k-1) if isempty(t) || i < t[1]]

"coefficient of e^{1..7} in alpha /\\ beta /\\ gamma, for a 2-, 2- and 3-form"
function top(alpha::Dict, beta::Dict, gamma::Dict{NTuple{3,Int},T}) where {T}
    tot = zero(T)
    for a in subsets(7, 2), b in subsets(7, 2), c in subsets(7, 3)
        length(union(a, b, c)) == 7 || continue
        ca = get(alpha, (a[1], a[2]), zero(T)); iszero(ca) && continue
        cb = get(beta,  (b[1], b[2]), zero(T)); iszero(cb) && continue
        cg = get(gamma, (c[1], c[2], c[3]), zero(T)); iszero(cg) && continue
        tot += permsign(vcat(a, b, c)) * ca * cb * cg
    end
    return tot
end

"the contraction i_a phi of a 3-form, a 2-form"
contract(phi::Dict{NTuple{3,Int},T}, a::Int) where {T} =
    Dict{NTuple{2,Int},T}((j, k) => alt(phi, a, j, k)
                          for j in 1:7 for k in (j+1):7 if !iszero(alt(phi, a, j, k)))

"all permutations of the vector x"
allperms(x) = length(x) <= 1 ? [x] :
    [vcat(x[i], p) for i in 1:length(x) for p in allperms(vcat(x[1:i-1], x[i+1:end]))]
