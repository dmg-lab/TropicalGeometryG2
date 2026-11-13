###############################################################################
# The curve C ∩ (B + t) of a translate B + t of the uniform Bergman fan B (rank
# 17 on 21 elements) with the cells C of a fan of the tropical G2-Grassmannian,
# exactly, cell by cell; a library of method_1_curve.jl, curve_rays.jl,
# core_trivalence.jl and bounded_cell_valences.jl.
#
# The piece C ∩ (B_T + t), |T| = 5, lies on the line of span(C) modulo the
# diagonal on which x - t is constant on T, so it is an interval of that line,
# cut by the 16 inequalities of B_T and the facets of C (`exact_piece`). A pass
# in floating point first keeps, with a margin, the subsets T whose piece may
# be more than a point (`candidate_subsets`). A piece of dimension at least 2
# means that the translation is not transverse to the cell. A piece of the
# curve has multiplicity m_C [ZZ^21 : L_C + L_T], with m_C the multiplicity of
# C over QQbar, L_C = span(C) ∩ ZZ^21 and L_T = {x_T constant}.
###############################################################################

include(joinpath(@__DIR__, "decorated_curve_graph.jl"))
import Combinatorics
import LinearAlgebra

"The census orbits of the cells whose part of the curve is a polytope for every generic translation (shifted_bergman_cells.jl)."
const BOUNDED_ORBITS = (16, 17, 19, 22)
"The range of the coordinates of random translations."
const SHIFT_RANGE = -(1 << 20):(1 << 20)
const FLOAT_MARGIN = 1e-6

# polymake, behind `exact_piece` for the subsets of deficient rank, is not thread-safe.
const POLYMAKE_LOCK = ReentrantLock()

"""
A basis of the span of the cone `C`, its lineality first, as the columns of a
21×6 matrix, and the matrix `A` of its facets, `A x ≤ 0`.
"""
function cell_span(C)
    rays = [collect(r) for r in rays_modulo_lineality(C).rays_modulo_lineality]
    columns = [[collect(r) for r in lineality_space(C)]; rays]
    basis = Int[]
    for i in eachindex(columns)
        rank(matrix(QQ, reduce(hcat, columns[[basis; i]]))) > length(basis) && push!(basis, i)
    end
    M = matrix(QQ, reduce(hcat, columns[basis]))
    A = matrix(QQ, linear_inequality_matrix(facets(C)))
    @req ncols(M) == 6 && all(<(0), A * matrix(QQ, 21, 1, sum(rays))) "the cone is not as expected"
    return M, A
end

"The direction `x` modulo the diagonal, with first coordinate 0 and largest coordinate ±1."
direction(x) = (v = normalized(x); v .// maximum(abs, v))

"""
The piece `C ∩ (B_T + t)` of the cell of span `M` and facets `A`, `AM = A * M`:
`nothing` when it is empty or a point, `(:segment, a, b)` with `a` and `b`
its ends, `(:leg, a, u)` with `a` its vertex and `u` its direction, or
`(:higher,)` when its dimension is at least 2, the points normalized. Only the
subsets of deficient rank reach polymake, under `POLYMAKE_LOCK`; they are few.
"""
function exact_piece(M, A, AM, t::Vector{QQFieldElem}, T)
    k, others = T[1], setdiff(1:21, T)
    E = vcat([M[i:i, :] - M[k:k, :] for i in T[2:end]]..., M[k:k, :])
    f = matrix(QQ, 5, 1, [[t[i] - t[k] for i in T[2:end]]; t[k]])
    solvable, y0, K = can_solve_with_solution_and_kernel(E, f; side=:right)
    solvable || return nothing
    G = vcat([M[k:k, :] - M[j:j, :] for j in others]..., AM)
    h = matrix(QQ, nrows(G), 1, [[t[k] - t[j] for j in others]; zeros(QQ, nrows(AM))])
    point(y) = normalized((M * y)[:, 1])
    if ncols(K) > 1
        return lock(POLYMAKE_LOCK) do
            P = polyhedron((G, h[:, 1]), (E, f[:, 1]))
            is_feasible(P) && dim(P) >= 1 || return nothing
            dim(P) >= 2 && return (:higher,)
            ends = [point(matrix(QQ, 6, 1, collect(v))) for v in vertices(P)]
            is_bounded(P) && return (:segment, ends[1], ends[2])
            return (:leg, only(ends), direction((M * matrix(QQ, 6, 1, collect(only(rays(P)))))[:, 1]))
        end
    end
    g, c = G * K, h - G * y0
    lo, hi = nothing, nothing
    for r in 1:nrows(G)
        if iszero(g[r, 1])
            c[r, 1] < 0 && return nothing
        elseif g[r, 1] > 0
            hi = isnothing(hi) ? c[r, 1] // g[r, 1] : min(hi, c[r, 1] // g[r, 1])
        else
            lo = isnothing(lo) ? c[r, 1] // g[r, 1] : max(lo, c[r, 1] // g[r, 1])
        end
    end
    @req !isnothing(lo) || !isnothing(hi) "a piece is a line"
    isnothing(lo) && return (:leg, point(y0 + hi * K), direction(-(M * K)[:, 1]))
    isnothing(hi) && return (:leg, point(y0 + lo * K), direction((M * K)[:, 1]))
    lo < hi || return nothing
    return (:segment, point(y0 + lo * K), point(y0 + hi * K))
end

"""
For every translation of `shifts`, the subsets of `subsets` whose piece of the
cell of span `Mf` and facets `AMf`, in floating point, may be more than a
point, with a margin relative to the size of the translation.
"""
function candidate_subsets(Mf, AMf, subsets, shifts)
    found = [Vector{Int}[] for _ in shifts]
    margins = [FLOAT_MARGIN * (1 + maximum(abs, t)) for t in shifts]
    for T in subsets
        k = T[1]
        E = vcat(Mf[T[2:end], :] .- Mf[k:k, :], Mf[k:k, :])
        F = LinearAlgebra.svd(E; full=true)
        r = count(>(1e-9 * F.S[1]), F.S)
        P = F.V[:, 1:r] * LinearAlgebra.Diagonal(1 ./ F.S[1:r]) * F.U[:, 1:r]'
        xd, ad = Mf * F.V[:, 6], AMf * F.V[:, 6]
        for (n, t) in enumerate(shifts)
            f = [t[T[2:end]] .- t[k]; t[k]]
            y0 = P * f
            if r < 5
                LinearAlgebra.norm(E * y0 - f) <= margins[n] && push!(found[n], T)
                continue
            end
            x0, a0 = Mf * y0, AMf * y0
            lo, hi, possible = -Inf, Inf, true
            # The constraint c + s g ≥ 0 on the parameter s of the line, relaxed by the margin.
            function constrain!(c, g)
                c += margins[n]
                if abs(g) < 1e-12
                    possible &= c >= 0
                elseif g > 0
                    lo = max(lo, -c / g)
                else
                    hi = min(hi, -c / g)
                end
            end
            for j in 1:21
                j in T || constrain!((x0[j] - t[j]) - (x0[k] - t[k]), xd[j] - xd[k])
            end
            foreach(i -> constrain!(-a0[i], -ad[i]), eachindex(a0))
            possible && lo <= hi && push!(found[n], T)
        end
    end
    return found
end

"""
The candidate subsets `found[c][n]` of `candidate_subsets` of every cell `c` of
`cells`, pairs of a span and facets, for every translation `shifts[n]`,
threaded over the cells; `approximations` are the cells in floating point, when
computed before.
"""
function cell_candidates(cells, shifts::Vector{Vector{QQFieldElem}}; approximations=nothing)
    subsets = collect(Combinatorics.combinations(1:21, 5))
    floats = [Float64.(t) for t in shifts]
    approximations = something(approximations, [(Float64.(Matrix(M)), Float64.(Matrix(A * M))) for (M, A) in cells])
    found = Vector{Vector{Vector{Vector{Int}}}}(undef, length(cells))
    Threads.@threads for c in eachindex(cells)
        found[c] = candidate_subsets(approximations[c]..., subsets, floats)
    end
    return found
end

"A basis of the lattice `span(M) ∩ ZZ^21`, as columns."
function saturated_lattice(M)
    N = Matrix(kernel(transpose(M); side=:right))
    return kernel(transpose(matrix(ZZ, ZZ.(N .* lcm(denominator.(N))))); side=:right)
end

"`[ZZ^21 : L + L_S]` with `L_S = {x_S constant}`: the index of the image of `L` under `x -> (x_s - x_S[1])`."
function bergman_lattice_index(L, S)
    D = matrix(ZZ, [j == S[i+1] ? 1 : j == S[1] ? -1 : 0 for i in 1:length(S)-1, j in 1:21])
    B = D * L
    @req rank(B) == length(S) - 1 "the cell is not transverse to B_S"
    s = snf(B)
    return Int(prod(s[i, i] for i in 1:length(S)-1))
end

"""
The pieces of the curve `C ∩ (B + t)` over the cells `cells`, pairs of a span
and facets, from the candidate subsets `found[c]` of every cell `c`, threaded
over the cells, by cell and then in the order of `found[c]`: named tuples of
the cell, the subset `T` and the piece of `exact_piece`, never empty or a
point, and with `weights`, the multiplicities of the cells, the multiplicity
`weights[c] [ZZ^21 : L_C + L_T]` of the piece, 0 for a piece of dimension at
least 2.
"""
function bergman_curve_pieces(cells, found, t::Vector{QQFieldElem}; weights=nothing)
    Piece = @NamedTuple{cell::Int, subset::Vector{Int}, piece::Tuple, multiplicity::Int}
    pieces = [Piece[] for _ in cells]
    Threads.@threads for c in eachindex(cells)
        M, A = cells[c]
        AM = A * M
        L = isnothing(weights) || isempty(found[c]) ? nothing : saturated_lattice(M)
        for T in found[c]
            piece = exact_piece(M, A, AM, t, T)
            isnothing(piece) && continue
            m = isnothing(weights) || first(piece) == :higher ? 0 : weights[c] * bergman_lattice_index(L, T)
            push!(pieces[c], (cell=c, subset=T, piece=piece, multiplicity=m))
        end
    end
    return reduce(vcat, pieces; init=Piece[])
end
