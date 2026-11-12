###############################################################################
# The valences of a genus-10 curve inside the cells of the orbits 16, 17, 19
# and 22 of the coarsest tropical G2-Grassmannian, whose part C ∩ (B + t) is a
# polytope for every generic translation t (shifted_bergman_cells.jl), for
# several kinds of translations:
#
#   uniform        uniformly random, as in mukai_curves/run_experiments.jl;
#   normal+fixed   a fixed part in the span W of the fan plus a random part
#                  orthogonal to W, the same fixed part for every trial;
#   normal         orthogonal to W only;
#   in W           in W only;
#   cell normal    orthogonal to the span of the cell only.
#
# Exact. The piece C ∩ (B_T + t), |T| = 5, lies on the line of span(C)
# modulo the diagonal on which x - t is constant on T, so it is an interval
# of that line, cut by the 16 inequalities of B_T and the facets of C. The
# vertices of the curve inside C lie on walls of B + t; those on the boundary
# of C lie on walls of the G2-Grassmannian. A piece of dimension at least 2
# means that the translation is not transverse to the cell.
#
#   julia --project=julia julia/studies/g2/curves/bounded_cell_valences.jl [--trials=N]
#   julia --project=julia julia/studies/g2/curves/bounded_cell_valences.jl --check=ARCHIVE,K,CELL
#
# The second form computes the part in CELL of the curve of experiment K of
# ARCHIVE, for comparison with decorated_curve_graph.jl.
###############################################################################

using Oscar
import Combinatorics
import Random
using TropicalCensus: data_path

const BOUNDED_ORBITS = (16, 17, 19, 22)
const SHIFT_RANGE = -(1 << 20):(1 << 20)
const SHIFT_KINDS = ("uniform", "normal+fixed", "normal", "in W", "cell normal")

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

"The point `x` modulo the diagonal, with first coordinate 0."
normalized(x) = [x[i] - x[1] for i in eachindex(x)]

"The direction `x` modulo the diagonal, with first coordinate 0 and largest coordinate ±1."
direction(x) = (v = normalized(x); v .// maximum(abs, v))

"""
The piece `C ∩ (B_T + t)` of the cell of span `M` and facets `A`, `AM = A * M`:
`nothing` when it is empty or a point, `(:segment, a, b)` with `a` and `b`
its ends, `(:leg, a, u)` with `a` its vertex and `u` its direction, or
`(:higher,)` when its dimension is at least 2, the points normalized.
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
        P = polyhedron((G, h[:, 1]), (E, f[:, 1]))
        is_feasible(P) && dim(P) >= 1 || return nothing
        dim(P) >= 2 && return (:higher,)
        ends = [point(matrix(QQ, 6, 1, collect(v))) for v in vertices(P)]
        is_bounded(P) && return (:segment, ends[1], ends[2])
        return (:leg, only(ends), direction((M * matrix(QQ, 6, 1, collect(only(rays(P)))))[:, 1]))
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
The part of `C ∩ (B + t)`: its segments, as pairs of normalized points, the
number of unbounded pieces and the number of pieces of dimension at least 2.
"""
function cell_curve(M, A, t::Vector{QQFieldElem})
    segments, unbounded, higher = Tuple{Vector{QQFieldElem},Vector{QQFieldElem}}[], 0, 0
    AM = A * M
    for T in Combinatorics.combinations(1:21, 5)
        piece = exact_piece(M, A, AM, t, T)
        isnothing(piece) && continue
        first(piece) == :segment && push!(segments, (piece[2], piece[3]))
        unbounded += first(piece) == :leg
        higher += first(piece) == :higher
    end
    return segments, unbounded, higher
end

"The valences of the vertices of `segments` inside the cell of facets `A`, and the number of vertices on its boundary."
function valences(segments, A)
    degree = Dict{Vector{QQFieldElem},Int}()
    foreach(s -> foreach(v -> degree[v] = get(degree, v, 0) + 1, s), segments)
    inside = [d for (v, d) in degree if all(<(0), A * matrix(QQ, 21, 1, v))]
    return inside, length(degree) - length(inside)
end

"The histogram `count` of `values` as `value^count`, in increasing order."
histogram(values) = join(("$v^$(count(==(v), values))" for v in sort!(unique(values))), " ")

"The translations of `kind` for the cell of span `M`; `normals` spans the complement of the span of the fan, `fixed` lies in it."
function shift(kind, M, normals, fixed, rng)
    random(n) = matrix(QQ, n, 1, rand(rng, SHIFT_RANGE, n))
    project(v) = v - normals * solve(transpose(normals) * normals, transpose(normals) * v; side=:right)
    t = kind == "uniform" ? random(21) :
        kind == "normal+fixed" ? fixed + normals * random(ncols(normals)) :
        kind == "normal" ? normals * random(ncols(normals)) :
        kind == "in W" ? project(random(21)) :
        kind == "cell normal" ? (K = kernel(transpose(M); side=:right); K * random(ncols(K))) :
        error("unknown kind of translation $kind")
    return t[:, 1]
end

function main(args)
    fan = load(data_path("g2", "grassmannian", "02_coarsest_fan.mrdi"))
    cones = collect(maximal_cones(fan))
    checks = [split(split(a, '=')[2], ',') for a in args if startswith(a, "--check=")]
    if !isempty(checks)
        archive, experiment, cell = only(checks)
        t = Vector{QQFieldElem}(QQ.(load(joinpath(archive, "translation_$(experiment).mrdi"))))
        M, A = cell_span(cones[parse(Int, cell)])
        segments, unbounded, higher = cell_curve(M, A, t)
        inside, boundary = valences(segments, A)
        println("cell $cell: $(length(segments)) segments, $unbounded unbounded, $higher of dimension ≥ 2; ",
                "valences inside: $(histogram(inside)); $boundary points on the boundary")
        return
    end
    trials = something(tryparse(Int, last([["3"]; [split(a, '=')[2] for a in args if startswith(a, "--trials=")]])), 3)
    rng = Random.Xoshiro(2026)
    span = vcat(matrix(QQ, reduce(vcat, permutedims.(collect.(rays_modulo_lineality(fan).rays_modulo_lineality)))),
                matrix(QQ, lineality_space(fan)))
    normals = kernel(span; side=:right)
    @req ncols(normals) == 6 "the fan does not span a 15-dimensional space"
    fixed = shift("in W", nothing, normals, nothing, rng)
    fixed = matrix(QQ, 21, 1, fixed)
    representatives = Dict(parse(Int, split(l, ',')[1]) => parse(Int, split(l, ',')[3])
                           for l in readlines(data_path("g2", "grassmannian", "05_symmetry_cell_orbits.csv"))[2:end])
    println("| orbit | cell | translation | segments | valences inside the cell | points on its boundary | unbounded | dimension ≥ 2 |")
    println("|---:|---:|---|---:|---|---:|---:|---:|")
    for orbit in BOUNDED_ORBITS
        cell = representatives[orbit]
        M, A = cell_span(cones[cell])
        for kind in SHIFT_KINDS, trial in 1:trials
            segments, unbounded, higher = cell_curve(M, A, shift(kind, M, normals, fixed, rng))
            inside, boundary = valences(segments, A)
            println("| $orbit | $cell | $kind | $(length(segments)) | $(histogram(inside)) | $boundary | $unbounded | $higher |")
            flush(stdout)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
