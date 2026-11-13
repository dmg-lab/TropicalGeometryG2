###############################################################################
# The valences of a genus-10 curve inside the cells of the orbits 16, 17, 19
# and 22 of the coarsest tropical G2-Grassmannian, whose part C ∩ (B + t) is a
# polytope for every generic translation t (shifted_bergman_cells.jl), for
# several kinds of translations:
#
#   uniform        uniformly random, as in run_experiments.jl;
#   normal+fixed   a fixed part in the span W of the fan plus a random part
#                  orthogonal to W, the same fixed part for every trial;
#   normal         orthogonal to W only;
#   in W           in W only;
#   cell normal    orthogonal to the span of the cell only.
#
# Exact: every piece C ∩ (B_T + t), |T| = 5, by `exact_piece` of
# bergman_curve.jl. The vertices of the curve inside C lie on walls of B + t;
# those on the boundary of C lie on walls of the G2-Grassmannian. A piece of
# dimension at least 2 means that the translation is not transverse to the
# cell.
#
#   julia --project=. studies/g2/curves/bounded_cell_valences.jl [--trials=N]
#   julia --project=. studies/g2/curves/bounded_cell_valences.jl --check=ARCHIVE,K,CELL
#
# The second form computes the part in CELL of the curve of experiment K of
# ARCHIVE, for comparison with decorated_curve_graph.jl.
###############################################################################

include(joinpath(@__DIR__, "bergman_curve.jl"))
import Random

const SHIFT_KINDS = ("uniform", "normal+fixed", "normal", "in W", "cell normal")

"""
The part of `C ∩ (B + t)`: its segments, as pairs of normalized points, the
number of unbounded pieces and the number of pieces of dimension at least 2.
"""
function cell_curve(M, A, t::Vector{QQFieldElem})
    pieces = [p.piece for p in bergman_curve_pieces([(M, A)], [collect(Combinatorics.combinations(1:21, 5))], t)]
    segments = Tuple{Vector{QQFieldElem},Vector{QQFieldElem}}[(p[2], p[3]) for p in pieces if first(p) == :segment]
    return segments, count(p -> first(p) == :leg, pieces), count(p -> first(p) == :higher, pieces)
end

"The valences of the vertices of `segments` inside the cell of facets `A`, and the number of vertices on its boundary."
function cell_valences(segments, A)
    degree = Dict{Vector{QQFieldElem},Int}()
    foreach(s -> foreach(v -> degree[v] = get(degree, v, 0) + 1, s), segments)
    inside = [d for (v, d) in degree if all(<(0), A * matrix(QQ, 21, 1, v))]
    return inside, length(degree) - length(inside)
end

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
        inside, boundary = cell_valences(segments, A)
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
            inside, boundary = cell_valences(segments, A)
            println("| $orbit | $cell | $kind | $(length(segments)) | $(histogram(inside)) | $boundary | $unbounded | $higher |")
            flush(stdout)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
