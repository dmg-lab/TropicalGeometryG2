###############################################################################
# The genus-10 curve of Method 1 of Section 10 of main.tex, trop(G2/P2) ∩
# (B + q) with B the uniform Bergman fan fan_uniform_linear_space_16_21.mrdi,
# exactly from the shift q alone, in the slice sum(x) = sum(q). Its pieces are
# those of the certified run of run_experiments.jl that found it
# (20 to 35 min on 24 cores of hydra), built by the same exact intersections in
# the same order, without the empty ones. The pairs of a maximal cell C of the raw fan
# and a maximal cone B_T of B, T the 5 coordinates where x - q is minimal, come
# from `bergman_curve_pieces` of bergman_curve.jl, cell by cell, instead of an
# LP for each of the 251 × 20349 pairs. The raw fan, not
# the coarsest: the curve crosses the diagonal of a quadrilateral of the
# coarsest fan, a wall of the raw fan, 5 times, and the stored pieces end there.
#
# The floating-point filter of `candidate_subsets` certifies nothing, so the
# curve is checked balanced, every piece with its multiplicity, and of degree
# 18 in every coordinate. The curve is balanced
# of degree 18, so a missing part would be balanced, hence have rays, and
# lower the degree.
#
#   julia --project=. -t 8 studies/g2/curves/method_1_curve.jl [SHIFT.mrdi [PIECES.mrdi]]
#
# SHIFT defaults to curves_g10/overnight_20261028/wpart_2722/translation_1.mrdi;
# the pieces are saved to PIECES when given.
###############################################################################

include(joinpath(@__DIR__, "bergman_curve.jl"))

module IntersectionLib
include(joinpath(@__DIR__, "intersection_lib.jl"))
end

const METHOD_1_SHIFT = project_path("curves_g10", "overnight_20261028", "wpart_2722", "translation_1.mrdi")
const METHOD_1_FILES = (
    raw_fan=data_path("g2", "grassmannian", "01_raw_fan.mrdi"),
    variety_complex=data_path("g2", "grassmannian", "01_variety_complex.mrdi"),
    bergman_fan=data_path("fan_uniform_linear_space_16_21.mrdi"),
    weights=study_path(:g2, "weights", "weights.csv"),
)
"The degree of G2/P2 in its Plücker embedding, and of its curve sections."
const G2_DEGREE = 18

"""
The multiplicities over QQbar of the maximal cones of the raw fan, in its order (weights.csv, kept in git). Each row's
weight is checked to lie in its cone, so that a raw fan computed again in another order cannot pass unnoticed.
"""
function raw_cell_weights(cones)
    rows = csv_rows(METHOD_1_FILES.weights)
    @req [parse(Int, row["cone"]) for row in rows] == eachindex(rows) "weights.csv does not list the raw cones in order"
    @req length(rows) == length(cones) "weights.csv has $(length(rows)) rows for $(length(cones)) raw cones"
    for (row, C) in zip(rows, cones)
        w = QQ.(parse.(Int, split(row["weight"], ';')))
        @req w in C "the weight of raw cone $(row["cone"]) in weights.csv is not in that cone of 01_raw_fan.mrdi"
    end
    return [parse(Int, row["mult_QQbar"]) for row in rows]
end

"""
The index of the maximal cone `B_T = {x : x_T constant and minimal}` of the
Bergman fan `bergman`, for every subset `T` of 5 coordinates.
"""
function bergman_cones(bergman)
    units = map(rays_modulo_lineality(bergman).rays_modulo_lineality) do r
        above = findall(>(minimum(r)), collect(r))
        @req length(above) == 1 "a ray of the Bergman fan is not some e_j"
        only(above)
    end
    incidence = maximal_cones(IncidenceMatrix, bergman)
    cone_of = Dict{Vector{Int},Int}()
    for i in 1:nrows(incidence)
        T = setdiff(1:21, units[collect(row(incidence, i))])
        @req length(T) == 5 && !haskey(cone_of, T) "the maximal cones of the Bergman fan are not the B_T, |T| = 5"
        cone_of[T] = i
    end
    return cone_of
end

"""
The curve `trop(G2/P2) ∩ (B + q)` in the slice `sum(x) = sum(q)`: its
one-dimensional `pieces`, the intersections of a maximal cone `B_T + q` with a
maximal cell of `01_variety_complex.mrdi`, in the order of the `pairs` of
their indices, with their `multiplicities`, and the graph of
`build_intersection_graph` of intersection_lib.jl on their segments. Errors
unless the curve is transverse, balanced and of degree 18 in every coordinate.
"""
function method_1_pieces(q::Vector{Int})
    @req length(q) == 21 "the shift is not in QQ^21"
    raw = load(METHOD_1_FILES.raw_fan)
    complex = load(METHOD_1_FILES.variety_complex)
    bergman = load(METHOD_1_FILES.bergman_fan)
    cones, polyhedra = collect(maximal_cones(raw)), collect(maximal_polyhedra(complex))
    same_cells = length(cones) == length(polyhedra) && all(((C, P),) -> polyhedron(C) == P, zip(cones, polyhedra))
    @req same_cells "the raw fan and the variety complex have other cells"
    weights = raw_cell_weights(cones)
    @req length(weights) == length(cones) "one weight per raw cone"
    cone_of = bergman_cones(bergman)
    cells = [cell_span(C) for C in cones]
    t = Vector{QQFieldElem}(QQ.(q))
    candidates = [only(f) for f in cell_candidates(cells, [t])]
    pieces = bergman_curve_pieces(cells, candidates, t; weights)
    for p in pieces
        first(p.piece) == :higher && error("a piece of dimension at least 2 in raw cell $(p.cell), T = $(p.subset): q is not transverse")
    end
    # (cone of B, raw cell, multiplicity) of the pieces of dimension 1, in the order of `hybrid_prefilter`: by cone of
    # B, then by cell.
    found = sort!([(cone_of[p.subset], p.cell, p.multiplicity) for p in pieces])
    pairs = [(i, c) for (i, c, _) in found]
    @req allunique(pairs) "a pair of a cone and a cell twice"

    # The exact intersections of `run_experiment` of intersection_lib.jl.
    bergman_cells = maximal_cones(bergman)
    cones_of_pairs = Dict(i => polyhedron(bergman_cells[i]) for i in unique(first.(pairs)))
    result = IntersectionLib.shifted_intersections(pairs, cones_of_pairs, polyhedra, q)
    @req all(P -> dim(P) == 1, result.list_P) "an exact intersection is not one-dimensional"

    multiplicities = [m for (_, _, m) in found]
    curve = balancing_and_degree(result.list_P, multiplicities)
    @req curve.balanced "the curve is not balanced: a piece is missing"
    @req all(==(G2_DEGREE), curve.degree) "the curve has degree $(curve.degree), not $G2_DEGREE in every coordinate: a piece is missing"
    return (pieces=result.list_P, multiplicities=multiplicities, pairs=pairs, graph=result.g,
            candidates=sum(length, candidates))
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 2 || error("usage: method_1_curve.jl [SHIFT.mrdi [PIECES.mrdi]]")
    shift_file = isempty(ARGS) ? METHOD_1_SHIFT : abspath(ARGS[1])
    seconds = @elapsed curve = method_1_pieces(Vector{Int}(load(shift_file)))
    segments = count(is_bounded, curve.pieces)
    println("shift $shift_file, $(Threads.nthreads()) threads, $(round(seconds; digits=1)) s")
    println("$(curve.candidates) candidate pairs of a raw cell and a subset T, $(length(curve.pieces)) pieces: ",
            "$segments segments and $(length(curve.pieces) - segments) rays, multiplicities $(histogram(curve.multiplicities))")
    println("graph of the segments: $(IntersectionLib.nv(curve.graph)) vertices, $(IntersectionLib.ne(curve.graph)) edges; ",
            "balanced, of degree $G2_DEGREE in every coordinate")
    if length(ARGS) == 2
        save(abspath(ARGS[2]), curve.pieces)
        println("wrote $(abspath(ARGS[2]))")
    end
end
