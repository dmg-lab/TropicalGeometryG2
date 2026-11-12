using Oscar
using ProgressMeter
using HiGHS
import LinearAlgebra
using TropicalCensus: coarsen_prevariety_fan, SupportLinkCoarsening

# The feasibility filters of intersection_lib.jl, kept in their own module so
# that its Graphs imports (vertices, degree, ...) do not shadow Oscar's.
module IntersectionLib
include(joinpath(@__DIR__, "..", "..", "..", "intersection_lib.jl"))
end
using .IntersectionLib: lp_prefilter, soplex_prefilter, hybrid_prefilter

"""
Exact H-representation of a cell: facets `F x ≤ g` and affine hull `E x = f`,
each row scaled to integer coefficients.
"""
struct CellHRep
    F :: Matrix{QQFieldElem}
    g :: Vector{QQFieldElem}
    E :: Matrix{QQFieldElem}
    f :: Vector{QQFieldElem}
end

function CellHRep(P::Polyhedron)
    function integral_rows(M)
        M = Matrix{QQFieldElem}(M)
        isempty(M) && return Matrix{QQFieldElem}(undef, 0, ambient_dim(P) + 1)
        return M .* [lcm(denominator.(row)) for row in eachrow(M)]
    end
    facets = integral_rows(Oscar.pm_object(P).FACETS)
    hull = integral_rows(Oscar.pm_object(P).AFFINE_HULL)
    return CellHRep(-facets[:, 2:end], facets[:, 1], hull[:, 2:end], -hull[:, 1])
end

"The rows `A x ≤ b` of a cell, equations as pairs of inequalities, for the feasibility filters."
inequalities(h::CellHRep) = (vcat(h.F, h.E, -h.E), vcat(h.g, h.f, -h.f))

"Rows of a ZZ-basis of the lattice points of the linear space parallel to the cell."
direction_lattice(h::CellHRep) =
    isempty(h.E) ? identity_matrix(ZZ, size(h.E, 2)) : transpose(kernel(matrix(ZZ, h.E); side=:right))

"`[ZZ^n : Λ_1 + Λ_2]`, or 0 when the lattices do not span `QQ^n`."
function lattice_index(lattice_1::ZZMatrix, lattice_2::ZZMatrix)
    divisors = snf(vcat(lattice_1, lattice_2))
    n = ncols(divisors)
    nrows(divisors) >= n && all(i -> !iszero(divisors[i, i]), 1:n) || return ZZ(0)
    return abs(prod(divisors[i, i] for i in 1:n))
end

"""
Exact LP in `(x, t)`: maximize `t ≤ 1` subject to `F x + t ≤ g` and `E x = f`
for both cells. `nothing` when the cells are disjoint, otherwise whether their
relative interiors meet (`t > 0`).
"""
function relative_interiors_meet(h1::CellHRep, h2::CellHRep)
    n = size(h1.E, 2)
    F = vcat(hcat(vcat(h1.F, h2.F), ones(QQFieldElem, size(h1.F, 1) + size(h2.F, 1))),
             permutedims([zeros(QQFieldElem, n); QQ(1)]))
    E = hcat(vcat(h1.E, h2.E), zeros(QQFieldElem, size(h1.E, 1) + size(h2.E, 1)))
    region = polyhedron((F, vcat(h1.g, h2.g, QQ(1))), (E, vcat(h1.f, h2.f)))
    t = optimal_value(linear_program(region, [zeros(Int, n); 1]))
    # t < 0 when the affine hulls meet but the cells do not.
    return isnothing(t) || t < 0 ? nothing : t > 0
end

"""
The nonempty intersection of a maximal cell `σ` of the reference complex with
the maximal cell `other_cell` of the other complex, kept in H-representation.
The cells meet transversely when their relative interiors meet and their
linear spans add up to the whole space; the multiplicity is then
`w_σ w_τ [ZZ^n : Λ_σ + Λ_τ]`, with `Λ` the lattice of the span, and 0 otherwise.
"""
struct CellIntersection
    other_cell   :: Int
    polyhedron   :: Polyhedron{QQFieldElem}
    dim          :: Int
    transverse   :: Bool
    multiplicity :: ZZRingElem
end

"""
Output of `intersect_complexes`. `intersections[i]` lists the maximal cells of
the other complex meeting `reference_cells[i]`; the keys are indices because
Oscar hashes every polyhedron of a given ambient dimension to the same value.
"""
struct ComplexIntersection
    complex         :: Union{PolyhedralComplex{QQFieldElem}, Nothing}
    reference_cells :: Vector{Polyhedron{QQFieldElem}}
    other_cells     :: Vector{Polyhedron{QQFieldElem}}
    intersections   :: Dict{Int, Vector{CellIntersection}}
end

function Base.show(io::IO, result::ComplexIntersection)
    meetings = reduce(vcat, values(result.intersections); init=CellIntersection[])
    print(io, "Intersection of $(length(result.reference_cells)) × $(length(result.other_cells)) maximal cells: ",
          "$(length(meetings)) meeting pairs, $(count(c -> c.transverse, meetings)) transverse, ",
          "total multiplicity $(sum(c -> c.multiplicity, meetings; init=ZZ(0)))")
end

"""
The cones over `cells` in `QQ^{1+n}`, spanned by `(1, v)` for the vertices and
`(0, r)` for the rays, reduced modulo the common lineality `(0, L)` so that
equal generators get equal rows. Returns the generators, the lineality rows
and the generators of each cone.
"""
function homogenized_cells(cells::Vector{Polyhedron{QQFieldElem}})
    n = ambient_dim(first(cells))
    L = lineality_dim(first(cells)) == 0 ? zero_matrix(QQ, 0, n) : matrix(QQ, lineality_space(first(cells)))
    lineality_rank, reduced = rref(L)
    reduced = Matrix{QQFieldElem}(reduced)[1:lineality_rank, :]
    pivots = [findfirst(!iszero, reduced[k, :]) for k in 1:lineality_rank]
    function canonical(v, is_ray)
        v = v - sum((v[pivots[k]] .* reduced[k, :] for k in 1:lineality_rank); init=zeros(QQFieldElem, n))
        is_ray && (v = v ./ abs(v[findfirst(!iszero, v)]))
        return [QQ(is_ray ? 0 : 1); v]
    end
    generators, index, generator_sets = Vector{QQFieldElem}[], Dict{Vector{QQFieldElem},Int}(), Vector{Int}[]
    for P in cells
        cone = [canonical(collect(v), false) for v in minimal_faces(P).base_points]
        append!(cone, canonical(collect(r), true) for r in rays_modulo_lineality(P).rays_modulo_lineality)
        push!(generator_sets, sort!([get!(() -> (push!(generators, g); length(generators)), index, g) for g in cone]))
    end
    return reduce(vcat, permutedims.(generators)), hcat(zeros(QQFieldElem, lineality_rank), reduced), generator_sets
end

"The polyhedral complex whose cells have the cones `generator_sets` over them."
function dehomogenized_complex(generator_sets, generators, lineality; non_redundant::Bool=false)
    points = [iszero(g[1]) ? g[2:end] : g[2:end] ./ g[1] for g in eachrow(generators)]
    return polyhedral_complex(QQ, incidence_matrix(generator_sets), matrix(QQ, reduce(vcat, permutedims.(points))),
                              findall(iszero, generators[:, 1]), matrix(QQ, lineality[:, 2:end]); non_redundant)
end

"The polyhedral complex with the given cells and their common lineality."
function complex_from_cells(cells::Vector{Polyhedron{QQFieldElem}}, n::Int)
    isempty(cells) && return Oscar._empty_complex(QQ, n)
    generators, lineality, generator_sets = homogenized_cells(cells)
    return dehomogenized_complex(generator_sets, generators, lineality)
end

"The maximal cells of `complex`, their H-representations and the lattices of their spans."
function cell_data(complex::PolyhedralComplex{QQFieldElem})
    cells = collect(maximal_polyhedra(complex))
    hreps = CellHRep.(cells)
    return (; cells, hreps, lattices=direction_lattice.(hreps))
end

"""
The nonzero rows `A x ≤ b` scaled to unit normals, in floating point, so that
slacks are distances. Polymake lists the far face `0 ≤ 1` among the facets of
an unbounded cell.
"""
function unit_rows(A::Matrix{QQFieldElem}, b::Vector{QQFieldElem})
    nonzero = [!all(iszero, row) for row in eachrow(A)]
    @req all(>=(0), b[.!nonzero]) "a row 0 ≤ b with b < 0: the cell is empty"
    A, b = Float64.(A[nonzero, :]), b[nonzero]
    norms = [LinearAlgebra.norm(row) for row in eachrow(A)]
    return A ./ norms, Float64.(b) ./ norms
end

"""
A reference cell `{p + B y : A y ≤ b}` in floating point for `slack_prefilter`:
`p` a point of its relative interior, `B` a basis of the directions of its
affine hull, modulo `(1,...,1)` when `modulo_ones`, and the facets in unit rows.
"""
struct SlackCell
    p :: Vector{Float64}
    B :: Matrix{Float64}
    A :: Matrix{Float64}
    b :: Vector{Float64}
end

function SlackCell(P::Polyhedron{QQFieldElem}, h::CellHRep, modulo_ones::Bool)
    n = size(h.E, 2)
    p = Vector{QQFieldElem}(relative_interior_point(P))
    E = modulo_ones ? vcat(h.E, permutedims(QQ.([1; zeros(Int, n - 1)]))) : h.E
    B = isempty(E) ? Matrix{QQFieldElem}(identity_matrix(QQ, n)) : Matrix{QQFieldElem}(kernel(matrix(QQ, E); side=:right))
    # The slacks of p, computed exactly before rounding.
    F, b = unit_rows(h.F, h.g - h.F * p)
    return SlackCell(Float64.(p), Float64.(B), F * Float64.(B), b)
end

"""
    slack_prefilter(reference_cells, reference_hreps, other_hreps; tolerance=1e-7, cap=1.0)

Floating-point filter: for every pair of a reference cell and a cell of the
other complex, HiGHS maximizes the common slack `t ≤ cap` of the facets of
both, in unit rows, within the affine hulls, written in the coordinates `y` of
the reference cell. The relative interiors meet when `t > 0`. Pairs with
`t > tolerance` are returned as meeting, pairs with `|t| ≤ tolerance`, which
touch on their boundaries, are dropped, and those HiGHS does not solve to
optimality are returned as undecided, for an exact LP. Nothing certifies the
pairs dropped. Threaded over the reference cells.
"""
function slack_prefilter(reference_cells::Vector{Polyhedron{QQFieldElem}}, reference_hreps::Vector{CellHRep},
                         other_hreps::Vector{CellHRep}; tolerance::Float64=1e-7, cap::Float64=1.0)
    n = isempty(reference_hreps) ? 0 : size(first(reference_hreps).E, 2)
    ones_n = ones(QQFieldElem, n)
    contains_ones(h) = all(iszero, h.F * ones_n) && all(iszero, h.E * ones_n)
    modulo_ones = all(contains_ones, reference_hreps) && all(contains_ones, other_hreps)
    # Polymake is not thread-safe: the relative interior points are computed here.
    slack_cells = [SlackCell(P, h, modulo_ones) for (P, h) in zip(reference_cells, reference_hreps)]
    others = [(unit_rows(h.F, h.g)..., unit_rows(h.E, h.f)...) for h in other_hreps]

    meeting = [Tuple{Int,Int}[] for _ in reference_cells]
    undecided = [Tuple{Int,Int}[] for _ in reference_cells]
    # Pairs by slack: below -tolerance, touching, up to 1e-4, above.
    bands = zeros(Int, 4, length(reference_cells))
    models = [HiGHS.Highs_create() for _ in 1:Threads.maxthreadid()]
    try
        for h in models
            HiGHS.Highs_setBoolOptionValue(h, "output_flag", false)
            HiGHS.Highs_setIntOptionValue(h, "threads", 1)
        end
        Threads.@threads :static for i in eachindex(slack_cells)
            h, c = models[Threads.threadid()], slack_cells[i]
            d = size(c.B, 2)
            HiGHS.Highs_clearModel(h)
            HiGHS.Highs_changeObjectiveSense(h, HiGHS.kHighsObjSenseMaximize)
            for _ in 1:d
                HiGHS.Highs_addCol(h, 0.0, -Inf, Inf, 0, C_NULL, C_NULL)
            end
            HiGHS.Highs_addCol(h, 1.0, -Inf, cap, 0, C_NULL, C_NULL)
            columns = Cint.(0:d)
            add_row!(lower, a, upper) = HiGHS.Highs_addRow(h, lower, upper, d + 1, columns, a)
            foreach(k -> add_row!(-Inf, [c.A[k, :]; 1.0], c.b[k]), axes(c.A, 1))
            n_base = size(c.A, 1)
            for (j, (F, g, E, f)) in enumerate(others)
                FB, EB = F * c.B, E * c.B
                values = f - E * c.p
                # An equation constant on the reference cell: the affine hulls are parallel.
                constant = [LinearAlgebra.norm(EB[k, :]) <= 1e-9 for k in axes(E, 1)]
                if any(constant)
                    any(k -> constant[k] && abs(values[k]) > 1e-6, axes(E, 1)) ? (bands[1, i] += 1) : push!(undecided[i], (i, j))
                    continue
                end
                for k in axes(F, 1)
                    add_row!(-Inf, [FB[k, :]; 1.0], g[k] - LinearAlgebra.dot(F[k, :], c.p))
                end
                foreach(k -> add_row!(values[k], [EB[k, :]; 0.0], values[k]), axes(E, 1))
                HiGHS.Highs_run(h)
                if HiGHS.Highs_getModelStatus(h) == HiGHS.kHighsModelStatusOptimal
                    t = HiGHS.Highs_getObjectiveValue(h)
                    bands[t < -tolerance ? 1 : t <= tolerance ? 2 : t <= 1e-4 ? 3 : 4, i] += 1
                    t > tolerance && push!(meeting[i], (i, j))
                elseif HiGHS.Highs_getModelStatus(h) == HiGHS.kHighsModelStatusInfeasible
                    # The affine hulls do not meet.
                    bands[1, i] += 1
                else
                    push!(undecided[i], (i, j))
                end
                HiGHS.Highs_deleteRowsByRange(h, Cint(n_base), Cint(n_base + size(F, 1) + size(E, 1) - 1))
            end
        end
    finally
        foreach(HiGHS.Highs_destroy, models)
    end
    totals = vec(sum(bands; dims=2))
    println("    slack filter: $(totals[1]) apart, $(totals[2]) touching, $(totals[3]) meeting with slack below 1e-4, ",
            "$(totals[4]) above, $(sum(length, undecided)) undecided")
    return reduce(vcat, meeting; init=Tuple{Int,Int}[]), reduce(vcat, undecided; init=Tuple{Int,Int}[])
end

"""
    intersect_complexes(reference, other; feasibility_solver=:highs, build_complex=true,
                        reference_weights, other_weights, reference_data=cell_data(reference))

Intersect the maximal cells of two polyhedral complexes, working with
H-representations throughout. The candidate pairs come from an LP filter,
threaded over the reference cells (start Julia with `-t`): `:highs` solves in
floating point, `:soplex` and `:hybrid` (HiGHS answers certified over QQ,
SoPlex otherwise) are exact, and `:certified` keeps the pairs HiGHS cannot
certify as candidates instead of asking SoPlex, so that it misses no pair.
Each candidate is then settled by one exact LP, which discards false positives. `:float`, the
`slack_prefilter`, decides in floating point which relative interiors meet,
with no exact LP but for the pairs HiGHS leaves undecided, and keeps no
touching pair: it is fast but certifies nothing. Only `build_complex` needs the
vertices of the intersections. Weights default to 1. Passing `reference_data`
saves recomputing it when the same reference is cut many times.
"""
function intersect_complexes(reference::PolyhedralComplex{QQFieldElem}, other::PolyhedralComplex{QQFieldElem};
                             feasibility_solver::Symbol=:highs, build_complex::Bool=true,
                             reference_weights::AbstractVector=ones(Int, n_maximal_polyhedra(reference)),
                             other_weights::AbstractVector=ones(Int, n_maximal_polyhedra(other)),
                             reference_data=nothing)
    n = ambient_dim(reference)
    @req ambient_dim(other) == n "the complexes live in different ambient spaces"
    @req feasibility_solver in (:highs, :soplex, :hybrid, :certified, :float) "feasibility_solver must be :highs, :soplex, :hybrid, :certified or :float"
    data = isnothing(reference_data) ? cell_data(reference) : reference_data
    reference_cells, reference_hreps, reference_lattices = data
    other_cells = collect(maximal_polyhedra(other))
    @req length(reference_weights) == length(reference_cells) "one weight per maximal cell of the reference"
    @req length(other_weights) == length(other_cells) "one weight per maximal cell of the other complex"
    other_hreps = CellHRep.(other_cells)

    reference_rows, other_rows = inequalities.(reference_hreps), inequalities.(other_hreps)
    meeting = Set{Tuple{Int,Int}}()
    pairs = if feasibility_solver === :float
        found, undecided = slack_prefilter(reference_cells, reference_hreps, other_hreps)
        union!(meeting, found)
        [found; undecided]
    elseif feasibility_solver === :highs
        float_rows(rows) = [(Float64.(A), Float64.(b)) for (A, b) in rows]
        lp_prefilter(float_rows(reference_rows), float_rows(other_rows), zeros(0, n), Float64[])
    elseif feasibility_solver === :soplex
        soplex_prefilter(reference_rows, other_rows, Matrix{QQFieldElem}(undef, 0, n), QQFieldElem[])
    else
        hybrid_prefilter(reference_rows, other_rows, Matrix{QQFieldElem}(undef, 0, n), QQFieldElem[];
                         soplex_fallback=feasibility_solver === :hybrid)
    end
    sort!(pairs)

    other_lattices = direction_lattice.(other_hreps)
    intersections = Dict(i => CellIntersection[] for i in eachindex(reference_cells))
    @showprogress for (i, j) in pairs
        h1, h2 = reference_hreps[i], other_hreps[j]
        meet = (i, j) in meeting || relative_interiors_meet(h1, h2)
        isnothing(meet) && continue
        P = polyhedron((vcat(h1.F, h2.F), vcat(h1.g, h2.g)), (vcat(h1.E, h2.E), vcat(h1.f, h2.f)))
        # Near a common relative interior point the intersection is the intersection of the affine hulls.
        d = meet ? n - rank(matrix(QQ, vcat(h1.E, h2.E))) : dim(P)
        index = meet ? lattice_index(reference_lattices[i], other_lattices[j]) : ZZ(0)
        multiplicity = reference_weights[i] * other_weights[j] * index
        push!(intersections[i], CellIntersection(j, P, d, !iszero(index), multiplicity))
    end

    cells = [c.polyhedron for i in sort!(collect(keys(intersections))) for c in intersections[i]]
    complex = build_complex ? complex_from_cells(cells, n) : nothing
    return ComplexIntersection(complex, reference_cells, other_cells, intersections)
end

"A cell of the coarsened intersection: the union of the pieces `σ_i ∩ τ_j`, listed as `(i, j)`."
struct CoarseCell
    polyhedron   :: Polyhedron{QQFieldElem}
    pieces       :: Vector{Tuple{Int,Int}}
    dim          :: Int
    transverse   :: Bool
    multiplicity :: ZZRingElem
end

"""
Output of `coarsen_intersection`. `by_reference[i]` lists the coarse cells
containing a piece of reference cell `i`; `coarsest_certified` is that of
`coarsen_prevariety_fan`.
"""
struct CoarseIntersection
    complex            :: PolyhedralComplex{QQFieldElem}
    cells              :: Vector{CoarseCell}
    by_reference       :: Dict{Int, Vector{Int}}
    coarsest_certified :: Bool
end

function Base.show(io::IO, coarse::CoarseIntersection)
    pieces = sum(c -> length(c.pieces), coarse.cells; init=0)
    print(io, "Coarsened intersection: $pieces pieces glued into $(length(coarse.cells)) cells, ",
          "coarsest certified: $(coarse.coarsest_certified)")
end

"""
    coarsen_intersection(result; strategy=SupportLinkCoarsening())

The coarsest polyhedral structure on the intersection with constant
multiplicities. Pieces are glued across the walls of both complexes when they
span the same affine space, have the same multiplicity and transversality, and
their union is convex and meets every other cell in a face. The exact
prevariety coarsener of TropicalCensus runs on the cone over the intersection.
Pieces that are faces of other pieces are dropped.
"""
function coarsen_intersection(result::ComplexIntersection; strategy=SupportLinkCoarsening())
    pieces = [(i, c) for i in sort!(collect(keys(result.intersections))) for c in result.intersections[i]]
    by_reference = Dict(i => Int[] for i in eachindex(result.reference_cells))
    isempty(pieces) && return CoarseIntersection(Oscar._empty_complex(QQ, ambient_dim(first(result.reference_cells))),
                                                 CoarseCell[], by_reference, true)
    generators, lineality, generator_sets = homogenized_cells([c.polyhedron for (_, c) in pieces])

    # A piece is a face of another exactly when its generators are among the other's.
    owners = Dict{Int, Vector{Int}}()
    for (k, generator_set) in enumerate(generator_sets), g in generator_set
        push!(get!(owners, g, Int[]), k)
    end
    is_face(k) = any(l -> l != k && issubset(generator_sets[k], generator_sets[l]) &&
                          (generator_sets[k] != generator_sets[l] || l < k), owners[first(generator_sets[k])])
    kept = filter(!is_face, eachindex(pieces))
    fan = polyhedral_fan(IncidenceMatrix(generator_sets[kept]), matrix(QQ, generators), matrix(QQ, lineality);
                         non_redundant=true)
    labels = [(c.multiplicity, c.transverse) for (_, c) in pieces[kept]]
    coarse = coarsen_prevariety_fan(fan; strategy, labels)

    incidence = maximal_cones(IncidenceMatrix, coarse.fan)
    coarse_generators = Matrix{QQFieldElem}(matrix(QQ, rays_modulo_lineality(coarse.fan).rays_modulo_lineality))
    complex = dehomogenized_complex([collect(row(incidence, k)) for k in 1:nrows(incidence)],
                                    coarse_generators, lineality; non_redundant=true)
    cells = CoarseCell[]
    for (k, P) in enumerate(maximal_polyhedra(complex))
        members = pieces[kept[coarse.original_cells[k]]]
        c = last(first(members))
        push!(cells, CoarseCell(P, [(i, piece.other_cell) for (i, piece) in members], c.dim, c.transverse, c.multiplicity))
        foreach(i -> push!(by_reference[i], k), unique(first.(members)))
    end
    return CoarseIntersection(complex, cells, by_reference, coarse.coarsest_certified)
end

"One row per coarse cell: dimension, transverse, multiplicity, and its pieces `i:j` separated by `;`."
function write_coarse_table(path, coarse::CoarseIntersection)
    open(path, "w") do io
        println(io, "cell,dim,transverse,multiplicity,pieces")
        for (k, c) in enumerate(coarse.cells)
            println(io, "$k,$(c.dim),$(c.transverse),$(c.multiplicity),", join(("$i:$j" for (i, j) in c.pieces), ';'))
        end
    end
end

"One row per meeting pair: reference cell, other cell, dimension, transverse, multiplicity."
function write_intersection_table(path, result::ComplexIntersection)
    open(path, "w") do io
        println(io, "reference_cell,other_cell,dim,transverse,multiplicity")
        for i in sort!(collect(keys(result.intersections))), c in result.intersections[i]
            println(io, "$i,$(c.other_cell),$(c.dim),$(c.transverse),$(c.multiplicity)")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (3, 4) || error("usage: intersect_complexes.jl REFERENCE.mrdi OTHER.mrdi OUTPUT_PREFIX [highs|soplex|hybrid|certified]")
    solver = length(ARGS) == 4 ? Symbol(ARGS[4]) : :highs
    result = intersect_complexes(load(ARGS[1]), load(ARGS[2]); feasibility_solver=solver)
    println(result)
    save(ARGS[3] * ".mrdi", result.complex)
    write_intersection_table(ARGS[3] * ".csv", result)
    coarse = coarsen_intersection(result)
    println(coarse)
    save(ARGS[3] * "_coarse.mrdi", coarse.complex)
    write_coarse_table(ARGS[3] * "_coarse.csv", coarse)
end
