###############################################################################
# hyperplane_intersection_census: curve sections of the tropical
# G2-Grassmannian by random shifted tropical hyperplanes.
#
# Experiment k draws 5 random shifted tropical hyperplanes of
# tropical_hyperplanes.jl with the random generator Xoshiro(k). Starting from
# the coarsest fan of data/g2/grassmannian, with the multiplicities of its
# maximal cells over QQbar (weights/weights.csv), it cuts the current complex
# by the next hyperplane (intersect_complexes) and coarsens the cut
# (coarsen_intersection) before cutting again. When the coarsening fails, the
# transverse pieces are kept as they are. After 4 cuts the complex is a
# tropical curve, of dimension 1 modulo its lineality (1,...,1). Its graph is
# decorated by the maximal cells of the G2-Grassmannian as in
# decorated_curve_graph.jl. The target is a curve of genus 10, the genus of
# the curve sections of G2/P2, whose core is trivalent. The fifth cut leaves
# points, 18 counted with multiplicity, the degree of G2/P2; the 18 points of
# a codimension-5 linear section of G2/P2 form a self-dual configuration.
#
#   julia --project=. -t 16 studies/g2/curves/hyperplane_intersection_census.jl \
#     [--first=K] [--experiments=N] [--shift=T] [--coordinates=I,J,...] [--support-sizes=A:B]
#     [--solver=certified|highs|hybrid|soplex|float] [--output=DIR]
#     [--from=K:S] [--candidates=C] [--point-candidates=P] [--shift-basis=default|FILE]
#
# Steering (--from, --candidates above 1 or --shift-basis): experiment k starts from stage S
# of experiment K of the same census, its files copied, and every later cut s
# is the best of C hyperplanes (P for the fifth, by default C) drawn with
# Xoshiro(1000k + s), the first curve on target ending the search. With
# --shift-basis their shifts are combinations of the rows of the basis of
# shift_basis.jl (default: data/mukai_curves/genus_10_shift_basis.mrdi) with
# coefficients in -T:T, by default 2^15. A candidate is discarded when the
# part of G2/P2 where all coordinates of the supports so far vanish has
# dimension at least 5 - s: it lies in the section, outside the torus. Before
# the curve the best cut has the smallest share of cells of multiplicity
# above 1, at the curve the target or else the largest genus, and at the end
# the most points. The candidates are stored in cut_s/, and steering.csv has a
# row for each.
#
# The supports of the first 4 hyperplanes have sizes drawn from A:B, by default
# 2:21. When a torus-invariant line of G2/P2, which has 2 nonzero coordinates,
# lies in all of them, it is a component of the section outside the torus, and
# the curve has degree below 18 (experiment 5: degree 17). The target is a
# balanced curve of degree 18 and genus 10 whose core is trivalent.
#
# With --coordinates the supports of the hyperplanes are drawn among these
# coordinates only. On 1,2,3,6,9,14,19,20,21, whose tie forms x_i - x_j are
# unimodular on every cell, the curve has multiplicity 1, but experiments 1
# and 2 gave a tree with 8 legs: the rest of the section lies outside the
# torus. The fifth hyperplane has full support among the coordinates, since a
# hyperplane misses the points where all its coordinates vanish.
#
# The default filter, certified, drops a pair of cells only on an exact Farkas
# certificate from HiGHS; every other pair is settled by an exact LP. HiGHS
# alone missed transverse pairs (80 of the 5420 of the first cut of experiment
# 1), even after a second solve without presolve (experiment 3), and left the
# curves unbalanced. The SoPlex fallback of hybrid and soplex crashed with a
# segmentation fault on experiment 1 (2026-10-27), even on one thread. The
# census checks that every curve is balanced. The floating-point filter float
# certifies nothing, but found the transverse pairs of certified on every cut of
# experiments 1 and 3 to 12 (2026-10-27); sample_genus_10.jl uses it. Run one experiment per process
# on the laptop, with a 7G cap: four in one process went past 5 GB, and so did
# experiment 4 alone.
#
# Experiments K to K+N-1 (default 1 to 1) go to DIR/experiment_k, by default
# in results/curves/hyperplane_census. Every object is stored in .mrdi:
#   hyperplanes.mrdi     the 5 hyperplanes: complexes, normal vectors, shifts
#   intersection_s.mrdi  stage s-1 cut by hyperplane s, and its pieces
#   stage_s.mrdi         its coarsening, the complex cut next, with multiplicities
#   curve.mrdi           stage 4 modulo lineality, its graphs and invariants
#   points.mrdi          stage 5 modulo lineality, with multiplicities and cells
# beside stage_s.csv (the cells and their pieces), curve_desubdivided.csv,
# curve_reduced.csv and curve_core.csv for build_curve_atlas.jl, and
# summary.csv. A stored object is loaded instead of recomputed, so an
# interrupted census resumes. DIR/census.csv collects the summaries. The
# shifts are drawn in -T:T, by default 2^20, so that the cuts are transverse.
###############################################################################

include(joinpath(@__DIR__, "tropical_hyperplanes.jl"))
include(joinpath(@__DIR__, "intersect_complexes.jl"))
include(joinpath(@__DIR__, "decorated_curve_graph.jl"))

using .IntersectionLib: load_shift_basis

# The first 4 hyperplanes cut the curve, the fifth its points.
const N_HYPERPLANES = 5
const CURVE_STAGE = 4
const TARGET_GENUS = 10
const DEGREE = 18
const CENSUS_HEADER = "experiment,status,supports,pieces,cells,certified,vertices,edges,legs,genus," *
                      "core_vertices,core_edges,core_valences,multiplicities,balanced," *
                      "points,point_multiplicities,degree,origin,candidates,target,seconds"
const STEERING_HEADER = "cut,candidate,support,base_locus,status,cells,share_above_1,mean_multiplicity," *
                        "genus,core_valences,points,degree,chosen"

"""
The coarsest fan of the tropical G2-Grassmannian as a polyhedral complex, and
the multiplicities over QQbar of its maximal cells, those of the raw cones they
merge.
"""
function g2_grassmannian()
    fan = load(data_path("g2", "grassmannian", "02_coarsest_fan.mrdi"))
    weights = ZZ.(coarsest_cell_weights())
    @req length(weights) == n_maximal_cones(fan) "02_coarsest_cells.csv does not list the cells of the fan"
    complex = polyhedral_complex(fan)
    @req all(((P, C),) -> P == polyhedron(C), zip(maximal_polyhedra(complex), maximal_cones(fan))) "the complex reorders the cones"
    return complex, weights
end

"""
The 5 hyperplanes of experiment `k` in `QQ^n`, drawn with `Xoshiro(k)` among
`coordinates`, the first 4 with supports of sizes in `support_sizes`. The last
one has full support: a hyperplane misses the points of the curve where all its
coordinates vanish, so that fewer than 18 points remain.
"""
function draw_hyperplanes(k, n; shift=1 << 20, coordinates=1:n, support_sizes=2:n)
    rng = Xoshiro(k)
    sizes(s) = s < N_HYPERPLANES ? support_sizes : (length(coordinates):length(coordinates))
    return [random_shifted_tropical_hyperplane(n; rng, shift, coordinates, support_sizes=sizes(s)) for s in 1:N_HYPERPLANES]
end

"""
The hyperplanes of experiment `k` of `draw_hyperplanes`, stored in `path`. A
stored file that differs from the draw after the curve stage, from an earlier
version of the census, is replaced, and the stages it cut are removed.
"""
function experiment_hyperplanes(path, k, n, shift, coordinates, support_sizes)
    drawn = draw_hyperplanes(k, n; shift, coordinates, support_sizes)
    if isfile(path)
        stored = load(path)
        drawn_among = haskey(stored, :coordinates) ? stored.coordinates : collect(1:n)
        @req drawn_among == collect(coordinates) "$path was drawn among the coordinates $drawn_among"
        same = [s <= length(stored.shifts) && stored.shifts[s] == drawn[s].shift && stored.normals[s] == drawn[s].normal
                for s in 1:N_HYPERPLANES]
        all(same) && return stored
        changed = findfirst(!, same)
        @req changed > CURVE_STAGE "$path is not drawn with Xoshiro($k), shift $shift and support sizes $support_sizes"
        stale = [["stage_$s.mrdi", "stage_$s.csv", "intersection_$s.mrdi"] for s in changed:N_HYPERPLANES]
        foreach(name -> rm(joinpath(dirname(path), name); force=true), [reduce(vcat, stale); "points.mrdi"; "summary.csv"])
    end
    save(path, (complexes=[h.complex for h in drawn], normals=[h.normal for h in drawn],
                shifts=[h.shift for h in drawn], seed=k, coordinates=collect(coordinates)))
    return load(path)
end

"Every transverse piece of `result` as a cell of its own, for a cut whose coarsening failed."
function uncoarsened_intersection(result::ComplexIntersection)
    pieces = [(i, c) for i in sort!(collect(keys(result.intersections))) for c in result.intersections[i] if c.transverse]
    by_reference = Dict(i => Int[] for i in eachindex(result.reference_cells))
    isempty(pieces) && return CoarseIntersection(Oscar._empty_complex(QQ, ambient_dim(first(result.reference_cells))),
                                                 CoarseCell[], by_reference, false)
    generators, lineality, generator_sets = homogenized_cells([c.polyhedron for (_, c) in pieces])
    complex = dehomogenized_complex(generator_sets, generators, lineality; non_redundant=true)
    cells = [CoarseCell(P, [(i, c.other_cell)], c.dim, c.transverse, c.multiplicity)
             for (P, (i, c)) in zip(maximal_polyhedra(complex), pieces)]
    foreach(((k, (i, _)),) -> push!(by_reference[i], k), enumerate(pieces))
    return CoarseIntersection(complex, cells, by_reference, false)
end

"""
The record of the cut `result` stored in `intersection_s.mrdi`, and that of its
coarsening, the complex cut next, stored in `stage_s.mrdi`. The cut is
transverse when all its cells are transverse of dimension `d`.
"""
function stage_records(result::ComplexIntersection, d::Int)
    meetings = [(i, c) for i in sort!(collect(keys(result.intersections))) for c in result.intersections[i]]
    intersection = (complex=result.complex, reference_cells=first.(meetings),
                    hyperplane_cells=[c.other_cell for (_, c) in meetings], dims=[c.dim for (_, c) in meetings],
                    transverse=[c.transverse for (_, c) in meetings], multiplicities=[c.multiplicity for (_, c) in meetings])
    coarse, coarsened = try
        (coarsen_intersection(result), true)
    catch error
        error isa InterruptException && rethrow()
        @warn "The coarsening failed; the transverse pieces are kept" exception=(error, catch_backtrace())
        (uncoarsened_intersection(result), false)
    end
    println(coarse)
    stage = (; complex=coarse.complex, multiplicities=[c.multiplicity for c in coarse.cells],
             reference_cells=[first.(c.pieces) for c in coarse.cells], hyperplane_cells=[last.(c.pieces) for c in coarse.cells],
             n_pieces=count(c -> c.transverse, last.(meetings)), coarsened, coarsest_certified=coarse.coarsest_certified,
             transverse=!isempty(coarse.cells) && all(c -> c.transverse && c.dim == d, coarse.cells))
    return intersection, stage, coarse
end

"The cells of a complex modulo its lineality `(1,...,1)`, in the hyperplane `x_1 = 0`: segments and rays of a curve, or points."
function modulo_lineality(complex)
    n = ambient_dim(complex)
    slice = polyhedron((Matrix{QQFieldElem}(undef, 0, n), QQFieldElem[]), (permutedims(QQ.([1; zeros(Int, n - 1)])), [QQ(0)]))
    return [intersect(P, slice) for P in maximal_polyhedra(complex)]
end

"Whether the multiplicities times the primitive directions of the edges and legs at every vertex add up to zero."
is_balanced(list, multiplicities) = all(iszero, values(balancing_sums(list, multiplicities).sums))

"A decorated graph as plain vectors, with the multiplicity of every edge and leg."
graph_record(G::DecoratedGraph, multiplicities) =
    (points=G.points, tails=[e.tail for e in G.edges], heads=[e.head for e in G.edges],
     cells=[e.cells for e in G.edges], pieces=[e.pieces for e in G.edges],
     multiplicities=[multiplicities[first(e.pieces)] for e in G.edges])

"""
The curve of the last stage modulo lineality, its graphs as in
`decorated_curve_graph.jl` and its invariants. The graph has the target genus
and its core is trivalent exactly when `target` holds.
"""
function curve_record(stage, cells)
    list = modulo_lineality(stage.complex)
    G = decorated_graph(list, cells)
    graphs = (subdivided=G, desubdivided=desubdivide(G), reduced=desubdivide(G; within_cells=true),
              core=desubdivide(core(G)))
    core_valences = graph_valences(graphs.core)
    g = genus(G)
    return (; polyhedra=list, multiplicities=stage.multiplicities,
            (name => graph_record(graphs[name], stage.multiplicities) for name in keys(graphs))...,
            genus=g, core_valences, balanced=is_balanced(list, stage.multiplicities),
            target=g == TARGET_GENUS && !isempty(core_valences) && all(==(3), core_valences)), graphs
end

"""
The points of the last stage modulo lineality, with their multiplicities, the
cells of the curve (`curve.mrdi`) they lie on, and the maximal cells of the
G2-Grassmannian containing them. Their total multiplicity is the degree.
"""
function point_record(stage, cells)
    list = modulo_lineality(stage.complex)
    @req all(P -> dim(P) == 0, list) "the last stage is not a set of points modulo lineality"
    return (; points=[collect(only(vertices(P))) for P in list], multiplicities=stage.multiplicities,
            curve_cells=stage.reference_cells, cells=[containing_cells(P, cells) for P in list])
end

"""
The `s`-th cut, of `complex` by `hyperplane`, and its coarsening: the records
of `stage_records` for cells of dimension `d`.
"""
function cut_by(complex, multiplicities, hyperplane, s, d; solver, label, reference_data=nothing)
    # Polymake objects live outside the Julia heap, so the GC does not see their weight.
    GC.gc(true)
    println("$label, cut $s: $(n_maximal_polyhedra(complex)) cells by ",
            "$(n_maximal_polyhedra(hyperplane)) cells of the hyperplane")
    seconds = @elapsed result = intersect_complexes(complex, hyperplane; reference_weights=multiplicities,
                                                    feasibility_solver=solver, reference_data)
    println(result, " ($(round(Int, seconds)) s)")
    seconds = @elapsed records = stage_records(result, d)
    println("  coarsened in $(round(Int, seconds)) s")
    flush(stdout)
    return records
end

"""
Run experiment `k` in `directory`, loading the objects already stored there,
and return its row of the census.
"""
function run_census_experiment(directory, k, grassmannian, weights, cells, labels; shift, solver, coordinates, support_sizes)
    mkpath(directory)
    file(name) = joinpath(directory, name)
    hyperplanes = experiment_hyperplanes(file("hyperplanes.mrdi"), k, ambient_dim(grassmannian), shift, coordinates,
                                         support_sizes)
    row = Dict{String,Any}("experiment" => k, "supports" => join((count(!iszero, a) for a in hyperplanes.normals), ';'))
    complex, multiplicities = grassmannian, weights
    stages = []
    for s in 1:N_HYPERPLANES
        if isfile(file("stage_$s.mrdi"))
            stage = load(file("stage_$s.mrdi"))
        else
            intersection, stage, coarse = cut_by(complex, multiplicities, hyperplanes.complexes[s], s,
                                                 dim(grassmannian) - s; solver, label="Experiment $k")
            save(file("intersection_$s.mrdi"), intersection)
            save(file("stage_$s.mrdi"), stage)
            write_coarse_table(file("stage_$s.csv"), coarse)
            flush(stdout)
        end
        @req n_maximal_polyhedra(stage.complex) == length(stage.multiplicities) "stage_$s.mrdi has one multiplicity per cell"
        push!(stages, stage)
        complex, multiplicities = stage.complex, stage.multiplicities
        stage.transverse || (row["status"] = "cut $s not transverse"; break)
    end
    return summarize_stages!(row, stages, file, cells, labels)
end

"""
Complete the census `row` of an experiment from its `stages`: the curve of
stage 4 and the points of stage 5, stored with `file`.
"""
function summarize_stages!(row, stages, file, cells, labels)
    row["pieces"] = join((stage.n_pieces for stage in stages), ';')
    row["cells"] = join((length(stage.multiplicities) for stage in stages), ';')
    row["certified"] = join((stage.coarsest_certified ? 1 : 0 for stage in stages), ';')
    haskey(row, "status") && return row
    for s in (CURVE_STAGE, N_HYPERPLANES)
        L = lineality_dim(stages[s].complex)
        L == 1 || (row["status"] = "lineality of dimension $L at stage $s"; return row)
    end

    curve, graphs = curve_record(stages[CURVE_STAGE], cells)
    save(file("curve.mrdi"), curve)
    for name in (:desubdivided, :reduced, :core)
        write_decorated_graph(file("curve_$name.csv"), graphs[name], labels...)
    end
    foreach(name -> describe(stdout, name, graphs[name]), keys(graphs))
    G = graphs.subdivided
    merge!(row, Dict("status" => "curve", "vertices" => length(G.points), "edges" => count(!is_leg, G.edges),
                     "legs" => count(is_leg, G.edges), "genus" => curve.genus,
                     "core_vertices" => length(graphs.core.points), "core_edges" => length(graphs.core.edges),
                     "core_valences" => histogram(curve.core_valences),
                     "multiplicities" => histogram(curve.multiplicities), "balanced" => curve.balanced))
    points = point_record(stages[N_HYPERPLANES], cells)
    save(file("points.mrdi"), points)
    degree = sum(points.multiplicities)
    merge!(row, Dict("points" => length(points.points), "point_multiplicities" => histogram(points.multiplicities),
                     "degree" => degree, "target" => curve.target && curve.balanced && degree == DEGREE))
    return row
end

"""
The projective dimension of the part of G2/P2, of ideal `I`, where all the
coordinates `covered` vanish; -1 when it is empty. Every hyperplane supported
on `covered` contains it.
"""
base_locus_dimension(I, covered) = dim(I + ideal(base_ring(I), gens(base_ring(I))[covered])) - 1

"""
The score of the cut `stage` as the `s`-th cut of a steered experiment, larger
being better, and its row of `steering.csv`. Before the curve it is the share
of cells of multiplicity 1, at the curve the target and then the genus, and at
the end the degree and the number of points.
"""
function candidate_score(stage, s, cells)
    m = stage.multiplicities
    share = isempty(m) ? 1.0 : count(>(1), m) / length(m)
    info = Dict{String,Any}("cells" => length(m), "share_above_1" => round(share; digits=3),
                            "mean_multiplicity" => isempty(m) ? "" : round(Float64(sum(m)) / length(m); digits=3))
    stage.transverse && lineality_dim(stage.complex) == 1 || return (false,), merge!(info, Dict("status" => "not transverse"))
    info["status"] = "cut"
    s < CURVE_STAGE && return (true, -share), info
    if s == CURVE_STAGE
        curve, _ = curve_record(stage, cells)
        merge!(info, Dict("genus" => curve.genus, "core_valences" => histogram(curve.core_valences)))
        return (true, curve.target && curve.balanced, curve.genus, -share), info
    end
    points = point_record(stage, cells)
    degree = sum(points.multiplicities)
    merge!(info, Dict("points" => length(points.points), "degree" => degree))
    return (true, degree == DEGREE, length(points.points)), info
end

"""
The `s`-th cut of steered experiment `k`: the best of `candidates` cuts of
`complex` by hyperplanes drawn with `Xoshiro(1000k + s)`, or `nothing` when
none is kept; at the curve, the first on target. Every candidate is stored in
`directory/cut_s/`, its row in `steering.csv`.
"""
function steered_cut(directory, k, s, complex, multiplicities, previous_normals, I, cells, d;
                     candidates, shift, shift_basis, solver)
    rng = Xoshiro(1000k + s)
    n = ambient_dim(complex)
    sizes = s < N_HYPERPLANES ? (2:n) : (n:n)
    covered = reduce(union, (findall(!iszero, a) for a in previous_normals); init=Int[])
    mkpath(joinpath(directory, "cut_$s"))
    reference_data = nothing  # the cells of `complex`, computed for the first candidate cut
    best, best_score, rows = nothing, nothing, Dict{String,Any}[]
    for j in 1:candidates
        h = random_shifted_tropical_hyperplane(n; rng, shift, shift_basis, support_sizes=sizes)
        support = findall(!iszero, h.normal)
        locus = base_locus_dimension(I, sort!(union(covered, support)))
        row = Dict{String,Any}("cut" => s, "candidate" => j, "support" => length(support), "base_locus" => locus)
        push!(rows, row)
        locus >= 5 - s && (row["status"] = "base locus"; continue)
        path = joinpath(directory, "cut_$s", "candidate_$j.mrdi")
        record = if isfile(path)
            load(path)
        else
            isnothing(reference_data) && (reference_data = cell_data(complex))
            intersection, stage, coarse = cut_by(complex, multiplicities, h.complex, s, d; solver, reference_data,
                                                 label="Experiment $k, candidate $j")
            write_coarse_table(joinpath(directory, "cut_$s", "candidate_$j.csv"), coarse)
            hyperplane = (complex=h.complex, normal=h.normal, shift=h.shift,
                          coefficients=something(h.coefficients, Int[]))
            save(path, (; hyperplane, intersection, stage))
            (; hyperplane, intersection, stage)
        end
        score, info = candidate_score(record.stage, s, cells)
        merge!(row, info)
        println("Experiment $k, cut $s, candidate $j: ", join(("$c=$(get(row, c, ""))" for c in split(STEERING_HEADER, ',')[3:end-1]), ' '))
        flush(stdout)
        first(score) && (isnothing(best_score) || score > best_score) && ((best, best_score) = ((j, record), score))
        # A curve on target ends the search.
        s == CURVE_STAGE && first(score) && score[2] && break
    end
    foreach(row -> row["chosen"] = !isnothing(best) && row["candidate"] == first(best), rows)
    path = joinpath(directory, "steering.csv")
    isfile(path) || write(path, STEERING_HEADER * '\n')
    open(io -> foreach(row -> println(io, join((get(row, c, "") for c in split(STEERING_HEADER, ',')), ',')), rows),
         path, "a")
    return best, count(row -> haskey(row, "cells"), rows)
end

"""
Run the steered experiment `k` in `directory` and return its row of the
census: the stages up to `from_stage` are those of the experiment in
`origin_directory`, numbered `origin`, and every later cut is a `steered_cut`.
"""
function run_steered_experiment(directory, k, grassmannian, weights, cells, labels, I; origin, origin_directory,
                                from_stage, candidates, point_candidates, shift, shift_basis, solver)
    mkpath(directory)
    file(name) = joinpath(directory, name)
    for s in 1:from_stage, name in ("stage_$s.mrdi", "stage_$s.csv", "intersection_$s.mrdi")
        isfile(file(name)) || cp(joinpath(origin_directory, name), file(name))
    end
    chosen = if isfile(file("hyperplanes.mrdi"))
        stored = load(file("hyperplanes.mrdi"))
        (complexes=collect(stored.complexes), normals=collect(stored.normals), shifts=collect(stored.shifts),
         coefficients=collect(stored.coefficients))
    elseif from_stage == 0
        (complexes=PolyhedralComplex{QQFieldElem}[], normals=Vector{QQFieldElem}[], shifts=Vector{Int}[],
         coefficients=Vector{Int}[])
    else
        inherited = load(joinpath(origin_directory, "hyperplanes.mrdi"))
        (complexes=collect(inherited.complexes[1:from_stage]), normals=collect(inherited.normals[1:from_stage]),
         shifts=collect(inherited.shifts[1:from_stage]), coefficients=[Int[] for _ in 1:from_stage])
    end
    save_chosen() = save(file("hyperplanes.mrdi"), (; chosen..., seed=k, origin, from_stage))
    row = Dict{String,Any}("experiment" => k, "origin" => from_stage > 0 ? "$origin:$from_stage" : "")
    stages = [load(file("stage_$s.mrdi")) for s in 1:from_stage]
    complex, multiplicities = isempty(stages) ? (grassmannian, weights) : (last(stages).complex, last(stages).multiplicities)
    evaluated = Int[]
    for s in from_stage+1:N_HYPERPLANES
        if isfile(file("stage_$s.mrdi")) && length(chosen.shifts) >= s
            stage = load(file("stage_$s.mrdi"))
        else
            best, n_evaluated = steered_cut(directory, k, s, complex, multiplicities, chosen.normals, I, cells,
                                            dim(grassmannian) - s; shift, shift_basis, solver,
                                            candidates=s == N_HYPERPLANES ? point_candidates : candidates)
            push!(evaluated, n_evaluated)
            isnothing(best) && (row["status"] = "no candidate kept at cut $s"; break)
            j, record = best
            stage = record.stage
            save(file("stage_$s.mrdi"), stage)
            save(file("intersection_$s.mrdi"), record.intersection)
            cp(joinpath(directory, "cut_$s", "candidate_$j.csv"), file("stage_$s.csv"); force=true)
            push!(chosen.complexes, record.hyperplane.complex)
            push!(chosen.normals, record.hyperplane.normal)
            push!(chosen.shifts, record.hyperplane.shift)
            push!(chosen.coefficients, record.hyperplane.coefficients)
            save_chosen()
        end
        push!(stages, stage)
        complex, multiplicities = stage.complex, stage.multiplicities
    end
    row["supports"] = join((count(!iszero, a) for a in chosen.normals), ';')
    row["candidates"] = join(evaluated, ';')
    return summarize_stages!(row, stages, file, cells, labels)
end

census_line(row) = join((get(row, column, "") for column in split(CENSUS_HEADER, ',')), ',')

is_summarized(directory) = isfile(joinpath(directory, "summary.csv"))

"The row of `summary.csv` in `directory`, by column, whatever version of the census wrote it."
function summary_row(directory)
    header, line = readlines(joinpath(directory, "summary.csv"))[1:2]
    return Dict(zip(split(header, ','), split(line, ',')))
end

"Collect the summaries of the experiments in `output` into `census.csv`."
function write_census(output)
    numbered = [(parse(Int, m[1]), joinpath(output, name)) for name in readdir(output)
                for m in (match(r"^experiment_(\d+)$", name),) if !isnothing(m)]
    lines = [census_line(summary_row(path)) for (_, path) in sort!(numbered) if is_summarized(path)]
    write(joinpath(output, "census.csv"), join([CENSUS_HEADER; lines], '\n') * '\n')
    return lines
end

function parse_options(args)
    options = Dict("first" => "1", "experiments" => "1", "shift" => "", "solver" => "certified",
                   "coordinates" => join(1:21, ','), "support-sizes" => "", "from" => "", "candidates" => "1", "point-candidates" => "",
                   "shift-basis" => "", "output" => project_path("results", "curves", "hyperplane_census"))
    for arg in args
        m = match(r"^--([a-z-]+)=(.+)$", arg)
        !isnothing(m) && haskey(options, m[1]) || error("unknown argument $arg; see the header of $(@__FILE__)")
        options[m[1]] = m[2]
    end
    coordinates = parse.(Int, split(options["coordinates"], ','))
    sizes = isempty(options["support-sizes"]) ? [2, length(coordinates)] : parse.(Int, split(options["support-sizes"], ':'))
    basis = options["shift-basis"] == "default" ? data_path("mukai_curves", "genus_10_shift_basis.mrdi") :
            isempty(options["shift-basis"]) ? nothing : abspath(options["shift-basis"])
    shift = isempty(options["shift"]) ? (isnothing(basis) ? 1 << 20 : 1 << 15) : parse(Int, options["shift"])
    origin, from_stage = isempty(options["from"]) ? (0, 0) : parse.(Int, split(options["from"], ':'))
    candidates = parse(Int, options["candidates"])
    point_candidates = isempty(options["point-candidates"]) ? candidates : parse(Int, options["point-candidates"])
    steered = from_stage > 0 || candidates > 1 || !isnothing(basis)
    steered && (coordinates == 1:21 && isempty(options["support-sizes"]) ||
                error("--coordinates and --support-sizes are not used by steered experiments"))
    output = abspath(options["output"])
    return (first=parse(Int, options["first"]), experiments=parse(Int, options["experiments"]), shift,
            solver=Symbol(options["solver"]), coordinates, support_sizes=first(sizes):last(sizes), output, steered,
            origin, from_stage, origin_directory=joinpath(output, "experiment_$origin"), candidates, point_candidates,
            shift_basis=isnothing(basis) ? nothing : load_shift_basis(basis))
end

if abspath(PROGRAM_FILE) == @__FILE__
    options = parse_options(ARGS)
    grassmannian, weights = g2_grassmannian()
    cells = collect(maximal_cones(load(data_path("g2", "grassmannian", "02_coarsest_fan.mrdi"))))
    labels = cell_labels()
    I = load(data_path("g2", "grassmannian", "00_ideal.mrdi"))
    for k in options.first:options.first+options.experiments-1
        directory = joinpath(options.output, "experiment_$k")
        is_summarized(directory) && (println("Experiment $k: done"); continue)
        start = time()
        row = try
            options.steered ?
                run_steered_experiment(directory, k, grassmannian, weights, cells, labels, I; options.origin,
                                       options.origin_directory, options.from_stage, options.candidates,
                                       options.point_candidates, options.shift,
                                       options.shift_basis, options.solver) :
                run_census_experiment(directory, k, grassmannian, weights, cells, labels; options.shift, options.solver,
                                      options.coordinates, options.support_sizes)
        catch error
            error isa InterruptException && rethrow()
            @error "Experiment $k failed" exception=(error, catch_backtrace())
            nothing
        end
        isnothing(row) && continue
        row["seconds"] = round(Int, time() - start)
        write(joinpath(directory, "summary.csv"), CENSUS_HEADER * '\n' * census_line(row) * '\n')
        println("Experiment $k: ", census_line(row))
        flush(stdout)
        GC.gc(true)
    end
    lines = write_census(options.output)
    targets = [line for line in lines if endswith(line, r",true,\d+$")]
    println("$(length(lines)) experiments in $(joinpath(options.output, "census.csv")), $(length(targets)) on target")
    foreach(println, targets)
end
