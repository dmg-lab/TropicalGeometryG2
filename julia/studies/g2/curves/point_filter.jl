###############################################################################
# The points of trop(G2/P2) ∩ (U(16,21) + t), as a filter in front of
# core_trivalence.jl. For a generic t, x - t attains its minimum exactly 6
# times at such a point, on a subset S, in the interior of a maximal cell C.
# Normalized by min(x - t) = 0, x is the point of span(C) with x = t on S,
# x = B B_S^{-1} t_S for a basis B of span(C), and it is a point of the
# intersection when x > t off S and A x < 0 for the facets A of C. With B a
# ZZ-basis of the lattice of span(C), its multiplicity is w_C |det B_S|, and
# the multiplicities add up to 18, the degree of G2/P2.
#
# These points are the 6-valent vertices of the curve
# trop(G2/P2) ∩ (U(17,21) + t). In the night search of 2026-10-27/28
# (41580 trials, wpart and uniform), 2.2% of the translations had 18 of them,
# 92% of those gave a curve of genus 10, and all 8 curves of genus 10 with a
# trivalent core were among them.
#
#   julia --project=julia -t 8 julia/studies/g2/curves/point_filter.jl KIND TRIALS [--first=K]
#         [--translation=FILE] [--batch=N] [--output=DIR]
#   julia --project=julia -t 8 julia/studies/g2/curves/point_filter.jl --check=TABLE
#         [--trials=K1;K2;...] [--sample=N] [--curves]
#
# The first form draws the translations of trials K..TRIALS of KIND as
# core_trivalence.jl does (Xoshiro(k)), computes their points, and the curve
# of those with 18 distinct points, into DIR/point_filter_KIND.csv (DIR
# defaults to results/curves/genus_10/point_filter), N trials at a time (1000
# by default: the factorizations of the 6×6 minors, once per batch, take most
# of the time); the trials already in the table are skipped. Curves of genus
# 10 with a trivalent core are saved as
# DIR/KIND_K/decorated_graph_K.mrdi with the CSV files of
# decorated_curve_graph.jl. The second form compares the points of trials of
# a table of core_trivalence.jl (the given ones, N random ones, by default
# those of genus 10 with a trivalent core) with the 6-valent vertices of
# their curves: their number from the table, and with --curves the points
# themselves, from the curves computed again.
###############################################################################

include("core_trivalence.jl")
import TropicalCensus: csv_rows
using TropicalCensus: study_path

const POINT_SUBSETS = collect(Combinatorics.combinations(1:21, 6))
const POINT_COLUMNS = (:points, :degree, :degenerate, :multiplicities, :point_orbits, :point_cells)

"The multiplicity over QQbar of every maximal cell of the coarsest fan, that of the raw cones it merges."
function cell_weights()
    raw = Dict(row["cone"] => parse(Int, row["mult_QQbar"]) for row in csv_rows(study_path(:g2, "weights", "weights.csv")))
    rows = csv_rows(data_path("g2", "grassmannian", "02_coarsest_cells.csv"))
    @req [parse(Int, row["coarsest_cell"]) for row in rows] == 1:length(rows) "02_coarsest_cells.csv does not list the cells in order"
    return [only(unique(raw[c] for c in split(row["original_cells"], ';'))) for row in rows]
end

"""
The cell of span `M` and facets `A` for `exact_point` and `candidate_points`:
an LLL-reduced ZZ-basis `B` of the lattice points of its span, as columns,
whose 6×6 minors are small integers and so exact in floating point, and its
facets in the coordinates of `B`, scaled to largest entry 1 in floating point.
"""
function point_cell(M, A)
    equations = transpose(kernel(transpose(M); side=:right))
    rows = [collect(equations[i, :]) for i in 1:nrows(equations)]
    integral = matrix(ZZ, reduce(vcat, [permutedims(numerator.(r .* lcm(denominator.(r)))) for r in rows]))
    B =transpose(lll(transpose(kernel(integral; side=:right))))
    @req ncols(B) == 6 "the span of the cell is not 6-dimensional"
    BQ = matrix(QQ, B)
    ABf = Float64.(Matrix(A * BQ))
    ABf ./= maximum(abs, ABf; dims=2)
    return (B=B, BQ=BQ, A=A, Bf=Float64.(Matrix(B)), ABf=ABf)
end

"""
For every column t of `T`, the indices in `POINT_SUBSETS` of the subsets S
whose point, the x in the span of `Bf` with x = t on S, may in floating point
lie in the cell of facets `ABf` with x ≥ t, up to the margins of the
translations. A subset of vanishing minor has no point for a generic
translation, and the minors are integers.
"""
function candidate_points(Bf, ABf, T, margins)
    N = size(T, 2)
    found = [Int[] for _ in 1:N]
    E, R = zeros(6, 6), zeros(6, N)
    for (s, S) in enumerate(POINT_SUBSETS)
        for a in 1:6, b in 1:6
            E[a, b] = Bf[S[a], b]
        end
        F = LinearAlgebra.lu!(E; check=false)
        abs(LinearAlgebra.det(F)) < 0.5 && continue
        for n in 1:N, a in 1:6
            R[a, n] = T[S[a], n]
        end
        LinearAlgebra.ldiv!(F, R)
        for n in 1:N
            m = margins[n]
            inside = all(1:size(ABf, 1)) do r
                sum(ABf[r, b] * R[b, n] for b in 1:6) <= m
            end
            inside || continue
            above = all(1:21) do j
                sum(Bf[j, b] * R[b, n] for b in 1:6) - T[j, n] >= -m
            end
            above && push!(found[n], s)
        end
    end
    return found
end

"""
The point of the subset `S` in the cell `cell` of `point_cell` for the
translation `t`, exactly: `nothing` when it is not a point of the
intersection, `:degenerate` when x - t also vanishes off S or x lies on the
boundary of the cell, that is, when t is not generic, and otherwise x.
"""
function exact_point(cell, t::Vector{QQFieldElem}, S)
    y = solve(cell.BQ[S, :], matrix(QQ, 6, 1, t[S]); side=:right)
    x = (cell.BQ * y)[:, 1]
    above = [x[j] - t[j] for j in 1:21 if !(j in S)]
    facets = (cell.A * matrix(QQ, 21, 1, x))[:, 1]
    all(>=(0), above) && all(<=(0), facets) || return nothing
    all(>(0), above) && all(<(0), facets) || return :degenerate
    return x
end

"""
The points of `C ∩ (U(16,21) + t)` over the cells `cells` of `point_cell`,
of multiplicities `weights`, for the translations `shifts`: for each, the
points as named tuples (cell, subset, point normalized, multiplicity) and the
number of degenerate ones.
"""
function intersection_points(cells, weights, shifts::Vector{Vector{QQFieldElem}})
    T = reduce(hcat, [Float64.(t) for t in shifts])
    margins = [FLOAT_MARGIN * (1 + maximum(abs, t)) for t in eachcol(T)]
    found = Vector{Vector{Vector{Int}}}(undef, length(cells))
    Threads.@threads for c in eachindex(cells)
        found[c] = candidate_points(cells[c].Bf, cells[c].ABf, T, margins)
    end
    results = Vector{Any}(undef, length(shifts))
    Threads.@threads for n in eachindex(shifts)
        points, degenerate = NamedTuple[], 0
        for c in eachindex(cells), s in found[c][n]
            S = POINT_SUBSETS[s]
            x = exact_point(cells[c], shifts[n], S)
            isnothing(x) && continue
            x === :degenerate && (degenerate += 1; continue)
            push!(points, (cell=c, subset=S, point=normalized(x), multiplicity=weights[c] * abs(det(cells[c].B[S, :]))))
        end
        results[n] = (points=points, degenerate=degenerate)
    end
    return results
end

"A row of the table for the points `result` of `intersection_points`."
function describe_points(result, orbit)
    points = result.points
    return (points=length(points), degree=sum(p -> p.multiplicity, points; init=ZZ(0)), degenerate=result.degenerate,
            multiplicities=histogram([p.multiplicity for p in points]), point_orbits=histogram([orbit[p.cell] for p in points]),
            point_cells=join(sort([p.cell for p in points]), ';'))
end

"Whether the translation of `result` passes the filter: 18 distinct points, so each of multiplicity 1."
passes(result) = result.degenerate == 0 && length(result.points) == 18

"The translation of a table of core_trivalence.jl, `;`-separated rationals `a` or `a//b`."
parse_translation(s) = [QQ(map(x -> parse(BigInt, x), split(q, "//"))...) for q in split(s, ';')]

"The number `k` of the term `6^k` of a histogram of valences, 0 without one."
six_valent(valences) = something(tryparse(Int, last([["0"]; [split(v, '^')[2] for v in split(valences) if startswith(v, "6^")]])), 0)

function check(table, cells, spans, weights, orbit, args, option)
    rows = csv_rows(table)
    by_trial = Dict(parse(Int, r["trial"]) => r for r in rows)
    trials = option("trials", "")
    sample = parse(Int, option("sample", "0"))
    chosen = !isempty(trials) ? parse.(Int, split(trials, ';')) :
             sample > 0 ? sort(Random.shuffle(Random.Xoshiro(0), collect(keys(by_trial)))[1:sample]) :
             sort([k for (k, r) in by_trial if r["genus"] == "10" && r["trivalent"] == "true"])
    shifts = [parse_translation(by_trial[k]["translation"]) for k in chosen]
    elapsed = @elapsed results = intersection_points(cells, weights, shifts)
    println("$(length(chosen)) trials of $table: points in $(round(elapsed; digits=1)) s on $(Threads.nthreads()) threads")
    agree = 0
    for (k, result) in zip(chosen, results)
        d = describe_points(result, orbit)
        expected = six_valent(by_trial[k]["curve_valences"])
        agree += d.points == expected && d.degree == 18
        (d.points == expected && d.degree == 18) ||
            println("trial $k: $(d.points) points, degree $(d.degree), $(d.degenerate) degenerate, but 6^$expected in the table")
    end
    println("$agree of $(length(chosen)) trials: as many points as 6-valent vertices, degree 18")
    "--curves" in args || return
    elapsed = @elapsed curve_results = curves(spans, shifts)
    same = 0
    for (k, result, curve) in zip(chosen, results, curve_results)
        G = curve.graph
        six = Set(G.points[findall(==(6), graph_valences(G))])
        same += six == Set(p.point for p in result.points)
        six == Set(p.point for p in result.points) || println("trial $k: the points differ from the 6-valent vertices")
    end
    println("$same of $(length(chosen)) trials: the points are the 6-valent vertices of the curve (curves in $(round(elapsed; digits=1)) s)")
end

function main(args)
    option(name, default) = last([[default]; [split(a, '='; limit=2)[2] for a in args if startswith(a, "--$name=")]])
    positional = filter(a -> !startswith(a, "--"), args)
    fan = load(data_path("g2", "grassmannian", "02_coarsest_fan.mrdi"))
    spans = [cell_span(C) for C in maximal_cones(fan)]
    cells = [point_cell(M, A) for (M, A) in spans]
    weights = cell_weights()
    @req length(weights) == length(cells) "one multiplicity per cell"
    orbit, number = cell_labels()

    table = option("check", "")
    isempty(table) || return check(table, cells, spans, weights, orbit, args, option)

    length(positional) == 2 && positional[1] in CORE_KINDS ||
        error("usage: point_filter.jl KIND TRIALS [--first=K] [--translation=FILE] [--batch=N] [--output=DIR], KIND in $CORE_KINDS")
    kind, trials = positional[1], parse(Int, positional[2])
    first_trial = parse(Int, option("first", "1"))
    batch = parse(Int, option("batch", "1000"))
    file = option("translation", project_path("curves_g10", "uniform_shifts", "translation_3.mrdi"))
    output = abspath(option("output", project_path("results", "curves", "genus_10", "point_filter")))
    mkpath(output)
    translation = trial_translations(fan, kind, load(String(file)))

    table = joinpath(output, "point_filter_$(kind).csv")
    done = isfile(table) ? Set(parse(Int, first(split(l, ','))) for l in readlines(table)[2:end]) : Set{Int}()
    isfile(table) || open(io -> println(io, join(["trial"; collect(string.(POINT_COLUMNS)); collect(string.(CORE_COLUMNS)); "translation"], ",")), table, "w")
    todo = [k for k in first_trial:trials if !(k in done)]
    println("$kind: $(length(todo)) of trials $first_trial-$trials to run, translation $file, $(Threads.nthreads()) threads, output $table")
    approximations = [(Float64.(Matrix(M)), Float64.(Matrix(A * M))) for (M, A) in spans]
    for range in Iterators.partition(todo, batch)
        shifts = translation.(range)
        elapsed = @elapsed results = intersection_points(cells, weights, shifts)
        selected = findall(passes, results)
        curve_elapsed = @elapsed found = isempty(selected) ? [] : curves(spans, shifts[selected]; approximations)
        curve_of = Dict(zip(selected, found))
        open(table, "a") do io
            for (n, (k, t, result)) in enumerate(zip(range, shifts, results))
                p = describe_points(result, orbit)
                p.degree == 18 || println("trial $k: degree $(p.degree), $(p.degenerate) degenerate points")
                d = haskey(curve_of, n) ? describe_curve(curve_of[n].graph, curve_of[n].higher, orbit) : nothing
                curve = isnothing(d) ? fill("", length(CORE_COLUMNS)) : [d[c] for c in CORE_COLUMNS]
                println(io, join([k; [p[c] for c in POINT_COLUMNS]; curve; join(string.(t), ";")], ","))
                isnothing(d) && continue
                target = d.trivalent && d.genus == 10
                target && save_curve(joinpath(output, "$(kind)_$(k)"), k, curve_of[n].graph, orbit, number)
                target && println("trial $k: genus 10, trivalent core $(d.core)")
            end
        end
        println("trials $(first(range))-$(last(range)): points $(round(elapsed; digits=1)) s, ",
                "$(length(selected)) with 18 points, curves $(round(curve_elapsed; digits=1)) s")
        flush(stdout)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
