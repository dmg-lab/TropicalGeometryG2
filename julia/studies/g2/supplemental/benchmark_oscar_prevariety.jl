#!/usr/bin/env julia

# Oscar's `tropical_prevariety` on the G2 flag Dressian and the G2-Dressian.
# In Oscar 1.8.2 it takes no keywords: it is reduce(common_refinement, ...) over
# the tropical hypersurfaces in the order given (TropicalGeometry/hypersurface.jl),
# and polymake's common_refinement intersects every pair of maximal cells, one
# convex hull per pair. The variants change what can be changed from outside:
# the generator order, the coordinates (the 72 generators in the 22 coordinates
# left by the six binomial linear forms, as for Gfan) and polymake's convex hull
# backend. One variant per process; nothing is written to data/.
#
#   julia --project=julia benchmark_oscar_prevariety.jl --case=flag78 --order=input --hull=default --mode=oscar
#
# --case   flag78 (78 generators, 28 coordinates), flag72 (the 72 generators left by the six
#          binomial linear forms, in 22 coordinates; 71 distinct, two Pluecker relations coincide),
#          dressian42 (42 generators of the G2 Grassmannian, 21 coordinates)
# --order  input (linears, Pluecker, incidence, quadric), colex (linears first, then by the
#          leaves 1..7 of the variables, largest leaf first, so leaf 7 enters last) or greedy
#          (linears first, then always the generator with the fewest new variables)
# --hull   default (polymake's preference: ppl), ppl, cdd, lrs, beneath_beyond
# --mode   oscar (the call tropical_prevariety(F, nu) itself) or steps (the same fold, timed per step)
# --limit  only the first N generators after ordering (tests)

using Oscar
using Printf

include(joinpath(@__DIR__, "..", "ideals", "G2GrassmannianIdeal.jl"))
include(joinpath(@__DIR__, "..", "ideals", "FlagG2Ideal.jl"))
using .G2GrassmannianIdeal
using .FlagG2Ideal

const NU = tropical_semiring_map(QQ, min)
const HULLS = Dict("ppl" => "ppl.convex_hull", "cdd" => "cdd.convex_hull", "lrs" => "lrs.convex_hull",
                   "beneath_beyond" => "beneath_beyond.convex_hull")
# Rays modulo lineality and maximal cells of the Gfan fans in data/g2/*/01_prevariety_fan.mrdi.
const EXPECTED = Dict("flag78" => (147, 2726), "flag72" => (147, 2726), "dressian42" => (51, 239))

function parse_options(arguments)
    options = Dict("case" => "flag78", "order" => "input", "hull" => "default", "mode" => "oscar", "limit" => "0")
    for argument in arguments
        m = match(r"^--([a-z]+)=(.*)$", argument)
        (isnothing(m) || !haskey(options, m[1])) && error("unknown argument $argument")
        options[m[1]] = m[2]
    end
    options["case"] in keys(EXPECTED) || error("unknown case $(options["case"])")
    options["order"] in ("input", "colex", "greedy") || error("unknown order $(options["order"])")
    options["hull"] == "default" || haskey(HULLS, options["hull"]) || error("unknown hull $(options["hull"])")
    options["mode"] in ("oscar", "steps") || error("unknown mode $(options["mode"])")
    return options
end

generators(case) = case == "flag78" ? build_full_flag_g2_ideal().generators :
                   case == "flag72" ? build_flag_g2_ideal().generators :
                   build_g2_grassmannian_ideal().generators

"Leaves 1..7 of a variable x_i or p_ij."
leaves(name) = Set(parse(Int, c) for c in name[2:end])

function colex_order(polynomials)
    names = string.(symbols(parent(first(polynomials))))
    key(f) = (total_degree(f) == 1 ? 0 : 1,
              sort!(collect(union((leaves(names[i]) for i in vars_indices(f))...)); rev=true))
    return sort(polynomials; by=key)
end

vars_indices(f) = [var_index(v) for v in vars(f)]

"Linear forms first, then always the generator with the fewest variables not met so far
(earliest on ties); two variables of a binomial linear form count as one."
function greedy_order(polynomials)
    class = collect(1:ngens(parent(first(polynomials))))
    root(i) = class[i] == i ? i : root(class[i])
    for f in polynomials
        total_degree(f) == 1 && length(f) == 2 || continue
        i, j = root.(vars_indices(f))
        class[max(i, j)] = min(i, j)
    end
    classes(f) = Set(root.(vars_indices(f)))
    ordered = filter(f -> total_degree(f) == 1, polynomials)
    remaining = filter(f -> total_degree(f) != 1, polynomials)
    met = union(Set{Int}(), classes.(ordered)...)
    while !isempty(remaining)
        k = argmin(f -> length(setdiff(classes(f), met)), remaining)
        k = findfirst(==(k), remaining)
        union!(met, classes(remaining[k]))
        push!(ordered, popat!(remaining, k))
    end
    return ordered
end

kind(f) = total_degree(f) == 1 ? "linear" : length(f) == 4 ? "quadric" :
          any(startswith(string(v), "x") for v in vars(f)) ? "incidence" : "pluecker"

function with_hull(f, hull)
    hull == "default" && return f()
    return Polymake.prefer(f, HULLS[hull]; application="polytope")
end

gib() = Sys.maxrss() / 2^30

function counts(complex)
    data = rays_modulo_lineality(complex)
    return (rays=length(data.rays_modulo_lineality), maximal_cells=n_maximal_polyhedra(complex),
            dimension=dim(complex), lineality_dimension=lineality_dim(complex))
end

"The fold of `tropical_prevariety`, timed per step. The cells of the previous
step are read before the next refinement, which would read them anyway."
function stepwise_prevariety(polynomials)
    complexes = [tropical_hypersurface(tropical_polynomial(f, NU)).polyhedralComplex for f in polynomials]
    current, start = complexes[1], time()
    for k in 2:length(complexes)
        refine = @elapsed current = common_refinement(current, complexes[k])
        read = @elapsed cells = n_maximal_polyhedra(current)
        @printf("step %2d/%d %-9s %8.2f s refine %7.2f s cells %7d maximal cells, dim %d, %6.0f s total, peak %.2f GiB\n",
                k, length(complexes), kind(polynomials[k]), refine, read, cells, dim(current), time() - start, gib())
        flush(stdout)
    end
    return current
end

function compute(polynomials, mode)
    mode == "oscar" && return tropical_prevariety(polynomials, NU)
    return stepwise_prevariety(polynomials)
end

"Compile the same code paths on the Pluecker ideal of Gr(2,5)."
function warm_up(options)
    ring, _ = polynomial_ring(QQ, ["p$i$j" for i in 1:4 for j in (i + 1):5])
    small = forget_grading.(gens(grassmann_pluecker_ideal(ring, 2, 5)))
    with_hull(() -> counts(compute(small, options["mode"])), options["hull"])
end

function main(arguments)
    options = parse_options(arguments)
    println("Oscar $(Base.pkgversion(Oscar)), Julia $VERSION, $(Threads.nthreads()) thread(s), ",
            "$(Sys.cpu_info()[1].model), host $(gethostname())")
    println("options: ", join(["$k=$v" for (k, v) in sort!(collect(options))], " "))
    warm = @elapsed redirect_stdout(devnull) do
        warm_up(options)
    end
    @printf("warm-up on Gr(2,5): %.1f s\n", warm)
    polynomials = generators(options["case"])
    options["order"] == "colex" && (polynomials = colex_order(polynomials))
    options["order"] == "greedy" && (polynomials = greedy_order(polynomials))
    options["mode"] == "steps" && println("order: ", join(kind.(polynomials), " "))
    limit = parse(Int, options["limit"])
    limit > 0 && (polynomials = polynomials[1:limit])
    println("$(length(polynomials)) generators in $(ngens(parent(first(polynomials)))) coordinates: ",
            join(["$(count(==(k), kind.(polynomials))) $k" for k in ("linear", "pluecker", "incidence", "quadric")], ", "))
    flush(stdout)
    stats = @timed with_hull(() -> compute(polynomials, options["mode"]), options["hull"])
    result = with_hull(() -> counts(stats.value), options["hull"])
    @printf("RESULT %s order=%s hull=%s mode=%s: %.1f s, %s, peak %.2f GiB\n", options["case"], options["order"],
            options["hull"], options["mode"], stats.time, result, gib())
    if limit == 0
        expected = EXPECTED[options["case"]]
        (result.rays, result.maximal_cells) == expected ||
            println("MISMATCH: Gfan gives $(expected[1]) rays and $(expected[2]) maximal cones")
    end
end

main(ARGS)
