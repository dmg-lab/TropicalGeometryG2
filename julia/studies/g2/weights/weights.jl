using Oscar
using ProgressMeter

# Multiplicities (weights) of the maximal cones of the tropical G2-Grassmannian trop(G) ⊂ R^21/R1.
#
# Oscar does not compute them: for an ideal that is not principal, binomial or affine linear, `tropical_variety`
# sets every multiplicity to 1 (`compute_weights_and_construct_tropical_variety`, "TODO: compute multiplicities"
# in src/TropicalGeometry/variety_prime.jl, Oscar 1.8.2). We compute them following Maclagan–Sturmfels,
# Definition 3.4.3: for w in the relative interior of a maximal cone σ,
#
#     mult(σ) = Σ_{P ∈ Ass_min(in_w(I))} mult(P, in_w(I)),     mult(P, in_w(I)) = length((S/Q_P)_P),
#
# where in_w(I) is taken in the Laurent polynomial ring S (i.e. saturated at the product of all variables) and Q_P
# is the P-primary component. Everything is homogeneous, so the length is deg(Q_P)/deg(P) (associativity formula).
#
# Two versions:
#   over QQ   (`mult_QQ`):   minimal primes and lengths of `primary_decomposition` over QQ;
#   over QQbar (`mult_QQbar`): each QQ-prime P splits into d_P conjugate geometric primes, each with the same length
#                             (char 0), so mult = Σ length_P · d_P, with d_P from `absolute_primary_decomposition`.
# MS15 work over an algebraically closed field, so `mult_QQbar` is the multiplicity in the sense of the book. The two
# differ exactly when some minimal prime over QQ is not geometrically prime (e.g. type 21 in main.tex).
#
# The cones are the 251 maximal cones of the raw fan data/g2/grassmannian/01_raw_fan.mrdi (constant initial ideal on
# each relative interior, min convention, coordinates p_12, ..., p_67 as in data/g2/grassmannian/00_ideal.mrdi).
#
# From the `julia` folder:
#     julia --project=. studies/g2/weights/weights.jl [--absolute] [--cones=a:b]
# computes mult_QQ for all maximal cones (with --absolute also mult_QQbar), with a progress bar, appending one row
# per cone to studies/g2/weights/weights.csv (an interrupted run resumes), and prints a summary.

const DATA = joinpath(@__DIR__, "..", "..", "..", "data", "g2", "grassmannian")
const IDEAL_PATH = joinpath(DATA, "00_ideal.mrdi")
const RAW_FAN_PATH = joinpath(DATA, "01_raw_fan.mrdi")
const WEIGHTS_PATH = joinpath(@__DIR__, "weights.csv")
const HEADER = "cone,weight,saturated,minimal_primes_QQ,lengths,mult_QQ,conjugates,mult_QQbar,seconds"

"Primitive integer vector on the ray of a rational vector."
function integer_ray(v)
    r = [BigInt(numerator(QQ(c))) // BigInt(denominator(QQ(c))) for c in v]
    d = reduce(lcm, denominator.(r); init=BigInt(1))
    n = [numerator(c * d) for c in r]
    g = reduce(gcd, n; init=BigInt(0))
    return Int.(iszero(g) ? n : n .÷ g)
end

"The maximal cones of the raw fan, each with an integer weight in its relative interior (the sum of its primitive rays)."
function maximal_cone_weights(path=RAW_FAN_PATH)
    F = load(path)
    return [sum(integer_ray(Vector(r)) for r in rays_modulo_lineality(C)[1]) for C in maximal_cones(F)]
end

"in_w(I) in the Laurent ring: the initial ideal (min convention) saturated at the product of all variables."
function laurent_initial_ideal(I, w)
    inw = initial(I, tropical_semiring_map(QQ, min), QQ.(w))
    R = base_ring(inw)
    sat = saturation(inw, ideal(R, [prod(gens(R))]))
    @assert !is_one(sat) "w = $w is not in the tropical variety (in_w(I) contains a monomial)"
    return (initial=inw, saturation=sat)
end

"Length of (S/Q)_P for a P-primary homogeneous ideal Q: deg(Q)/deg(P)."
function primary_length(Q, P)
    Q == P && return 1
    dQ, dP = degree(Q), degree(P)
    @assert dQ % dP == 0 "degree of the primary component is not a multiple of the degree of its prime"
    return Int(dQ ÷ dP)
end

"Indices of the minimal primes in a list of primes."
minimal_indices(primes) = [i for i in eachindex(primes)
                           if !any(j != i && is_subset(primes[j], primes[i]) for j in eachindex(primes))]

"""
    cone_multiplicity(I, w; absolute=false)

The multiplicity of the maximal cone containing `w` in its relative interior (MS15, Definition 3.4.3): over QQ, and
with `absolute=true` also over QQbar. Returns a named tuple with the minimal primes' lengths and conjugate counts.
"""
function cone_multiplicity(I, w; absolute::Bool=false)
    L = laurent_initial_ideal(I, w)
    S = L.saturation
    if absolute
        components = absolute_primary_decomposition(S)      # (Q, P over QQ, a geometric prime, number of conjugates)
        Ps = [c[2] for c in components]
        keep = minimal_indices(Ps)
        lengths = [primary_length(components[i][1], components[i][2]) for i in keep]
        conjugates = [Int(components[i][4]) for i in keep]
        mult_QQbar = sum(lengths .* conjugates)
    else
        components = primary_decomposition(S)               # (Q, P) over QQ
        Ps = [c[2] for c in components]
        keep = minimal_indices(Ps)
        lengths = [primary_length(components[i][1], components[i][2]) for i in keep]
        conjugates = missing
        mult_QQbar = missing
    end
    return (saturated=L.initial == S, minimal_primes=length(keep), lengths=lengths, mult_QQ=sum(lengths),
            conjugates=conjugates, mult_QQbar=mult_QQbar)
end

"Cones already in the results file."
function completed(path)
    isfile(path) || return Set{Int}()
    return Set(parse(Int, split(l, ",")[1]) for l in Iterators.drop(eachline(path), 1) if !isempty(l))
end

"""
    compute_weights(; absolute=false, cones=nothing, path=WEIGHTS_PATH)

Compute the multiplicities of the maximal cones (all, or the indices in `cones`), skipping cones already in `path`,
and append one row per cone. Shows a progress bar.
"""
function compute_weights(; absolute::Bool=false, cones=nothing, path=WEIGHTS_PATH)
    I = load(IDEAL_PATH)
    weights = maximal_cone_weights()
    indices = cones === nothing ? collect(eachindex(weights)) : collect(cones)
    done = completed(path)
    todo = filter(k -> !(k in done), indices)
    println("$(length(weights)) maximal cones; computing $(length(todo)) (", absolute ? "over QQ and QQbar" : "over QQ", ")")
    isfile(path) || open(io -> println(io, HEADER), path, "w")
    progress = Progress(length(todo); desc="maximal cones: ", showspeed=true)
    for k in todo
        w = weights[k]
        t = @elapsed m = cone_multiplicity(I, w; absolute)
        row = join((k, join(w, ";"), m.saturated, m.minimal_primes, join(m.lengths, ";"), m.mult_QQ,
                    m.conjugates === missing ? "" : join(m.conjugates, ";"),
                    m.mult_QQbar === missing ? "" : m.mult_QQbar, round(t; digits=2)), ",")
        open(io -> println(io, row), path, "a")
        next!(progress; showvalues=[(:cone, k), (:mult_QQ, m.mult_QQ), (:mult_QQbar, something(m.mult_QQbar, "-"))])
    end
    finish!(progress)
    return path
end

"Summary of the results file: how many cones have which multiplicity."
function summarize(path=WEIGHTS_PATH)
    rows = [split(l, ",") for l in Iterators.drop(eachline(path), 1) if !isempty(l)]
    count_by(col) = sort(collect(Dict(v => count(r -> r[col] == v, rows) for v in unique(r[col] for r in rows))))
    println("\n$(length(rows)) maximal cones computed")
    println("mult over QQ:    ", join(("$v: $n cones" for (v, n) in count_by(6)), ", "))
    any(r -> !isempty(r[8]), rows) &&
        println("mult over QQbar: ", join(("$v: $n cones" for (v, n) in count_by(8) if !isempty(v)), ", "))
    println("in_w(I) already saturated: ", count(r -> r[3] == "true", rows), " of $(length(rows))")
    differ = [r[1] for r in rows if !isempty(r[8]) && r[6] != r[8]]
    isempty(differ) || println("cones where QQ and QQbar differ: ", join(differ, ", "))
end

function main(args)
    absolute = "--absolute" in args
    cones = nothing
    for a in args
        if startswith(a, "--cones=")
            a1, b1 = parse.(Int, split(split(a, "=")[2], ":"))
            cones = a1:b1
        end
    end
    summarize(compute_weights(; absolute, cones))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
