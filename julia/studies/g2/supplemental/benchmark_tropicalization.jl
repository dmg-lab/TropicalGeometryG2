#!/usr/bin/env julia

# Timings of the tropicalization of the G2 ideals over QQ outside the census:
# the prevarieties with Gfan's threaded intersection and with Oscar's
# `tropical_prevariety`, and the Grassmannian with Oscar's `tropical_variety`.
# A small example first compiles the same code paths, each repetition starts
# from a new ideal, and nothing is written to data/.
#
#   julia --project=. studies/g2/supplemental/benchmark_tropicalization.jl <case> [repetitions] [gfan threads]
#
# Cases:
#   dressian_gfan                    Gfan on the 42 Grassmannian generators
#   flag_dressian_gfan               Gfan on the 78 flag generators
#   dressian                         tropical_prevariety of the 42 Grassmannian generators
#   flag_dressian                    tropical_prevariety of the 78 flag generators
#   grassmannian                     tropical_variety with its default keywords
#   grassmannian_skip_decomposition  tropical_variety as in the census, without primary decomposition

using Oscar
using TropicalGeometryG2

include(joinpath(@__DIR__, "..", "ideals", "G2GrassmannianIdeal.jl"))
include(joinpath(@__DIR__, "..", "ideals", "FlagG2Ideal.jl"))
using .G2GrassmannianIdeal
using .FlagG2Ideal

const NU = tropical_semiring_map(QQ, min)
const GFAN = something(Sys.which("gfan_tropicalprevariety"), "")
const GFAN_BITS = 64

fan_counts(fan) = (rays=length(rays_modulo_lineality(fan).rays_modulo_lineality), maximal_cells=n_maximal_cones(fan),
                   dimension=dim(fan), lineality_dimension=lineality_dim(fan))
complex_counts(complex) = fan_counts(TropicalGeometryG2.as_polyhedral_fan(complex))

"Seconds of `gfan_tropicalprevariety` alone, without reading its output."
function gfan_seconds(polynomials, threads)
    directory = mktempdir()
    input = joinpath(directory, "input.txt")
    TropicalGeometryG2.write_basis_input(input, string.(symbols(parent(first(polynomials)))), string.(polynomials))
    command = `$GFAN -j $threads --bits $GFAN_BITS`
    return @elapsed run(pipeline(command, stdin=input, stdout=joinpath(directory, "prevariety.gfan"),
                                 stderr=joinpath(directory, "gfan.log")))
end

# Each case returns the computation on a new ideal and how to count its result;
# the Gfan cases also return their generators.
const CASES = Dict(
    "dressian_gfan" => threads -> begin
        generators = build_g2_grassmannian_ideal().generators
        (computation=() -> gfan_prevariety(generators; threads=threads, bits=GFAN_BITS, gfan=GFAN),
         counts=fan_counts, generators=generators)
    end,
    "flag_dressian_gfan" => threads -> begin
        generators = build_full_flag_g2_ideal().generators
        (computation=() -> gfan_prevariety(generators; threads=threads, bits=GFAN_BITS, gfan=GFAN),
         counts=fan_counts, generators=generators)
    end,
    "dressian" => threads -> begin
        generators = build_g2_grassmannian_ideal().generators
        (computation=() -> tropical_prevariety(generators, NU), counts=complex_counts)
    end,
    "flag_dressian" => threads -> begin
        generators = build_full_flag_g2_ideal().generators
        (computation=() -> tropical_prevariety(generators, NU), counts=complex_counts)
    end,
    "grassmannian" => threads -> begin
        ideal_value = build_g2_grassmannian_ideal().ideal
        (computation=() -> tropical_variety(ideal_value, NU),
         counts=components -> [complex_counts(polyhedral_complex(component)) for component in components])
    end,
    "grassmannian_skip_decomposition" => threads -> begin
        ideal_value = build_g2_grassmannian_ideal().ideal
        (computation=() -> tropical_variety(ideal_value, NU; skip_decomposition=true),
         counts=variety -> complex_counts(polyhedral_complex(variety)))
    end,
)

"Compile the prevariety and variety code on the Pluecker ideal of Gr(2,5)."
function warm_up()
    ring, _ = polynomial_ring(QQ, ["p$i$j" for i in 1:4 for j in (i + 1):5])
    generators = forget_grading.(gens(grassmann_pluecker_ideal(ring, 2, 5)))
    gfan_seconds(generators, 1)
    fan_counts(gfan_prevariety(generators; threads=1, bits=GFAN_BITS, gfan=GFAN))
    complex_counts(tropical_prevariety(generators, NU))
    complex_counts(polyhedral_complex(only(tropical_variety(ideal(ring, generators), NU))))
    complex_counts(polyhedral_complex(tropical_variety(ideal(ring, generators), NU; skip_decomposition=true)))
end

function main(arguments)
    name = first(arguments)
    haskey(CASES, name) || error("unknown case $name; choose one of $(join(sort!(collect(keys(CASES))), ", "))")
    repetitions = length(arguments) > 1 ? parse(Int, arguments[2]) : 1
    threads = length(arguments) > 2 ? parse(Int, arguments[3]) : 1
    isempty(GFAN) && error("gfan_tropicalprevariety was not found")
    println("Oscar $(Base.pkgversion(Oscar)), Julia $VERSION with $(Threads.nthreads()) thread, $(Sys.cpu_info()[1].model)")
    warm_up_seconds = @elapsed warm_up()
    println("warm-up on Gr(2,5): $(round(warm_up_seconds; digits=1)) s")
    flush(stdout)
    for repetition in 1:repetitions
        case = CASES[name](threads)
        if haskey(case, :generators)
            println("$name, run $repetition, $threads Gfan thread(s), $GFAN_BITS bits: ",
                    "gfan_tropicalprevariety alone $(round(gfan_seconds(case.generators, threads); digits=2)) s")
        end
        stats = @timed case.computation()
        println("$name, run $repetition: $(round(stats.time; digits=2)) s; result $(case.counts(stats.value)); ",
                "peak memory of the process so far $(round(Sys.maxrss() / 2^30; digits=2)) GiB")
        flush(stdout)
    end
end

main(ARGS)
