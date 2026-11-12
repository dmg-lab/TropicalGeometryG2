#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "intersection_lib.jl"))

const JULIA_DIR = normpath(joinpath(@__DIR__, ".."))
const DATA_DIR = joinpath(JULIA_DIR, "data")

const CONFIGS = Dict(
    7 => (
        ambient_dim=16,
        default_count=10,
        translation_range=-20_000:20_000,
        default_solver=:highs,
        grassmannian=joinpath(DATA_DIR, "mukai_curves", "genus_7_fan.mrdi"),
        linear_space=joinpath(DATA_DIR, "mukai_curves", "genus_7_linear_space.mrdi"),
    ),
    8 => (
        ambient_dim=15,
        default_count=5,
        translation_range=-2_000:2_000,
        default_solver=:highs,
        grassmannian=joinpath(DATA_DIR, "mukai_curves", "genus_8_fan.mrdi"),
        linear_space=joinpath(DATA_DIR, "mukai_curves", "genus_8_linear_space.mrdi"),
    ),
    9 => (
        ambient_dim=14,
        default_count=10,
        translation_range=-200_000:200_000,
        default_solver=:highs,
        grassmannian=joinpath(DATA_DIR, "mukai_curves", "genus_9_fan.mrdi"),
        linear_space=joinpath(DATA_DIR, "mukai_curves", "genus_9_linear_space.mrdi"),
    ),
    10 => (
        ambient_dim=21,
        default_count=2,
        translation_range=-(1 << 20):(1 << 20),
        default_solver=:hybrid,
        grassmannian=joinpath(DATA_DIR, "g2", "grassmannian", "01_variety_complex.mrdi"),
        linear_space=joinpath(DATA_DIR, "fan_uniform_linear_space_16_21.mrdi"),
        # With --shift-basis, a translation is a random integral combination of
        # the rows of this basis, written by studies/g2/curves/shift_basis.jl.
        shift_basis=joinpath(DATA_DIR, "mukai_curves", "genus_10_shift_basis.mrdi"),
        coefficient_range=-(1 << 15):(1 << 15),
    ),
)

function usage(io::IO=stdout)
    println(io, "Usage: julia --project=julia julia/mukai_curves/run_experiments.jl GENUS [COUNT] [SOLVER] [OUTPUT_DIR]")
    println(io, "         [--shift-basis[=FILE]] [--coefficients=LOW:HIGH]")
    println(io, "  GENUS:  7, 8, 9, or 10")
    println(io, "  SOLVER: highs, soplex, hybrid, or certified (hybrid without SoPlex; defaults are genus-specific)")
    println(io, "  --shift-basis: translations are random integral combinations of the rows of FILE,")
    println(io, "                 by default the genus-10 basis, with coefficients in LOW:HIGH")
end

"""
The basis of the translations and the range of their coefficients, or
`nothing` for uniformly random translations, from the options `flags`.
"""
function parse_shifts(flags, config)
    basis_files = [flag == "--shift-basis" ? get(config, :shift_basis, nothing) : split(flag, '='; limit=2)[2]
                   for flag in flags if flag == "--shift-basis" || startswith(flag, "--shift-basis=")]
    ranges = [split(flag, '='; limit=2)[2] for flag in flags if startswith(flag, "--coefficients=")]
    all(flag -> flag == "--shift-basis" || startswith(flag, "--shift-basis=") || startswith(flag, "--coefficients="), flags) ||
        (usage(stderr); error("unknown option among $flags"))
    length(basis_files) <= 1 && length(ranges) <= 1 || error("pass --shift-basis and --coefficients at most once")
    isempty(basis_files) && (isempty(ranges) || error("--coefficients needs --shift-basis"); return nothing)
    file = only(basis_files)
    isnothing(file) && error("there is no default shift basis for this genus; pass --shift-basis=FILE")
    rows = Matrix{QQFieldElem}(load(abspath(file)).basis)
    all(isone ∘ denominator, rows) || error("the shift basis $file is not integral")
    coefficients = get(config, :coefficient_range, -(1 << 15):(1 << 15))
    if !isempty(ranges)
        low, high = parse.(Int, split(only(ranges), ':'; limit=2))
        coefficients = low:high
    end
    return (basis=Int.(numerator.(rows)), coefficients=coefficients, file=abspath(file))
end

function parse_arguments(args)
    flags = filter(arg -> startswith(arg, "--"), args)
    args = filter(arg -> !startswith(arg, "--"), args)
    isempty(args) && (usage(stderr); exit(2))
    genus = tryparse(Int, args[1])
    genus in keys(CONFIGS) || error("GENUS must be one of 7, 8, 9, or 10")
    config = CONFIGS[genus]

    count = length(args) >= 2 ? something(tryparse(Int, args[2]), 0) : config.default_count
    count > 0 || error("COUNT must be a positive integer")

    solver = length(args) >= 3 ? Symbol(lowercase(args[3])) : config.default_solver
    solver in (:highs, :soplex, :hybrid, :certified) || error(
        "SOLVER must be highs, soplex, hybrid, or certified",
    )

    output_dir = length(args) >= 4 ? abspath(args[4]) :
        joinpath(JULIA_DIR, "results", "curves", "genus_$genus")
    length(args) <= 4 || (usage(stderr); error("too many arguments"))
    shifts = parse_shifts(flags, config)
    isnothing(shifts) || size(shifts.basis) == (config.ambient_dim, config.ambient_dim) ||
        error("the shift basis is not a square matrix of size $(config.ambient_dim)")
    return genus, count, solver, output_dir, config, shifts
end

function as_polyhedra(object)
    if applicable(maximal_polyhedra, object)
        return collect(maximal_polyhedra(object))
    end
    return polyhedron.(collect(maximal_cones(object)))
end

function reclaim_experiment_memory!()
    GC.gc(true)
    GC.gc(true)
    if Sys.islinux()
        try
            ccall(:malloc_trim, Cint, (Csize_t,), 0)
        catch error
            @warn "Could not trim the native heap" exception=(error, catch_backtrace())
        end
    end
end

function run_all(genus, count, solver, output_dir, config, shifts=nothing)
    mkpath(output_dir)
    result_path(k) = joinpath(output_dir, "result_$(k).mrdi")
    translation_path(k) = joinpath(output_dir, "translation_$(k).mrdi")
    coefficients_path(k) = joinpath(output_dir, "coefficients_$(k).mrdi")
    graph_path(k) = joinpath(output_dir, "graph_$(k).mrdi")
    polyhedra_path(k) = joinpath(output_dir, "list_polyhedra_$(k).mrdi")

    completed = Dict{Int,ExperimentResult}()
    for k in 1:count
        old = filter(isfile, [joinpath(output_dir, "$(name)_$(k).jls") for name in ("result", "translation", "graph")])
        isempty(old) || isfile(result_path(k)) ||
            error("$(join(old, ", ")) is a checkpoint of an earlier version, written with Julia's Serialization; " *
                  "these runs write .mrdi files and do not resume it: use another OUTPUT_DIR")
        isfile(result_path(k)) || continue
        if !isfile(graph_path(k)) || !isfile(polyhedra_path(k))
            @warn "Ignoring incomplete experiment checkpoint" experiment=k
            continue
        end
        completed[k] = load_experiment_result(result_path(k))
    end

    grassmannian = load(config.grassmannian)
    linear_space = load(config.linear_space)
    tg_polys = as_polyhedra(grassmannian)
    tl_polys = as_polyhedra(linear_space)
    all(ambient_dim(polyhedron) == config.ambient_dim for polyhedron in tl_polys) ||
        error("linear-space input has the wrong ambient dimension")

    isnothing(shifts) || println("Translations: random combinations of the rows of $(shifts.file), ",
                                 "coefficients in $(shifts.coefficients)")
    exact_filter = solver in (:soplex, :hybrid, :certified)
    println("Extracting $(length(tl_polys)) linear-space H-representations...")
    tl_hreps = exact_filter ? [get_hrep_qq(p) for p in tl_polys] :
        [get_hrep_f64(p) for p in tl_polys]
    println("Extracting $(length(tg_polys)) Mukai-fan H-representations...")
    tg_hreps = exact_filter ? [get_hrep_qq(p) for p in tg_polys] :
        [get_hrep_f64(p) for p in tg_polys]

    for k in 1:count
        if haskey(completed, k)
            println("Experiment $k / $count [checkpoint complete; skipping]")
            continue
        end

        translation = if isfile(translation_path(k))
            load(translation_path(k))
        elseif isnothing(shifts)
            value = rand(config.translation_range, config.ambient_dim)
            atomic_oscar_save(translation_path(k), value)
            value
        else
            coefficients = rand(shifts.coefficients, size(shifts.basis, 1))
            value = vec(permutedims(coefficients) * shifts.basis)
            atomic_oscar_save(coefficients_path(k), coefficients)
            atomic_oscar_save(translation_path(k), value)
            value
        end
        translation isa Vector{Int} || error("translation $k has the wrong type")
        length(translation) == config.ambient_dim || error(
            "translation $k has the wrong dimension",
        )

        println("Experiment $k / $count (t = $(translation[1:4])...)")
        try
            completed[k] = run_experiment(
                k,
                translation,
                tl_polys,
                tg_polys,
                tl_hreps,
                tg_hreps;
                output_dir=output_dir,
                feasibility_solver=solver,
            )
        finally
            reclaim_experiment_memory!()
        end
        ordered = [completed[i] for i in sort!(collect(keys(completed)))]
        atomic_oscar_save(joinpath(output_dir, "all_results.mrdi"), experiment_record.(ordered))
    end

    println("\nSummary")
    @printf("%-4s  %-8s  %-8s  %-8s  %-6s\n", "exp", "nodes", "edges", "genus", "target")
    for k in sort!(collect(keys(completed)))
        result = completed[k]
        observed_genus = result.n_edges - result.n_nodes + result.n_components
        @printf(
            "%-4d  %-8d  %-8d  %-8d  %-6s\n",
            k,
            result.n_nodes,
            result.n_edges,
            observed_genus,
            observed_genus == genus ? "yes" : "no",
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_all(parse_arguments(ARGS)...)
end
