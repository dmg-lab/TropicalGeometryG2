@enum CensusStage begin
    tropicalize
    coarsen
    normalize
    trees
    symmetry
    degenerations
    atlas
    projections
    summarize
    zenodo
end

const CENSUS_STAGES = instances(CensusStage)
const CENSUS_STAGE_NAMES = Dict(string(stage) => stage for stage in CENSUS_STAGES)

Base.@kwdef mutable struct CensusOptions
    stages::Vector{CensusStage} = CensusStage[]
    target::String = "all"
    engine::Symbol = :auto
    resume::Bool = false
    force::Bool = false
    no_pdf::Bool = false
    workers::Int = min(2, max(1, Sys.CPU_THREADS - 1))
    threads::Int = min(8, Sys.CPU_THREADS)
    bits::Int = 64
    gfan_chunksize::Int = 0
    fresh_processes::Bool = false
    verify_final::Bool = false
    status_only::Bool = true
end

Base.@kwdef struct StudyTarget
    name::String
    stages::Set{CensusStage}
    callbacks::Dict{CensusStage, Function}
    artifacts::Dict{CensusStage, Vector{String}} = Dict{CensusStage, Vector{String}}()
    # `() -> problem` of the parallel tropicalization that --verify-final rechecks.
    tropical_problem::Union{Nothing, Function} = nothing
end

Base.@kwdef struct StudyCensus
    name::String
    targets::Dict{String, StudyTarget}
end

# Every target keeps all it produces in one data directory: the numbered stage
# artifacts, the coarsening report, and the stage provenance in `provenance/`.

"The provenance record of `stage`, in `provenance/` beside the artifacts in `directory`."
stage_provenance_path(directory::AbstractString, stage) =
    joinpath(directory, "provenance", "$(string(stage)).toml")

"The coarsening report, beside the artifacts in `directory`."
coarsest_report_path(directory::AbstractString) = joinpath(directory, "02_coarsest_report.md")

function parse_stage_list(text::AbstractString)
    names = filter(!isempty, strip.(split(text, ',')))
    isempty(names) && throw(ArgumentError("the stage list is empty"))
    return [get(CENSUS_STAGE_NAMES, name) do
        throw(ArgumentError("unknown census stage: $name"))
    end for name in names]
end

function parse_census_options(arguments)
    options = CensusOptions()
    from_stage = nothing
    through_stage = nothing
    for argument in arguments
        if startswith(argument, "--stage=")
            options.stages = parse_stage_list(split(argument, "=", limit=2)[2])
            options.status_only = false
        elseif startswith(argument, "--from=")
            from_stage = only(parse_stage_list(split(argument, "=", limit=2)[2]))
            options.status_only = false
        elseif startswith(argument, "--through=")
            through_stage = only(parse_stage_list(split(argument, "=", limit=2)[2]))
            options.status_only = false
        elseif startswith(argument, "--target=")
            options.target = split(argument, "=", limit=2)[2]
        elseif startswith(argument, "--engine=")
            options.engine = Symbol(split(argument, "=", limit=2)[2])
        elseif startswith(argument, "--workers=")
            options.workers = parse(Int, split(argument, "=", limit=2)[2])
        elseif startswith(argument, "--threads=")
            options.threads = parse(Int, split(argument, "=", limit=2)[2])
        elseif startswith(argument, "--bits=")
            options.bits = parse(Int, split(argument, "=", limit=2)[2])
        elseif startswith(argument, "--gfan-chunksize=")
            options.gfan_chunksize = parse(Int, split(argument, "=", limit=2)[2])
        elseif argument == "--fresh-processes"
            options.fresh_processes = true
        elseif argument == "--verify-final"
            options.verify_final = true
            options.status_only = false
        elseif argument == "--resume"
            options.resume = true
        elseif argument == "--force"
            options.force = true
        elseif argument == "--no-pdf"
            options.no_pdf = true
        elseif argument in ("--help", "-h")
            options.status_only = true
        else
            throw(ArgumentError("unknown census option: $argument"))
        end
    end
    isempty(options.stages) && !isnothing(from_stage) && begin
        first_index = findfirst(==(from_stage), CENSUS_STAGES)
        last_index = isnothing(through_stage) ? length(CENSUS_STAGES) :
                     findfirst(==(through_stage), CENSUS_STAGES)
        first_index <= last_index || throw(ArgumentError("--from occurs after --through"))
        options.stages = collect(CENSUS_STAGES[first_index:last_index])
    end
    isnothing(from_stage) && !isnothing(through_stage) &&
        throw(ArgumentError("--through requires --from"))
    options.engine in (:auto, :oscar, :parallel) ||
        throw(ArgumentError("engine must be auto, oscar, or parallel"))
    options.bits in (0, 32, 64) || throw(ArgumentError("bits must be 0, 32, or 64"))
    options.gfan_chunksize >= 0 || throw(ArgumentError("gfan chunk size cannot be negative"))
    # Each new process would start the tropicalization again.
    options.fresh_processes && options.force &&
        throw(ArgumentError("--fresh-processes continues with --resume and cannot be combined with --force"))
    if options.verify_final
        # The recheck reads the finished run and writes beside it in work/ only.
        isempty(options.stages) || throw(ArgumentError("--verify-final runs without --stage, --from or --through"))
        options.target == "all" && throw(ArgumentError("--verify-final rechecks one target: give --target=NAME"))
        options.force && throw(ArgumentError("--verify-final reads the finished run and cannot be combined with --force"))
        options.fresh_processes && throw(ArgumentError("--verify-final cannot be combined with --fresh-processes"))
        any(startswith(argument, "--bits=") for argument in arguments) && options.bits != 0 &&
            throw(ArgumentError("--verify-final computes with exact integers, --bits=0"))
        options.bits = 0
    end
    return options
end

function artifact_digest(path)
    return isfile(path) ? bytes2hex(sha256(read(path))) : "missing"
end

function census_configuration(options::CensusOptions)
    return Dict{String, Any}(
        "engine" => string(options.engine),
        "resume" => options.resume,
        "force" => options.force,
        "no_pdf" => options.no_pdf,
        "workers" => options.workers,
        "threads" => options.threads,
        "integer_bits" => options.bits,
        "gfan_chunksize" => options.gfan_chunksize,
    )
end

"Whether a stage keeps its `outputs`: they all exist and `--force` is not given."
function skip_current(outputs, options::CensusOptions)
    current = all(isfile, outputs) && !options.force
    current && println("  using existing canonical artifact(s)")
    return current
end

"""
    fan_metadata(fan; with_f_vector=false, with_quotient_f_vector=false)

The counts and dimensions of `fan` that a stage records in its provenance:
rays modulo lineality, maximal cells, ambient dimension, dimension, lineality
dimension and purity. Optionally also the f-vector (`f_vector`) and the
f-vector modulo lineality from the zero cone on
(`f_vector_modulo_lineality_including_zero`).
"""
function fan_metadata(fan; with_f_vector=false, with_quotient_f_vector=false)
    metadata = Dict{String, Any}(
        "rays" => length(rays_modulo_lineality(fan).rays_modulo_lineality),
        "maximal_cells" => n_maximal_cones(fan),
        "ambient_dimension" => ambient_dim(fan),
        "dimension" => dim(fan),
        "lineality_dimension" => lineality_dim(fan),
        "pure" => is_pure(fan),
    )
    if with_f_vector || with_quotient_f_vector
        counts, lineality = Int.(f_vector(fan)), lineality_dim(fan)
        with_f_vector && (metadata["f_vector"] = counts)
        with_quotient_f_vector && (metadata["f_vector_modulo_lineality_including_zero"] =
            iszero(lineality) ? vcat([1], counts) : counts[lineality:end])
    end
    return metadata
end

function stage_dependencies(stage)
    stage == tropicalize && return String[]
    stage == coarsen && return ["tropicalize"]
    stage == normalize && return ["coarsen"]
    stage == trees && return ["normalize"]
    stage == symmetry && return ["normalize", "trees"]
    stage == degenerations && return ["normalize", "symmetry"]
    stage == atlas && return ["normalize", "trees", "symmetry"]
    stage == projections && return ["normalize", "symmetry"]
    stage == summarize && return ["normalize", "symmetry", "degenerations"]
    stage == zenodo && return ["tropicalize", "coarsen", "normalize", "symmetry"]
    return String[]
end

function software_versions()
    return Dict{String, Any}(
        "julia" => string(VERSION),
        "oscar" => string(Base.pkgversion(Oscar)),
        "tropical_census" => something(
            isnothing(Base.pkgversion(@__MODULE__)) ? nothing : string(Base.pkgversion(@__MODULE__)),
            "development",
        ),
    )
end

function write_provenance(
    path::AbstractString;
    study,
    target,
    stage,
    inputs=String[],
    outputs=String[],
    configuration=Dict{String, Any}(),
    metadata=Dict{String, Any}(),
)
    document = Dict{String, Any}(
        "study" => string(study),
        "target" => string(target),
        "stage" => string(stage),
        "generated_at" => string(now()),
        "software" => software_versions(),
        "stage_dependencies" => stage_dependencies(stage),
        "configuration" => Dict(string(key) => value for (key, value) in configuration),
        "inputs" => [Dict("path" => relpath(file, PACKAGE_ROOT), "sha256" => artifact_digest(file)) for file in inputs],
        "outputs" => [Dict("path" => relpath(file, PACKAGE_ROOT), "sha256" => artifact_digest(file)) for file in outputs],
        "metadata" => Dict(string(key) => value for (key, value) in metadata),
    )
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, document; sorted=true)
    end
    return path
end

"The record of `--verify-final` for `problem`, beside its exact recomputation."
verify_final_record_path(problem::TropicalProblem) =
    joinpath(problem.default_output, "verification_bits_0", "comparison_with_run.toml")

"The fan the finished run in `output` certifies, from its marker `certified_fan.txt`, or nothing."
function certified_run_fan(output)
    marker = joinpath(output, "certified_fan.txt")
    return isfile(marker) ? normpath(joinpath(output, strip(read(marker, String)))) : nothing
end

"""
Whether `--verify-final` has found the fan of the current run of `problem` again:
the record must name the fan the run certifies now, and both fans must be
unchanged since.
"""
function verify_final_state(problem::TropicalProblem)
    record_path = verify_final_record_path(problem)
    isfile(record_path) || return "not run"
    record = TOML.parsefile(record_path)
    output = problem.default_output
    run_fan = certified_run_fan(output)
    unchanged(name) = isfile(joinpath(output, record[name])) &&
        artifact_digest(joinpath(output, record[name])) == record[name * "_sha256"]
    !isnothing(run_fan) && relpath(run_fan, output) == record["run_fan"] && unchanged("run_fan") &&
        unchanged("exact_fan") || return "stale, the run or its recheck has changed"
    return record["same_fan"] ? "same fan as the $(record["run_bits"])-bit run" :
        "DIFFERENT fan from the $(record["run_bits"])-bit run"
end

"""
Whether the gfan fans `a` and `b` of `parse_gfan_fan` are the same fan: the same
lineality space, and the same maximal cones as sets of rays modulo it. Unlike
their texts, this depends neither on the order of rays and cones nor on the
basis Gfan prints for the lineality space, which changes with `--bits`.
"""
function same_gfan_fan(a, b)
    (a.ambient_dimension, a.dimension, a.lineality_dimension) ==
        (b.ambient_dimension, b.dimension, b.lineality_dimension) || return false
    length(a.orth_lineality) == a.ambient_dimension - a.lineality_dimension &&
        length(b.lineality) == b.lineality_dimension ||
        throw(ArgumentError("a gfan fan without a full basis of its lineality space or of its orthogonal complement"))
    # The lineality space of `a` is where its orthogonal rows vanish; that of `b` has the same dimension.
    all(iszero(sum(row .* v)) for v in b.lineality, row in a.orth_lineality) || return false
    # The rays of both fans in the coordinates of `a` modulo the lineality space.
    cones(fan) = Set(Set(gfan_ray_key(a, fan.rays[i + 1]) for i in cone) for cone in fan.maximal_cones)
    return cones(a) == cones(b)
end

"""
    verify_final_tropicalization(problem, options)

The exact recheck of the finished parallel tropicalization of `problem`
(`--verify-final`). Gfan computes the prevariety of the final basis, the
checkpoint `basis_polynomials.txt` of the run, again with exact integers, and
every cone is certified as in the run, all in `verification_bits_0/` beside the
checkpoint. The certified fan must be the one the run certified, up to the
order of its rays and cones and the basis of its lineality space
(`same_gfan_fan`); the result is recorded in
`verify_final_record_path(problem)`. The census artifacts and their provenance
are not touched.
"""
function verify_final_tropicalization(problem::TropicalProblem, options::CensusOptions)
    engine = options.engine == :auto ? problem.default_engine : options.engine
    engine == :parallel || throw(ArgumentError(
        "--verify-final rechecks the parallel engine; $(problem.name) is tropicalized by $engine"))
    output = problem.default_output
    run_fan = certified_run_fan(output)
    isnothing(run_fan) && error("--verify-final needs the finished run in $output, with its marker certified_fan.txt")
    policy = read(joinpath(output, "certified_policy.txt"), String)
    run_bits = parse(Int, match(r"(?m)^GFAN_BITS (\d+)$", policy).captures[1])
    # A recheck that fails leaves no record, rather than the result of an earlier one.
    rm(verify_final_record_path(problem); force=true)
    exact_fan = run_tropicalization(problem, TropicalizationConfig(
        engine=:parallel, output=output, resume=true, verify_final=true, bits=0,
        workers=options.workers, threads=options.threads, gfan_chunksize=options.gfan_chunksize))
    same = read(run_fan) == read(exact_fan) || same_gfan_fan(parse_gfan_fan(run_fan), parse_gfan_fan(exact_fan))
    record = Dict{String, Any}(
        "run_fan" => relpath(run_fan, output), "run_fan_sha256" => artifact_digest(run_fan), "run_bits" => run_bits,
        "exact_fan" => relpath(exact_fan, output), "exact_fan_sha256" => artifact_digest(exact_fan),
        "same_fan" => same, "generated_at" => string(now()), "software" => software_versions(),
    )
    open(verify_final_record_path(problem), "w") do io
        TOML.print(io, record; sorted=true)
    end
    same || error("the exact fan $exact_fan differs from the fan $run_fan of the $run_bits-bit run")
    println("The exact fan is the fan of the $run_bits-bit run: $run_fan")
    return exact_fan
end

function print_census_status(census::StudyCensus, io=stdout)
    println(io, "$(census.name) tropical census")
    println(io, "Choose --stage=<name> or --from=<name> --through=<name>; no computation runs implicitly.")
    println(io, "Stages: ", join(string.(CENSUS_STAGES), ", "))
    for name in sort!(collect(keys(census.targets)))
        target = census.targets[name]
        println(io, "  $name: ", join(string.(sort!(collect(target.stages); by=Int)), ", "))
        for stage in sort!(collect(keys(target.artifacts)); by=Int)
            files = target.artifacts[stage]
            present = count(isfile, files)
            state = present == length(files) ? "ready" :
                    iszero(present) ? "missing" : "partial ($present/$(length(files)))"
            println(io, "    $(string(stage)): $state")
        end
        isnothing(target.tropical_problem) ||
            println(io, "    verify-final: ", verify_final_state(target.tropical_problem()))
    end
end

function run_census(census::StudyCensus, arguments=ARGS)
    options = arguments isa CensusOptions ? arguments : parse_census_options(arguments)
    if options.status_only
        print_census_status(census)
        return nothing
    end
    target_names = options.target == "all" ? sort!(collect(keys(census.targets))) : [options.target]
    all(haskey(census.targets, target) for target in target_names) ||
        throw(ArgumentError("unknown target $(options.target)"))
    if options.verify_final
        target = census.targets[only(target_names)]
        isnothing(target.tropical_problem) && throw(ArgumentError(
            "--verify-final: the target $(target.name) has no parallel tropicalization to recheck"))
        println("[$(census.name)/$(target.name)] verify-final")
        verify_final_tropicalization(target.tropical_problem(), options)
        return nothing
    end
    for target_name in target_names
        target = census.targets[target_name]
        for stage in options.stages
            stage in target.stages || continue
            callback = get(target.callbacks, stage) do
                error("target $target_name advertises $stage but provides no implementation")
            end
            println("[$(census.name)/$target_name] $stage")
            callback(target, options)
        end
    end
    return nothing
end
