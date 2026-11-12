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
    status_only::Bool = true
end

Base.@kwdef struct StudyTarget
    name::String
    stages::Set{CensusStage}
    callbacks::Dict{CensusStage, Function}
    artifacts::Dict{CensusStage, Vector{String}} = Dict{CensusStage, Vector{String}}()
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
