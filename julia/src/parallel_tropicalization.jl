#!/usr/bin/env julia

using Distributed
using Random
using SHA

const PARALLEL_SOURCE_DIR = @__DIR__
const PROJECT_DIR = normpath(joinpath(PARALLEL_SOURCE_DIR, ".."))
const TROPICAL_CERTIFICATION_POLICY = "weighted_initial_torus_saturation_v1"

Base.@kwdef struct TropicalProblem
    name::String
    slug::String
    provider_file::String
    provider_module::Symbol
    builder::Symbol
    output_prefix::String
    default_output::String
    published_prevariety::Union{Nothing, String} = nothing
    published_fan::Union{Nothing, String} = nothing
    published_complex::Union{Nothing, String} = nothing
    published_coarsest_fan::Union{Nothing, String} = nothing
    coarsest_provenance::Union{Nothing, String} = nothing
    coarsest_report::Union{Nothing, String} = nothing
    default_engine::Symbol = :parallel
    # Provider names of a symmetry group and of its signed coordinate images.
    # `:taxa`: the whole group of taxon permutations, acting on Pluecker
    # coordinates. `:variables`: gfan-style generators `p` or `p => signs`,
    # sending each variable `k` to `signs[k]` times variable `p[k]`; no images.
    symmetry_candidates::Union{Nothing, Symbol} = nothing
    symmetry_coordinate_images::Union{Nothing, Symbol} = nothing
    symmetry_action::Symbol = :taxa
    certification_strategy::Symbol = :groebner
    purity_verifier::Union{Nothing, Symbol} = nothing
    orbit_restricted_intersection::Bool = false
    whole_interior_pruning::Bool = false
    # Merging the prevariety's cones before the orbit checks only orders the
    # checks under :groebner certification. On a large prevariety, false spares
    # its time and memory by checking the orbits of gfan's cones directly.
    coarsen_prevariety::Bool = true
    # With verified symmetry, gfan computes one cone per orbit (--symmetry),
    # which is expanded before the checks. On the path A-B-C of the hexagon
    # study this cut gfan's prevariety from 71 to 15 minutes.
    gfan_symmetry::Bool = false
    # Gfan's threads list the rays and cones of a prevariety in a different
    # order on every run, and every later numbering follows them. true saves
    # each prevariety in the order of `canonical_gfan_text`.
    canonical_gfan_order::Bool = false
end

"Reject unsupported assumptions before using the balanced-support theorem."
function validate_support_purity_certificate(certificate, ideal_dimension)
    isnothing(certificate) && error("balanced support requires an exact purity certificate")
    for flag in (:verified, :exact, :pure, :torus_saturated)
        hasproperty(certificate, flag) && getproperty(certificate, flag) === true ||
            error("balanced support requires certified $flag=true")
    end
    hasproperty(certificate, :dimension) && certificate.dimension == ideal_dimension ||
        error("purity certificate dimension does not match the defining ideal")
    hasproperty(certificate, :method) && !isempty(string(certificate.method)) ||
        error("purity certificate must identify its exact verification method")
    return true
end

function purity_certificate_signature(certificate)
    isnothing(certificate) && return "none"
    fields = [string(flag) * "=" * string(getproperty(certificate, flag))
              for flag in (:verified, :exact, :pure, :torus_saturated, :dimension, :method)]
    if hasproperty(certificate, :certificate_path)
        path = certificate.certificate_path
        isfile(path) || error("missing exact purity certificate $path")
        push!(fields, "certificate_sha256=" * bytes2hex(sha256(read(path))))
    end
    return bytes2hex(sha256(join(fields, '\n')))
end

"Require a complete finite taxon group before using representative-only checks."
function validate_tropical_permutation_group(permutations)
    isempty(permutations) && error("the symmetry group is empty")
    ntaxa = length(first(permutations))
    keys = Set(Tuple(permutation) for permutation in permutations)
    length(keys) == length(permutations) || error("duplicate symmetry permutations")
    Tuple(1:ntaxa) in keys || error("the symmetry group must contain the identity")
    all(sort(collect(permutation)) == collect(1:ntaxa) for permutation in permutations) ||
        error("invalid taxon permutation")
    for left in permutations, right in permutations
        Tuple(left[right]) in keys || error("symmetry candidates are not a complete group")
    end
    return true
end

"Largest-dimensional representative cones are processed first."
prioritize_tropical_jobs(jobs) = sort(jobs; by=job -> (-job.dimension, job.representative))

Base.@kwdef mutable struct TropicalizationConfig
    engine::Symbol = :parallel
    workers::Int = min(2, max(1, Sys.CPU_THREADS - 1))
    threads::Int = min(8, Sys.CPU_THREADS)
    witness_limit::Int = 8
    max_rounds::Int = typemax(Int)
    output::String = ""
    resume::Bool = false
    prepare_only::Bool = false
    verify_final::Bool = false
    bits::Int = 64
    # Restriction cones per pass of the orbit-restricted gfan; 0 keeps gfan's 1000.
    gfan_chunksize::Int = 0
    # Leave with FRESH_PROCESS_STATUS once the restrictions of a round are
    # written and once the round is checked, for the caller to run the command
    # again with --resume. A process that has expanded and checked a fan of 51
    # million cones keeps 150 GB, and then stalls and leaves gfan too little.
    fresh_processes::Bool = false
    gfan::String = something(Sys.which("gfan_tropicalprevariety"), "")
    force::Bool = false
end

function tropicalization_usage(problem::TropicalProblem, io::IO=stdout)
    println(io, """
    Usage:
      julia --project=.. <problem-runner>.jl [options]

    Problem:
      $(problem.name)

    Options:
      --workers=N          Oscar/Singular worker processes (default: 2)
      --threads=N          gfan prevariety threads (default: up to 8)
      --witness-limit=N    New witnesses collected per round (default: 8)
      --max-rounds=N       Stop safely after N refinement rounds
      --output=PATH        Checkpoint/output directory
      --resume             Continue from output/basis_polynomials.txt
      --prepare-only       Build and export the reduced ideal, then stop
      --bits=N             gfan integer arithmetic: 64, 32, or 0 (exact)
      --gfan-chunksize=N   Restriction cones per orbit-restricted gfan pass (default: 1000)
      --fresh-processes    Leave with status 75 after each phase of a round; run again with --resume
      --verify-final       Independently recompute and certify the final basis
      --gfan=PATH          Path to gfan_tropicalprevariety
      --help               Show this help
    """)
end

function parse_tropicalization_config(problem::TropicalProblem, arguments)
    config = TropicalizationConfig(engine=problem.default_engine, output=problem.default_output)
    for argument in arguments
        if argument == "--help"
            tropicalization_usage(problem)
            exit(0)
        elseif argument == "--resume"
            config.resume = true
        elseif argument == "--prepare-only"
            config.prepare_only = true
        elseif argument == "--verify-final"
            config.verify_final = true
        elseif argument == "--force"
            config.force = true
        elseif startswith(argument, "--engine=")
            config.engine = Symbol(split(argument, "=", limit=2)[2])
        elseif startswith(argument, "--workers=")
            config.workers = parse(Int, split(argument, "=", limit=2)[2])
        elseif startswith(argument, "--threads=")
            config.threads = parse(Int, split(argument, "=", limit=2)[2])
        elseif startswith(argument, "--witness-limit=")
            config.witness_limit = parse(Int, split(argument, "=", limit=2)[2])
        elseif startswith(argument, "--max-rounds=")
            config.max_rounds = parse(Int, split(argument, "=", limit=2)[2])
        elseif startswith(argument, "--output=")
            config.output = abspath(split(argument, "=", limit=2)[2])
        elseif startswith(argument, "--bits=")
            config.bits = parse(Int, split(argument, "=", limit=2)[2])
        elseif startswith(argument, "--gfan-chunksize=")
            config.gfan_chunksize = parse(Int, split(argument, "=", limit=2)[2])
        elseif argument == "--fresh-processes"
            config.fresh_processes = true
        elseif startswith(argument, "--gfan=")
            config.gfan = abspath(split(argument, "=", limit=2)[2])
        else
            throw(ArgumentError("unknown option: $argument"))
        end
    end
    config.workers >= 1 || throw(ArgumentError("--workers must be positive"))
    config.threads >= 1 || throw(ArgumentError("--threads must be positive"))
    config.witness_limit >= 1 || throw(ArgumentError("--witness-limit must be positive"))
    config.max_rounds >= 0 || throw(ArgumentError("--max-rounds cannot be negative"))
    config.bits in (0, 32, 64) || throw(ArgumentError("--bits must be 0, 32, or 64"))
    config.gfan_chunksize >= 0 || throw(ArgumentError("--gfan-chunksize cannot be negative"))
    config.engine in (:auto, :oscar, :parallel) ||
        throw(ArgumentError("--engine must be auto, oscar, or parallel"))
    if config.verify_final
        config.resume || throw(ArgumentError("--verify-final requires --resume"))
        config.bits == 0 || throw(ArgumentError("--verify-final requires --bits=0"))
        config.prepare_only && throw(ArgumentError("--verify-final conflicts with --prepare-only"))
    end
    engine = config.engine == :auto ? problem.default_engine : config.engine
    if engine == :parallel
        isempty(config.gfan) && throw(ArgumentError("gfan_tropicalprevariety was not found"))
        isfile(config.gfan) || throw(ArgumentError("gfan executable does not exist: $(config.gfan)"))
    end
    isfile(problem.provider_file) ||
        throw(ArgumentError("ideal provider does not exist: $(problem.provider_file)"))
    return config
end

function atomic_write(path::String, contents::AbstractString)
    mkpath(dirname(path))
    temporary = path * ".tmp.$(getpid())"
    open(temporary, "w") do io
        write(io, contents)
        flush(io)
    end
    mv(temporary, path; force=true)
end

function section_lines(text::String, section::String, next_section::String)
    lines = split(text, '\n')
    stripped = strip.(lines)
    first_line = findfirst(==(section), stripped)
    isnothing(first_line) && throw(ArgumentError("missing $section section in gfan output"))
    last_line = findnext(==(next_section), stripped, first_line + 1)
    isnothing(last_line) && throw(ArgumentError("missing $next_section section in gfan output"))
    return lines[(first_line + 1):(last_line - 1)]
end

function section_lines(text::String, section::String)
    lines = split(text, '\n')
    stripped = strip.(lines)
    first_line = findfirst(==(section), stripped)
    isnothing(first_line) && throw(ArgumentError("missing $section section in gfan output"))
    last_line = findnext(
        line -> !isnothing(match(r"^[A-Z][A-Z_]*$", line)),
        stripped,
        first_line + 1,
    )
    stop = isnothing(last_line) ? length(lines) : last_line - 1
    return lines[(first_line + 1):stop]
end

function parse_integer_row(line::AbstractString)
    data = strip(first(split(line, '#', limit=2)))
    isempty(data) && return nothing
    return parse.(Int, split(data))
end

function parse_gfan_fan(path::String)
    text = read(path, String)
    occursin("MAXIMAL_CONES", text) || throw(ArgumentError("gfan produced no fan at $path"))

    ray_rows = [row for row in parse_integer_row.(section_lines(text, "RAYS", "N_RAYS")) if !isnothing(row)]
    lineality_rows = [row for row in parse_integer_row.(section_lines(text, "LINEALITY_SPACE", "ORTH_LINEALITY_SPACE")) if !isnothing(row)]
    orth_lineality_next_section = occursin(r"(?:^|\n)F_VECTOR\n", text) ?
        "F_VECTOR" : "MAXIMAL_CONES"
    orth_lineality_rows = [
        row
        for row in parse_integer_row.(section_lines(
            text,
            "ORTH_LINEALITY_SPACE",
            orth_lineality_next_section,
        ))
        if !isnothing(row)
    ]

    maximal_section = if occursin(r"(?:^|\n)MAXIMAL_CONES\n", text)
        "MAXIMAL_CONES"
    elseif occursin(r"(?:^|\n)MAXIMAL_CONES_OF_CLOSURE\n", text)
        # gfan 0.8's threaded prevariety engine uses this name because its
        # internal traversal is through half-open cones.
        "MAXIMAL_CONES_OF_CLOSURE"
    else
        throw(ArgumentError("malformed maximal-cones section"))
    end
    maximal_lines = section_lines(text, maximal_section)

    maximal_cones = Vector{Vector{Int}}()
    for line in maximal_lines
        data = strip(first(split(line, '#', limit=2)))
        isempty(data) && continue
        startswith(data, "{") || continue
        inside = strip(replace(data, r"[{}]" => ""))
        push!(maximal_cones, isempty(inside) ? Int[] : parse.(Int, split(inside)))
    end
    isempty(maximal_cones) && throw(ArgumentError("fan contains no maximal cones"))

    dimension_match = match(r"(?:^|\n)DIM\n(\d+)", text)
    lineality_match = match(r"(?:^|\n)LINEALITY_DIM\n(\d+)", text)
    ambient_match = match(r"(?:^|\n)AMBIENT_DIM\n(\d+)", text)
    any(isnothing, (lineality_match, ambient_match)) &&
        throw(ArgumentError("fan dimension header is incomplete"))

    lineality_dimension = parse(Int, lineality_match.captures[1])
    dimension = if !isnothing(dimension_match)
        parse(Int, dimension_match.captures[1])
    else
        f_vector_data = filter(
            !isempty,
            strip.(first.(split.(section_lines(text, "F_VECTOR"), '#', limit=2))),
        )
        isempty(f_vector_data) && throw(ArgumentError("fan has neither DIM nor F_VECTOR"))
        # F_VECTOR starts at the lineality space and ends at the fan dimension.
        lineality_dimension + length(parse.(Int, split(first(f_vector_data)))) - 1
    end

    return (
        text=text,
        ambient_dimension=parse(Int, ambient_match.captures[1]),
        dimension=dimension,
        lineality_dimension=lineality_dimension,
        rays=ray_rows,
        lineality=lineality_rows,
        orth_lineality=orth_lineality_rows,
        maximal_cones=maximal_cones,
    )
end

"""
Two rays of a gfan fan agree modulo its lineality space iff their pairings with
a basis of its orthogonal complement agree; dividing by the gcd also identifies
positive multiples.
"""
function gfan_ray_key(fan, ray)
    pairing = [sum(row .* ray) for row in fan.orth_lineality]
    return pairing .÷ foldl(gcd, pairing; init=0)
end

"""
The text of the gfan fan `fan` of `parse_gfan_fan` with its rays sorted by
`gfan_ray_key`, and the cones of each dimension in each cone section by their
sorted ray indices. Every other section is kept as it is, so a fan in this
order is unchanged.
"""
function canonical_gfan_text(fan)
    ray_keys = [gfan_ray_key(fan, ray) for ray in fan.rays]
    allunique(ray_keys) || throw(ArgumentError("two gfan rays agree modulo the lineality space"))
    order = sortperm(ray_keys)
    # Gfan numbers the rays from 0.
    new_index = invperm(order) .- 1
    is_header(line) = !isnothing(match(r"^[A-Z][A-Z_]*$", strip(line)))
    data(line) = strip(first(split(line, '#', limit=2)))
    lines = split(fan.text, '\n')
    output = String[]
    i = 1
    while i <= length(lines)
        push!(output, lines[i])
        i += 1
        is_header(lines[i - 1]) || continue
        name = strip(lines[i - 1])
        stop = something(findnext(is_header, lines, i), length(lines) + 1)
        body = lines[i:(stop - 1)]
        i = stop
        if name == "RAYS"
            is_row = line -> !isnothing(parse_integer_row(line))
            rows = [join(fan.rays[old], ' ') * "\t# $(new - 1)" for (new, old) in enumerate(order)]
        elseif name in ("CONES", "MAXIMAL_CONES", "MAXIMAL_CONES_OF_CLOSURE")
            # Gfan lists the cones by dimension, with a comment at the first of each.
            is_row = line -> startswith(data(line), "{")
            dimension = -1
            cones = Tuple{Int, Vector{Int}}[]
            for line in filter(is_row, body)
                comment = match(r"#\s*Dimension\s+(\d+)", line)
                isnothing(comment) || (dimension = parse(Int, only(comment.captures)))
                inside = strip(replace(data(line), r"[{}]" => ""))
                push!(cones, (dimension, sort!([new_index[parse(Int, ray) + 1] for ray in split(inside)])))
            end
            sort!(cones)
            rows = ["{" * join(cone, ' ') * "}" *
                    (d >= 0 && (j == 1 || first(cones[j - 1]) != d) ? "\t# Dimension $d" : "")
                    for (j, (d, cone)) in enumerate(cones)]
        elseif name in ("AMBIENT_DIM", "DIM", "LINEALITY_DIM", "N_RAYS", "LINEALITY_SPACE",
                        "ORTH_LINEALITY_SPACE", "F_VECTOR", "SIMPLICIAL", "PURE")
            append!(output, body)
            continue
        else
            throw(ArgumentError("no canonical order for the gfan section $name"))
        end
        count(is_row, body) == length(rows) || throw(ArgumentError("malformed gfan section $name"))
        # The rows come first, then the blank line before the next section.
        append!(output, rows, filter(!is_row, body))
    end
    return join(output, '\n')
end

"Rewrite the gfan fan at `path` in the order of `canonical_gfan_text`, and parse it."
function canonicalize_gfan_fan(path::String)
    fan = parse_gfan_fan(path)
    text = canonical_gfan_text(fan)
    text == fan.text && return fan
    atomic_write(path, text)
    return parse_gfan_fan(path)
end

function cone_jobs(fan)
    jobs = NamedTuple[]
    for cone in fan.maximal_cones
        cone_rays = [fan.rays[zero_based_ray + 1] for zero_based_ray in cone]
        weight = zeros(Int, fan.ambient_dimension)
        for ray in cone_rays
            weight .+= ray
        end
        push!(jobs, (weight=weight, rays=cone_rays, lineality=fan.lineality))
    end
    return jobs
end

basis_input_text(variable_names, basis_polynomials) = "Q[" * join(variable_names, ",") * "]\n" *
    "{\n" * join(("  " * f for f in basis_polynomials), ",\n") * "\n}\n"

write_basis_input(path, variable_names, basis_polynomials) =
    atomic_write(path, basis_input_text(variable_names, basis_polynomials))

function next_round_number(output_directory::String, resume::Bool)
    resume || return 1
    isdir(output_directory) || return 1
    completed_rounds = Int[]
    for entry in readdir(output_directory)
        match_result = match(r"^round_(\d+)$", entry)
        # A directory is committed only when its witness file has been written
        # atomically after all cone jobs for that round have returned.
        witness_path = joinpath(output_directory, entry, "new_witnesses.txt")
        if !isnothing(match_result) && isfile(witness_path)
            push!(completed_rounds, parse(Int, match_result.captures[1]))
        end
    end
    return isempty(completed_rounds) ? 1 : maximum(completed_rounds) + 1
end

"Identify the original round of an uncertified completion checkpoint."
function checkpoint_round_number(output_directory::String, fan_path::String)
    directory = dirname(abspath(fan_path))
    basename(directory) in ("orbit_restricted", "interior_pruning") && (directory = dirname(directory))
    dirname(directory) == abspath(output_directory) || return nothing
    basename(fan_path) in ("prevariety.gfan", "certified_support.gfan", "expanded_prevariety.gfan") || return nothing
    matched = match(r"^round_(\d+)$", basename(directory))
    return isnothing(matched) ? nothing : parse(Int, only(matched.captures))
end

"""
Run Gfan's threaded prevariety on `input_path` and save the fan at `fan_path`.
Its threads list the rays and cones in a different order on every run; with
`canonical=true` the fan is saved in the order of `canonical_gfan_text`.
"""
function run_prevariety(config, input_path, fan_path, log_path; symmetry::Bool=false,
                        canonical::Bool=false)
    # A fan becomes a checkpoint only after Gfan succeeds and its output
    # parses. Remove any older target before pairing a new basis with it.
    partial_path = fan_path * ".partial"
    if isfile(fan_path)
        suffix = prevariety_checkpoint_is_parseable(fan_path) ? ".previous" : ".invalid"
        mv(fan_path, fan_path * suffix; force=true)
    end
    # With --symmetry, Gfan reads permutation generators after the polynomials
    # and outputs only one maximal cone per orbit.
    symmetry_flag = symmetry ? ["--symmetry"] : String[]
    command = `$(config.gfan) -j $(config.threads) --bits $(config.bits) $symmetry_flag`
    println("  threaded prevariety: ", command)
    flush(stdout)
    exit_code = open(input_path, "r") do input
        open(partial_path, "w") do output
            open(log_path, "w") do errors
                process = run_reporting(pipeline(ignorestatus(command), stdin=input, stdout=output, stderr=errors), "gfan")
                process.exitcode
            end
        end
    end
    if exit_code != 0 && config.bits != 0
        println(
            "  $(config.bits)-bit intersection failed; " *
            "retrying with arbitrary-precision integers",
        )
        command = `$(config.gfan) -j $(config.threads) --bits 0 $symmetry_flag`
        exit_code = open(input_path, "r") do input
            open(partial_path, "w") do output
                open(log_path, "a") do errors
                    process = run_reporting(pipeline(ignorestatus(command), stdin=input, stdout=output, stderr=errors), "gfan")
                    process.exitcode
                end
            end
        end
    end
    exit_code == 0 || error("gfan failed; see $log_path and $partial_path")
    fan = canonical ? canonicalize_gfan_fan(partial_path) : parse_gfan_fan(partial_path)
    prevariety_checkpoint_is_parseable(partial_path) || error("gfan produced an invalid fan; retained $partial_path")
    mv(partial_path, fan_path; force=true)
    return fan
end

"Permutation vectors of `1:nvars` from a `PermGroup` or a list of generators."
function gfan_symmetry_generators(symmetry, nvars::Int)
    symmetry isa PermGroup && (symmetry = gens(symmetry))
    generators = [
        generator isa PermGroupElem ? [i^generator for i in 1:nvars] : Int.(collect(generator))
        for generator in symmetry
    ]
    for generator in generators
        sort(generator) == 1:nvars ||
            throw(ArgumentError("$generator is not a permutation of the $nvars variables"))
    end
    return generators
end

"""
Expand the orbit representatives output by `gfan --symmetry` to all maximal
cones, keeping Gfan's 0-based cone indices.
"""
function expand_gfan_orbits(fan, generators)
    rays = Vector{Int}[]
    index = Dict{Vector{Int}, Int}()
    ray_index(ray) = get!(index, gfan_ray_key(fan, ray)) do
        push!(rays, ray)
        length(rays) - 1
    end
    queue = [sort!([ray_index(fan.rays[ray + 1]) for ray in cone]) for cone in fan.maximal_cones]
    cones = Set(queue)
    # The group is finite, so closing the orbits under the generators suffices.
    while !isempty(queue)
        cone = pop!(queue)
        for permutation in generators
            image = sort!([ray_index(rays[ray + 1][permutation]) for ray in cone])
            image in cones && continue
            push!(cones, image)
            push!(queue, image)
        end
    end
    # Gfan calls a representative maximal when no other representative
    # contains it, but an image of another one may. In a fan, a cone lies in
    # another exactly when its rays do.
    cones = sort!(collect(cones))
    containing = Dict{Int, Vector{Int}}()
    for (i, cone) in enumerate(cones), ray in cone
        push!(get!(containing, ray, Int[]), i)
    end
    maximal = [cone for cone in cones if isempty(cone) ? length(cones) == 1 :
               !any(j -> length(cones[j]) > length(cone) && issubset(cone, cones[j]), containing[first(cone)])]
    return merge(fan, (rays=rays, maximal_cones=maximal))
end

"""
Run gfan with `--symmetry` for the verified group of `export_data`, then expand
the orbit representatives and save the whole fan at `fan_path`. Gfan's orbit
output goes to a separate file, so an interrupted run never leaves a partial
fan under the checkpoint name. With `canonical=true` the whole fan is saved in
the order of `canonical_gfan_text`.
"""
function run_symmetric_prevariety(config, input_path, fan_path, log_path, export_data;
                                  canonical::Bool=false)
    n = length(export_data.variable_names)
    generators = [symmetry_coordinate_permutation(permutation, n, export_data.symmetry_action)
                  for permutation in export_data.symmetry_permutations]
    # The basis file stays a plain list of polynomials for the checkpoint checks.
    symmetric_input = input_path * ".symmetric"
    open(symmetric_input, "w") do io
        write(io, read(input_path, String))
        # Gfan numbers the variables from 0.
        println(io, "{", join(("(" * join(p .- 1, ",") * ")" for p in generators), ","), "}")
    end
    orbit_path = joinpath(dirname(fan_path), "prevariety_orbits.gfan")
    orbits = run_prevariety(config, symmetric_input, orbit_path, log_path; symmetry=true)
    fan = expand_gfan_orbits(orbits, generators)
    row_matrix(rows) = isempty(rows) ? zeros(Int, 0, fan.ambient_dimension) : permutedims(reduce(hcat, rows))
    expanded = polyhedral_fan(IncidenceMatrix([cone .+ 1 for cone in fan.maximal_cones]),
                              row_matrix(fan.rays), row_matrix(fan.lineality); non_redundant=true)
    atomic_write(fan_path, _gfan_fan_text(expanded, fan))
    println("  expanded $(length(orbits.maximal_cones)) orbit representatives to $(length(fan.maximal_cones)) maximal cones")
    flush(stdout)
    return canonical ? canonicalize_gfan_fan(fan_path) : parse_gfan_fan(fan_path)
end

"""
    gfan_prevariety(polynomials; threads=Sys.CPU_THREADS, bits=64, symmetry=nothing, gfan=...,
                    non_redundant=true)

Tropical prevariety of `polynomials` over QQ with trivial valuation, in the min
convention of `tropical_semiring_map(QQ, min)`, as an Oscar polyhedral fan.
Same result as `tropical_prevariety`, but computed by Gfan's threaded
intersection, which is orders of magnitude faster than Oscar's pairwise
common refinement.

`symmetry` is a `PermGroup` or a list of generators, each a `PermGroupElem` or
a permutation vector of `1:nvars`, of a group permuting the variables. Gfan
then computes one maximal cone per orbit (`--symmetry`), and the orbits are
expanded so that the full fan is returned. The group must permute the tropical
hypersurfaces of `polynomials`; this is not checked.

With `non_redundant=true`, polymake takes Gfan's rays and maximal cones as
they are. Recomputing them took 15 minutes on a prevariety of 37,986 cones.
This is only valid when the rays and the incidence matrix really describe a
fan: the rays distinct modulo the lineality space, each cone listed by exactly
its extreme rays, no cone inside another, and any two cones meeting in a
common face. Gfan's output, expanded or not, satisfies this; polymake does not
check it. `non_redundant=false` lets polymake remove redundant rays and cones.
"""
function gfan_prevariety(
    polynomials::AbstractVector{<:MPolyRingElem};
    threads::Int=Sys.CPU_THREADS,
    bits::Int=64,
    symmetry=nothing,
    gfan::String=something(Sys.which("gfan_tropicalprevariety"), ""),
    non_redundant::Bool=true,
)
    isempty(gfan) && throw(ArgumentError("gfan_tropicalprevariety was not found"))
    isempty(polynomials) && throw(ArgumentError("empty list of polynomials"))
    variable_names = string.(symbols(parent(first(polynomials))))
    generators = isnothing(symmetry) ? Vector{Int}[] :
        gfan_symmetry_generators(symmetry, length(variable_names))
    # Removed when Julia exits, so the log named in a Gfan error stays readable.
    directory = mktempdir()
    input_path = joinpath(directory, "input.txt")
    write_basis_input(input_path, variable_names, string.(polynomials))
    # Gfan outputs nothing for an empty list of generators.
    if !isempty(generators)
        # Gfan numbers the variables from 0.
        open(input_path, "a") do io
            println(io, "{", join(("(" * join(p .- 1, ",") * ")" for p in generators), ","), "}")
        end
    end
    config = (gfan=gfan, threads=threads, bits=bits)
    fan = run_prevariety(config, input_path, joinpath(directory, "prevariety.gfan"), joinpath(directory, "gfan.log");
                         symmetry=!isempty(generators))
    isempty(generators) || (fan = expand_gfan_orbits(fan, generators))
    # An empty Vector{Union{}} of rows would dispatch to Oscar's PointVector
    # method, so give empty ray and lineality lists their ambient width.
    row_matrix(rows) = isempty(rows) ? zeros(Int, 0, fan.ambient_dimension) : permutedims(reduce(hcat, rows))
    # Gfan uses the max convention; negate rays to match the min convention.
    return polyhedral_fan(
        IncidenceMatrix([cone .+ 1 for cone in fan.maximal_cones]),
        -row_matrix(fan.rays),
        row_matrix(fan.lineality);
        non_redundant=non_redundant,
    )
end

function initialize_workers(count::Int, problem::TropicalProblem)
    # Distributed treats a String as one argument. Keep the project path and
    # startup flag as separate arguments, including when the path has spaces.
    exeflags = `--project=$(PROJECT_DIR) --startup-file=no`
    worker_ids = addprocs(count; exeflags=exeflags)
    worker_file = joinpath(PARALLEL_SOURCE_DIR, "parallel_worker.jl")
    remote_include(worker, path) = remotecall_fetch(
        Core.eval,
        worker,
        Main,
        :(Base.include(Main, $path); nothing),
    )
    @sync for worker in worker_ids
        @async begin
            # Register the same package identity as the parent before sending
            # closures defined in it. Including the module source directly
            # creates an unregistered Main module that cannot deserialize
            # package-owned closures.
            remotecall_fetch(Core.eval, worker, Main, :(using TropicalCensus; nothing))
            remote_include(worker, worker_file)
            # Provider loading is a separate remote call so Julia 1.12 sees
            # its bindings in the worker's current world age.
            remote_include(worker, problem.provider_file)
            remotecall_fetch(
                (provider_file, provider_module, builder, candidates, coordinate_images, strategy, purity_verifier,
                 symmetry_action) ->
                    Main.initialize_parallel_trop_worker(
                        provider_file,
                        provider_module,
                        builder,
                        candidates,
                        coordinate_images,
                        strategy,
                        purity_verifier,
                        symmetry_action,
                    ),
                worker,
                problem.provider_file,
                problem.provider_module,
                problem.builder,
                problem.symmetry_candidates,
                problem.symmetry_coordinate_images,
                problem.certification_strategy,
                problem.purity_verifier,
                problem.symmetry_action,
            )
        end
    end
    return worker_ids
end

"The workers of a round, started again when its gfan ran without them."
function restored_workers(worker_ids, config::TropicalizationConfig, problem::TropicalProblem)
    issubset(worker_ids, workers()) && return worker_ids
    rmprocs(intersect(worker_ids, workers()))
    return timed(() -> initialize_workers(config.workers, problem), "Starting $(config.workers) workers")
end

"Exit status of a run with `fresh_processes` that has more to do: run it again with --resume."
const FRESH_PROCESS_STATUS = 75

"End this process after a phase whose results are on disk; its workers end with it."
function leave_for_fresh_process(done)
    println("  ", done, " Leaving with status $FRESH_PROCESS_STATUS for a fresh process: run again with --resume.")
    flush(stdout)
    exit(FRESH_PROCESS_STATUS)
end

function collect_witnesses(jobs, basis_set, config; membership_only::Bool=false, checks_path=nothing)
    new_witnesses = String[]
    checked = 0
    valid = 0
    failures = 0
    chunk_size = max(config.workers, 2 * config.workers)
    progress = TaskProgress(membership_only ? "support-membership checks" : "whole-cone checks", length(jobs))

    for first_index in 1:chunk_size:length(jobs)
        last_index = min(length(jobs), first_index + chunk_size - 1)
        chunk = jobs[first_index:last_index]
        # The worker code is deliberately not loaded in the coordinator (that
        # would keep another large Oscar instance in memory).  Resolve the
        # function in each worker's Main module.
        results = pmap(
            job -> membership_only ? Main.parallel_trop_support_probe(job) :
                Main.parallel_trop_cone_certificate(job),
            chunk;
            batch_size=1,
        )
        for result in results
            checked += 1
            job = jobs[checked]
            if result.valid
                valid += 1
            else
                failures += 1
                orbit = hasproperty(result, :witness_orbit) ? result.witness_orbit : [result.witness]
                for witness in orbit
                    if !(witness in basis_set) && !(witness in new_witnesses)
                        push!(new_witnesses, witness)
                    end
                end
            end
            if !isnothing(checks_path)
                open(checks_path, "a") do io
                    println(io, join([
                        membership_only ? "support_membership" : "whole_cone",
                        hasproperty(job, :representative) ? job.representative : checked,
                        hasproperty(job, :orbit_size) ? job.orbit_size : 1,
                        hasproperty(job, :dimension) ? job.dimension : "",
                        result.valid, result.reason, result.basis_size, result.seconds,
                        hasproperty(result, :sample_index) ? result.sample_index : "",
                        hasproperty(result, :probe_weights) ?
                            join((join(weight, ' ') for weight in result.probe_weights), ';') : join(job.weight, ' '),
                    ], ','))
                end
            end
            # A large intermediate prevariety can contain thousands of cones.
            # Report every failure; the progress lines count the others.
            if !result.valid
                note!(progress, "cone $checked/$(length(jobs)): refine ($(result.reason)), " *
                      "GB=$(result.basis_size), $(round(result.seconds; digits=2))s")
            end
            advance!(progress)
        end
        !membership_only && length(new_witnesses) >= config.witness_limit && break
    end
    finish!(progress)
    return (witnesses=new_witnesses, checked=checked, valid=valid, failures=failures)
end

"One job per cell orbit of a group of taxon (`action=:taxa`) or variable permutations."
function tree_orbit_jobs(fan, permutations; action::Symbol=:taxa)
    validate_tropical_permutation_group(permutations)
    ntaxa = length(first(permutations))
    ncoordinates = ambient_dim(fan)
    symmetry_coordinate_permutation(first(permutations), ncoordinates, action)
    normalized = timed("Normalizing the rays of $(_count(n_maximal_cones(fan))) cones") do
        if action == :variables || ncoordinates == ntaxa + binomial(ntaxa, 2)
            # Flag rays have an independent point block, and variable symmetries
            # have no taxa. Primitive integral rows preserve the complete cones
            # and actual lineality without choosing a tree gauge.
            data = rays_modulo_lineality(fan)
            rows = [primitive_integer(collect(row)) for row in data.rays_modulo_lineality]
            rays = isempty(rows) ? zero_matrix(QQ, 0, ncoordinates) : rows_matrix(rows; ring=QQ)
            lineality = matrix(QQ, data.lineality_basis)
            integral_fan = polyhedral_fan(maximal_cones(IncidenceMatrix, fan), rays, lineality;
                                           non_redundant=true)
            (fan=integral_fan, rays=rays)
        else
            normalize_fan(fan, TreeGauge(ntaxa), collect(1:ncoordinates))
        end
    end
    lineality = Matrix{QQFieldElem}(matrix(QQ, rays_modulo_lineality(normalized.fan).lineality_basis))
    return cell_orbit_jobs(merge(normalized, (cells=fan_cells(normalized.fan), lineality=lineality)),
                           permutations; action=action)
end

"""
Whether the rounds read the cells of a prevariety from its gfan file, without
a Polymake fan: variable symmetries need no tree gauge, and only coarsening
needs the fan.
"""
reads_gfan_cells(export_data, coarsen::Bool) = !coarsen && export_data.symmetry_action == :variables

"""
One job per cell orbit of the fan of a gfan file under variable symmetries, in
the provider's unreduced coordinates. The same as `tree_orbit_jobs` of its
`save_lifted_fan`, whose rays, lineality and cells are those of the file:
gfan's rays are primitive. Polymake keeps 0.5 kB per cone for that fan and
takes 3 kB more each time its cells are read back, which it never returns to
the system: 13 GB and 80 GB for the 26 million cones of the path A-B-C.
"""
function gfan_orbit_jobs(fan_path, permutations)
    validate_tropical_permutation_group(permutations)
    normalized = timed("Reading the cells of the fan") do
        fan = parse_gfan_fan(fan_path)
        symmetry_coordinate_permutation(first(permutations), fan.ambient_dimension, :variables)
        row_matrix(rows) = isempty(rows) ? zero_matrix(QQ, 0, fan.ambient_dimension) : rows_matrix(rows; ring=QQ)
        # gfan uses max; the surrounding project uses Oscar's min convention.
        (rays=row_matrix([primitive_integer(-ray) for ray in fan.rays]),
         cells=[sort!(cone .+ 1) for cone in fan.maximal_cones],
         lineality=Matrix{QQFieldElem}(row_matrix(fan.lineality)))
    end
    return cell_orbit_jobs(normalized, permutations; action=:variables)
end

"""
One job per orbit of the cells of `normalized`: the ray indices `cells` of its
maximal cones, its normalized `rays` and a basis of its `lineality` space.
"""
function cell_orbit_jobs(normalized, permutations; action::Symbol)
    cells, lineality = normalized.cells, normalized.lineality
    symmetry = taxon_symmetry(normalized.rays, cells; candidates=permutations, lineality=lineality, action=action)
    validate_tropical_permutation_group(symmetry.permutations)
    rays = Matrix{QQFieldElem}(normalized.rays)
    lineality_rows = [Int.(primitive_integer(row)) for row in eachrow(lineality)]
    progress = TaskProgress("orbit representatives", length(symmetry.cell_orbits))
    jobs = NamedTuple[]
    for orbit in symmetry.cell_orbits
        representative = first(orbit)
        cone_rays = [Int.(-rays[index, :]) for index in cells[representative]]
        weight = isempty(cone_rays) ? zeros(Int, size(rays, 2)) : reduce(+, cone_rays)
        dimension = rank(matrix(QQ, vcat(lineality, rays[cells[representative], :])))
        push!(jobs, (weight=weight, rays=cone_rays, lineality=lineality_rows,
                     representative=representative, orbit_size=length(orbit), dimension=dimension))
        advance!(progress)
    end
    return (jobs=prioritize_tropical_jobs(jobs), normalized=normalized, symmetry=symmetry)
end

function prepare_symmetric_round(fan_path, directory, export_data; strategy::Symbol=:groebner,
                                 coarsen::Bool=true)
    export_data.lift_indices == collect(1:length(export_data.variable_names)) ||
        error("symmetry requires the provider's unreduced coordinates")
    coarsen || strategy == :groebner || error("balanced support needs the coarsened prevariety")
    plain = reads_gfan_cells(export_data, coarsen)
    fine_path = joinpath(directory, "fine_tree_prevariety.mrdi")
    fine_fan = plain ? nothing : timed("Saving the candidate fan") do
        save_lifted_fan(fan_path, fine_path, export_data.lift_indices).fan
    end
    coarse_path = joinpath(directory, "coarsened_prevariety.mrdi")
    source_marker = coarse_path * ".source_sha256"
    source_sha = bytes2hex(sha256(read(fan_path)))
    coarse_fan = if !coarsen
        fine_fan
    elseif isfile(coarse_path) && isfile(source_marker) && strip(read(source_marker, String)) == source_sha
        load(coarse_path)
    else
        println("  Coarsening the candidate prevariety support before tree-orbit checks...")
        flush(stdout)
        result = coarsen_prevariety_fan(fine_fan; strategy=EqualSpanCoarsening())
        result.support_preserved || error("prevariety coarsening did not certify support preservation")
        save(coarse_path, result.fan)
        open(joinpath(directory, "coarsening_summary.toml"), "w") do io
            TOML.print(io, Dict(
                "source_fan_sha256" => source_sha,
                "support_preserved" => result.support_preserved,
                "coarsest_certified" => result.coarsest_certified,
                "open_equal_span_walls" => result.open_equal_span_walls,
                "input_cones" => result.input_cones,
                "output_cones" => result.output_cones,
                "rolled_back_merges" => result.rolled_back_merges,
                "removed_redundant_cones" => result.removed_redundant_cones,
                "strata" => [Dict("dimension" => stratum.dimension,
                                  "input_cones" => stratum.input_cones,
                                  "output_cones" => stratum.output_cones)
                             for stratum in result.strata],
            ))
        end
        open(joinpath(directory, "coarsened_prevariety_cells.csv"), "w") do io
            println(io, "coarsened_cell,original_fine_cells")
            for (cell, originals) in enumerate(result.original_cells)
                println(io, "$cell,$(join(originals, ';'))")
            end
        end
        atomic_write(source_marker, source_sha * "\n")
        result.fan
    end
    coarsen || (println("  Computing the orbits of the uncoarsened prevariety's cones..."); flush(stdout))
    coarse = plain ? gfan_orbit_jobs(fan_path, export_data.symmetry_permutations) :
        tree_orbit_jobs(coarse_fan, export_data.symmetry_permutations; action=export_data.symmetry_action)
    timed("Saving the orbits") do
        plain || save(joinpath(directory, "coarsened_normalized_tree_fan.mrdi"), coarse.normalized.fan)
        save(joinpath(directory, "coarsened_normalized_tree_rays.mrdi"), coarse.normalized.rays)
        save_symmetry_artifacts(joinpath(directory, "coarsened_symmetry"), coarse.symmetry)
    end
    high_jobs = filter(job -> job.dimension > export_data.ideal_dimension, coarse.jobs)
    atomic_write(joinpath(directory, "orbit_checks.csv"),
        "phase,representative,orbit_size,dimension,valid,reason,groebner_generators,seconds,sample_index,tested_weights\n")
    open(joinpath(directory, "symmetry_summary.toml"), "w") do io
        TOML.print(io, Dict(
            "source_fan_sha256" => source_sha,
            "ideal_dimension" => export_data.ideal_dimension,
            "fine_maximal_cones" => plain ? length(coarse.normalized.cells) : n_maximal_cones(fine_fan),
            "prevariety_coarsened" => coarsen,
            "coarsened_maximal_cones" => length(coarse.normalized.cells),
            "coarsened_group_order" => length(coarse.symmetry.permutations),
            "coarsened_cell_orbits" => length(coarse.jobs),
            "high_dimensional_orbits" => length(high_jobs),
            "high_dimensional_cells" => sum(job.orbit_size for job in high_jobs; init=0),
            "completion_policy" => strategy == :balanced_support ?
                "exact interior membership of expected-dimensional support-cell orbits, gated by exact purity" :
                "whole-cone certification of every fine-fan orbit representative",
        ))
    end
    println("  $(length(coarse.jobs)) coarsened tree-cell orbits; $(length(high_jobs)) above ideal dimension $(export_data.ideal_dimension)")
    flush(stdout)
    # Without coarsening, the fine orbits are already computed.
    return (fine_fan=fine_fan, coarse=coarse, high_jobs=high_jobs, fine=coarsen ? nothing : coarse)
end

"Remove only proven empty relative interiors, retaining every proper face."
function collect_whole_interior_pruning(prepared, export_data, basis_set, config, directory,
                                        source_fan_path)
    validate_support_purity_certificate(export_data.purity_certificate, export_data.ideal_dimension)
    pruning_directory = joinpath(directory, "interior_pruning")
    mkpath(pruning_directory)
    jobs = prepared.high_jobs
    checks_path = joinpath(pruning_directory, "whole_interior_checks.csv")
    atomic_write(checks_path, "representative,orbit_size,dimension,excluded,reason,groebner_generators,seconds,weight\n")
    results = NamedTuple[]
    chunk_size = max(config.workers, 2 * config.workers)
    for first_index in 1:chunk_size:length(jobs)
        chunk = jobs[first_index:min(length(jobs), first_index + chunk_size - 1)]
        chunk_results = pmap(job -> Main.parallel_trop_whole_interior_exclusion(job), chunk; batch_size=1)
        append!(results, chunk_results)
        for (job, result) in zip(chunk, chunk_results)
            open(checks_path, "a") do io
                println(io, join([job.representative, job.orbit_size, job.dimension,
                                  result.excluded, result.reason, result.basis_size, result.seconds,
                                  join(job.weight, ' ')], ','))
            end
            println("  whole-interior orbit $(job.representative): ",
                    result.excluded ? "excluded by exact term dominance" : "no exclusion proof",
                    ", GB=$(result.basis_size), ", round(result.seconds; digits=2), "s")
            flush(stdout)
        end
    end
    all(result.excluded for result in results) || return nothing
    # Codimension one is essential: preserving all proper faces then gives
    # a dimension-at-most-d candidate still containing the tropical support.
    all(job.dimension == export_data.ideal_dimension + 1 for job in jobs) ||
        error("whole-interior pruning requires codimension-one excess cells")
    constraints = sort!(unique(vcat(collect(basis_set),
        reduce(vcat, [result.witness_orbit for result in results]))))
    constraints = only(pmap(polynomials -> Main.parallel_trop_symmetrize_polynomials(polynomials),
                             [constraints]; batch_size=1))
    constraints_path = joinpath(pruning_directory, "exclusion_constraints.gfan")
    write_basis_input(constraints_path, export_data.variable_names, constraints)
    skeleton = prevariety_dimension_skeleton(prepared.coarse.normalized.fan, export_data.ideal_dimension)
    skeleton.skeleton_certified && skeleton.support_included || error("dimension skeleton was not certified")
    skeleton_min_path = joinpath(pruning_directory, "dimension_skeleton_min.mrdi")
    save(skeleton_min_path, skeleton.fan)
    coarse = tree_orbit_jobs(skeleton.fan, export_data.symmetry_permutations; action=export_data.symmetry_action)
    all(job.dimension <= export_data.ideal_dimension for job in coarse.jobs) ||
        error("dimension skeleton contains an oversized cell")
    save(joinpath(pruning_directory, "coarsened_prevariety.mrdi"), skeleton.fan)
    save(joinpath(pruning_directory, "coarsened_normalized_tree_fan.mrdi"), coarse.normalized.fan)
    save(joinpath(pruning_directory, "coarsened_normalized_tree_rays.mrdi"), coarse.normalized.rays)
    save_symmetry_artifacts(joinpath(pruning_directory, "coarsened_symmetry"), coarse.symmetry)
    # Serialize the full candidate with the same MIN/MAX change as selected
    # support serialization. Every candidate maximal cell is retained here.
    skeleton_path = save_certified_support_fan(coarse.normalized,
        collect(1:n_maximal_cones(coarse.normalized.fan)), source_fan_path, pruning_directory)
    mv(skeleton_path, joinpath(pruning_directory, "dimension_skeleton.gfan"); force=true)
    skeleton_path = joinpath(pruning_directory, "dimension_skeleton.gfan")
    mv(joinpath(pruning_directory, "certified_support_min.mrdi"),
       joinpath(pruning_directory, "dimension_skeleton_normalized_min.mrdi"); force=true)
    mv(joinpath(pruning_directory, "certified_support_max.mrdi"),
       joinpath(pruning_directory, "dimension_skeleton_max.mrdi"); force=true)
    proof_records = [Dict(
        "representative" => job.representative, "orbit_size" => job.orbit_size,
        "dimension" => job.dimension, "weight" => job.weight,
        "rays_max" => job.rays, "lineality" => job.lineality,
        "witness" => result.witness, "witness_in_defining_ideal" => true,
        "dominant_term_index" => result.proof.dominant_index,
        "dominant_exponent" => string.(result.proof.dominant_exponent),
        "comparisons" => [Dict("competitor_index" => comparison.competitor_index,
            "competitor_exponent" => string.(comparison.competitor_exponent),
            "ray_differences" => string.(comparison.ray_differences),
            "lineality_differences" => string.(comparison.lineality_differences))
            for comparison in result.proof.comparisons])
        for (job, result) in zip(jobs, results)]
    proof_path = joinpath(pruning_directory, "whole_interior_exclusion.toml")
    open(proof_path, "w") do io
        TOML.print(io, Dict(
            "policy" => "verified-whole-relative-interior-dominance-and-dimension-skeleton-v1",
            "source_fan_sha256" => bytes2hex(sha256(read(source_fan_path))),
            "exclusion_constraints_sha256" => bytes2hex(sha256(read(constraints_path))),
            "skeleton_fan_sha256" => bytes2hex(sha256(read(skeleton_path))),
            "purity_certificate_sha256" => purity_certificate_signature(export_data.purity_certificate),
            "ideal_dimension" => export_data.ideal_dimension,
            "signed_actions_preserve_defining_ideal" => true,
            "group_preserves_candidate_fan" => true,
            "group_order" => length(prepared.coarse.symmetry.permutations),
            (export_data.symmetry_action == :taxa ? "taxon_permutations" : "variable_permutations") =>
                prepared.coarse.symmetry.permutations,
            "excluded_high_dimensional_orbits" => length(jobs),
            "excluded_high_dimensional_cells" => sum(job.orbit_size for job in jobs),
            "all_proper_faces_retained" => true,
            "tropical_support_contained_in_skeleton" => true,
            "skeleton_counts" => Dict(string(key) => getproperty(skeleton.counts, key)
                                      for key in propertynames(skeleton.counts)),
            "skeleton_provenance" => [Dict(string(key) => getproperty(record, key)
                                           for key in propertynames(record)) for record in skeleton.provenance],
            "checks" => proof_records))
    end
    atomic_write(joinpath(pruning_directory, "coarsened_prevariety.mrdi.source_sha256"),
                 bytes2hex(sha256(read(skeleton_path))) * "\n")
    open(joinpath(pruning_directory, "symmetry_summary.toml"), "w") do io
        TOML.print(io, Dict("source_fan_sha256" => bytes2hex(sha256(read(skeleton_path))),
            "ideal_dimension" => export_data.ideal_dimension,
            "coarsened_maximal_cones" => n_maximal_cones(skeleton.fan),
            "coarsened_cell_orbits" => length(coarse.jobs),
            "coarsened_group_order" => length(coarse.symmetry.permutations),
            "high_dimensional_orbits" => 0,
            "completion_policy" => "exact balanced-support membership after verified whole-interior exclusions"))
    end
    println("  Retained dimension skeleton: $(n_maximal_cones(skeleton.fan)) cells in $(length(coarse.jobs)) tree orbits.")
    flush(stdout)
    selected = collect_balanced_support((fine_fan=skeleton.fan, coarse=coarse, high_jobs=NamedTuple[]),
        export_data, config, pruning_directory, skeleton_path)
    support_certificate_path = joinpath(pruning_directory, "support_certificate.toml")
    certificate = TOML.parsefile(support_certificate_path)
    certificate["candidate_containment_proof_sha256"] = bytes2hex(sha256(read(proof_path)))
    certificate["original_prevariety_sha256"] = bytes2hex(sha256(read(source_fan_path)))
    open(support_certificate_path, "w") do io
        TOML.print(io, certificate)
    end
    return merge(selected, (phase=:whole_interior_pruning, containment_proof_path=proof_path))
end

function collect_symmetric_witnesses(prepared, export_data, basis_set, config, directory;
                                     strategy::Symbol=:groebner, source_fan_path=nothing,
                                     whole_interior_pruning::Bool=false)
    if whole_interior_pruning && strategy == :balanced_support && !isempty(prepared.high_jobs) &&
       all(job.dimension == export_data.ideal_dimension + 1 for job in prepared.high_jobs)
        pruned = collect_whole_interior_pruning(prepared, export_data, basis_set, config,
                                                directory, source_fan_path)
        isnothing(pruned) || return pruned
    end
    checks_path = joinpath(directory, "orbit_checks.csv")
    support = collect_witnesses(prepared.high_jobs, basis_set, config;
                                membership_only=true, checks_path=checks_path)
    if !isempty(support.witnesses)
        return merge(support, (complete=false, phase=:support_membership))
    end
    if strategy == :balanced_support && isempty(prepared.high_jobs)
        return collect_balanced_support(prepared, export_data, config, directory, source_fan_path)
    end
    # A tropical interior does not certify a complete coarsened cone. Fine
    # cones retain the tropical-basis walls needed by the original proof.
    if hasproperty(prepared, :fine) && !isnothing(prepared.fine)
        # Already saved as coarsened_symmetry: a copy takes 3 GB for 17 million cells.
        fine = prepared.fine
    else
        fine = tree_orbit_jobs(prepared.fine_fan, export_data.symmetry_permutations; action=export_data.symmetry_action)
        save_symmetry_artifacts(joinpath(directory, "fine_symmetry"), fine.symmetry)
    end
    println("  Whole-cone certification: $(length(fine.jobs)) fine tree-cell orbit representatives")
    flush(stdout)
    result = collect_witnesses(fine.jobs, basis_set, config; checks_path=checks_path)
    complete = isempty(result.witnesses) && result.failures == 0 && result.checked == length(fine.jobs)
    return merge(result, (complete=complete, phase=:whole_cone))
end

"Select entire expected-dimensional cone orbits by the balanced-support theorem."
function collect_balanced_support(prepared, export_data, config, directory, source_fan_path)
    validate_support_purity_certificate(export_data.purity_certificate, export_data.ideal_dimension)
    all(job.dimension <= export_data.ideal_dimension for job in prepared.coarse.jobs) ||
        error("balanced-support selection cannot certify higher-dimensional candidate cells")
    jobs = filter(job -> job.dimension == export_data.ideal_dimension, prepared.coarse.jobs)
    isempty(jobs) && error("candidate support has no cells of the certified tropical dimension")
    selected_cells = Int[]
    selected_orbits = 0
    checked = 0
    chunk_size = max(config.workers, 2 * config.workers)
    checks_path = joinpath(directory, "support_membership.csv")
    orbits_by_representative = Dict(first(orbit) => orbit for orbit in prepared.coarse.symmetry.cell_orbits)
    atomic_write(checks_path, "representative,orbit_size,dimension,tropical,reason,groebner_generators,seconds,weight,cells\n")
    progress = TaskProgress("balanced-support membership", length(jobs))
    for first_index in 1:chunk_size:length(jobs)
        chunk = jobs[first_index:min(length(jobs), first_index + chunk_size - 1)]
        results = pmap(job -> Main.parallel_trop_support_membership(job), chunk; batch_size=1)
        for (job, result) in zip(chunk, results)
            checked += 1
            orbit = orbits_by_representative[job.representative]
            if result.valid
                append!(selected_cells, orbit)
                selected_orbits += 1
            end
            open(checks_path, "a") do io
                println(io, join([job.representative, job.orbit_size, job.dimension,
                                  result.valid, result.reason, result.basis_size, result.seconds,
                                  join(job.weight, ' '), join(orbit, ';')], ','))
            end
            advance!(progress)
        end
    end
    isempty(selected_cells) && error("exact membership rejected every cell of a nonempty pure tropical variety")
    sort!(selected_cells)
    fan_path = save_certified_support_fan(prepared.coarse.normalized, selected_cells,
                                          source_fan_path, directory)
    lower_cells = sum(job.orbit_size for job in prepared.coarse.jobs
                      if job.dimension < export_data.ideal_dimension; init=0)
    certificate = export_data.purity_certificate
    open(joinpath(directory, "support_certificate.toml"), "w") do io
        TOML.print(io, Dict(
            "certification_strategy" => "balanced_support",
            "proof" => "positive balanced pure d-cycle contained in a candidate fan of dimension at most d",
            "ideal_dimension" => export_data.ideal_dimension,
            "purity_certificate_sha256" => purity_certificate_signature(certificate),
            "purity_method" => string(certificate.method),
            "source_fan_sha256" => bytes2hex(sha256(read(source_fan_path))),
            "certified_fan_sha256" => bytes2hex(sha256(read(fan_path))),
            "checked_orbits" => checked,
            "selected_orbits" => selected_orbits,
            "selected_cells" => length(selected_cells),
            "excluded_expected_dimension_cells" => sum(job.orbit_size for job in jobs; init=0) - length(selected_cells),
            "excluded_lower_dimension_maximal_cells" => lower_cells,
            "initial_ideal_constancy_certified" => false,
        ))
    end
    println("  Balanced support certified: $(length(selected_cells)) cells in $selected_orbits orbits; $lower_cells lower-dimensional maximal cells excluded by purity.")
    flush(stdout)
    return (witnesses=String[], checked=checked, valid=selected_orbits, failures=0,
            complete=true, phase=:balanced_support, certified_fan_path=fan_path)
end

function save_certified_support_fan(normalized, selected_cells, source_fan_path, directory)
    cells = fan_cells(normalized.fan)[selected_cells]
    used_rays = sort!(unique(reduce(vcat, cells)))
    ray_index = Dict(original => index for (index, original) in enumerate(used_rays))
    incidence = IncidenceMatrix([[ray_index[index] for index in cell] for cell in cells])
    rays = Matrix{QQFieldElem}(normalized.rays)[used_rays, :]
    lineality = matrix(QQ, rays_modulo_lineality(normalized.fan).lineality_basis)
    min_fan = polyhedral_fan(incidence, matrix(QQ, rays), lineality; non_redundant=true)
    max_fan = polyhedral_fan(incidence, matrix(QQ, -rays), -lineality; non_redundant=true)
    n_maximal_cones(min_fan) == length(selected_cells) || error("selected cone incidences changed during serialization")
    save(joinpath(directory, "certified_support_min.mrdi"), min_fan)
    save(joinpath(directory, "certified_support_max.mrdi"), max_fan)
    path = joinpath(directory, "certified_support.gfan")
    atomic_write(path, _gfan_fan_text(max_fan, parse_gfan_fan(source_fan_path)))
    return path
end

function gfan_basis_polynomials(path, variable_names)
    lines = readlines(path)
    isempty(lines) && return nothing
    first(lines) == "Q[" * join(variable_names, ',') * "]" || return nothing
    return [String(strip(replace(strip(line), r",$" => ""))) for line in lines[2:end]
            if !isempty(strip(line)) && !(strip(line) in ("{", "}"))]
end

function reusable_symmetric_prevariety(input_path, fan_path, basis, export_data, worker)
    isfile(input_path) && isfile(fan_path) || return false
    prevariety_checkpoint_is_parseable(fan_path) || return false
    old_basis = gfan_basis_polynomials(input_path, export_data.variable_names)
    isnothing(old_basis) && return false
    signatures = remotecall_fetch(
        (old, current) -> (Main.parallel_trop_basis_support_signature(old),
                           Main.parallel_trop_basis_support_signature(current)),
        worker, old_basis, basis,
    )
    return signatures[1] == signatures[2]
end

"Reject empty, truncated-metadata, or structurally invalid cached fan files."
function prevariety_checkpoint_is_parseable(fan_path)
    isfile(fan_path) && filesize(fan_path) > 0 || return false
    try
        fan = parse_gfan_fan(fan_path)
        0 <= fan.lineality_dimension <= fan.dimension <= fan.ambient_dimension || return false
        all(length(row) == fan.ambient_dimension for row in
            Iterators.flatten((fan.rays, fan.lineality, fan.orth_lineality))) || return false
        count_match = match(r"(?:^|\n)N_RAYS\n\s*(\d+)", fan.text)
        isnothing(count_match) && return false
        parse(Int, only(count_match.captures)) == length(fan.rays) || return false
        all(0 <= index < length(fan.rays) for cone in fan.maximal_cones for index in cone) || return false
        section = occursin(r"(?:^|\n)MAXIMAL_CONES\n", fan.text) ?
            "MAXIMAL_CONES" : "MAXIMAL_CONES_OF_CLOSURE"
        for line in section_lines(fan.text, section)
            data = strip(first(split(line, '#', limit=2)))
            startswith(data, "{") && !endswith(data, "}") && return false
        end
        return true
    catch
        return false
    end
end

function verification_marker_is_current(marker_path, basis_input_path, fan_path, bits;
                                         strategy::Symbol=:groebner, purity_signature::String="none",
                                         canonical::Bool=false)
    isfile(marker_path) || return false
    isfile(basis_input_path) || return false
    isfile(fan_path) || return false
    expected = join([
        "CERTIFICATION_POLICY " * TROPICAL_CERTIFICATION_POLICY,
        "CERTIFICATION_STRATEGY " * string(strategy),
        "PURITY_CERTIFICATE " * purity_signature,
        "BASIS_SHA256 " * bytes2hex(sha256(read(basis_input_path))),
        "FAN_SHA256 " * bytes2hex(sha256(read(fan_path))),
        "GFAN_BITS $bits",
        # A fan certified before its prevariety was in canonical order is checked again.
        (canonical ? ["GFAN_ORDER canonical"] : String[])...,
        "",
    ], '\n')
    return read(marker_path, String) == expected
end

function write_verification_marker(marker_path, basis_input_path, fan_path, bits;
                                    strategy::Symbol=:groebner, purity_signature::String="none",
                                    canonical::Bool=false)
    atomic_write(marker_path, join([
        "CERTIFICATION_POLICY " * TROPICAL_CERTIFICATION_POLICY,
        "CERTIFICATION_STRATEGY " * string(strategy),
        "PURITY_CERTIFICATE " * purity_signature,
        "BASIS_SHA256 " * bytes2hex(sha256(read(basis_input_path))),
        "FAN_SHA256 " * bytes2hex(sha256(read(fan_path))),
        "GFAN_BITS $bits",
        (canonical ? ["GFAN_ORDER canonical"] : String[])...,
        "",
    ], '\n'))
end

function verify_final_basis(config, export_data, basis_polynomials; strategy::Symbol=:groebner,
                            whole_interior_pruning::Bool=false, canonical::Bool=false)
    verification_directory = joinpath(config.output, "verification_bits_$(config.bits)")
    mkpath(verification_directory)
    input_path = joinpath(verification_directory, "basis.gfan")
    fan_path = joinpath(verification_directory, "prevariety.gfan")
    log_path = joinpath(verification_directory, "gfan.log")
    witness_path = joinpath(verification_directory, "new_witnesses.txt")
    marker_path = joinpath(verification_directory, "certified_fan.txt")
    result_path = strategy == :balanced_support ?
        joinpath(verification_directory, "certified_support.gfan") : fan_path
    pruned_result_path = joinpath(verification_directory, "interior_pruning", "certified_support.gfan")
    whole_interior_pruning && isfile(pruned_result_path) && (result_path = pruned_result_path)
    purity_signature = purity_certificate_signature(export_data.purity_certificate)
    write_basis_input(input_path, export_data.variable_names, basis_polynomials)

    if verification_marker_is_current(
        marker_path,
        input_path,
        result_path,
        config.bits,
        strategy=strategy,
        purity_signature=purity_signature,
        canonical=canonical,
    )
        fan = parse_gfan_fan(result_path)
        println(
            "Using exact final-verification checkpoint $result_path " *
            "($(length(fan.maximal_cones)) maximal cones).",
        )
        return result_path
    end

    println(
        "Verifying the final $(length(basis_polynomials))-polynomial basis " *
        "with gfan --bits $(config.bits)...",
    )
    fan = run_prevariety(config, input_path, fan_path, log_path; canonical=canonical)
    if strategy == :balanced_support
        prepared = prepare_symmetric_round(fan_path, verification_directory, export_data; strategy=strategy)
        certification = collect_symmetric_witnesses(prepared, export_data, Set(basis_polynomials),
            config, verification_directory; strategy=strategy, source_fan_path=fan_path,
            whole_interior_pruning=whole_interior_pruning)
        certification.complete || error("exact support verification requires further tropical-basis refinement")
        result_path = hasproperty(certification, :certified_fan_path) ? certification.certified_fan_path : fan_path
        write_verification_marker(marker_path, input_path, result_path, config.bits;
                                  strategy=strategy, purity_signature=purity_signature, canonical=canonical)
        println("Arbitrary-precision tropical support certified exactly: $result_path")
        return result_path
    end
    jobs = if isnothing(export_data.symmetry_permutations)
        result = cone_jobs(fan)
        Random.seed!(0x510e527f)
        shuffle!(result)
        result
    else
        min_fan = save_lifted_fan(fan_path, joinpath(verification_directory, "fine_tree_fan.mrdi"),
                                  export_data.lift_indices).fan
        result = tree_orbit_jobs(min_fan, export_data.symmetry_permutations; action=export_data.symmetry_action)
        save_symmetry_artifacts(joinpath(verification_directory, "fine_symmetry"), result.symmetry)
        result.jobs
    end
    println(
        "  exact candidate fan: dim=$(fan.dimension), " *
        "lineality=$(fan.lineality_dimension), $(length(jobs)) maximal cones",
    )

    certification = collect_witnesses(jobs, Set(basis_polynomials), config)
    atomic_write(
        witness_path,
        join(certification.witnesses, "\n") *
        (isempty(certification.witnesses) ? "" : "\n"),
    )
    if certification.failures != 0
        error(
            "arbitrary-precision final verification found " *
            "$(certification.failures) non-tropical cone(s); artifacts retained in " *
            verification_directory,
        )
    end
    certification.checked == length(jobs) || error(
        "arbitrary-precision final verification stopped before checking every cone",
    )
    write_verification_marker(marker_path, input_path, fan_path, config.bits;
                              strategy=strategy, purity_signature=purity_signature, canonical=canonical)
    println(
        "All $(length(jobs)) arbitrary-precision cone representatives certified exactly. " *
        "Verification marker: $marker_path",
    )
    return fan_path
end

function collect_initial_ideal_groups(jobs, config)
    group_ids = Int[]
    signature_to_group = Dict{String, Int}()
    checked = 0
    chunk_size = max(config.workers, 2 * config.workers)
    progress = TaskProgress("initial ideals", length(jobs))

    for first_index in 1:chunk_size:length(jobs)
        last_index = min(length(jobs), first_index + chunk_size - 1)
        chunk = jobs[first_index:last_index]
        results = pmap(
            job -> Main.parallel_trop_initial_ideal_signature(job.weight),
            chunk;
            batch_size=1,
        )
        for result in results
            group = get!(
                signature_to_group,
                result.signature,
                length(signature_to_group) + 1,
            )
            push!(group_ids, group)
            checked += 1
        end
        advance!(progress, length(results))
    end
    println("  $(length(signature_to_group)) Gröbner cones")
    return group_ids
end

function write_coarsening_checkpoint(path, fan_path, group_ids)
    digest = bytes2hex(sha256(read(fan_path)))
    number_of_groups = isempty(group_ids) ? 0 : maximum(group_ids)
    contents = join([
        "FAN_SHA256 $digest",
        "N_CONES $(length(group_ids))",
        "N_GROUPS $number_of_groups",
        "GROUP_IDS",
        join(group_ids, '\n'),
        "",
    ], '\n')
    atomic_write(path, contents)
end

function load_coarsening_checkpoint(path, fan_path, number_of_cones)
    isfile(path) || return nothing
    lines = readlines(path)
    length(lines) >= 4 || return nothing
    lines[1] == "FAN_SHA256 " * bytes2hex(sha256(read(fan_path))) || return nothing
    lines[2] == "N_CONES $number_of_cones" || return nothing
    startswith(lines[3], "N_GROUPS ") || return nothing
    lines[4] == "GROUP_IDS" || return nothing
    group_ids = parse.(Int, lines[5:end])
    length(group_ids) == number_of_cones || return nothing
    number_of_groups = parse(Int, split(lines[3])[2])
    isempty(group_ids) || maximum(group_ids) == number_of_groups || return nothing
    return group_ids
end

function coarsening_outputs_are_current(marker_path, fan_path, group_ids, output_paths)
    isfile(marker_path) || return false
    all(isfile, output_paths) || return false
    expected = join([
        "FAN_SHA256 " * bytes2hex(sha256(read(fan_path))),
        "N_GROUPS $(maximum(group_ids))",
        "",
    ], '\n')
    return read(marker_path, String) == expected
end

function write_coarsening_completion(marker_path, fan_path, group_ids)
    atomic_write(marker_path, join([
        "FAN_SHA256 " * bytes2hex(sha256(read(fan_path))),
        "N_GROUPS $(maximum(group_ids))",
        "",
    ], '\n'))
end

"Preserve trusted native scalar facts under an injective coordinate lift."
function _seed_gfan_scalar_metadata(lifted_fan, source_fan, lift_indices)
    # Dropping a source coordinate can change dimensions, purity and ray
    # independence. Duplication and sign reversal preserve all three facts.
    Set(lift_indices) == Set(1:source_fan.ambient_dimension) || return lifted_fan
    metadata = Dict{Symbol, Any}()
    if occursin(r"(?:^|\n)DIM\n", source_fan.text)
        metadata[:FAN_DIM] = source_fan.dimension
    end
    for (header, property) in (("PURE", :PURE), ("SIMPLICIAL", :SIMPLICIAL))
        matched = match(Regex("(?:^|\\n)" * header * "\\n([^\\n]*)"), source_fan.text)
        isnothing(matched) && continue
        value = strip(first(split(only(matched.captures), '#'; limit=2)))
        value in ("0", "1") || throw(ArgumentError("invalid Gfan $header scalar header"))
        metadata[property] = value == "1"
    end
    isempty(metadata) && return lifted_fan
    # These are immutable Polymake facts. Supply them with the original
    # incidence data in a constructor transaction, before it is committed.
    object = Oscar.pm_object(lifted_fan)
    seeded = Oscar.Polymake.bigobject(Oscar.Polymake.bigobject_qualifiedname(object);
        RAYS=object.RAYS, LINEALITY_SPACE=object.LINEALITY_SPACE,
        MAXIMAL_CONES=object.MAXIMAL_CONES, metadata...)
    return Oscar.polyhedral_fan(seeded)
end

function _save_lifted_fan_oscar(fan_path, output_path, lift_indices)
    fan = parse_gfan_fan(fan_path)
    original_ambient_dimension = length(lift_indices)
    all(index -> 1 <= index <= fan.ambient_dimension, lift_indices) ||
        throw(ArgumentError("lift index lies outside the reduced ambient space"))
    # gfan uses max; the surrounding project uses Oscar's min convention.
    lifted_rays = [[-row[i] for i in lift_indices] for row in fan.rays]
    lifted_lineality = [[row[i] for i in lift_indices] for row in fan.lineality]
    cone_incidence = [cone .+ 1 for cone in fan.maximal_cones]

    ray_matrix = isempty(lifted_rays) ?
        Oscar.zero_matrix(Oscar.QQ, 0, original_ambient_dimension) :
        Oscar.matrix(Oscar.QQ, lifted_rays)
    lineality_matrix = isempty(lifted_lineality) ?
        Oscar.zero_matrix(Oscar.QQ, 0, original_ambient_dimension) :
        Oscar.matrix(Oscar.QQ, lifted_lineality)
    lifted_fan = Oscar.polyhedral_fan(
        Oscar.IncidenceMatrix(cone_incidence),
        ray_matrix,
        lineality_matrix;
        non_redundant=true,
    )
    lifted_fan = _seed_gfan_scalar_metadata(lifted_fan, fan, lift_indices)
    mkpath(dirname(output_path))
    Oscar.save(output_path, lifted_fan)
    return (fan=lifted_fan, maximal_cones=Oscar.n_maximal_cones(lifted_fan))
end

function save_lifted_fan(fan_path, output_path, lift_indices)
    # Loading Oscar only for final serialization saves one large Oscar/Singular
    # instance throughout the parallel certification.  invokelatest crosses
    # the world-age boundary introduced by this deliberately late import.
    @eval import Oscar
    return Base.invokelatest(_save_lifted_fan_oscar, fan_path, output_path, lift_indices)
end

function _primitive_integer_row(row)
    denominators = BigInt[BigInt(denominator(entry)) for entry in row]
    common_denominator = foldl(lcm, denominators; init=big(1))
    integer_row = BigInt[
        BigInt(numerator(entry)) * div(common_denominator, BigInt(denominator(entry)))
        for entry in row
    ]
    common_divisor = foldl(gcd, abs.(integer_row); init=big(0))
    iszero(common_divisor) && return integer_row
    return div.(integer_row, common_divisor)
end

function _gfan_fan_text(polyhedral_fan, source_fan)
    rml = Oscar.rays_modulo_lineality(polyhedral_fan)
    ray_rows = _primitive_integer_row.(collect(rml.rays_modulo_lineality))
    lineality_rows = _primitive_integer_row.(collect(rml.lineality_basis))
    incidence = Oscar.maximal_cones(Oscar.IncidenceMatrix, polyhedral_fan)
    maximal_cones = [
        sort!([index - 1 for index in Oscar.row(incidence, row_index)])
        for row_index in axes(incidence, 1)
    ]
    lineality_dimension = Oscar.lineality_dim(polyhedral_fan)

    row_text(rows) = join((join(row, " ") for row in rows), '\n')
    cone_text = join(
        ("{" * join(cone, " ") * "}" for cone in maximal_cones),
        '\n',
    )
    return join([
        "_application PolyhedralFan",
        "_version 2.2",
        "_type PolyhedralFan",
        "",
        "AMBIENT_DIM",
        string(Oscar.ambient_dim(polyhedral_fan)),
        "",
        "DIM",
        string(Oscar.dim(polyhedral_fan)),
        "",
        "LINEALITY_DIM",
        string(lineality_dimension),
        "",
        "RAYS",
        row_text(ray_rows),
        "",
        "N_RAYS",
        string(length(ray_rows)),
        "",
        "LINEALITY_SPACE",
        row_text(lineality_rows),
        "",
        "ORTH_LINEALITY_SPACE",
        row_text(source_fan.orth_lineality),
        "",
        "MAXIMAL_CONES",
        cone_text,
        "",
        "SIMPLICIAL",
        Oscar.is_simplicial(polyhedral_fan) ? "1" : "0",
        "",
        "PURE",
        Oscar.is_pure(polyhedral_fan) ? "1" : "0",
        "",
    ], '\n')
end

function _coarsen_and_save_oscar(
    fine_fan_path,
    group_ids,
    reduced_gfan_path,
    reduced_oscar_path,
    lifted_oscar_path,
    lift_indices,
)
    fine_fan = parse_gfan_fan(fine_fan_path)
    length(group_ids) == length(fine_fan.maximal_cones) ||
        throw(ArgumentError("one coarsening group is required per fine maximal cone"))
    number_of_groups = maximum(group_ids)
    sort!(unique(group_ids)) == collect(1:number_of_groups) ||
        throw(ArgumentError("coarsening group identifiers must be consecutive"))

    grouped_ray_indices = [Set{Int}() for _ in 1:number_of_groups]
    for (cone, group) in zip(fine_fan.maximal_cones, group_ids)
        union!(grouped_ray_indices[group], cone)
    end
    lineality_matrix = Oscar.matrix(Oscar.QQ, fine_fan.lineality)
    coarse_cones = [
        Oscar.positive_hull(
            Oscar.matrix(
                Oscar.QQ,
                [fine_fan.rays[index + 1] for index in sort!(collect(indices))],
            ),
            lineality_matrix;
            non_redundant=false,
        )
        for indices in grouped_ray_indices
    ]
    coarse_fan = Oscar.polyhedral_fan(coarse_cones; non_redundant=true)
    Oscar.n_maximal_cones(coarse_fan) == number_of_groups || error(
        "coarsening produced $(Oscar.n_maximal_cones(coarse_fan)) maximal cones " *
        "from $number_of_groups initial ideals",
    )

    Oscar.save(reduced_oscar_path, coarse_fan)

    rml = Oscar.rays_modulo_lineality(coarse_fan)
    reduced_rays = collect(rml.rays_modulo_lineality)
    reduced_lineality = collect(rml.lineality_basis)
    lifted_rays = [[-ray[index] for index in lift_indices] for ray in reduced_rays]
    lifted_lineality = [
        [direction[index] for index in lift_indices]
        for direction in reduced_lineality
    ]
    incidence = Oscar.maximal_cones(Oscar.IncidenceMatrix, coarse_fan)
    lifted_fan = Oscar.polyhedral_fan(
        incidence,
        Oscar.matrix(Oscar.QQ, lifted_rays),
        Oscar.matrix(Oscar.QQ, lifted_lineality);
        non_redundant=true,
    )
    Oscar.save(lifted_oscar_path, lifted_fan)
    # The Oscar artifacts are the primary results. Write the auxiliary gfan
    # representation last; in particular, do not compute a full face vector.
    atomic_write(reduced_gfan_path, _gfan_fan_text(coarse_fan, fine_fan))
    return (
        fan=lifted_fan,
        maximal_cones=Oscar.n_maximal_cones(lifted_fan),
        rays=length(reduced_rays),
    )
end

function coarsen_and_save(
    fine_fan_path,
    group_ids,
    reduced_gfan_path,
    reduced_oscar_path,
    lifted_oscar_path,
    lift_indices,
)
    @eval import Oscar
    return Base.invokelatest(
        _coarsen_and_save_oscar,
        fine_fan_path,
        group_ids,
        reduced_gfan_path,
        reduced_oscar_path,
        lifted_oscar_path,
        lift_indices,
    )
end

"Publish a certified support fan without claiming any initial-ideal constancy."
function publish_certified_support(problem, config, fan_path, lift_indices)
    max_path = joinpath(config.output, problem.output_prefix * "_reduced_max.gfan")
    min_path = joinpath(config.output, problem.output_prefix * "_Fan.mrdi")
    cp(fan_path, max_path; force=true)
    result = save_lifted_fan(fan_path, min_path, lift_indices)
    source_max = joinpath(dirname(fan_path), "certified_support_max.mrdi")
    isfile(source_max) && cp(source_max,
        joinpath(config.output, problem.output_prefix * "_reduced_max.mrdi"); force=true)
    cp(fan_path, joinpath(config.output, problem.output_prefix * "_refined_reduced_max.gfan"); force=true)
    cp(min_path, joinpath(config.output, problem.output_prefix * "_refined_Fan.mrdi"); force=true)
    if !isnothing(problem.published_fan)
        mkpath(dirname(problem.published_fan))
        cp(min_path, problem.published_fan; force=true)
    end
    println("Published certified tropical support ($(result.maximal_cones) maximal cones): $min_path")
    println("Initial ideals may vary inside these support cells; Gröbner grouping was skipped.")
    if !isnothing(problem.published_coarsest_fan)
        input_fan = something(problem.published_fan, min_path)
        coarsen_fan_file(input_fan, problem.published_coarsest_fan;
            label="Trop($(problem.name))",
            provenance_path=something(problem.coarsest_provenance,
                joinpath(config.output, problem.output_prefix * "_coarsest_cells.csv")),
            report_path=something(problem.coarsest_report,
                joinpath(config.output, problem.output_prefix * "_COARSEST.md")))
    end
    return result
end

function run_parallel_tropicalization(
    problem::TropicalProblem,
    config::TropicalizationConfig,
)
    isempty(config.output) && (config.output = problem.default_output)
    mkpath(config.output)
    println("Parallel tropicalization of $(problem.name)")
    println(
        "  workers=$(config.workers), gfan threads=$(config.threads), " *
        "gfan bits=$(config.bits)",
    )
    println("  output=$(config.output)")

    worker_ids = initialize_workers(config.workers, problem)
    export_data = remotecall_fetch(
        () -> Main.parallel_trop_export_data(),
        first(worker_ids),
    )
    strategy = problem.certification_strategy
    strategy in (:groebner, :balanced_support) || error("unknown tropical certification strategy $strategy")
    strategy == :balanced_support &&
        validate_support_purity_certificate(export_data.purity_certificate, export_data.ideal_dimension)
    purity_signature = purity_certificate_signature(export_data.purity_certificate)
    canonical = problem.canonical_gfan_order
    canonical && problem.orbit_restricted_intersection &&
        error("the orbit-restricted intersection has no canonical gfan order")
    println(
        "  torus coordinates: $(length(export_data.original_variable_names)) -> " *
        "$(length(export_data.variable_names)) variables; ",
        "$(length(export_data.generators)) distinct generators",
    )

    basis_path = joinpath(config.output, "basis_polynomials.txt")
    if config.verify_final && !isfile(basis_path)
        error("final verification needs the existing basis checkpoint $basis_path")
    end
    if !config.resume && isfile(basis_path)
        error(
            "output already contains a checkpoint at $basis_path; " *
            "use --resume or choose a different --output directory",
        )
    end
    basis_polynomials = if config.resume && isfile(basis_path)
        filter(!isempty, String.(strip.(readlines(basis_path))))
    else
        copy(export_data.generators)
    end
    symmetric = !isnothing(export_data.symmetry_permutations)
    if symmetric
        before = length(basis_polynomials)
        basis_polynomials = remotecall_fetch(
            polynomials -> Main.parallel_trop_symmetrize_polynomials(polynomials),
            first(worker_ids), basis_polynomials,
        )
        println("  verified signed ideal symmetry: $(length(export_data.symmetry_permutations)) " *
                "$(export_data.symmetry_action == :taxa ? "taxon" : "variable") actions; ideal dimension=$(export_data.ideal_dimension)" *
                (export_data.declared_dimension ? " (declared by the provider)" : ""))
        println("  symmetry-closed polynomial supports: $before -> $(length(basis_polynomials))")
        flush(stdout)
    end

    if config.verify_final
        try
            verify_final_basis(config, export_data, basis_polynomials; strategy=strategy,
                               whole_interior_pruning=problem.whole_interior_pruning, canonical=canonical)
        finally
            rmprocs(worker_ids)
        end
        return
    end

    atomic_write(basis_path, join(basis_polynomials, "\n") * "\n")

    ideal_input_path = joinpath(config.output, problem.slug * "_reduced.gfan")
    write_basis_input(ideal_input_path, export_data.variable_names, basis_polynomials)
    if config.prepare_only
        println("Prepared $ideal_input_path")
        rmprocs(worker_ids)
        return
    end

    certificate_path = joinpath(config.output, "certified_fan.txt")
    certificate_policy_path = joinpath(config.output, "certified_policy.txt")
    completed = false
    final_fan_path = ""
    if config.resume && isfile(certificate_path)
        relative_fan_path = strip(read(certificate_path, String))
        final_fan_path = normpath(joinpath(config.output, relative_fan_path))
        isfile(final_fan_path) || error(
            "certification marker points to a missing fan: $final_fan_path",
        )
        completed = verification_marker_is_current(
            certificate_policy_path, ideal_input_path, final_fan_path, config.bits,
            strategy=strategy, purity_signature=purity_signature, canonical=canonical,
        )
        println(completed ? "Using certified fan checkpoint $final_fan_path" :
                "Existing fan checkpoint needs certification under the current exact membership policy.")
    end

    round_number = next_round_number(config.output, config.resume)
    if config.resume && !completed && !isempty(final_fan_path)
        checkpoint_round = checkpoint_round_number(config.output, final_fan_path)
        checkpoint_directory = dirname(final_fan_path)
        basename(checkpoint_directory) in ("interior_pruning", "orbit_restricted") &&
            (checkpoint_directory = dirname(checkpoint_directory))
        old_input = joinpath(checkpoint_directory, "basis.gfan")
        if !isnothing(checkpoint_round) && reusable_symmetric_prevariety(
            old_input, final_fan_path, basis_polynomials, export_data, first(worker_ids),
        )
            # An old predicate may have falsely certified this fan. Recheck
            # its existing round rather than discarding the costly intersection.
            round_number = checkpoint_round
        end
    end
    rounds_run = 0
    while !completed && rounds_run < config.max_rounds
        rounds_run += 1
        round_directory = joinpath(config.output, "round_" * lpad(round_number, 3, '0'))
        mkpath(round_directory)
        input_path = joinpath(round_directory, "basis.gfan")
        fan_path = joinpath(round_directory, "prevariety.gfan")
        orbit_fan_path = joinpath(round_directory, "orbit_restricted", "expanded_prevariety.gfan")
        problem.orbit_restricted_intersection && isfile(orbit_fan_path) && (fan_path = orbit_fan_path)
        log_path = joinpath(round_directory, "gfan.log")
        reuse_fan = (symmetric || config.resume) && reusable_symmetric_prevariety(
            input_path, fan_path, basis_polynomials, export_data, first(worker_ids),
        )
        write_basis_input(input_path, export_data.variable_names, basis_polynomials)

        println("Round $round_number: $(length(basis_polynomials)) basis polynomials")
        flush(stdout)
        fan = if reuse_fan
            println("  Reusing the existing gfan prevariety with identical polynomial supports.")
            # It may have been saved before its order was canonical.
            timed(() -> canonical ? canonicalize_gfan_fan(fan_path) : parse_gfan_fan(fan_path),
                  "Reading the candidate fan")
        elseif problem.orbit_restricted_intersection && round_number > 1
            symmetric || error("orbit-restricted intersection requires verified signed symmetry")
            intersection = run_orbit_restricted_prevariety(problem, config, export_data,
                basis_polynomials; round_number=round_number, worker=first(worker_ids),
                release_workers=let released = worker_ids; () -> rmprocs(released) end)
            fan_path = intersection.fan_path
            timed(() -> parse_gfan_fan(fan_path), "Reading the candidate fan")
        elseif symmetric && problem.gfan_symmetry
            run_symmetric_prevariety(config, input_path, fan_path, log_path, export_data; canonical=canonical)
        else
            run_prevariety(config, input_path, fan_path, log_path; canonical=canonical)
        end
        # Symmetric rounds check orbit representatives rather than every cone.
        jobs = symmetric ? nothing : cone_jobs(fan)
        Random.seed!(0x6a09e667 + round_number)
        symmetric || shuffle!(jobs)
        println(
            "  candidate fan: dim=$(fan.dimension), lineality=$(fan.lineality_dimension), ",
            "$(length(fan.maximal_cones)) maximal cones",
        )

        basis_set = Set(basis_polynomials)
        certification = if symmetric
            prepared = prepare_symmetric_round(fan_path, round_directory, export_data; strategy=strategy,
                                               coarsen=problem.coarsen_prevariety)
            worker_ids = restored_workers(worker_ids, config, problem)
            collect_symmetric_witnesses(prepared, export_data, basis_set, config, round_directory;
                                        strategy=strategy, source_fan_path=fan_path,
                                        whole_interior_pruning=problem.whole_interior_pruning)
        else
            collect_witnesses(jobs, basis_set, config)
        end
        witness_path = joinpath(round_directory, "new_witnesses.txt")
        atomic_write(witness_path, join(certification.witnesses, "\n") *
                                   (isempty(certification.witnesses) ? "" : "\n"))

        if isempty(certification.witnesses)
            certification.failures == 0 || error(
                "non-tropical cones were found but no new witness was produced; " *
                "checkpoint retained in $round_directory",
            )
            if symmetric
                certification.complete || error("fine symmetry-orbit certification stopped early")
            else
                certification.checked == length(jobs) || error("internal certification stopped early")
            end
            completed = true
            final_fan_path = hasproperty(certification, :certified_fan_path) ?
                certification.certified_fan_path : fan_path
            atomic_write(certificate_path, relpath(final_fan_path, config.output) * "\n")
            write_verification_marker(certificate_policy_path, ideal_input_path, final_fan_path, config.bits;
                                      strategy=strategy, purity_signature=purity_signature, canonical=canonical)
            println("Tropicalization complete under $strategy certification: $final_fan_path")
            break
        end

        append!(basis_polynomials, certification.witnesses)
        symmetric && sort!(basis_polynomials)
        atomic_write(basis_path, join(basis_polynomials, "\n") * "\n")
        write_basis_input(ideal_input_path, export_data.variable_names, basis_polynomials)
        println("  added $(length(certification.witnesses)) Gröbner witnesses")
        config.fresh_processes && leave_for_fresh_process("Round $round_number is checked.")
        round_number += 1
    end

    if !completed
        rmprocs(worker_ids)
        println("Stopped after $rounds_run rounds; resume with --resume.")
        return
    end

    fine_fan = parse_gfan_fan(final_fan_path)
    fine_jobs = cone_jobs(fine_fan)

    if !isnothing(problem.published_prevariety)
        first_prevariety = joinpath(config.output, "round_001", "prevariety.gfan")
        isfile(first_prevariety) || error(
            "the initial-generator prevariety checkpoint is missing: $first_prevariety",
        )
        save_lifted_fan(
            first_prevariety,
            problem.published_prevariety,
            export_data.lift_indices,
        )
        println("Published initial-generator prevariety to $(problem.published_prevariety).")
    end

    if strategy == :balanced_support
        rmprocs(worker_ids)
        publish_certified_support(problem, config, final_fan_path, export_data.lift_indices)
        return
    end

    coarsening_checkpoint = joinpath(config.output, "groebner_coarsening.txt")
    group_ids = load_coarsening_checkpoint(
        coarsening_checkpoint,
        final_fan_path,
        length(fine_jobs),
    )
    if isnothing(group_ids)
        println(
            "Coarsening $(length(fine_jobs)) tropical-basis cones by weight-initial ideal...",
        )
        group_ids = collect_initial_ideal_groups(fine_jobs, config)
        write_coarsening_checkpoint(coarsening_checkpoint, final_fan_path, group_ids)
    else
        println(
            "Using Gröbner coarsening checkpoint with $(maximum(group_ids)) cones.",
        )
    end

    rmprocs(worker_ids)

    refined_gfan_path = joinpath(
        config.output,
        problem.output_prefix * "_refined_reduced_max.gfan",
    )
    cp(final_fan_path, refined_gfan_path; force=true)
    refined_lifted_path = joinpath(
        config.output,
        problem.output_prefix * "_refined_Fan.mrdi",
    )
    refined = save_lifted_fan(
        final_fan_path,
        refined_lifted_path,
        export_data.lift_indices,
    )
    println(
        "Saved the fine $(length(export_data.original_variable_names))-coordinate " *
        "min-convention fan to $refined_lifted_path ",
        "($(refined.maximal_cones) maximal cones).",
    )

    coarse_gfan_path = joinpath(
        config.output,
        problem.output_prefix * "_reduced_max.gfan",
    )
    coarse_reduced_path = joinpath(
        config.output,
        problem.output_prefix * "_reduced_max.mrdi",
    )
    coarse_lifted_path = joinpath(
        config.output,
        problem.output_prefix * "_Fan.mrdi",
    )
    coarse_marker_path = joinpath(config.output, "groebner_coarsened.txt")
    coarse_output_paths = [coarse_gfan_path, coarse_reduced_path, coarse_lifted_path]
    if coarsening_outputs_are_current(
        coarse_marker_path,
        final_fan_path,
        group_ids,
        coarse_output_paths,
    )
        coarse_fan = parse_gfan_fan(coarse_gfan_path)
        println(
            "Using Gröbner-coarsened fan checkpoint at $coarse_lifted_path ",
            "($(length(coarse_fan.maximal_cones)) maximal cones, ",
            "$(length(coarse_fan.rays)) rays).",
        )
    else
        coarse = coarsen_and_save(
            final_fan_path,
            group_ids,
            coarse_gfan_path,
            coarse_reduced_path,
            coarse_lifted_path,
            export_data.lift_indices,
        )
        write_coarsening_completion(coarse_marker_path, final_fan_path, group_ids)
        println(
            "Saved the Gröbner-coarsened " *
            "$(length(export_data.original_variable_names))-coordinate " *
            "min-convention fan to ",
            "$coarse_lifted_path ($(coarse.maximal_cones) maximal cones, ",
            "$(coarse.rays) rays).",
        )
    end

    if !isnothing(problem.published_fan)
        mkpath(dirname(problem.published_fan))
        publication_is_current = isfile(problem.published_fan) &&
            sha256(read(coarse_lifted_path)) == sha256(read(problem.published_fan))
        if publication_is_current
            println("Published fan is already current at $(problem.published_fan).")
        else
            cp(coarse_lifted_path, problem.published_fan; force=true)
            println("Published fan to $(problem.published_fan).")
        end
    end

    if !isnothing(problem.published_coarsest_fan)
        input_fan = something(problem.published_fan, coarse_lifted_path)
        provenance = something(
            problem.coarsest_provenance,
            joinpath(config.output, problem.output_prefix * "_coarsest_cells.csv"),
        )
        report = something(
            problem.coarsest_report,
            joinpath(config.output, problem.output_prefix * "_COARSEST.md"),
        )
        artifacts = [problem.published_coarsest_fan, provenance, report]
        if all(isfile, artifacts) && all(mtime(path) >= mtime(input_fan) for path in artifacts)
            println("Using current coarsest support-fan checkpoint $(problem.published_coarsest_fan).")
        else
            result = Base.invokelatest(
                coarsen_fan_file,
                input_fan,
                problem.published_coarsest_fan;
                label="Trop($(problem.name))",
                provenance_path=provenance,
                report_path=report,
            )
            println(
                "Published coarsest support fan to $(problem.published_coarsest_fan) " *
                "($(result.output_cones) maximal cones, $(result.output_rays) rays).",
            )
        end
    end
end

include(joinpath(PARALLEL_SOURCE_DIR, "orbit_restricted_gfan.jl"))

function as_polyhedral_fan(complex)
    incidence = ray_indices(maximal_polyhedra(complex))
    data = rays_modulo_lineality(complex)
    return polyhedral_fan(
        incidence,
        data.rays_modulo_lineality,
        data.lineality_basis;
        non_redundant=true,
    )
end

function run_oscar_tropicalization(
    problem::TropicalProblem,
    config::TropicalizationConfig,
)
    output = problem.published_fan
    isnothing(output) && throw(ArgumentError("Oscar tropicalization requires published_fan"))
    prevariety_output = problem.published_prevariety
    outputs = isnothing(prevariety_output) ? [output] : [prevariety_output, output]
    !isnothing(problem.published_complex) && push!(outputs, problem.published_complex)
    all(isfile, outputs) && !config.force && return (
        prevariety=isnothing(prevariety_output) ? nothing : load(prevariety_output),
        fan=load(output),
        path=output,
        cached=true,
    )
    isfile(problem.provider_file) || error("missing ideal provider $(problem.provider_file)")
    Base.include(@__MODULE__, problem.provider_file)
    provider = getproperty(@__MODULE__, problem.provider_module)
    builder = getproperty(provider, problem.builder)
    data = Base.invokelatest(builder)
    nu = tropical_semiring_map(coefficient_ring(data.ideal), min)
    prevariety = if isnothing(prevariety_output)
        nothing
    else
        precomplex = tropical_prevariety(gens(data.ideal), nu)
        prefan = as_polyhedral_fan(precomplex)
        mkpath(dirname(prevariety_output))
        save(prevariety_output, prefan)
        prefan
    end
    tropical = tropical_variety(
        data.ideal,
        nu;
        skip_saturation=false,
        skip_decomposition=true,
    )
    complex = polyhedral_complex(tropical)
    fan = as_polyhedral_fan(complex)
    mkpath(dirname(output))
    !isnothing(problem.published_complex) && save(problem.published_complex, complex)
    save(output, fan)
    return (
        prevariety=prevariety,
        fan=fan,
        complex=complex,
        ideal=data.ideal,
        path=output,
        cached=false,
    )
end

"Run a study-independent tropicalization with the selected engine."
function run_tropicalization(
    problem::TropicalProblem,
    config::TropicalizationConfig=TropicalizationConfig(
        engine=problem.default_engine,
        output=problem.default_output,
    ),
)
    engine = config.engine == :auto ? problem.default_engine : config.engine
    config.workers >= 1 || throw(ArgumentError("workers must be positive"))
    config.threads >= 1 || throw(ArgumentError("threads must be positive"))
    config.bits in (0, 32, 64) || throw(ArgumentError("bits must be 0, 32, or 64"))
    isfile(problem.provider_file) ||
        throw(ArgumentError("ideal provider does not exist: $(problem.provider_file)"))
    engine == :oscar && return run_oscar_tropicalization(problem, config)
    if engine == :parallel
        isempty(config.gfan) && throw(ArgumentError("gfan_tropicalprevariety was not found"))
        isfile(config.gfan) ||
            throw(ArgumentError("gfan executable does not exist: $(config.gfan)"))
        return run_parallel_tropicalization(problem, config)
    end
    throw(ArgumentError("unsupported tropicalization engine: $engine"))
end
