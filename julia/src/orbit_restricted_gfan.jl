"Hash large native table files without holding their contents in memory."
orbit_stream_sha256(path) = bytes2hex(open(sha256, path))

"The native table grammar is triangular: polytope P has P-1 edge matrices."
function native_intersection_table_is_complete(path, edges)
    isfile(path) && !isempty(edges) || return false
    depth = 0
    p = q = e = f = 0
    opened = closed = false
    buffer = Vector{UInt8}(undef, 1 << 20)
    open(path, "r") do io
        while !eof(io)
            count = readbytes!(io, buffer, length(buffer))
            for index in 1:count
                byte = buffer[index]
                byte in (0x09, 0x0a, 0x0d, 0x20) && continue
                closed && return false
                if byte == UInt8('{')
                    if depth == 0
                        opened && return false
                        opened = true
                    elseif depth == 1
                        p += 1
                        p <= length(edges) || return false
                        q = 0
                    elseif depth == 2
                        q += 1
                        q < p || return false
                        e = 0
                    else
                        return false
                    end
                    depth += 1
                elseif byte == UInt8('(')
                    depth == 3 || return false
                    depth = 4
                    f = 0
                elseif byte == UInt8(')')
                    depth == 4 && f == edges[q] || return false
                    e += 1
                    depth = 3
                elseif byte == UInt8('}')
                    if depth == 3
                        e == edges[p] || return false
                    elseif depth == 2
                        q == p-1 || return false
                    elseif depth == 1
                        p == length(edges) || return false
                        closed = true
                    else
                        return false
                    end
                    depth -= 1
                elseif byte in (UInt8('0'), UInt8('1'))
                    depth == 4 || return false
                    f += 1
                elseif byte == UInt8(',')
                    depth > 0 || return false
                else
                    return false
                end
            end
        end
        return opened && closed && depth == 0
    end
end

function orbit_table_cache_key(config, basis_input_path, bits, number_of_polynomials)
    return Dict{String, Any}(
        "policy" => "gfan-ordered-basis-intersection-tables-v1",
        "basis_input_sha256" => orbit_stream_sha256(basis_input_path),
        "gfan_sha256" => orbit_stream_sha256(config.gfan),
        "integer_bits" => bits,
        "number_of_polynomials" => number_of_polynomials,
    )
end

function orbit_write_toml(path, document)
    io = IOBuffer()
    TOML.print(io, document; sorted=true)
    atomic_write(path, String(take!(io)))
end

function orbit_table_cache_paths(directory, bits)
    prefix = joinpath(directory, "intersection_tables_bits_$(bits)")
    return (table=prefix * ".dat", marker=prefix * ".toml",
            partial=prefix * ".dat.partial", pending=prefix * ".pending.toml")
end

function promote_orbit_table_cache(paths, key, log_path)
    isfile(paths.partial) && isfile(paths.pending) && isfile(log_path) || return false
    pending = TOML.parsefile(paths.pending)
    all(get(pending, field, nothing) == value for (field, value) in key) || return false
    edges = [parse(Int, only(matched.captures)) for matched in
             eachmatch(r"Number of edges in polytope:(\d+)", read(log_path, String))]
    length(edges) == key["number_of_polynomials"] || return false
    native_intersection_table_is_complete(paths.partial, edges) || return false
    mv(paths.partial, paths.table; force=true)
    marker = merge(key, Dict("status" => "complete", "edge_counts" => edges,
        "table_sha256" => orbit_stream_sha256(paths.table), "table_bytes" => filesize(paths.table)))
    orbit_write_toml(paths.marker, marker)
    return true
end

function reusable_orbit_table_cache(paths, key)
    isfile(paths.table) && isfile(paths.marker) || return false
    marker = TOML.parsefile(paths.marker)
    get(marker, "status", "") == "complete" || return false
    all(get(marker, field, nothing) == value for (field, value) in key) || return false
    get(marker, "table_bytes", -1) == filesize(paths.table) || return false
    return get(marker, "table_sha256", "") == orbit_stream_sha256(paths.table)
end

"Run one native restricted intersection, reusing completed matching pair tables."
function run_native_orbit_intersection(config, basis_path, restriction_input_path,
                                       directory, number_of_polynomials)
    fan_path = joinpath(directory, "restricted_prevariety.gfan")
    partial_path = fan_path * ".partial"
    isfile(fan_path) && mv(fan_path, fan_path * ".previous"; force=true)
    for bits in unique([config.bits, 0])
        log_path = joinpath(directory, "gfan_bits_$(bits).log")
        key = orbit_table_cache_key(config, basis_path, bits, number_of_polynomials)
        paths = orbit_table_cache_paths(directory, bits)
        promote_orbit_table_cache(paths, key, log_path)
        load_tables = reusable_orbit_table_cache(paths, key)
        table_option = if load_tables
            ["--loadintersectiondatafrom", paths.table]
        else
            # Never relabel an older partial table with a new basis key.
            # Gfan only opens this file after preprocessing, so an early
            # interruption must leave no apparently current older table.
            isfile(paths.partial) && mv(paths.partial, paths.partial * ".invalid"; force=true)
            orbit_write_toml(paths.pending, key)
            ["--saveintersectiondataas", paths.partial]
        end
        # Gfan keeps the refined cones of all chunks; larger chunks trade memory for parallelism.
        chunk_option = config.gfan_chunksize > 0 ? ["--chunksize", string(config.gfan_chunksize)] : String[]
        command = `$(config.gfan) -j $(config.threads) --bits $bits --halfopenrestrictions $chunk_option $table_option`
        println("  Orbit-restricted Gfan: $command")
        println(load_tables ? "  Reusing matching native intersection tables." :
                "  Computing native intersection tables; a completed table survives interrupted traversal.")
        flush(stdout)
        exit_code = open(restriction_input_path, "r") do input
            open(partial_path, "w") do output
                open(log_path, "w") do errors
                    process = run_reporting(pipeline(ignorestatus(command), stdin=input, stdout=output, stderr=errors), "gfan")
                    process.exitcode
                end
            end
        end
        load_tables || promote_orbit_table_cache(paths, key, log_path)
        if exit_code == 0
            prevariety_checkpoint_is_parseable(partial_path) ||
                error("restricted Gfan output is invalid; retained $partial_path")
            mv(partial_path, fan_path; force=true)
            return (fan_path=fan_path, integer_bits=bits,
                    table_path=paths.table, reused_tables=load_tables, log_path=log_path)
        end
        bits == 0 && error("restricted Gfan failed; see $log_path and $partial_path")
        println("  $bits-bit restricted intersection failed; retrying in arbitrary precision.")
        flush(stdout)
    end
    error("restricted Gfan produced no result")
end

"Cells, rays and lineality of a normalized fan, read once for all its restrictions."
function closed_restriction_data(normalized)
    rays = Matrix{QQFieldElem}(normalized.rays)
    # The orbit jobs keep the cells they read: Polymake does not return the memory of another reading.
    hasproperty(normalized, :cells) && return (cells=normalized.cells, rays=rays, lineality=normalized.lineality)
    return (cells=fan_cells(normalized.fan), rays=rays,
            lineality=Matrix{QQFieldElem}(matrix(QQ, rays_modulo_lineality(normalized.fan).lineality_basis)))
end

"Exact closed H-representation of one normalized min-convention fan cone."
closed_orbit_restriction(normalized, representative) =
    _closed_orbit_restriction(closed_restriction_data(normalized), representative)

_closed_orbit_restriction(data, representative) =
    closed_cone_restriction(-data.rays[data.cells[representative], :], data.lineality)

"Exact closed H-representation of the cone spanned by the rows of `rays` and the space of `lineality`."
function closed_cone_restriction(rays, lineality)
    cone = positive_hull(matrix(QQ, rays), matrix(QQ, lineality); non_redundant=true)
    inequalities = [_primitive_integer_row(-collect(normal_vector(facet)))
                    for facet in facets(Halfspace, cone)]
    equations = [_primitive_integer_row(collect(row)) for row in eachrow(
        Matrix{QQFieldElem}(transpose(kernel(matrix(QQ, vcat(lineality, rays)); side=:right))))]
    all(all(dot(ray, normal) >= 0 for ray in eachrow(rays)) &&
        all(iszero(dot(direction, normal)) for direction in eachrow(lineality))
        for normal in inequalities) || error("closed restriction inequality orientation failed")
    all(all(iszero(dot(row, normal)) for row in eachrow(vcat(lineality, rays)))
        for normal in equations) || error("closed restriction equations failed")
    return (ambient_dimension=size(rays, 2), inequalities=inequalities,
            equations=equations, strict=Vector{Vector{BigInt}}())
end

function closed_orbit_restriction_text(restriction)
    matrix_text(rows) = "{" * join(("(" * join(row, ',') * ")" for row in rows), ',') * "}"
    return "(" * string(restriction.ambient_dimension) * "," * matrix_text(restriction.inequalities) *
           "," * matrix_text(restriction.equations) * ",{})"
end

closed_orbit_restrictions_text(restrictions) =
    "{\n" * join((closed_orbit_restriction_text(restriction) for restriction in restrictions), ",\n") * "\n}\n"

"The restriction texts of a batch of cones, each given by its rays, as a worker computes them."
function closed_cone_restriction_texts(cones, lineality)
    lineality = QQ.(lineality)
    texts = [closed_orbit_restriction_text(closed_cone_restriction(QQ.(rays), lineality)) for rays in cones]
    # Julia does not count the memory of the Polymake cones and collects them
    # too rarely: uncollected, they held 30 kB per restriction.
    GC.gc()
    return texts
end

"""
Write a basis followed by any chosen exact closed previous-cone restrictions.

The workers compute them, `batch` cones at a time, and the restrictions of
`batches` batches are written before the next ones are computed.
"""
function write_closed_orbit_restricted_input(path, variable_names, basis, normalized,
                                             representatives;
                                             basis_path=path * ".basis.gfan", batch=500, batches=40)
    # Reading the cells of the fan for each representative was quadratic.
    data = closed_restriction_data(normalized)
    # Machine integers, as in the orbit jobs: as exact rationals the rays of
    # 20,000 cones cost this process half a gigabyte of numbers to collect.
    rays = Int.(data.rays)
    lineality = Rational{BigInt}.(data.lineality)
    cone_rays(representative) = -rays[data.cells[representative], :]
    write_basis_input(basis_path, variable_names, basis)
    progress = TaskProgress("closed orbit restrictions", length(representatives))
    # An interrupted run leaves one partial file, which the next one overwrites.
    partial = path * ".partial"
    # The workers are idle here, and one process took three hours for 880,000
    # cones. Written as they come: kept as integer rows, these took 10 GB.
    open(partial, "w") do io
        write(io, read(basis_path, String), "{\n")
        written = 0
        for chunk in Iterators.partition(representatives, batch * batches)
            cones = [cone_rays.(part) for part in Iterators.partition(chunk, batch)]
            for texts in pmap(part -> closed_cone_restriction_texts(part, lineality), cones), text in texts
                written > 0 && write(io, ",\n")
                write(io, text)
                written += 1
            end
            advance!(progress, length(chunk))
        end
        write(io, "\n}\n")
    end
    mv(partial, path; force=true)
    return (input_path=path, basis_path=basis_path, representatives=Int.(collect(representatives)))
end

function previous_orbit_round_fan(config, round_number)
    directory = joinpath(config.output, "round_" * lpad(round_number-1, 3, '0'))
    candidates = [joinpath(directory, "orbit_restricted", "expanded_prevariety.gfan"),
                  joinpath(directory, "prevariety.gfan")]
    for path in candidates
        prevariety_checkpoint_is_parseable(path) && return path
    end
    error("orbit restriction requires a completed previous prevariety in $directory")
end

"Collect in this process and in its workers: a Julia process that waits keeps its garbage."
collect_all_processes() = foreach(process -> remotecall_wait(GC.gc, process), procs())

"Certificate entries that a finished gfan run adds to those of its restriction input."
const GFAN_RESULT_KEYS = ("integer_bits", "restricted_fan_sha256", "native_intersection_table_sha256")

"""
Whether the certificate `earlier` was written for the restriction input on
disk and for these `inputs`: the restrictions are a function of the previous
fan, the group and the basis, so they need not be computed again.
"""
function certifies_restriction_input(earlier, inputs, basis_path, restriction_path)
    isfile(basis_path) && isfile(restriction_path) || return false
    all(get(earlier, key, nothing) == value for (key, value) in inputs) || return false
    haskey(earlier, "number_of_coarsened_previous_cells") && haskey(earlier, "number_of_representatives") ||
        return false
    return inputs["basis_input_sha256"] == orbit_stream_sha256(basis_path) &&
           get(earlier, "restriction_input_sha256", "") == orbit_stream_sha256(restriction_path)
end

"Compute the cell orbits of the previous fan, write one closed restriction per orbit and certify them."
function orbit_restriction_certificate(problem, config, export_data, basis, inputs, directory,
                                       previous_fan_path, previous_directory, basis_path, restriction_path)
    plain = reads_gfan_cells(export_data, problem.coarsen_prevariety)
    orbits = if plain
        gfan_orbit_jobs(previous_fan_path, export_data.symmetry_permutations)
    else
        previous_coarse = joinpath(previous_directory, "coarsened_prevariety.mrdi")
        source_marker = previous_coarse * ".source_sha256"
        coarse = if problem.coarsen_prevariety && isfile(previous_coarse) && isfile(source_marker) &&
                    strip(read(source_marker, String)) == inputs["source_fan_sha256"]
            load(previous_coarse)
        else
            min_fan = timed("Saving the previous fan") do
                save_lifted_fan(previous_fan_path, joinpath(directory, "previous_prevariety_min.mrdi"),
                                export_data.lift_indices).fan
            end
            # Restricting to the uncoarsened cones only gives gfan more restrictions.
            if problem.coarsen_prevariety
                result = coarsen_prevariety_fan(min_fan; strategy=EqualSpanCoarsening())
                result.support_preserved || error("previous support coarsening was not certified")
                result.fan
            else
                min_fan
            end
        end
        save(joinpath(directory, "previous_coarsened_prevariety.mrdi"), coarse)
        tree_orbit_jobs(coarse, export_data.symmetry_permutations; action=export_data.symmetry_action)
    end
    restricted_cells = length(orbits.normalized.cells)
    covered = sort!(reduce(vcat, orbits.symmetry.cell_orbits))
    covered == collect(1:restricted_cells) || error("previous cell orbit coverage is incomplete")
    save_symmetry_artifacts(joinpath(directory, "restriction_symmetry"), orbits.symmetry)
    plain || save(joinpath(directory, "restriction_normalized_fan.mrdi"), orbits.normalized.fan)
    save(joinpath(directory, "restriction_normalized_rays.mrdi"), orbits.normalized.rays)
    written = write_closed_orbit_restricted_input(restriction_path, export_data.variable_names,
        basis, orbits.normalized, [job.representative for job in orbits.jobs]; basis_path=basis_path)
    return merge(inputs, Dict{String, Any}(
        "restriction_input_sha256" => orbit_stream_sha256(restriction_path),
        "support_preserved" => true, "all_previous_cells_covered" => true,
        "previous_basis_subset" => true, "invariant_basis_closed" => true,
        "closed_restrictions" => true, "group_closed" => true,
        # The cells read from the file are its cones.
        "number_of_previous_cells" => plain ? restricted_cells :
            length(parse_gfan_fan(previous_fan_path).maximal_cones),
        "number_of_coarsened_previous_cells" => restricted_cells,
        "number_of_representatives" => length(written.representatives),
        "restriction_group_order" => length(orbits.symmetry.permutations),
    ))
end

"Refine every previous cell orbit with the same invariant new basis, then expand."
function run_orbit_restricted_prevariety(problem, config, export_data, basis;
                                         round_number::Int, worker::Int,
                                         previous_fan_path=nothing, release_workers=() -> nothing)
    round_number > 1 || error("orbit restriction needs a previous refinement round")
    isnothing(export_data.symmetry_permutations) && error("orbit restriction requires verified signed symmetry")
    export_data.lift_indices == collect(1:length(export_data.variable_names)) ||
        error("orbit restriction requires the provider's unreduced coordinates")
    # The fans and the checks of the previous round are no longer needed.
    collect_all_processes()
    directory = joinpath(config.output, "round_" * lpad(round_number, 3, '0'), "orbit_restricted")
    mkpath(directory)
    output_path = joinpath(directory, "expanded_prevariety.gfan")
    isfile(output_path) && mv(output_path, output_path * ".previous"; force=true)
    previous_fan_path = isnothing(previous_fan_path) ?
        previous_orbit_round_fan(config, round_number) : previous_fan_path
    previous_directory = joinpath(config.output, "round_" * lpad(round_number-1, 3, '0'))
    previous_basis_path = joinpath(previous_directory, "basis.gfan")
    previous_basis = gfan_basis_polynomials(previous_basis_path, export_data.variable_names)
    isnothing(previous_basis) && error("previous basis variables do not match orbit-restricted basis")
    signatures = remotecall_fetch((previous, current) -> (
        Main.parallel_trop_basis_support_signature(previous),
        Main.parallel_trop_basis_support_signature(current),
        Main.parallel_trop_basis_support_signature(Main.parallel_trop_symmetrize_polynomials(current))),
        worker, previous_basis, basis)
    issubset(Set(split(signatures[1], '\n')), Set(split(signatures[2], '\n'))) ||
        error("new orbit-restricted basis does not contain every previous hypersurface")
    signatures[2] == signatures[3] || error("orbit-restricted basis supports are not symmetry closed")
    # A cell of an uncoarsened previous fan lies in one cone of every earlier
    # hypersurface: only the polynomials added since can cut it. Gfan's
    # traversal costs in proportion to the polynomials it is given: by the 24
    # new of 133 polynomials, 608 cells took 6 CPU seconds instead of 34, for
    # the same cones. Coarsened cells cross earlier hypersurfaces: all are needed.
    refining = problem.coarsen_prevariety ? basis : setdiff(basis, previous_basis)
    isempty(refining) && error("orbit restriction found no polynomial added since the previous round")
    println("  Gfan refines the previous cells by $(length(refining)) of the $(length(basis)) polynomials.")
    flush(stdout)
    source_sha = orbit_stream_sha256(previous_fan_path)
    basis_path = joinpath(directory, "basis.gfan")
    restriction_path = joinpath(directory, "restricted_input.gfan")
    coordinate_actions = [symmetry_coordinate_permutation(permutation, length(export_data.variable_names),
                                                          export_data.symmetry_action)
                          for permutation in export_data.symmetry_permutations]
    coordinate_group_sha = bytes2hex(sha256(join((join(permutation, ",") for permutation in
        sort!(unique(Tuple.(coordinate_actions)))), '\n') * "\n"))
    certificate_path = joinpath(directory, "restriction_certificate.toml")
    earlier = isfile(certificate_path) ? TOML.parsefile(certificate_path) : Dict{String, Any}()
    # What the restrictions are computed from: the rest of the certificate describes them.
    inputs = Dict{String, Any}(
        "policy" => problem.coarsen_prevariety ? "all-previous-cell-orbits-closed-common-refinement-v1" :
                    "all-previous-cell-orbits-closed-refinement-by-new-polynomials-v1",
        "source_fan_sha256" => source_sha,
        "basis_input_sha256" => bytes2hex(sha256(basis_input_text(export_data.variable_names, refining))),
        "group_order" => length(export_data.symmetry_permutations),
        "coordinate_group_sha256" => coordinate_group_sha,
        "previous_basis_support_sha256" => bytes2hex(sha256(signatures[1])),
        "current_basis_support_sha256" => bytes2hex(sha256(signatures[2])),
    )
    # The orbits of 51 million cells and their restrictions take an hour and
    # 150 GB: on --resume, keep the input that an earlier run certified.
    prepared = config.resume && certifies_restriction_input(earlier, inputs, basis_path, restriction_path)
    certificate = if prepared
        println("  Reusing the restriction input certified for these inputs.")
        flush(stdout)
        Dict{String, Any}(key => value for (key, value) in earlier if !(key in GFAN_RESULT_KEYS))
    else
        orbit_restriction_certificate(problem, config, export_data, refining, inputs, directory,
                                      previous_fan_path, previous_directory, basis_path, restriction_path)
    end
    # Gfan can take hours here: on --resume, keep its output when the
    # certificate of the run that produced it agrees on every input.
    earlier_fan_path = joinpath(directory, "restricted_prevariety.gfan")
    reuse = config.resume && isfile(earlier_fan_path) &&
        all(get(earlier, key, nothing) == value for (key, value) in certificate) &&
        get(earlier, "restricted_fan_sha256", "") == orbit_stream_sha256(earlier_fan_path)
    reuse || orbit_write_toml(certificate_path, certificate)
    println("  Restricted intersection covers all $(certificate["number_of_coarsened_previous_cells"]) previous maximal cells through $(certificate["number_of_representatives"]) closed orbit representatives.")
    flush(stdout)
    # Only for a new input: the next process must find it and run gfan, not write it again.
    config.fresh_processes && !prepared && !reuse &&
        get(earlier, "restriction_input_sha256", "") != certificate["restriction_input_sha256"] &&
        leave_for_fresh_process("The restriction input of round $round_number is written.")
    # Gfan runs for hours beside this process, which then expands its output
    # alone: forty waiting workers held 100 GB that gfan needed.
    release_workers()
    GC.gc()
    native = if reuse
        println("  Reusing the restricted Gfan output certified for these inputs.")
        flush(stdout)
        (fan_path=earlier_fan_path, integer_bits=earlier["integer_bits"],
         table_path=orbit_table_cache_paths(directory, earlier["integer_bits"]).table)
    else
        run_native_orbit_intersection(config, basis_path, restriction_path, directory, length(refining))
    end
    certificate["integer_bits"] = native.integer_bits
    certificate["restricted_fan_sha256"] = orbit_stream_sha256(native.fan_path)
    certificate["native_intersection_table_sha256"] = isfile(native.table_path) ? orbit_stream_sha256(native.table_path) : "unavailable"
    orbit_write_toml(certificate_path, certificate)
    result = expand_restricted_prevariety(native.fan_path, output_path, coordinate_actions;
        previous_fan_path=previous_fan_path, provenance_path=joinpath(directory, "expansion_provenance.toml"),
        restriction_certificate_path=certificate_path,
        restriction_input_path=restriction_path, basis_input_path=basis_path, progress=progress_steps)
    result.support_checked || error("orbit expansion did not validate its support proof")
    return result
end
