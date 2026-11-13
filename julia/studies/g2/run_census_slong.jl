#!/usr/bin/env julia

"""
The tropical G2-Grassmannian over QQ in the coordinates of the long tripod
S_long of Section 8: the 18 Pluecker coordinates p_ij, ij not in {14, 25, 36},
and u = p25 - p36, v = p36 - p14, w = p14 - p25 in the slots of p14, p25, p36
(`--target=slong`, artifacts in data/g2/grassmannian_Slong). The target
`short_long` tropicalizes in the 24 coordinates of S_short and S_long together
(data/g2/grassmannian_Slong/short_long), and the stage `projections` of `slong`
compares both tropicalizations through it.

Both are tropicalized by `run_tropicalization` with the parallel engine: Gfan
prevarieties, Groebner certification of one cone per orbit of the Weyl group,
given as signed variable permutations. Points of trop_S_long are not tree
metrics, so rays are normalized by orthogonal projection, and there is no tree
stage.
"""

using Oscar
using TropicalGeometryG2
using TOML

include(joinpath(@__DIR__, "ideals", "G2GrassmannianSlongIdeal.jl"))
using .G2GrassmannianSlongIdeal
include(joinpath(@__DIR__, "ideals", "G2Symmetry.jl"))
using .G2Symmetry: G2_GROUP_CANDIDATES
include(joinpath(@__DIR__, "census", "g2_chevalley.jl"))

const PROVIDER = study_path(:g2, "ideals", "G2GrassmannianSlongIdeal.jl")

const TARGETS = Dict(
    "slong" => (
        directory=data_path("g2", "grassmannian_Slong"),
        work=work_path("g2", "grassmannian_Slong", "tropicalize"),
        name="G2/P2 in S_long", slug="g2_p2_slong", prefix="TGrass_G2_Slong",
        builder=:build_g2_grassmannian_slong_ideal, generators=:SLONG_WEYL_GENERATORS,
        label="tropical G2-Grassmannian over QQ in the coordinates of S_long",
    ),
    "short_long" => (
        directory=data_path("g2", "grassmannian_Slong", "short_long"),
        work=work_path("g2", "grassmannian_Slong", "short_long_tropicalize"),
        name="G2/P2 in S_short and S_long", slug="g2_p2_short_long", prefix="TGrass_G2_ShortLong",
        builder=:build_g2_grassmannian_short_long_ideal, generators=:SHORT_LONG_WEYL_GENERATORS,
        label="tropical G2-Grassmannian over QQ in the coordinates of S_short and S_long",
    ),
)

spec(target) = TARGETS[target.name]
file(target, name) = joinpath(spec(target).directory, name)
symmetry_files(target) = collect(symmetry_artifact_paths(file(target, "05_symmetry")))
manifest_file(target, stage) = stage_provenance_path(spec(target).directory, stage)
ideal_builder(target) = getproperty(G2GrassmannianSlongIdeal, spec(target).builder)

function record(target, stage, options; inputs, outputs, metadata=Dict{String, Any}())
    write_provenance(
        manifest_file(target, stage);
        study="G2", target=basename(spec(target).directory), stage=stage,
        inputs=inputs, outputs=outputs, configuration=census_configuration(options),
        metadata=metadata,
    )
end

fan_counts(fan) = (length(rays_modulo_lineality(fan).rays_modulo_lineality), n_maximal_cones(fan))

function slong_problem(name)
    s = TARGETS[name]
    return TropicalProblem(
        name=s.name, slug=s.slug,
        provider_file=PROVIDER, provider_module=:G2GrassmannianSlongIdeal, builder=s.builder,
        output_prefix=s.prefix, default_output=s.work, published_fan=joinpath(s.directory, "01_raw_fan.mrdi"),
        symmetry_candidates=s.generators, symmetry_action=:variables,
        default_engine=:parallel,
    )
end

function stage_tropicalize(target, options)
    s = spec(target)
    ideal_path, output = file(target, "00_ideal.mrdi"), file(target, "01_raw_fan.mrdi")
    mkpath(s.directory)
    (!isfile(ideal_path) || options.force) && save(ideal_path, ideal_builder(target)().ideal)
    problem = slong_problem(target.name)
    engine = options.engine == :auto ? :parallel : options.engine
    if !skip_current([output], options)
        run_tropicalization(problem, TropicalizationConfig(
            engine=engine, output=s.work, resume=options.resume, force=options.force,
            workers=options.workers, threads=options.threads, bits=options.bits,
            gfan_chunksize=options.gfan_chunksize, fresh_processes=options.fresh_processes,
        ))
    end
    isfile(output) || error("the tropicalization did not publish $output; run again with --resume")
    fan = load(output)
    rays, cells = fan_counts(fan)
    println("  $rays rays, $cells maximal cones; dimension $(dim(fan)), lineality $(lineality_dim(fan))")
    record(target, tropicalize, options; inputs=[ideal_path], outputs=[output],
           metadata=Dict("effective_engine" => string(engine), "rays" => rays, "maximal_cells" => cells,
                         "symmetry" => "Weyl group, signed variable permutations $(s.generators)"))
end

function stage_coarsen(target, options)
    input, output = file(target, "01_raw_fan.mrdi"), file(target, "02_coarsest_fan.mrdi")
    cells_path = file(target, "02_coarsest_cells.csv")
    isfile(input) || error("stage 1 artifact is missing: $input")
    if !skip_current([output, cells_path], options)
        coarsen_fan_file(input, output; strategy=SupportLinkCoarsening(), label="Trop($(spec(target).name))",
                         provenance_path=cells_path, report_path=coarsest_report_path(spec(target).directory))
    end
    rays, cells = fan_counts(load(output))
    println("  coarsest: $rays rays, $cells maximal cells")
    record(target, coarsen, options; inputs=[input], outputs=[output, cells_path],
           metadata=Dict("rays" => rays, "maximal_cells" => cells))
end

function stage_normalize(target, options)
    input = file(target, "02_coarsest_fan.mrdi")
    outputs = [file(target, "03_normalized_coarsest_fan.mrdi"), file(target, "03_normalized_rays.mrdi")]
    isfile(input) || error("stage 2 artifact is missing: $input")
    if !skip_current(outputs, options)
        # No tree gauge: each ray is the primitive integer vector of its
        # orthogonal projection onto the complement of the lineality space.
        normalized = normalize_fan(load(input))
        save(outputs[1], normalized.fan)
        save(outputs[2], normalized.rays)
    end
    rays, cells = fan_counts(load(outputs[1]))
    record(target, normalize, options; inputs=[input], outputs=outputs,
           metadata=Dict("rays" => rays, "maximal_cells" => cells,
                         "gauge" => "orthogonal projection modulo lineality, primitive integer"))
end

"The 12 Weyl group elements as unsigned permutations of the coordinates of `target`."
function weyl_coordinate_permutations(target)
    generators = getproperty(G2GrassmannianSlongIdeal, spec(target).generators)
    n = length(first(first(generators)))
    group = TropicalGeometryG2.signed_permutation_group(generators, n)
    return group.permutations
end

function stage_symmetry(target, options)
    fan_path, rays_path = file(target, "03_normalized_coarsest_fan.mrdi"), file(target, "03_normalized_rays.mrdi")
    outputs = symmetry_files(target)
    if !skip_current(outputs, options)
        fan, rays = load(fan_path), load(rays_path)
        lineality = rays_modulo_lineality(fan).lineality_basis
        # In the S_long slots, u, v, w are permuted as p14, p25, p36, so the
        # taxon permutations act as on Pluecker coordinates. The 24 combined
        # coordinates have no taxa: there the Weyl group acts on variables.
        result = target.name == "slong" ?
            taxon_symmetry(rays, TropicalGeometryG2.fan_cells(fan); candidates=G2_GROUP_CANDIDATES, lineality=lineality) :
            taxon_symmetry(rays, TropicalGeometryG2.fan_cells(fan); candidates=weyl_coordinate_permutations(target),
                           lineality=lineality, action=:variables)
        length(result.permutations) >= 12 || error("the Weyl group does not preserve the fan")
        TropicalGeometryG2.save_symmetry_artifacts(file(target, "05_symmetry"), result)
    end
    result = TropicalGeometryG2.load_symmetry_artifacts(file(target, "05_symmetry"))
    println("  group order $(length(result.permutations)), $(length(result.ray_orbits)) ray orbits, ",
            "$(length(result.cell_orbits)) maximal-cell orbits")
    record(target, symmetry, options; inputs=[fan_path, rays_path], outputs=outputs,
           metadata=Dict("group_order" => length(result.permutations), "cell_orbits" => length(result.cell_orbits)))
end

function stage_degenerations(target, options)
    target.name == "slong" || error("stage 6 is defined for the target slong")
    fan_path, rays_path = file(target, "03_normalized_coarsest_fan.mrdi"), file(target, "03_normalized_rays.mrdi")
    ideal_path, prefix = file(target, "00_ideal.mrdi"), file(target, "05_symmetry")
    output = file(target, "06_canonical_interior_degenerations.csv")
    symmetry_result = TropicalGeometryG2.load_symmetry_artifacts(prefix)
    expected_orbits = length(symmetry_result.cell_orbits)
    recorded_orbits() = Set(parse(Int, first(split(line, ',', limit=2)))
                            for line in readlines(output)[2:end] if !isempty(line))
    continuing = options.resume && !options.force && isfile(output)
    completed = continuing ? recorded_orbits() : Set{Int}()
    if options.force || !isfile(output) || (options.resume && length(completed) < expected_orbits)
        config = DegenerationConfig(
            expected_orbits=expected_orbits,
            coordinate_images=slong_coordinate_images,
            chevalley_test=g2_chevalley_simplex,
        )
        if options.force
            objects = TropicalGeometryG2.degeneration_object_directory(output)
            isdir(objects) && rm(objects; recursive=true)
        end
        open(output, continuing ? "a" : "w") do io
            continuing || TropicalGeometryG2.write_degenerations_header(io)
            classify_degenerations(
                load(ideal_path), load(fan_path), load(rays_path), symmetry_result;
                config=config, skip_orbits=completed,
                on_record=record -> TropicalGeometryG2.write_degenerations_record(io, output, record),
            )
        end
        recorded_orbits() == Set(1:expected_orbits) || error("degeneration census is incomplete")
    else
        skip_current([output], options)
    end
    record(target, degenerations, options;
           inputs=vcat([ideal_path, fan_path, rays_path], symmetry_files(target)),
           outputs=degeneration_artifact_paths(output),
           metadata=Dict("weight_policy" => "equal sum of normalized coarsest-cell extreme rays",
                         "cell_orbits" => expected_orbits,
                         "monomial_freeness" => "verified by torus saturation"))
end

function stage_summarize(target, options)
    directory = spec(target).directory
    inputs = vcat([file(target, "01_raw_fan.mrdi"), file(target, "03_normalized_coarsest_fan.mrdi")],
                  symmetry_files(target), filter(isfile, [file(target, "06_canonical_interior_degenerations.csv")]))
    all(isfile, inputs) || error("the summary requires the artifacts of stages 1, 3 and 5")
    output = summary_path(directory)
    skip_current([output], options) ||
        write_census_summary(directory; label=spec(target).label, stage1_fan=file(target, "01_raw_fan.mrdi"))
    record(target, summarize, options; inputs=inputs, outputs=[output])
end

"""
The fan, rays, lineality, cells, cell orbits and toric cells of a census
directory, and the polyhedral cone of each maximal cell.
"""
function census_fan(directory)
    fan = load(joinpath(directory, "03_normalized_coarsest_fan.mrdi"))
    rays = Matrix{QQFieldElem}(load(joinpath(directory, "03_normalized_rays.mrdi")))
    lineality = Matrix{QQFieldElem}(matrix(QQ, rays_modulo_lineality(fan).lineality_basis))
    cells = TropicalGeometryG2.fan_cells(fan)
    orbits = TropicalGeometryG2.load_symmetry_artifacts(joinpath(directory, "05_symmetry")).cell_orbits
    orbit_of = zeros(Int, length(cells))
    for (index, orbit) in enumerate(orbits)
        orbit_of[orbit] .= index
    end
    census = joinpath(directory, "06_canonical_interior_degenerations.csv")
    toric = Set{Int}()
    if isfile(census)
        rows = TropicalGeometryG2.csv_rows(census)
        sort([parse(Int, row["orbit"]) for row in rows]) == collect(eachindex(orbits)) ||
            error("the stage-6 census in $directory is incomplete")
        toric = Set(TropicalGeometryG2.degeneration_cells(rows, "toric"))
    end
    cones = [positive_hull(rays[cell, :], lineality) for cell in cells]
    return (fan=fan, rays=rays, lineality=lineality, cells=cells, orbits=orbits, orbit_of=orbit_of,
            toric=toric, cones=cones, dimension=dim(fan) - size(lineality, 1))
end

"""
The image of the cell with ray rows `rays` under the coordinate projection
`columns`, in the census fan `target`: its dimension modulo lineality and, if
it is full, the maximal cell containing it and whether it is that whole cell.
"""
function projected_cell(rays, columns, target)
    image = rays[:, columns]
    lineality_rank = rank(matrix(QQ, target.lineality))
    dimension = rank(matrix(QQ, vcat(image, target.lineality))) - lineality_rank
    dimension == target.dimension || return (dimension=dimension, cell=0, orbit=0, equal=false)
    point = vec(sum(image; dims=1))
    containing = [index for (index, cone) in enumerate(target.cones)
                  if point in cone && all(collect(row) in cone for row in eachrow(image))]
    length(containing) == 1 || error("the image of a cell lies in $(length(containing)) target cells")
    cell = only(containing)
    image_cone = positive_hull(image, target.lineality)
    equal = all(collect(row) in image_cone for row in eachrow(target.rays[target.cells[cell], :]))
    return (dimension=dimension, cell=cell, orbit=target.orbit_of[cell], equal=equal)
end

"""
Whether the cones with ray rows `pieces` (each a matrix) cover the maximal
cell `cell` of `target` with disjoint interiors: they lie in it, meet in lower
dimension, and their exact volumes, truncated by a linear form positive on
the cell, add up to that of the cell. Computed in coordinates on the span of
the cell modulo lineality.
"""
function covers_cell(pieces, cell, target)
    annihilator = Matrix{QQFieldElem}(kernel(matrix(QQ, target.lineality); side=:right))
    cell_rays = target.rays[target.cells[cell], :] * annihilator
    _, reduced = rref(matrix(QQ, cell_rays))
    basis = Matrix{QQFieldElem}(reduced)[1:target.dimension, :]
    coordinates(rows) = Matrix{QQFieldElem}(solve(matrix(QQ, basis), matrix(QQ, rows * annihilator); side=:left))
    whole = positive_hull(coordinates(target.rays[target.cells[cell], :]))
    parts = [positive_hull(coordinates(piece)) for piece in pieces]
    all(issubset(part, whole) for part in parts) || return false
    all(dim(intersect(parts[i], parts[j])) < target.dimension for i in eachindex(parts) for j in 1:(i - 1)) ||
        return false
    # Oscar's facets of a cone are the halfspaces a x <= 0, so minus their sum is positive on the cell.
    positive = -sum(normal_vector(facet) for facet in facets(whole))
    all(sum(positive .* collect(ray)) > 0 for ray in rays(whole)) || error("the truncating form is not positive")
    truncated(cone) = volume(intersect(polyhedron(cone), polyhedron(reshape(collect(positive), 1, :), [1])))
    return sum(truncated, parts) == truncated(whole)
end

"""
Compare trop(G2/P2) in S_short (the Pluecker census in data/g2/grassmannian)
with trop in S_long through the tropicalization in both: one maximal cell per
orbit of the latter is projected to each. Projections commute with the Weyl
group, so the cells of an orbit O over one cell of a target orbit T number
|O|/|T|.
"""
function stage_projections(target, options)
    target.name == "slong" || error("the projections stage belongs to the target slong")
    both_directory, long_directory = TARGETS["short_long"].directory, TARGETS["slong"].directory
    short_directory = data_path("g2", "grassmannian")
    outputs = [file(target, "07_projections.csv"), file(target, "07_projections.toml")]
    inputs = vcat(
        [joinpath(directory, name) for directory in (both_directory, short_directory, long_directory)
         for name in ("03_normalized_coarsest_fan.mrdi", "03_normalized_rays.mrdi",
                      symmetry_artifact_paths("05_symmetry")...)],
        # the toric columns come from the stage-6 census of both targets
        [joinpath(directory, "06_canonical_interior_degenerations.csv") for directory in (short_directory, long_directory)],
    )
    all(isfile, inputs) || error("the projections need stages 3 and 5 of short_long, and stages 3, 5 and 6 " *
                                 "of slong and of the Pluecker census")
    if !skip_current(outputs, options)
        both, short, long = census_fan.((both_directory, short_directory, long_directory))
        short_columns = collect(1:21)
        long_columns = [k in SLONG_POSITIONS ? 21 + findfirst(==(k), SLONG_POSITIONS) : k for k in 1:21]
        for (columns, fan) in ((short_columns, short), (long_columns, long))
            rank(matrix(QQ, vcat(both.lineality[:, columns], fan.lineality))) == rank(matrix(QQ, fan.lineality)) ||
                error("a projection does not map lineality into lineality")
        end
        targets = Dict("short" => (short_columns, short), "long" => (long_columns, long))
        rows = Dict{String, Any}[]
        for (index, orbit) in enumerate(both.orbits)
            representative = first(orbit)
            row = Dict{String, Any}("orbit" => index, "size" => length(orbit), "representative" => representative)
            for (name, (columns, fan)) in targets
                image = projected_cell(both.rays[both.cells[representative], :], columns, fan)
                row["$(name)_dimension"] = image.dimension
                row["$(name)_cell"] = image.cell
                row["$(name)_orbit"] = image.orbit
                row["$(name)_equal"] = image.equal
                row["$(name)_toric"] = image.cell > 0 && image.cell in fan.toric
            end
            push!(rows, row)
        end
        header = ["orbit", "size", "representative",
                  [name * suffix for name in ("short", "long") for suffix in ("_dimension", "_cell", "_orbit", "_equal", "_toric")]...]
        open(outputs[1], "w") do io
            println(io, join(header, ','))
            for row in rows
                println(io, join((row[column] for column in header), ','))
            end
        end
        summary = Dict{String, Any}("both" => Dict{String, Any}(
            "maximal_cells" => length(both.cells), "maximal_cell_orbits" => length(both.orbits),
            "rays" => size(both.rays, 1)))
        for (name, (columns, fan)) in targets
            vertical = filter(row -> row["$(name)_orbit"] == 0, rows)
            lying = Dict{Int, Vector{Any}}()
            for row in rows
                row["$(name)_orbit"] > 0 && push!(get!(lying, row["$(name)_orbit"], Any[]), row)
            end
            sort(collect(keys(lying))) == collect(eachindex(fan.orbits)) ||
                error("the projections miss maximal-cell orbits of the $name fan")
            per_orbit = Dict{String, Any}[]
            for (orbit, members) in sort!(collect(lying); by=first)
                over = sum(row["size"] for row in members) // length(fan.orbits[orbit])
                denominator(over) == 1 || error("orbit sizes are not equivariant")
                if over > 1
                    # The pieces over one cell of the orbit must cover it.
                    representative = first(fan.orbits[orbit])
                    pieces = [both.rays[both.cells[cell], columns] for row in members for cell in both.orbits[row["orbit"]]
                              if projected_cell(both.rays[both.cells[cell], :], columns, fan).cell == representative]
                    length(pieces) == over || error("$(length(pieces)) cells lie over cell $representative, expected $over")
                    covers_cell(pieces, representative, fan) ||
                        error("the cells over cell $representative of the $name fan do not cover it")
                end
                push!(per_orbit, Dict{String, Any}(
                    "orbit" => orbit, "orbit_size" => length(fan.orbits[orbit]),
                    "cells_over_each" => Int(over), "toric" => first(fan.orbits[orbit]) in fan.toric,
                    "lifted_bijectively" => over == 1 && only(members)["$(name)_equal"],
                    "source_orbits" => [row["orbit"] for row in members],
                ))
            end
            summary[name] = Dict{String, Any}(
                "maximal_cells" => length(fan.cells), "maximal_cell_orbits" => length(fan.orbits),
                "vertical_source_cells" => sum((row["size"] for row in vertical); init=0),
                "vertical_source_orbits" => length(vertical),
                "cells_lifted_bijectively" => sum((o["orbit_size"] for o in per_orbit if o["lifted_bijectively"]); init=0),
                "orbits_lifted_bijectively" => count(o -> o["lifted_bijectively"], per_orbit),
                "orbits" => per_orbit,
            )
        end
        open(io -> TOML.print(io, summary; sorted=true), outputs[2], "w")
    end
    summary = TOML.parsefile(outputs[2])
    for name in ("short", "long")
        s = summary[name]
        println("  $name: $(s["maximal_cells"]) cells; $(s["cells_lifted_bijectively"]) lifted bijectively; ",
                "$(s["vertical_source_cells"]) cells of the 24-coordinate fan are vertical")
    end
    record(target, projections, options; inputs=inputs, outputs=outputs,
           metadata=Dict("source" => "short_long/03_normalized_coarsest_fan.mrdi",
                         "short_columns" => "1:21", "long_columns" => "p_ij, with u, v, w = 22, 23, 24 in the slots of 14, 25, 36"))
end

function target_definition(name)
    callbacks = Dict{CensusStage, Function}(
        tropicalize => stage_tropicalize, coarsen => stage_coarsen, normalize => stage_normalize,
        symmetry => stage_symmetry, summarize => stage_summarize,
    )
    name == "slong" && (callbacks[degenerations] = stage_degenerations)
    name == "slong" && (callbacks[projections] = stage_projections)
    path(n) = joinpath(TARGETS[name].directory, n)
    artifacts = Dict(
        tropicalize => [path("01_raw_fan.mrdi")],
        coarsen => [path("02_coarsest_fan.mrdi"), path("02_coarsest_cells.csv")],
        normalize => [path("03_normalized_coarsest_fan.mrdi"), path("03_normalized_rays.mrdi")],
        symmetry => collect(symmetry_artifact_paths(path("05_symmetry"))),
        summarize => [summary_path(TARGETS[name].directory)],
    )
    name == "slong" && (artifacts[degenerations] = [path("06_canonical_interior_degenerations.csv")])
    return StudyTarget(name=name, stages=Set(keys(callbacks)), callbacks=callbacks, artifacts=artifacts,
                       tropical_problem=() -> slong_problem(name))
end

const SLONG_CENSUS = StudyCensus(
    name="G2 S_long",
    targets=Dict(name => target_definition(name) for name in keys(TARGETS)),
)

run_census(SLONG_CENSUS)
