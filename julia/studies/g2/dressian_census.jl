# Census of the two Dressians of G2: the tropical prevarieties of the defining
# generators of the Grassmannian (`dressian`) and of the flag variety
# (`flag_dressian`). The generators have coefficients 1 and -1, so both fans
# are the same over every field. They are computed and stored once, in
# data/g2/dressian and data/g2/flag_dressian, for the censuses over QQ and
# GF(2). Neither fan is pure, so the coarsest fan comes from
# `coarsen_prevariety_fan`; the later stages are those of the tropical varieties.

include(joinpath(@__DIR__, "supplemental", "compute_flag_prevariety.jl"))

const DRESSIAN_SYMMETRY_SUFFIXES = (
    "_taxon_permutations.mrdi", "_ray_actions.mrdi", "_cell_actions.mrdi", "_cell_orbits.csv",
)
# Rays and maximal cones of the coarsest fan, and the orbits of its faces by dimension.
const DRESSIAN_TARGETS = Dict(
    "dressian" => (label="G2-Dressian", flag=false, coarsest=(51, 239), orbit_f_vector=[11, 33, 31, 1]),
    "flag_dressian" => (label="G2 flag Dressian", flag=true, coarsest=(141, 2627),
                        orbit_f_vector=[27, 173, 388, 274, 6]),
)

dressian_path(target, name) = joinpath(target_directory(target.name), name)

function record_dressian_stage(target, stage, options; inputs, outputs, metadata)
    write_provenance(
        manifest_file(target.name, stage);
        study="G2", target=target.name, stage=stage,
        inputs=inputs, outputs=outputs,
        configuration=census_configuration(options),
        metadata=metadata,
    )
end

"Every field gives the same prevariety when the coefficients are 1 and -1."
function require_unit_coefficients(polynomials)
    all(isone(c) || isone(-c) for f in polynomials for c in coefficients(f)) ||
        error("a defining generator has a coefficient other than 1 and -1")
end

"The f-vector modulo lineality, starting with the zero cone."
function quotient_f_vector(fan)
    counts, lineality = Int.(f_vector(fan)), lineality_dim(fan)
    return iszero(lineality) ? vcat([1], counts) : counts[lineality:end]
end

dressian_fan_counts(fan) = Dict{String, Any}(
    "rays" => length(rays_modulo_lineality(fan).rays_modulo_lineality),
    "maximal_cells" => n_maximal_cones(fan),
    "dimension" => dim(fan),
    "lineality_dimension" => lineality_dim(fan),
    "pure" => is_pure(fan),
)

function stage_dressian_tropicalize(target, options)
    output = dressian_path(target, "01_prevariety_fan.mrdi")
    if DRESSIAN_TARGETS[target.name].flag
        require_unit_coefficients(build_full_flag_g2_ideal().generators)
        # The Gfan computation records its own provenance.
        skip_current([output], options) || FlagDressianPrevariety.compute_flag_prevariety(options)
        return
    end
    data = build_g2_grassmannian_ideal()
    require_unit_coefficients(gens(data.ideal))
    computed = !skip_current([output], options)
    if computed
        prevariety = tropical_prevariety(gens(data.ideal), tropical_semiring_map(QQ, min))
        mkpath(dirname(output))
        save(output, TropicalCensus.as_polyhedral_fan(prevariety))
    end
    metadata = dressian_fan_counts(load(output))
    metadata["defining_generators"] = ngens(data.ideal)
    metadata["coefficients"] = "1 and -1 only: the prevariety is the same over every field"
    metadata["reused_published_artifacts"] = !computed
    record_dressian_stage(target, tropicalize, options;
        inputs=filter(isfile, [study_path(:g2, "ideals", "G2GrassmannianIdeal.jl"), base_file("00_ideal.mrdi")]),
        outputs=[output], metadata=metadata)
end

function write_dressian_coarsening_report(path, label, input, output, source, result)
    line(name, fan) = "| $name | $(length(rays_modulo_lineality(fan).rays_modulo_lineality)) | " *
        "$(n_maximal_cones(fan)) | `($(join(quotient_f_vector(fan), ", ")))` |"
    open(path, "w") do io
        println(io, "# Coarsest fan structure of the $label")
        println(io)
        println(io, "Input: `$(relpath(input, dirname(path)))`")
        println(io)
        println(io, "Output: `$(relpath(output, dirname(path)))`")
        println(io)
        println(io, "The f-vectors count cones modulo the common $(result.lineality_dimension)-dimensional lineality space; the first entry is the zero cone.")
        println(io)
        println(io, "| fan | rays | maximal cones | quotient f-vector |")
        println(io, "|---|---:|---:|---|")
        println(io, line("generator prevariety", source))
        println(io, line("coarsest fan", result.fan))
        println(io)
        println(io, "| affine dimension | input maximal cones | output maximal cones |")
        println(io, "|---:|---:|---:|")
        for stratum in result.strata
            println(io, "| $(stratum.dimension) | $(stratum.input_cones) | $(stratum.output_cones) |")
        end
        println(io)
        println(io, "Maximal cones are grouped by their linear span, and each convex group is merged (`coarsen_prevariety_fan`, equal-span strategy). ",
                "The support is preserved, $(result.rolled_back_merges) merges were rolled back, and the output has $(result.open_equal_span_walls) open equal-span walls, ",
                "so no two of its maximal cones can be merged.")
    end
end

function stage_dressian_coarsen(target, options)
    spec = DRESSIAN_TARGETS[target.name]
    input = dressian_path(target, "01_prevariety_fan.mrdi")
    output, cells_path = dressian_path(target, "02_coarsest_fan.mrdi"), dressian_path(target, "02_coarsest_cells.csv")
    report = coarsest_report_path(target_directory(target.name))
    isfile(input) || error("the generator prevariety is missing: $input")
    if !skip_current([output, cells_path, report], options)
        source = load(input)
        result = coarsen_prevariety_fan(source; strategy=EqualSpanCoarsening())
        result.support_preserved && result.coarsest_certified ||
            error("the coarsening of the $(spec.label) is not certified")
        # Rebuild the fan, so that the saved file holds no cached polymake properties.
        data = rays_modulo_lineality(result.fan)
        coarsest = polyhedral_fan(
            maximal_cones(IncidenceMatrix, result.fan),
            matrix(QQ, data.rays_modulo_lineality), matrix(QQ, data.lineality_basis);
            non_redundant=true,
        )
        save(output, coarsest)
        dimensions = dim.(collect(maximal_cones(result.fan)))
        open(cells_path, "w") do io
            println(io, "coarsest_cell,dimension,original_cells,original_cell_count,coarsest_rays,ray_count")
            for (cell, rays) in enumerate(TropicalCensus.fan_cells(result.fan))
                originals = result.original_cells[cell]
                println(io, "$cell,$(dimensions[cell]),\"$(join(originals, ';'))\",$(length(originals)),",
                        "\"$(join(rays, ';'))\",$(length(rays))")
            end
        end
        write_dressian_coarsening_report(report, spec.label, input, output, source, result)
    end
    # The f-vectors are in the report: polymake needs minutes for those of the flag.
    fan = load(output)
    metadata = dressian_fan_counts(fan)
    metadata["open_equal_span_walls"] = open_equal_span_walls(fan)
    counts = (metadata["rays"], metadata["maximal_cells"])
    counts == spec.coarsest || error("$(spec.label) count regression: expected $(spec.coarsest), got $counts")
    record_dressian_stage(target, coarsen, options;
        inputs=[input], outputs=[output, cells_path, report], metadata=metadata)
end

function stage_dressian_normalize(target, options)
    input = dressian_path(target, "02_coarsest_fan.mrdi")
    outputs = [dressian_path(target, "03_normalized_coarsest_fan.mrdi"), dressian_path(target, "03_normalized_rays.mrdi")]
    isfile(input) || error("stage 2 artifact is missing: $input")
    if !skip_current(outputs, options)
        normalized = DRESSIAN_TARGETS[target.name].flag ?
            normalize_fan(load(input), G2_GAUGE, 8:28; anchor_coordinate=7) :
            normalize_fan(load(input), G2_GAUGE)
        save(outputs[1], normalized.fan)
        save(outputs[2], normalized.rays)
    end
    source, normalized = dressian_fan_counts(load(input)), dressian_fan_counts(load(outputs[1]))
    source == normalized || error("tree gauge normalization changed the Dressian fan counts or dimensions")
    record_dressian_stage(target, normalize, options; inputs=[input], outputs=outputs, metadata=normalized)
end

function stage_dressian_trees(target, options)
    inputs = [dressian_path(target, "03_normalized_coarsest_fan.mrdi"), dressian_path(target, "03_normalized_rays.mrdi")]
    outputs = [dressian_path(target, "04_canonical_cell_weights.mrdi"),
               dressian_path(target, "04_normalized_cell_trees.mrdi"), dressian_path(target, "04_trees.csv")]
    all(isfile, inputs) || error("stage 3 artifacts are missing")
    if !skip_current(outputs, options)
        fan, rays = load(inputs[1]), load(inputs[2])
        weights = canonical_cell_weights(fan, rays)
        cell_trees = DRESSIAN_TARGETS[target.name].flag ?
            flag_normalized_cell_trees(fan, weights) : normalized_cell_trees(weights, G2_GAUGE)
        TropicalCensus.save_unless_current(outputs[1], TropicalCensus.rows_matrix(weights; ring=QQ))
        TropicalCensus.save_unless_current(outputs[2], TropicalCensus.rows_matrix(cell_trees; ring=QQ))
        DRESSIAN_TARGETS[target.name].flag ?
            write_flag_tree_census(outputs[3], rays, cell_trees) :
            TropicalCensus.write_tree_census(outputs[3], rays, cell_trees)
    end
    record_dressian_stage(target, trees, options;
        inputs=inputs, outputs=outputs, metadata=dressian_fan_counts(load(inputs[1])))
end

"Faces and their orbits under the ray actions, in each dimension modulo lineality."
function dressian_face_orbits(fan, ray_actions)
    lineality = lineality_dim(fan)
    return map((lineality + 1):dim(fan)) do dimension
        incidence = IncidenceMatrix(cones(fan, dimension))
        faces = [sort!(collect(row(incidence, index))) for index in 1:nrows(incidence)]
        face_index = Dict{Vector{Int}, Int}(face => index for (index, face) in enumerate(faces))
        actions = [TropicalCensus.induced_cell_action(action, faces, face_index) for action in ray_actions]
        any(isnothing, actions) && error("the taxon group does not permute the faces of dimension $dimension")
        (dimension=dimension - lineality, faces=length(faces), orbits=length(action_orbits(actions, length(faces))))
    end
end

function stage_dressian_symmetry(target, options)
    spec = DRESSIAN_TARGETS[target.name]
    inputs = [dressian_path(target, "03_normalized_coarsest_fan.mrdi"), dressian_path(target, "03_normalized_rays.mrdi")]
    all(isfile, inputs) || error("stage 3 artifacts are missing")
    prefix, face_orbits = dressian_path(target, "05_symmetry"), dressian_path(target, "05_symmetry_face_orbits.csv")
    outputs = vcat([prefix * suffix for suffix in DRESSIAN_SYMMETRY_SUFFIXES], [face_orbits])
    if !skip_current(outputs, options)
        fan, rays = load(inputs[1]), load(inputs[2])
        action = spec.flag ? flag_taxon_symmetry(fan, rays, G2_GROUP_CANDIDATES) :
            taxon_symmetry(rays, TropicalCensus.fan_cells(fan);
                candidates=G2_GROUP_CANDIDATES, lineality=rays_modulo_lineality(fan).lineality_basis)
        length(action.permutations) == 12 ||
            error("expected the order-12 G2 taxon group, found $(length(action.permutations))")
        TropicalCensus.save_symmetry_artifacts(prefix, action)
        open(face_orbits, "w") do io
            println(io, "dimension_modulo_lineality,faces,orbits")
            for record in dressian_face_orbits(fan, action.ray_actions)
                println(io, "$(record.dimension),$(record.faces),$(record.orbits)")
            end
        end
    end
    action = TropicalCensus.load_symmetry_artifacts(prefix)
    counts = [parse.(Int, split(line, ',')) for line in readlines(face_orbits)[2:end]]
    orbit_f_vector = [record[3] for record in counts]
    orbit_f_vector == spec.orbit_f_vector ||
        error("$(spec.label) orbit regression: expected $(spec.orbit_f_vector), got $orbit_f_vector")
    record_dressian_stage(target, symmetry, options; inputs=inputs, outputs=outputs,
        metadata=Dict{String, Any}(
            "group_order" => length(action.permutations),
            "ray_orbits" => length(action.ray_orbits),
            "cell_orbits" => length(action.cell_orbits),
            "f_vector_modulo_lineality" => [record[2] for record in counts],
            "orbit_f_vector_modulo_lineality" => orbit_f_vector,
        ))
end

function dressian_target_definition(name)
    file(artifact) = joinpath(target_directory(name), artifact)
    callbacks = Dict(
        tropicalize => stage_dressian_tropicalize,
        coarsen => stage_dressian_coarsen,
        normalize => stage_dressian_normalize,
        trees => stage_dressian_trees,
        symmetry => stage_dressian_symmetry,
        summarize => stage_summarize,
    )
    artifacts = Dict(
        tropicalize => [file("01_prevariety_fan.mrdi")],
        coarsen => [file("02_coarsest_fan.mrdi"), file("02_coarsest_cells.csv")],
        normalize => [file("03_normalized_coarsest_fan.mrdi"), file("03_normalized_rays.mrdi")],
        trees => [file("04_canonical_cell_weights.mrdi"), file("04_normalized_cell_trees.mrdi"), file("04_trees.csv")],
        symmetry => vcat([file("05_symmetry" * suffix) for suffix in DRESSIAN_SYMMETRY_SUFFIXES],
                         [file("05_symmetry_face_orbits.csv")]),
        summarize => [summary_path(target_directory(name))],
    )
    return StudyTarget(name=name, stages=Set(keys(callbacks)), callbacks=callbacks, artifacts=artifacts)
end
