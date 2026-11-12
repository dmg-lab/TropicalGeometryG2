#!/usr/bin/env julia

using Oscar
using TropicalCensus
using TOML

include(joinpath(@__DIR__, "ideals", "G2GrassmannianIdeal.jl"))
include(joinpath(@__DIR__, "ideals", "FlagG2Ideal.jl"))
using .G2GrassmannianIdeal
using .FlagG2Ideal

const G2_ROOT_DIRECTIONS = [
     1  -1   0  -1   1   0  0
    -1   0   1   1   0  -1  0
     0   1  -1   0  -1   1  0
]
const G2_GAUGE = TreeGauge(7, G2_ROOT_DIRECTIONS)
const G2_GROUP_CANDIDATES = [vcat(permutation, [7]) for permutation in all_permutations(collect(1:6))]
const G2_FORM_SUPPORT = Set([
    (1, 2, 3), (4, 5, 6), (1, 4, 7), (2, 5, 7), (3, 6, 7),
])

include(joinpath(@__DIR__, "flag_census_helpers.jl"))
include(joinpath(@__DIR__, "dressian_census.jl"))
include(joinpath(@__DIR__, "zenodo_census.jl"))
using .G2Zenodo: zenodo_target_definition

# Each target keeps its artifacts, coarsening report and provenance in one directory.
const G2_TARGET_DIRECTORIES = Dict(
    "base" => "grassmannian", "flag" => "flag",
    "dressian" => "dressian", "flag_dressian" => "flag_dressian",
)
target_directory(target) = data_path("g2", G2_TARGET_DIRECTORIES[target])
base_file(name) = data_path("g2", "grassmannian", name)
flag_file(name) = data_path("g2", "flag", name)
manifest_file(target, stage) = stage_provenance_path(target_directory(target), stage)

function skip_current(outputs, options)
    current = all(isfile, outputs)
    current && !options.force && println("  using existing canonical artifact(s)")
    return current && !options.force
end

function save_ideal_if_needed(path, builder, options)
    if !isfile(path) || options.force
        mkpath(dirname(path))
        save(path, builder().ideal)
    end
end

function g2_problem(target)
    if target == "base"
        return TropicalProblem(
            name="G2/P2",
            slug="g2_p2",
            provider_file=study_path(:g2, "ideals", "G2GrassmannianIdeal.jl"),
            provider_module=:G2GrassmannianIdeal,
            builder=:build_g2_grassmannian_ideal,
            output_prefix="TGrass_G2",
            default_output=work_path("g2", "grassmannian", "tropicalize"),
            published_fan=base_file("01_raw_fan.mrdi"),
            published_complex=base_file("01_variety_complex.mrdi"),
            default_engine=:oscar,
        )
    end
    return TropicalProblem(
        name="G2/B",
        slug="g2_flag",
        provider_file=study_path(:g2, "ideals", "FlagG2Ideal.jl"),
        provider_module=:FlagG2Ideal,
        builder=:build_full_flag_g2_ideal,
        output_prefix="TFlag_G2",
        default_output=work_path("g2", "flag", "tropicalize_symmetry_full28"),
        published_fan=flag_file("01_raw_fan.mrdi"),
        symmetry_candidates=:G2_FLAG_PERMUTATIONS,
        symmetry_coordinate_images=:flag_coordinate_images,
        certification_strategy=:balanced_support,
        purity_verifier=:verify_flag_tropical_purity,
        whole_interior_pruning=true,
        # The restricted Gfan expansion currently requires one-dimensional
        # lineality; the full flag fan has four-dimensional lineality.
        orbit_restricted_intersection=false,
        default_engine=:parallel,
        # The same ray and cell numbering on every run.
        canonical_gfan_order=true,
    )
end

function stage_tropicalize(target, options)
    output = target.name == "base" ? base_file("01_raw_fan.mrdi") : flag_file("01_raw_fan.mrdi")
    ideal_path = target.name == "base" ? base_file("00_ideal.mrdi") : flag_file("00_ideal.mrdi")
    builder = target.name == "base" ? build_g2_grassmannian_ideal : build_full_flag_g2_ideal
    problem = g2_problem(target.name)
    engine = options.engine == :auto ? problem.default_engine : options.engine
    # The generator prevarieties are the `dressian` and `flag_dressian` targets.
    outputs = target.name == "base" ? [output, base_file("01_variety_complex.mrdi")] : [output]
    save_ideal_if_needed(ideal_path, builder, options)
    unless_cached = !skip_current(outputs, options)
    if unless_cached
        run_tropicalization(
            problem,
            TropicalizationConfig(
                engine=engine,
                output=problem.default_output,
                resume=options.resume,
                force=options.force,
                workers=options.workers,
                threads=options.threads,
                bits=options.bits,
                gfan_chunksize=options.gfan_chunksize,
                fresh_processes=options.fresh_processes,
            ),
        )
    end
    fan = load(output)
    counts = (length(rays_modulo_lineality(fan).rays_modulo_lineality), n_maximal_cones(fan))
    target.name == "base" && counts != (51, 251) &&
        error("G2 raw base count regression: expected (51, 251), got $counts")
    metadata = Dict{String, Any}(
        "effective_engine" => string(engine),
        "rays" => counts[1],
        "maximal_cells" => counts[2],
        "reused_published_artifacts" => !unless_cached,
    )
    if target.name == "flag"
        (ambient_dim(fan), dim(fan), lineality_dim(fan)) == (28, 8, 4) ||
            error("flag tropical variety must have dimensions (ambient, affine, lineality) = (28,8,4)")
        metadata["configured_symmetry_group_order"] = length(G2_FLAG_PERMUTATIONS)
        metadata["configured_certification_strategy"] = string(problem.certification_strategy)
        metadata["configured_whole_interior_pruning"] = problem.whole_interior_pruning
        metadata["cone_priority"] = "dimensions above ideal dimension eight first"
        purity = verify_flag_tropical_purity(build_full_flag_g2_ideal(); reuse_only=true)
        if engine == :parallel && purity.verified
            certificate = TOML.parsefile(purity.certificate_path)
            recorded = certificate["initial_ideal_path"]
            evidence = isabspath(recorded) ? recorded : joinpath(dirname(purity.certificate_path), recorded)
            published_evidence = flag_file("00_verified_prime_initial_ideal.mrdi")
            published_certificate = flag_file("00_ideal_purity.toml")
            cp(evidence, published_evidence; force=true)
            certificate["initial_ideal_path"] = basename(published_evidence)
            open(published_certificate, "w") do io
                TOML.print(io, certificate; sorted=true)
            end
            append!(outputs, [published_certificate, published_evidence])
            metadata["purity_defining_ideal_sha256"] = certificate["defining_ideal_sha256"]
            metadata["purity_verification_policy"] = certificate["verification_policy"]
        end
    end
    write_provenance(
        manifest_file(target.name, tropicalize);
        study="G2",
        target=target.name,
        stage=tropicalize,
        inputs=[ideal_path],
        outputs=outputs,
        configuration=census_configuration(options),
        metadata=metadata,
    )
end

function stage_coarsen(target, options)
    input = target.name == "base" ? base_file("01_raw_fan.mrdi") : flag_file("01_raw_fan.mrdi")
    output = target.name == "base" ? base_file("02_coarsest_fan.mrdi") : flag_file("02_coarsest_fan.mrdi")
    provenance = target.name == "base" ? base_file("02_coarsest_cells.csv") : flag_file("02_coarsest_cells.csv")
    report = coarsest_report_path(target_directory(target.name))
    isfile(input) || error("stage 1 artifact is missing: $input")
    if !skip_current([output, provenance], options)
        coarsen_fan_file(
            input,
            output;
            strategy=target.name == "base" ? SupportLinkCoarsening() : EqualSpanCoarsening(),
            label=target.name == "base" ? "Trop(G2/P2)" : "Trop(G2/B)",
            provenance_path=provenance,
            report_path=report,
        )
    end
    fan = load(output)
    counts = (length(rays_modulo_lineality(fan).rays_modulo_lineality), n_maximal_cones(fan))
    expected = target.name == "base" ? (51, 248) : (141, 2714)
    counts == expected || error("G2 coarsest count regression: expected $expected, got $counts")
    write_provenance(
        manifest_file(target.name, coarsen);
        study="G2", target=target.name, stage=coarsen,
        inputs=[input], outputs=[output, provenance],
        configuration=census_configuration(options),
        metadata=Dict("rays" => counts[1], "maximal_cells" => counts[2]),
    )
end

function stage_normalize(target, options)
    input = target.name == "base" ? base_file("02_coarsest_fan.mrdi") : flag_file("02_coarsest_fan.mrdi")
    output_fan = target.name == "base" ? base_file("03_normalized_coarsest_fan.mrdi") : flag_file("03_normalized_coarsest_fan.mrdi")
    output_rays = target.name == "base" ? base_file("03_normalized_rays.mrdi") : flag_file("03_normalized_rays.mrdi")
    isfile(input) || error("stage 2 artifact is missing: $input")
    if !skip_current([output_fan, output_rays], options)
        fan = load(input)
        normalized = target.name == "base" ?
            normalize_fan(fan, G2_GAUGE) :
            normalize_fan(fan, G2_GAUGE, 8:28; anchor_coordinate=7)
        save(output_fan, normalized.fan)
        save(output_rays, normalized.rays)
    end
    fan = load(output_fan)
    expected = target.name == "base" ? (51, 248) : (141, 2714)
    actual = (length(rays_modulo_lineality(fan).rays_modulo_lineality), n_maximal_cones(fan))
    actual == expected || error("normalized fan count regression: expected $expected, got $actual")
    write_provenance(
        manifest_file(target.name, normalize);
        study="G2", target=target.name, stage=normalize,
        inputs=[input], outputs=[output_fan, output_rays],
        configuration=census_configuration(options),
        metadata=Dict("rays" => actual[1], "maximal_cells" => actual[2]),
    )
end

function stage_trees(target, options)
    fan_path = target.name == "base" ? base_file("03_normalized_coarsest_fan.mrdi") : flag_file("03_normalized_coarsest_fan.mrdi")
    rays_path = target.name == "base" ? base_file("03_normalized_rays.mrdi") : flag_file("03_normalized_rays.mrdi")
    weights_path = target.name == "base" ? base_file("04_canonical_cell_weights.mrdi") : flag_file("04_canonical_cell_weights.mrdi")
    cell_trees_path = target.name == "base" ? base_file("04_normalized_cell_trees.mrdi") : flag_file("04_normalized_cell_trees.mrdi")
    census_path = target.name == "base" ? base_file("04_trees.csv") : flag_file("04_trees.csv")
    all(isfile, [fan_path, rays_path]) || error("stage 3 artifacts are missing")
    if !skip_current([weights_path, cell_trees_path, census_path], options)
        fan, rays = load(fan_path), load(rays_path)
        weights = canonical_cell_weights(fan, rays)
        cell_trees = target.name == "base" ?
            normalized_cell_trees(weights, G2_GAUGE) : flag_normalized_cell_trees(fan, weights)
        TropicalCensus.save_unless_current(weights_path, TropicalCensus.rows_matrix(weights; ring=QQ))
        TropicalCensus.save_unless_current(cell_trees_path, TropicalCensus.rows_matrix(cell_trees; ring=QQ))
        target.name == "base" ?
            TropicalCensus.write_tree_census(census_path, rays, cell_trees) :
            write_flag_tree_census(census_path, rays, cell_trees)
    end
    write_provenance(
        manifest_file(target.name, trees);
        study="G2", target=target.name, stage=trees,
        inputs=[fan_path, rays_path], outputs=[weights_path, cell_trees_path, census_path],
        configuration=census_configuration(options),
        metadata=Dict("rays" => size(load(rays_path), 1), "maximal_cells" => n_maximal_cones(load(fan_path))),
    )
end

function stage_symmetry(target, options)
    fan_path = target.name == "base" ? base_file("03_normalized_coarsest_fan.mrdi") : flag_file("03_normalized_coarsest_fan.mrdi")
    rays_path = target.name == "base" ? base_file("03_normalized_rays.mrdi") : flag_file("03_normalized_rays.mrdi")
    prefix = target.name == "base" ? base_file("05_symmetry") : flag_file("05_symmetry")
    outputs = [prefix * suffix for suffix in ("_taxon_permutations.mrdi", "_ray_actions.mrdi", "_cell_actions.mrdi", "_cell_orbits.csv")]
    if !skip_current(outputs, options)
        fan, rays = load(fan_path), load(rays_path)
        cells = TropicalCensus.fan_cells(fan)
        symmetry_result = target.name == "base" ?
            taxon_symmetry(
                rays,
                cells;
                candidates=G2_GROUP_CANDIDATES,
                lineality=rays_modulo_lineality(fan).lineality_basis,
            ) :
            flag_taxon_symmetry(fan, rays, G2_GROUP_CANDIDATES)
        length(symmetry_result.permutations) == 12 ||
            error("expected the order-12 G2 taxon group, found $(length(symmetry_result.permutations))")
        target.name == "base" && length(symmetry_result.cell_orbits) != 31 &&
            error("expected 31 coarsest-cell orbits")
        TropicalCensus.save_symmetry_artifacts(prefix, symmetry_result)
    end
    symmetry_result = TropicalCensus.load_symmetry_artifacts(prefix)
    write_provenance(
        manifest_file(target.name, symmetry);
        study="G2", target=target.name, stage=symmetry,
        inputs=[fan_path, rays_path], outputs=outputs,
        configuration=census_configuration(options),
        metadata=Dict("group_order" => length(symmetry_result.permutations), "cell_orbits" => length(symmetry_result.cell_orbits)),
    )
end

function permutation_orientation(values)
    inversions = count(i < j && values[i] > values[j] for i in eachindex(values) for j in eachindex(values))
    return isodd(inversions) ? -1 : 1
end

function g2_leaf_signs(permutation)
    for mask in 0:(2^7 - 1)
        signs = [iszero(mask & (1 << (index - 1))) ? 1 : -1 for index in 1:7]
        coefficients = Dict{NTuple{3, Int}, Int}()
        for triple in G2_FORM_SUPPORT
            image = [permutation[index] for index in triple]
            key = Tuple(sort(image))
            coefficients[key] = permutation_orientation(image) * prod(signs[index] for index in triple)
        end
        Set(keys(coefficients)) == G2_FORM_SUPPORT || continue
        length(unique(values(coefficients))) == 1 && return signs
    end
    error("taxon permutation has no signed lift preserving the G2 three-form")
end

g2_coordinate_images(permutation, variables) = TropicalCensus.pluecker_coordinate_images(
    permutation,
    variables;
    leaf_signs=g2_leaf_signs(permutation),
)

function g2_flag_coordinate_images(permutation, variables)
    length(variables) == 28 || throw(ArgumentError("expected 28 flag coordinates"))
    signs = g2_leaf_signs(permutation)
    return vcat(
        [signs[i] * variables[permutation[i]] for i in 1:7],
        TropicalCensus.pluecker_coordinate_images(
            permutation, variables[8:28]; leaf_signs=signs,
        ),
    )
end

include(joinpath(@__DIR__, "g2_chevalley.jl"))

function stage_degenerations(target, options)
    file = target.name == "base" ? base_file : flag_file
    fan_path, rays_path = file("03_normalized_coarsest_fan.mrdi"), file("03_normalized_rays.mrdi")
    ideal_path, prefix = file("00_ideal.mrdi"), file("05_symmetry")
    output = file("06_canonical_interior_degenerations.csv")
    expected_orbits = target.name == "base" ? 31 : 270
    continuing = options.resume && !options.force && isfile(output)
    completed = continuing ? Set(
        parse(Int, first(split(line, ',', limit=2)))
        for line in readlines(output)[2:end] if !isempty(line)
    ) : Set{Int}()
    all(orbit -> 1 <= orbit <= expected_orbits, completed) ||
        error("degeneration checkpoint contains an invalid orbit")
    if options.force || !isfile(output) || (options.resume && length(completed) < expected_orbits)
        config = target.name == "base" ?
            DegenerationConfig(
                expected_orbits=expected_orbits,
                coordinate_images=g2_coordinate_images,
                chevalley_test=g2_chevalley_simplex,
            ) :
            DegenerationConfig(
                expected_orbits=expected_orbits,
                # The vector coordinates (omega_1) and the Pluecker coordinates
                # (omega_2): the body is that of rho, the sum of both slices.
                body_gradings=[collect(1:7), collect(8:28)],
                compute_saturation=false,
                compute_degree=false,
                certified_tropical=true,
                coordinate_images=g2_flag_coordinate_images,
                weight_coordinate_indices=g2_flag_weight_coordinate_indices,
            )
        if options.force
            object_directory = TropicalCensus.degeneration_object_directory(output)
            isdir(object_directory) && rm(object_directory; recursive=true)
        end
        mkpath(dirname(output))
        open(output, continuing ? "a" : "w") do io
            continuing || TropicalCensus.write_degenerations_header(io)
            classify_degenerations(
                load(ideal_path), load(fan_path), load(rays_path),
                TropicalCensus.load_symmetry_artifacts(prefix);
                config=config,
                skip_orbits=completed,
                on_record=record -> TropicalCensus.write_degenerations_record(io, output, record),
            )
        end
        recorded = Set(
            parse(Int, first(split(line, ',', limit=2)))
            for line in readlines(output)[2:end] if !isempty(line)
        )
        recorded == Set(1:expected_orbits) ||
            error("degeneration census is incomplete: found $(length(recorded)) of $expected_orbits orbits")
    else
        skip_current([output], options)
    end
    outputs = degeneration_artifact_paths(output)
    write_provenance(
        manifest_file(target.name, degenerations);
        study="G2", target=target.name, stage=degenerations,
        inputs=vcat(
            [ideal_path, fan_path, rays_path],
            [prefix * suffix for suffix in ("_taxon_permutations.mrdi", "_ray_actions.mrdi", "_cell_actions.mrdi", "_cell_orbits.csv")],
        ), outputs=outputs,
        configuration=census_configuration(options),
        metadata=Dict(
            "weight_policy" => "equal sum of normalized coarsest-cell extreme rays",
            "cell_orbits" => expected_orbits,
            "newton_okounkov_body" => target.name == "base" ?
                "degree-one values, for prime initial ideals" :
                "Minkowski sum of the omega_1 (coordinates 1-7) and omega_2 (8-28) slices, " *
                "for prime initial ideals",
            "monomial_freeness" => target.name == "base" ?
                "verified by torus saturation" : "certified tropical fan membership",
        ),
    )
end

function stage_projections(target, options)
    prefix = flag_file("07_projection")
    outputs = [prefix * suffix for suffix in ("_certificate.toml", "_orbits.csv", "_cells.csv")]
    symmetry_suffixes = ("_taxon_permutations.mrdi", "_ray_actions.mrdi", "_cell_actions.mrdi", "_cell_orbits.csv")
    inputs = vcat(
        [file(name) for file in (flag_file, base_file)
         for name in ("03_normalized_coarsest_fan.mrdi", "03_normalized_rays.mrdi")],
        [file("05_symmetry" * suffix) for file in (flag_file, base_file) for suffix in symmetry_suffixes],
    )
    all(isfile, inputs) || error("projection certificates require normalized fans and symmetry artifacts for both G2 targets")
    source_fan, source_rays = load(inputs[1]), load(inputs[2])
    target_fan, target_rays = load(inputs[3]), load(inputs[4])
    source_symmetry = TropicalCensus.load_symmetry_artifacts(flag_file("05_symmetry"))
    target_symmetry = TropicalCensus.load_symmetry_artifacts(base_file("05_symmetry"))
    arguments = (source_fan, source_rays, target_fan, target_rays, source_symmetry, target_symmetry)
    keywords = (
        projection_columns=8:28,
        source_coordinate_indices=g2_flag_weight_coordinate_indices,
        target_coordinate_indices=coordinate_permutation,
    )
    if skip_current(outputs, options)
        certificate = TOML.parsefile(first(outputs))
        verify_ray_projection_certificate(certificate, arguments...; keywords...)
        orbit_io, cell_io = IOBuffer(), IOBuffer()
        TropicalCensus.write_projection_tables(orbit_io, cell_io, certificate)
        read(outputs[2], String) == String(take!(orbit_io)) &&
            read(outputs[3], String) == String(take!(cell_io)) ||
            error("projection tables disagree with the verified certificate; rerun with --force")
    else
        certificate = certify_ray_projection(arguments...; keywords...)
        write_ray_projection_certificate(prefix, certificate)
    end
    length(certificate["orbits"]) == 270 && length(certificate["cells"]) == 2714 ||
        error("G2 projection orbit/cell count regression")
    targets = unique(record["target_cell"] for record in certificate["cells"])
    length(targets) == 248 || error("flag maximal-cell projections do not cover all Grassmannian maximal cells")
    println("  certified 270 orbit representatives and all 2714 cell transports onto 248 Grassmannian maximal cells")
    write_provenance(
        manifest_file("flag", projections);
        study="G2", target="flag", stage=projections,
        inputs=inputs, outputs=outputs,
        configuration=census_configuration(options),
        metadata=Dict(
            "projection" => "forget point coordinates x1,...,x7; retain p12,...,p67 (coordinates 8:28)",
            "proof" => "positive rational ray identities modulo target lineality; equality of generator sets; equivariant cell transport",
            "source_maximal_cells" => 2714, "source_cell_orbits" => 270,
            "target_maximal_cells" => 248, "group_order" => length(source_symmetry.permutations),
        ),
    )
end

# The paper also quotes the f-vector of the stage-1 fan: for the Dressians, the generator prevariety.
const G2_SUMMARY_TARGETS = Dict(
    "base" => (label="tropical G2-Grassmannian over QQ", stage1="01_raw_fan.mrdi"),
    "flag" => (label="tropical G2 flag variety over QQ", stage1="01_raw_fan.mrdi"),
    "dressian" => (label="G2-Dressian", stage1="01_prevariety_fan.mrdi"),
    "flag_dressian" => (label="G2 flag Dressian", stage1="01_prevariety_fan.mrdi"),
)

function stage_summarize(target, options)
    directory, spec = target_directory(target.name), G2_SUMMARY_TARGETS[target.name]
    file(name) = joinpath(directory, name)
    inputs = vcat(
        [file(spec.stage1), file("03_normalized_coarsest_fan.mrdi")],
        [file("05_symmetry" * suffix) for suffix in ("_taxon_permutations.mrdi", "_ray_actions.mrdi", "_cell_actions.mrdi", "_cell_orbits.csv")],
        filter(isfile, [file("06_canonical_interior_degenerations.csv")]),
    )
    all(isfile, inputs) || error("the summary requires the artifacts of stages 1, 3 and 5")
    output = summary_path(directory)
    skip_current([output], options) || write_census_summary(directory; label=spec.label, stage1_fan=file(spec.stage1))
    write_provenance(
        manifest_file(target.name, summarize);
        study="G2", target=target.name, stage=summarize,
        inputs=inputs, outputs=[output],
        configuration=census_configuration(options),
    )
end

function target_definition(name; degeneration_stage=false)
    stages = Set([tropicalize, coarsen, normalize, trees, symmetry, summarize])
    degeneration_stage && push!(stages, degenerations)
    callbacks = Dict(
        tropicalize => stage_tropicalize,
        coarsen => stage_coarsen,
        normalize => stage_normalize,
        trees => stage_trees,
        symmetry => stage_symmetry,
        summarize => stage_summarize,
    )
    degeneration_stage && (callbacks[degenerations] = stage_degenerations)
    file = name == "base" ? base_file : flag_file
    artifacts = Dict(
        tropicalize => name == "base" ?
            [file("01_raw_fan.mrdi"), file("01_variety_complex.mrdi")] : [file("01_raw_fan.mrdi")],
        coarsen => [file("02_coarsest_fan.mrdi"), file("02_coarsest_cells.csv")],
        normalize => [file("03_normalized_coarsest_fan.mrdi"), file("03_normalized_rays.mrdi")],
        trees => [file("04_canonical_cell_weights.mrdi"), file("04_normalized_cell_trees.mrdi"), file("04_trees.csv")],
        symmetry => [
            file("05_symmetry_taxon_permutations.mrdi"),
            file("05_symmetry_ray_actions.mrdi"),
            file("05_symmetry_cell_actions.mrdi"),
            file("05_symmetry_cell_orbits.csv"),
        ],
        summarize => [summary_path(target_directory(name))],
    )
    degeneration_stage && (artifacts[degenerations] = [file("06_canonical_interior_degenerations.csv")])
    if name == "flag"
        push!(stages, projections)
        callbacks[projections] = stage_projections
        artifacts[projections] = [file("07_projection" * suffix)
                                for suffix in ("_certificate.toml", "_orbits.csv", "_cells.csv")]
    end
    return StudyTarget(name=name, stages=stages, callbacks=callbacks, artifacts=artifacts)
end

const G2_CENSUS = StudyCensus(
    name="G2",
    targets=Dict(
        "base" => target_definition("base"; degeneration_stage=true),
        "flag" => target_definition("flag"; degeneration_stage=true),
        "dressian" => dressian_target_definition("dressian"),
        "flag_dressian" => dressian_target_definition("flag_dressian"),
        # Last in --target=all; it also needs run_census_gf2.jl.
        "zenodo" => zenodo_target_definition(),
    ),
)

run_census(G2_CENSUS)
