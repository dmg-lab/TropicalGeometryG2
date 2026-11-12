#!/usr/bin/env julia

"""
The G2/P2 (`--target=base`) and G2/B (`--target=flag`) censuses with
coefficients in GF(2), using the trivial valuation and the min convention.
Canonical artifacts, reports, and provenance are isolated from the
characteristic-zero census in data/g2/grassmannian_gf2 and data/g2/flag_gf2.
"""
module G2GF2Census

using Oscar
using TropicalCensus
using SHA
using TOML

include(joinpath(@__DIR__, "ideals", "G2GrassmannianGF2Ideal.jl"))
using .G2GrassmannianGF2Ideal
include(joinpath(@__DIR__, "ideals", "FlagG2GF2Ideal.jl"))
using .FlagG2GF2Ideal
include(joinpath(@__DIR__, "supplemental", "compute_gf2_tropicalization.jl"))
include(joinpath(@__DIR__, "flag_census_helpers.jl"))

const G2_ROOT_DIRECTIONS = [
     1  -1   0  -1   1   0  0
    -1   0   1   1   0  -1  0
     0   1  -1   0  -1   1  0
]
const G2_GAUGE = TreeGauge(7, G2_ROOT_DIRECTIONS)
const G2_GROUP_CANDIDATES = G2_GF2_PERMUTATIONS
const SYMMETRY_SUFFIXES = (
    "_taxon_permutations.mrdi", "_ray_actions.mrdi", "_cell_actions.mrdi", "_cell_orbits.csv",
)

gf2_file(parts...) = data_path("g2", "grassmannian_gf2", parts...)
gf2_work(parts...) = work_path("g2", "grassmannian_gf2", parts...)
flag_gf2_file(parts...) = data_path("g2", "flag_gf2", parts...)
flag_gf2_work(parts...) = work_path("g2", "flag_gf2", parts...)
is_flag(target) = target.name == "flag"
target_file(target) = is_flag(target) ? flag_gf2_file : gf2_file
manifest_file(target, stage) = stage_provenance_path(target_file(target)(), stage)
# The generator prevariety of the Grassmannian is the G2-Dressian. Its
# generators have coefficients 1 and -1, so it is stored once for all fields.
dressian_prevariety() = data_path("g2", "dressian", "01_prevariety_fan.mrdi")
symmetry_files(target) = [target_file(target)("05_symmetry" * suffix) for suffix in SYMMETRY_SUFFIXES]

function skip_current(outputs, options)
    current = all(isfile, outputs)
    current && !options.force && println("  using existing canonical artifact(s)")
    return current && !options.force
end

function fan_metadata(fan)
    data = rays_modulo_lineality(fan)
    lineality_dimension = rank(matrix(QQ, data.lineality_basis))
    full_f_vector = Int.(f_vector(fan))
    quotient_f_vector = iszero(lineality_dimension) ? vcat([1], full_f_vector) :
        full_f_vector[lineality_dimension:end]
    return Dict{String, Any}(
        "rays" => length(data.rays_modulo_lineality),
        "maximal_cells" => n_maximal_cones(fan),
        "ambient_dimension" => ambient_dim(fan),
        "dimension" => dim(fan),
        "lineality_dimension" => lineality_dimension,
        "pure" => is_pure(fan),
        "f_vector_modulo_lineality_including_zero" => quotient_f_vector,
    )
end

"""
Counts and dimensions of a saved fan. The f-vector of a flag fan takes minutes
and every stage reports each fan, so the result is kept by file hash.
"""
function saved_fan_metadata(path)
    cache = work_path("g2", "fan_metadata", bytes2hex(sha256(read(path))) * ".toml")
    if !isfile(cache)
        mkpath(dirname(cache))
        metadata = fan_metadata(load(path))
        open(io -> TOML.print(io, metadata; sorted=true), cache, "w")
    end
    return TOML.parsefile(cache)
end

function record_stage(target, stage, options; inputs, outputs, metadata=Dict{String, Any}())
    write_provenance(
        manifest_file(target, stage);
        study="G2", target=is_flag(target) ? "flag_gf2" : "grassmannian_gf2", stage=stage,
        inputs=inputs, outputs=outputs, configuration=census_configuration(options),
        metadata=merge(Dict{String, Any}("coefficient_field" => "GF(2)", "convention" => "min"), metadata),
    )
    write_summary(target)
end

"Write a compact report from the artifacts currently available."
function write_summary(target)
    gf2_file = target_file(target)
    mkpath(gf2_file())
    open(gf2_file("census.md"), "w") do io
        println(io, is_flag(target) ? "# Tropical G2 flag variety over GF(2)" :
                    "# Tropical G2-Grassmannian over GF(2)")
        println(io)
        println(io, "The polynomial ring has coefficients in GF(2), with the trivial valuation. All fan coordinates and weights use the min convention and ",
                is_flag(target) ? "the 28 coordinates x1, ..., x7, p12, ..., p67." :
                    "the 21 Pluecker coordinates p12, ..., p67.")
        println(io)
        println(io, "| Fan | Rays modulo lineality | Maximal cells | Dimension | Lineality dimension | Pure | Quotient f-vector including zero |")
        println(io, "|---|---:|---:|---:|---:|---|---|")
        for (label, path) in (
            is_flag(target) ?
                ("Candidate prevariety (defining generators and Grassmannian tropical basis)",
                 gf2_file("01_prevariety_fan.mrdi")) :
                ("Generator prevariety (the G2-Dressian, in `../dressian/`)", dressian_prevariety()),
            ("Certified tropical fan", gf2_file("01_raw_fan.mrdi")),
            ("Coarsest support fan", gf2_file("02_coarsest_fan.mrdi")),
            ("Normalized coarsest fan", gf2_file("03_normalized_coarsest_fan.mrdi")),
        )
            isfile(path) || continue
            metadata = saved_fan_metadata(path)
            values = Any[metadata[key] for key in ("rays", "maximal_cells", "dimension", "lineality_dimension", "pure")]
            push!(values, "(" * join(metadata["f_vector_modulo_lineality_including_zero"], ", ") * ")")
            println(io, "| $label | ", join(values, " | "), " |")
        end
        if all(isfile, symmetry_files(target))
            action = TropicalCensus.load_symmetry_artifacts(gf2_file("05_symmetry"))
            println(io)
            println(io, "Taxon symmetry group order: $(length(action.permutations)). Ray orbits: $(length(action.ray_orbits)). Maximal-cell orbits: $(length(action.cell_orbits)).")
        end
        path = gf2_file("06_canonical_interior_degenerations.csv")
        if isfile(path)
            lines = readlines(path)
            columns = split(first(lines), ',')
            rows = filter(!isempty, lines[2:end])
            println(io)
            println(io, "Canonical interior degenerations recorded: $(length(rows)). Each weight is the exact equal positive sum of the normalized extreme rays of one coarsest maximal-cell representative. The resulting initial ideal is computed over GF(2); it need not be constant on the entire coarsest cell.")
            counts = String[]
            for name in ("prime", "binomial", "monomial_free", "toric")
                column = findfirst(==(name), columns)
                value = count(row -> split(row, ',')[column] == "true", rows)
                push!(counts, "$name: $value")
            end
            println(io)
            println(io, "Degeneration orbit counts: ", join(counts, "; "), ".")
        end
        println(io)
        println(io, "Stage provenance, including input and output SHA-256 hashes, is stored in `provenance/`. Initial ideals",
                is_flag(target) ? ", Newton-Okounkov bodies, and their two slices" : " and torus saturations",
                " are stored in `06_canonical_interior_degenerations_objects/`.")
        if isfile(gf2_file("07_projection_certificate.toml"))
            certificate = TOML.parsefile(gf2_file("07_projection_certificate.toml"))
            images = length(unique(record["target_cell"] for record in certificate["cells"]))
            println(io)
            println(io, "Projection certificate: forgetting x1, ..., x7 maps each of the $(length(certificate["cells"])) maximal cells onto a maximal cell of the tropical G2-Grassmannian over GF(2); the images are all $images of its maximal cells. See [`07_projection_certificate.toml`](07_projection_certificate.toml).")
        end
        if isfile(gf2_file("07_verification.toml"))
            println(io)
            println(io, "Independent saved-artifact verification passed; see [`07_verification.toml`](07_verification.toml). This checks tropical-basis ideal membership, whole-cone certificate coverage, symmetry and canonical weights, and exact product saturations of all saved initial ideals.")
        end
    end
end

function flag_problem()
    return TropicalProblem(
        name="G2/B over GF(2)", slug="g2_flag_gf2",
        provider_file=study_path(:g2, "ideals", "FlagG2GF2Ideal.jl"),
        provider_module=:FlagG2GF2Ideal,
        builder=:build_full_flag_g2_gf2_ideal,
        output_prefix="TFlag_G2_GF2",
        default_output=flag_gf2_work("tropicalize"),
        published_prevariety=flag_gf2_file("01_prevariety_fan.mrdi"),
        published_fan=flag_gf2_file("01_raw_fan.mrdi"),
        symmetry_candidates=:G2_FLAG_GF2_PERMUTATIONS,
        symmetry_coordinate_images=:flag_gf2_coordinate_images,
        certification_strategy=:balanced_support,
        purity_verifier=:verify_flag_tropical_purity,
        whole_interior_pruning=true,
        default_engine=:parallel,
        # The same ray and cell numbering on every run.
        canonical_gfan_order=true,
    )
end

"""
Tropicalize the flag ideal over GF(2) by balanced support, as over QQ. Purity
in dimension eight comes from an exact prime initial ideal over GF(2). The
candidate fan is the prevariety of the 78 generators and of the tropical basis
of the Grassmannian over GF(2), which lies in the flag ideal.
"""
function stage_flag_tropicalize(target, options)
    options.engine in (:auto, :parallel) ||
        throw(ArgumentError("GF(2) uses the exact parallel prevariety/Groebner engine"))
    problem = flag_problem()
    ideal_path = flag_gf2_file("00_ideal.mrdi")
    outputs = [problem.published_prevariety, problem.published_fan]
    data = build_full_flag_g2_gf2_ideal()
    characteristic(coefficient_ring(data.ideal)) == 2 || error("the flag ideal is not over GF(2)")
    mkpath(flag_gf2_file())
    (!isfile(ideal_path) || options.force) && save(ideal_path, data.ideal)
    computed = !skip_current(outputs, options)
    if computed
        run_tropicalization(problem, TropicalizationConfig(
            engine=:parallel, output=problem.default_output,
            resume=options.resume, force=options.force,
            workers=options.workers, threads=options.threads, bits=options.bits,
            gfan_chunksize=options.gfan_chunksize,
        ))
        isfile(problem.published_fan) || error("GF(2) flag tropicalization did not finish")
    end
    fan = load(problem.published_fan)
    (ambient_dim(fan), dim(fan), lineality_dim(fan)) == (28, 8, 4) ||
        error("flag tropical variety must have dimensions (ambient, affine, lineality) = (28,8,4)")

    # Retain the purity proof, the candidate polynomials, and the exact
    # membership checks with the census.
    purity = verify_flag_tropical_purity(data; reuse_only=true)
    purity.verified || error("the GF(2) flag purity certificate is missing; rerun with --force")
    certificate = TOML.parsefile(purity.certificate_path)
    recorded = certificate["initial_ideal_path"]
    evidence = isabspath(recorded) ? recorded : joinpath(dirname(purity.certificate_path), recorded)
    published_evidence = flag_gf2_file("00_verified_prime_initial_ideal.mrdi")
    published_certificate = flag_gf2_file("00_ideal_purity.toml")
    cp(evidence, published_evidence; force=true)
    certificate["initial_ideal_path"] = basename(published_evidence)
    open(published_certificate, "w") do io
        TOML.print(io, certificate; sorted=true)
    end
    append!(outputs, [published_certificate, published_evidence])
    marker = joinpath(problem.default_output, "certified_fan.txt")
    if isfile(marker)
        final = normpath(joinpath(problem.default_output, strip(read(marker, String))))
        for (source, name) in (
            (joinpath(problem.default_output, "basis_polynomials.txt"), "01_candidate_basis.txt"),
            (joinpath(dirname(final), "support_certificate.toml"), "01_support_certificate.toml"),
            (joinpath(dirname(final), "support_membership.csv"), "01_support_membership.csv"),
        )
            isfile(source) || error("missing tropicalization evidence: $source")
            cp(source, flag_gf2_file(name); force=true)
            push!(outputs, flag_gf2_file(name))
        end
    end
    prevariety = load(problem.published_prevariety)
    record_stage(target, tropicalize, options;
        inputs=[problem.provider_file, ideal_path, FlagG2GF2Ideal.GF2_GRASSMANNIAN_TROPICAL_BASIS],
        outputs=outputs,
        metadata=merge(saved_fan_metadata(problem.published_fan), Dict{String, Any}(
            "valuation" => "trivial",
            "defining_generators" => length(data.defining_generators),
            "candidate_polynomials" => length(data.generators),
            "candidate_seed" => "tropical basis of the G2 Grassmannian over GF(2)",
            "prevariety_rays" => length(rays_modulo_lineality(prevariety).rays_modulo_lineality),
            "prevariety_maximal_cells" => n_maximal_cones(prevariety),
            "certification" => "balanced support with exact purity and interior membership",
            "purity_verification_policy" => certificate["verification_policy"],
            "reused_published_artifacts" => !computed,
        )))
end

function stage_tropicalize(target, options)
    is_flag(target) && return stage_flag_tropicalize(target, options)
    outputs = vcat([gf2_file(name) for name in (
        "00_ideal.mrdi", "01_raw_fan.mrdi",
        "01_tropical_basis.txt", "01_tropical_basis_supports.gfan",
        "01_certified_fan_max.gfan", "01_cone_certificates.csv",
        "01_certification_policy.toml",
    )], [manifest_file(target, tropicalize)])
    if !skip_current(outputs, options)
        G2GF2Tropicalization.compute_gf2_tropicalization(options)
    end
    ideal_value = load(gf2_file("00_ideal.mrdi"))
    characteristic(coefficient_ring(ideal_value)) == 2 || error("the saved ideal is not over GF(2)")
    write_summary(target)
end

function stage_coarsen(target, options)
    file = target_file(target)
    input = file("01_raw_fan.mrdi")
    output = file("02_coarsest_fan.mrdi")
    certificate = file("02_coarsest_cells.csv")
    report = file("02_coarsest_report.md")
    isfile(input) || error("the GF(2) tropical fan is missing: $input")
    if !skip_current([output, certificate, report], options)
        coarsen_fan_file(input, output;
            strategy=is_flag(target) ? EqualSpanCoarsening() : SupportLinkCoarsening(),
            label=is_flag(target) ? "Trop(G2/B) over GF(2)" : "Trop(G2/P2) over GF(2)",
            provenance_path=certificate, report_path=report)
    end
    record_stage(target, coarsen, options; inputs=[input], outputs=[output, certificate, report],
        metadata=saved_fan_metadata(output))
end

function stage_normalize(target, options)
    file = target_file(target)
    input = file("02_coarsest_fan.mrdi")
    outputs = [file("03_normalized_coarsest_fan.mrdi"), file("03_normalized_rays.mrdi")]
    isfile(input) || error("the GF(2) coarsest fan is missing: $input")
    if !skip_current(outputs, options)
        normalized = is_flag(target) ?
            normalize_fan(load(input), G2_GAUGE, 8:28; anchor_coordinate=7) :
            normalize_fan(load(input), G2_GAUGE)
        save(outputs[1], normalized.fan)
        save(outputs[2], normalized.rays)
    end
    source, normalized = saved_fan_metadata(input), saved_fan_metadata(outputs[1])
    source == normalized || error("tree gauge normalization changed the GF(2) fan counts or dimensions")
    record_stage(target, normalize, options; inputs=[input], outputs=outputs, metadata=normalized)
end

function stage_trees(target, options)
    file = target_file(target)
    inputs = [file("03_normalized_coarsest_fan.mrdi"), file("03_normalized_rays.mrdi")]
    outputs = [file("04_canonical_cell_weights.mrdi"), file("04_normalized_cell_trees.mrdi"), file("04_trees.csv")]
    all(isfile, inputs) || error("the normalized GF(2) fan and rays are missing")
    if !skip_current(outputs, options)
        fan, rays = load(inputs[1]), load(inputs[2])
        weights = canonical_cell_weights(fan, rays)
        cell_trees = is_flag(target) ?
            flag_normalized_cell_trees(fan, weights) : normalized_cell_trees(weights, G2_GAUGE)
        TropicalCensus.save_unless_current(outputs[1], TropicalCensus.rows_matrix(weights; ring=QQ))
        TropicalCensus.save_unless_current(outputs[2], TropicalCensus.rows_matrix(cell_trees; ring=QQ))
        is_flag(target) ? write_flag_tree_census(outputs[3], rays, cell_trees) :
            TropicalCensus.write_tree_census(outputs[3], rays, cell_trees)
    end
    record_stage(target, trees, options; inputs=inputs, outputs=outputs,
        metadata=saved_fan_metadata(inputs[1]))
end

function stage_symmetry(target, options)
    file = target_file(target)
    inputs = [file("03_normalized_coarsest_fan.mrdi"), file("03_normalized_rays.mrdi")]
    all(isfile, inputs) || error("the normalized GF(2) fan and rays are missing")
    outputs = symmetry_files(target)
    if !skip_current(outputs, options)
        fan, rays = load(inputs[1]), load(inputs[2])
        action = is_flag(target) ? flag_taxon_symmetry(fan, rays, G2_GROUP_CANDIDATES) :
            taxon_symmetry(rays, TropicalCensus.fan_cells(fan);
                candidates=G2_GROUP_CANDIDATES, lineality=rays_modulo_lineality(fan).lineality_basis)
        length(action.permutations) == 12 ||
            error("expected the order-12 G2 taxon group, found $(length(action.permutations))")
        TropicalCensus.save_symmetry_artifacts(file("05_symmetry"), action)
    end
    action = TropicalCensus.load_symmetry_artifacts(file("05_symmetry"))
    record_stage(target, symmetry, options; inputs=inputs, outputs=outputs,
        metadata=Dict("group_order" => length(action.permutations),
            "ray_orbits" => length(action.ray_orbits), "cell_orbits" => length(action.cell_orbits)))
end

include(joinpath(@__DIR__, "g2_chevalley.jl"))

function stage_degenerations(target, options)
    file = target_file(target)
    inputs = vcat([file(name) for name in
        ("00_ideal.mrdi", "03_normalized_coarsest_fan.mrdi", "03_normalized_rays.mrdi")],
        symmetry_files(target))
    all(isfile, inputs) || error("GF(2) degenerations require the ideal, normalized fan, and symmetry artifacts")
    action = TropicalCensus.load_symmetry_artifacts(file("05_symmetry"))
    expected_orbits = length(action.cell_orbits)
    output = file("06_canonical_interior_degenerations.csv")
    continuing = options.resume && !options.force && isfile(output)
    completed = continuing ? Set(parse(Int, first(split(line, ',', limit=2)))
        for line in readlines(output)[2:end] if !isempty(line)) : Set{Int}()
    all(orbit -> 1 <= orbit <= expected_orbits, completed) ||
        error("GF(2) degeneration checkpoint contains an invalid orbit")
    if options.force || !isfile(output) || (options.resume && length(completed) < expected_orbits)
        ideal_value = load(inputs[1])
        characteristic(coefficient_ring(ideal_value)) == 2 || error("the saved ideal is not over GF(2)")
        config = is_flag(target) ?
            DegenerationConfig(expected_orbits=expected_orbits,
                # The point coordinates (omega_1) and the Pluecker coordinates
                # (omega_2): the body is that of rho, the sum of both slices.
                body_gradings=[collect(1:7), collect(8:28)],
                compute_saturation=false, compute_degree=false, certified_tropical=true,
                coordinate_images=flag_gf2_coordinate_images,
                weight_coordinate_indices=g2_flag_weight_coordinate_indices) :
            DegenerationConfig(expected_orbits=expected_orbits,
                coordinate_images=g2_gf2_coordinate_images, chevalley_test=g2_chevalley_simplex)
        if options.force
            directory = TropicalCensus.degeneration_object_directory(output)
            isdir(directory) && rm(directory; recursive=true)
        end
        mkpath(dirname(output))
        open(output, continuing ? "a" : "w") do io
            continuing || TropicalCensus.write_degenerations_header(io)
            classify_degenerations(ideal_value, load(inputs[2]), load(inputs[3]), action;
                config=config, skip_orbits=completed,
                on_record=record -> TropicalCensus.write_degenerations_record(io, output, record))
        end
    else
        skip_current([output], options)
    end
    recorded = Set(parse(Int, first(split(line, ',', limit=2)))
        for line in readlines(output)[2:end] if !isempty(line))
    recorded == Set(1:expected_orbits) ||
        error("GF(2) degeneration census is incomplete: found $(length(recorded)) of $expected_orbits orbits; use --resume")
    record_stage(target, degenerations, options; inputs=inputs, outputs=degeneration_artifact_paths(output),
        metadata=Dict("cell_orbits" => expected_orbits,
            "weight_policy" => "equal sum of normalized coarsest-cell extreme rays",
            "newton_okounkov_body" => is_flag(target) ?
                "Minkowski sum of the omega_1 (coordinates 1-7) and omega_2 (8-28) slices, " *
                "for prime initial ideals" : "degree-one values, for prime initial ideals",
            "monomial_freeness" => is_flag(target) ?
                "certified tropical fan membership" : "verified by torus saturation"))
end

"Certify that forgetting the point maps every flag cell onto a Grassmannian cell."
function stage_projections(target, options)
    prefix = flag_gf2_file("07_projection")
    outputs = [prefix * suffix for suffix in ("_certificate.toml", "_orbits.csv", "_cells.csv")]
    inputs = vcat(
        [file(name) for file in (flag_gf2_file, gf2_file)
         for name in ("03_normalized_coarsest_fan.mrdi", "03_normalized_rays.mrdi")],
        [file("05_symmetry" * suffix) for file in (flag_gf2_file, gf2_file) for suffix in SYMMETRY_SUFFIXES],
    )
    all(isfile, inputs) ||
        error("projection certificates require normalized fans and symmetry artifacts for both GF(2) targets")
    source_fan, source_rays = load(inputs[1]), load(inputs[2])
    target_fan, target_rays = load(inputs[3]), load(inputs[4])
    source_symmetry = TropicalCensus.load_symmetry_artifacts(flag_gf2_file("05_symmetry"))
    target_symmetry = TropicalCensus.load_symmetry_artifacts(gf2_file("05_symmetry"))
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
    length(certificate["orbits"]) == length(source_symmetry.cell_orbits) &&
        length(certificate["cells"]) == n_maximal_cones(source_fan) ||
        error("GF(2) projection certificate does not cover every flag cell and orbit")
    images = unique(record["target_cell"] for record in certificate["cells"])
    length(images) == n_maximal_cones(target_fan) ||
        error("GF(2) flag maximal-cell projections do not cover all Grassmannian maximal cells")
    println("  certified $(length(certificate["orbits"])) orbit representatives and all ",
            "$(length(certificate["cells"])) cell transports onto $(length(images)) Grassmannian maximal cells")
    record_stage(target, projections, options; inputs=inputs, outputs=outputs,
        metadata=Dict(
            "projection" => "forget point coordinates x1,...,x7; retain p12,...,p67 (coordinates 8:28)",
            "proof" => "positive rational ray identities modulo target lineality; equality of generator sets; equivariant cell transport",
            "source_maximal_cells" => n_maximal_cones(source_fan),
            "source_cell_orbits" => length(source_symmetry.cell_orbits),
            "target_maximal_cells" => n_maximal_cones(target_fan),
            "group_order" => length(source_symmetry.permutations),
        ))
end

"The statements of the paper, with the maximal cells and toric degenerations compared with those over QQ."
function stage_summarize(target, options)
    file = target_file(target)
    reference = data_path("g2", is_flag(target) ? "flag" : "grassmannian")
    inputs = vcat(
        [file("01_raw_fan.mrdi"), file("03_normalized_coarsest_fan.mrdi")], symmetry_files(target),
        [file("06_canonical_interior_degenerations.csv")],
        [joinpath(reference, name) for name in ("03_normalized_coarsest_fan.mrdi", "06_canonical_interior_degenerations.csv")],
    )
    all(isfile, inputs) || error("the summary requires stages 1, 3, 5 and 6 over GF(2) and stages 3 and 6 over QQ")
    output = summary_path(file())
    skip_current([output], options) || write_census_summary(file();
        label=is_flag(target) ? "tropical G2 flag variety over GF(2)" : "tropical G2-Grassmannian over GF(2)",
        stage1_fan=file("01_raw_fan.mrdi"), reference=reference)
    record_stage(target, summarize, options; inputs=inputs, outputs=[output])
end

function target_definition(name; tropicalization_callback=stage_tropicalize)
    file = name == "flag" ? flag_gf2_file : gf2_file
    stages = Set([tropicalize, coarsen, normalize, trees, symmetry, degenerations, summarize])
    callbacks = Dict(tropicalize => tropicalization_callback, coarsen => stage_coarsen,
        normalize => stage_normalize, trees => stage_trees, symmetry => stage_symmetry,
        degenerations => stage_degenerations, summarize => stage_summarize)
    artifacts = Dict(
        tropicalize => [file(artifact) for artifact in (name == "flag" ?
            ("00_ideal.mrdi", "01_prevariety_fan.mrdi", "01_raw_fan.mrdi") : ("00_ideal.mrdi", "01_raw_fan.mrdi"))],
        coarsen => [file("02_coarsest_fan.mrdi"), file("02_coarsest_cells.csv")],
        normalize => [file("03_normalized_coarsest_fan.mrdi"), file("03_normalized_rays.mrdi")],
        trees => [file("04_canonical_cell_weights.mrdi"), file("04_normalized_cell_trees.mrdi"), file("04_trees.csv")],
        symmetry => [file("05_symmetry" * suffix) for suffix in SYMMETRY_SUFFIXES],
        degenerations => [file("06_canonical_interior_degenerations.csv")],
        summarize => [summary_path(file())],
    )
    if name == "flag"
        push!(stages, projections)
        callbacks[projections] = stage_projections
        artifacts[projections] = [file("07_projection" * suffix)
                                  for suffix in ("_certificate.toml", "_orbits.csv", "_cells.csv")]
    end
    return StudyTarget(name=name, stages=stages, callbacks=callbacks, artifacts=artifacts)
end

const G2_GF2_CENSUS = StudyCensus(name="G2 over GF(2)",
    targets=Dict("base" => target_definition("base"), "flag" => target_definition("flag")))

end # module G2GF2Census

if abspath(PROGRAM_FILE) == @__FILE__
    G2GF2Census.run_census(G2GF2Census.G2_GF2_CENSUS)
end
