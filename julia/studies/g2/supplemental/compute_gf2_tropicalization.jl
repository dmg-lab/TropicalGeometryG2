module G2GF2Tropicalization

using Oscar
using TropicalCensus
using SHA
using TOML

include(joinpath(@__DIR__, "..", "ideals", "G2GrassmannianGF2Ideal.jl"))
using .G2GrassmannianGF2Ideal

gf2_file(parts...) = data_path("g2", "grassmannian_gf2", parts...)
gf2_work(parts...) = work_path("g2", "grassmannian_gf2", parts...)

function gf2_problem()
    return TropicalProblem(
        name="G2/P2 over GF(2)", slug="g2_p2_gf2",
        provider_file=study_path(:g2, "ideals", "G2GrassmannianGF2Ideal.jl"),
        provider_module=:G2GrassmannianGF2Ideal,
        builder=:build_g2_grassmannian_gf2_ideal,
        output_prefix="TGrass_G2_GF2",
        default_output=gf2_work("tropicalize"),
        published_fan=gf2_file("01_raw_fan.mrdi"),
        default_engine=:parallel,
        symmetry_candidates=:G2_GF2_PERMUTATIONS,
        symmetry_coordinate_images=:g2_gf2_coordinate_images,
    )
end

"Reuse the G2-Dressian, the saved generator prevariety, only after checking every support."
function seed_prevariety(data, problem)
    source = data_path("g2", "dressian", "01_prevariety_fan.mrdi")
    source_ideal = load(data_path("g2", "grassmannian", "00_ideal.mrdi"))
    support(f) = Tuple(sort!(Tuple.(collect(exponents(f)))))
    symbols(base_ring(source_ideal)) == symbols(data.ring) || error("source coordinate order differs")
    sort!(support.(gens(source_ideal))) == sort!(support.(data.generators)) ||
        error("reduction to GF(2) changed the defining polynomial supports")
    round = joinpath(problem.default_output, "round_001")
    mkpath(round)
    fan_path = joinpath(round, "prevariety.gfan")
    # A populated checkpoint is validated by the parallel engine itself.
    if !isfile(fan_path)
        fan = load(source)
        rml = rays_modulo_lineality(fan)
        rays = matrix(QQ, rml.rays_modulo_lineality)
        lineality = matrix(QQ, rml.lineality_basis)
        max_fan = polyhedral_fan(maximal_cones(IncidenceMatrix, fan), -rays, -lineality;
                                 non_redundant=true)
        orth = transpose(kernel(lineality; side=:right))
        template = (orth_lineality=TropicalCensus._primitive_integer_row.(eachrow(orth)),)
        TropicalCensus.atomic_write(fan_path, TropicalCensus._gfan_fan_text(max_fan, template))
        TropicalCensus.write_basis_input(joinpath(round, "basis.gfan"),
                                         string.(symbols(data.ring)), string.(data.generators))
        println("Seeded the GF(2) computation from the saved $(n_maximal_cones(fan))-cone generator prevariety.")
    end
    return source
end

function compute_gf2_tropicalization(options=CensusOptions())
    options.engine in (:auto, :parallel) ||
        throw(ArgumentError("GF(2) uses the exact parallel prevariety/Groebner engine"))
    problem = gf2_problem()
    mkpath(gf2_file())
    data = build_g2_grassmannian_gf2_ideal()
    save(gf2_file("00_ideal.mrdi"), data.ideal)
    source = seed_prevariety(data, problem)
    run_tropicalization(problem, TropicalizationConfig(
        engine=:parallel, output=problem.default_output,
        resume=options.resume, force=options.force,
        workers=options.workers, threads=options.threads, bits=options.bits,
        witness_limit=24,
    ))
    isfile(problem.published_fan) || error("GF(2) tropicalization did not finish")
    # Retain the tropical basis and exact cone-check evidence with the census.
    basis_path = joinpath(problem.default_output, "basis_polynomials.txt")
    cp(basis_path, gf2_file("01_tropical_basis.txt"); force=true)
    cp(joinpath(problem.default_output, problem.slug * "_reduced.gfan"),
       gf2_file("01_tropical_basis_supports.gfan"); force=true)
    marker = joinpath(problem.default_output, "certified_fan.txt")
    final = normpath(joinpath(problem.default_output, strip(read(marker, String))))
    cp(final, gf2_file("01_certified_fan_max.gfan"); force=true)
    checks = joinpath(dirname(final), "orbit_checks.csv")
    cp(checks, gf2_file("01_cone_certificates.csv"); force=true)
    policy = Dict{String, Any}(lowercase(key)=>value for (key, value) in
        split.(filter(!isempty, readlines(joinpath(problem.default_output, "certified_policy.txt"))), ' '; limit=2))
    policy["gfan_bits"] = parse(Int, policy["gfan_bits"])
    policy["coefficient_field"] = "GF(2)"
    policy["valuation"] = "trivial"
    policy["gfan_input_role"] = "integer lifts for supports only; all algebra over GF(2)"
    policy["ideal_sha256"] = bytes2hex(sha256(read(gf2_file("00_ideal.mrdi"))))
    policy["source_prevariety_sha256"] = bytes2hex(sha256(read(source)))
    open(gf2_file("01_certification_policy.toml"), "w") do io
        TOML.print(io, policy; sorted=true)
    end
    fan = load(problem.published_fan)
    write_provenance(stage_provenance_path(gf2_file(), tropicalize);
        study="G2", target="base_gf2", stage=tropicalize,
        inputs=[source, problem.provider_file, gf2_file("00_ideal.mrdi")],
        outputs=[problem.published_fan,
                 gf2_file("01_tropical_basis.txt"), gf2_file("01_tropical_basis_supports.gfan"),
                 gf2_file("01_certified_fan_max.gfan"),
                 gf2_file("01_cone_certificates.csv"), gf2_file("01_certification_policy.toml")],
        configuration=census_configuration(options),
        metadata=Dict("coefficient_field"=>"GF(2)", "valuation"=>"trivial", "convention"=>"min",
                      "ideal_dimension"=>dim(data.ideal), "ideal_degree"=>Int(degree(data.ideal)),
                      "rays"=>length(rays_modulo_lineality(fan).rays_modulo_lineality),
                      "maximal_cells"=>n_maximal_cones(fan),
                      "lineality_dimension"=>lineality_dim(fan),
                      "certification"=>"whole_cone_groebner_and_torus_saturation",
                      "gfan_input_role"=>"integer lifts for polynomial supports only; algebra over GF(2)"))
    return fan
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    using .G2GF2Tropicalization
    G2GF2Tropicalization.compute_gf2_tropicalization(G2GF2Tropicalization.parse_census_options(ARGS))
end
