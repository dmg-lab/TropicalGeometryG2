#!/usr/bin/env julia

"""Check the completed GF(2) census using its saved exact artifacts.

This checks certificate coverage and the saved algebra; it does not repeat
the whole-cone Gröbner computation. Its only output is 07_verification.toml.
"""
module G2GF2CensusVerification

using Oscar
using TropicalCensus
using SHA
using TOML
using Dates

include(joinpath(@__DIR__, "..", "ideals", "G2GrassmannianGF2Ideal.jl"))
using .G2GrassmannianGF2Ideal
include(project_path("src", "parallel_worker.jl"))

gf2_file(parts...) = data_path("g2", "grassmannian_gf2", parts...)
file_digest(path) = bytes2hex(sha256(read(path)))
require(condition, message) = condition || error(message)
support(polynomial) = Tuple(sort!(Tuple.(collect(exponents(polynomial)))))

function csv_fields(line)
    characters = collect(line)
    fields = String[]
    buffer = IOBuffer()
    quoted = false
    index = 1
    while index <= length(characters)
        character = characters[index]
        if character == '"'
            if quoted && index < length(characters) && characters[index + 1] == '"'
                print(buffer, '"')
                index += 1
            else
                quoted = !quoted
            end
        elseif character == ',' && !quoted
            push!(fields, String(take!(buffer)))
        else
            print(buffer, character)
        end
        index += 1
    end
    require(!quoted, "unterminated CSV quote")
    push!(fields, String(take!(buffer)))
    return fields
end

function csv_records(path)
    lines = filter(!isempty, readlines(path))
    require(!isempty(lines), "empty CSV: $path")
    header = csv_fields(first(lines))
    return [begin
        fields = csv_fields(line)
        require(length(fields) == length(header), "CSV column count differs: $path")
        Dict(zip(header, fields))
    end for line in lines[2:end]]
end

integer_list(text) = isempty(text) ? Int[] : parse.(Int, split(text, ';'))
rational(text) = begin
    pieces = split(text, '/')
    length(pieces) == 1 ? QQ(parse(BigInt, only(pieces))) :
        QQ(parse(BigInt, pieces[1]), parse(BigInt, pieces[2]))
end

function fan_data(fan)
    data = rays_modulo_lineality(fan)
    return (rays=Matrix{QQFieldElem}(matrix(QQ, data.rays_modulo_lineality)),
            lineality=Matrix{QQFieldElem}(matrix(QQ, data.lineality_basis)),
            cells=TropicalCensus.fan_cells(fan))
end

function fan_counts(fan)
    data = fan_data(fan)
    return Dict("ambient_dimension" => ambient_dim(fan), "dimension" => dim(fan),
                "lineality_dimension" => rank(matrix(QQ, data.lineality)),
                "rays" => size(data.rays, 1), "maximal_cells" => length(data.cells))
end

"Exact cone keys after projection through one common lineality quotient."
function quotient_cone_keys(data, quotient)
    keys = [Tuple(TropicalCensus._primitive_integer_row(
        vec(transpose(collect(ray)) * quotient))) for ray in eachrow(data.rays)]
    require(all(any(!iszero, key) for key in keys), "a fan ray lies in lineality")
    require(length(unique(keys)) == length(keys), "duplicate quotient rays")
    return Set(Tuple(sort!(keys[cell])) for cell in data.cells)
end

function require_equal_fans_modulo_lineality(left, right, label)
    require(ambient_dim(left) == ambient_dim(right), "$label: ambient dimensions differ")
    a, b = fan_data(left), fan_data(right)
    arank, brank = rank(matrix(QQ, a.lineality)), rank(matrix(QQ, b.lineality))
    require(arank == brank == rank(matrix(QQ, vcat(a.lineality, b.lineality))),
            "$label: lineality spans differ")
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, a.lineality); side=:right))
    require(quotient_cone_keys(a, quotient) == quotient_cone_keys(b, quotient),
            "$label: maximal cones differ in the lineality quotient")
    return true
end

function require_gf2(ideal_value, label)
    field = coefficient_ring(ideal_value)
    require(field isa FinField && characteristic(field) == 2 && degree(field) == 1,
            "$label is not over GF(2)")
end

function min_fan(max_fan)
    rays = matrix(QQ, max_fan.rays)
    lineality = matrix(QQ, max_fan.lineality)
    return polyhedral_fan(IncidenceMatrix([cell .+ 1 for cell in max_fan.maximal_cones]),
                         -rays, lineality; non_redundant=true)
end

function verify_gf2_census()
    names = ["00_ideal.mrdi", "01_raw_fan.mrdi",
             "01_tropical_basis.txt", "01_tropical_basis_supports.gfan",
             "01_certified_fan_max.gfan", "01_cone_certificates.csv",
             "01_certification_policy.toml", "02_coarsest_fan.mrdi",
             "03_normalized_coarsest_fan.mrdi", "03_normalized_rays.mrdi",
             "04_canonical_cell_weights.mrdi", "04_normalized_cell_trees.mrdi", "04_trees.csv",
             "05_symmetry_taxon_permutations.mrdi", "05_symmetry_ray_actions.mrdi",
             "05_symmetry_cell_actions.mrdi", "05_symmetry_cell_orbits.csv",
             "06_canonical_interior_degenerations.csv"]
    missing_files = filter(name -> !isfile(gf2_file(name)), names)
    require(isempty(missing_files), "the census is incomplete: $(join(missing_files, ", "))")

    policy = TOML.parsefile(gf2_file("01_certification_policy.toml"))
    require(policy["certification_policy"] == TropicalCensus.TROPICAL_CERTIFICATION_POLICY,
            "unknown exact-certification policy")
    require(policy["certification_strategy"] == "groebner", "expected whole-cone certification")
    require(policy["coefficient_field"] == "GF(2)" && policy["valuation"] == "trivial",
            "coefficient field or valuation differs")
    for (key, path) in (
        ("basis_sha256", gf2_file("01_tropical_basis_supports.gfan")),
        ("fan_sha256", gf2_file("01_certified_fan_max.gfan")),
        ("ideal_sha256", gf2_file("00_ideal.mrdi")),
        ("source_prevariety_sha256", data_path("g2", "dressian", "01_prevariety_fan.mrdi")),
    )
        require(policy[key] == file_digest(path), "certificate SHA-256 mismatch: $key")
    end
    println("Verified coefficient policy and artifact hashes.")
    flush(stdout)

    defining_ideal = load(gf2_file("00_ideal.mrdi"))
    require_gf2(defining_ideal, "saved defining ideal")
    ring = base_ring(defining_ideal)
    source_ideal = load(data_path("g2", "grassmannian", "00_ideal.mrdi"))
    require(symbols(base_ring(source_ideal)) == symbols(ring), "source coordinate order differs")
    require(sort!(support.(gens(source_ideal))) == sort!(support.(gens(defining_ideal))),
            "defining supports differ from the source prevariety ideal")
    provider = build_g2_grassmannian_gf2_ideal()
    provider_polynomials = [_parallel_trop_parse_polynomial(parallel_trop_gfan_polynomial(f), ring)
                            for f in provider.generators]
    require(ideal(ring, provider_polynomials) == defining_ideal,
            "saved defining ideal differs from the GF(2) provider")
    basis_text = filter(!isempty, strip.(readlines(gf2_file("01_tropical_basis.txt"))))
    support_text = TropicalCensus.gfan_basis_polynomials(
        gf2_file("01_tropical_basis_supports.gfan"), string.(symbols(ring)))
    require(!isnothing(support_text) && sort(basis_text) == sort(support_text),
            "exported basis differs from its hashed support input")
    basis = [_parallel_trop_parse_polynomial(text, ring) for text in basis_text]
    require(all(!iszero, basis), "the tropical basis contains a zero polynomial")
    require(all(iszero(normal_form(f, defining_ideal)) for f in basis),
            "a tropical-basis polynomial does not belong to the GF(2) defining ideal")
    require(ideal(ring, basis) == defining_ideal, "the tropical basis generates a different ideal")
    println("Verified $(length(basis)) tropical-basis polynomials in the GF(2) ideal.")
    flush(stdout)

    # The generator prevariety is the G2-Dressian, stored once for all fields.
    prevariety = load(data_path("g2", "dressian", "01_prevariety_fan.mrdi"))
    certified_max = TropicalCensus.parse_gfan_fan(gf2_file("01_certified_fan_max.gfan"))
    certified_min = min_fan(certified_max)
    fine = TropicalCensus.tree_orbit_jobs(certified_min, G2_GF2_PERMUTATIONS)
    checks = csv_records(gf2_file("01_cone_certificates.csv"))
    require(all(record["valid"] == "true" for record in checks), "a final cone check failed")
    whole = filter(record -> record["phase"] == "whole_cone", checks)
    require(length(whole) == length(fine.jobs), "whole-cone certificate count differs")
    jobs = Dict(job.representative => job for job in fine.jobs)
    require(length(unique(record["representative"] for record in whole)) == length(whole),
            "duplicate whole-cone certificate representatives")
    for record in whole
        representative = parse(Int, record["representative"])
        require(haskey(jobs, representative), "unknown whole-cone representative")
        job = jobs[representative]
        require(record["reason"] == "contained", "a cone has no containment certificate")
        require(parse(Int, record["orbit_size"]) == job.orbit_size &&
                parse(Int, record["dimension"]) == job.dimension,
                "whole-cone certificate geometry differs")
        require(parse.(Int, split(record["tested_weights"])) == job.weight,
                "whole-cone certificate weight differs")
    end
    require(sum(parse(Int, record["orbit_size"]) for record in whole) ==
            length(certified_max.maximal_cones), "whole-cone certificates omit fine cones")
    println("Verified complete coverage of $(length(certified_max.maximal_cones)) fine cones.")
    flush(stdout)

    raw = load(gf2_file("01_raw_fan.mrdi"))
    coarse = load(gf2_file("02_coarsest_fan.mrdi"))
    normalized = load(gf2_file("03_normalized_coarsest_fan.mrdi"))
    require(dim(raw) == dim(coarse) == dim(normalized) == dim(defining_ideal),
            "the tropical support dimension differs from the defining ideal")
    require(ambient_dim(raw) == ambient_dim(coarse) == ambient_dim(normalized) == 21,
            "the fan is not in the 21 Pluecker coordinates")
    require(fan_counts(coarse) == fan_counts(normalized), "normalization changed fan counts")
    require_equal_fans_modulo_lineality(coarse, normalized, "normalized support fan")
    rays = load(gf2_file("03_normalized_rays.mrdi"))
    normalized_data = fan_data(normalized)
    saved_ray_rows = Matrix{QQFieldElem}(rays)
    require(size(saved_ray_rows) == size(normalized_data.rays) &&
            all(TropicalCensus.quotient_positive_equivalent(
                vec(saved_ray_rows[index, :]), vec(normalized_data.rays[index, :]),
                normalized_data.lineality) for index in axes(saved_ray_rows, 1)),
            "saved normalized rays differ from the normalized fan modulo lineality")
    action = TropicalCensus.load_symmetry_artifacts(gf2_file("05_symmetry"))
    require(length(action.permutations) == 12, "expected the order-12 taxon group")
    TropicalCensus.validate_tropical_permutation_group(action.permutations)
    recomputed = taxon_symmetry(rays, normalized_data.cells;
        candidates=G2_GF2_PERMUTATIONS, lineality=normalized_data.lineality)
    for field in (:permutations, :ray_actions, :cell_actions, :ray_orbits, :cell_orbits)
        require(getproperty(action, field) == getproperty(recomputed, field),
                "saved symmetry differs from exact geometry: $field")
    end
    require(TropicalCensus.verify_ideal_preservation(
        defining_ideal, action.permutations, g2_gf2_coordinate_images),
        "the saved symmetry does not preserve the GF(2) defining ideal")
    weights = canonical_cell_weights(normalized, rays)
    saved_weights = load(gf2_file("04_canonical_cell_weights.mrdi"))
    require(saved_weights == TropicalCensus.rows_matrix(weights; ring=QQ),
            "saved canonical weights differ from exact positive ray sums")
    require(TropicalCensus.verify_canonical_weight_action(
        weights, action, normalized_data.lineality), "canonical weights are not equivariant")
    cell_trees = normalized_cell_trees(weights, TreeGauge(7, [1 -1 0 -1 1 0 0; -1 0 1 1 0 -1 0]))
    require(load(gf2_file("04_normalized_cell_trees.mrdi")) == TropicalCensus.rows_matrix(cell_trees; ring=QQ),
            "saved normalized cell trees differ from the normalized ray sums")
    orbit_records = csv_records(gf2_file("05_symmetry_cell_orbits.csv"))
    require(length(orbit_records) == length(action.cell_orbits), "orbit CSV count differs")
    for (index, (record, orbit)) in enumerate(zip(orbit_records, action.cell_orbits))
        require(integer_list(record["cells"]) == orbit &&
                parse(Int, record["orbit"]) == index &&
                parse(Int, record["size"]) == length(orbit) &&
                parse(Int, record["representative"]) == first(orbit), "orbit CSV differs")
    end
    trees = csv_records(gf2_file("04_trees.csv"))
    require(all(record["kind"] in ("ray", "cell") for record in trees), "unknown tree record kind")
    require(sort([parse(Int, r["index"]) for r in trees if r["kind"] == "ray"]) ==
            collect(1:size(rays, 1)), "ray tree census is incomplete")
    require(sort([parse(Int, r["index"]) for r in trees if r["kind"] == "cell"]) ==
            collect(eachindex(weights)), "cell tree census is incomplete")
    println("Verified normalized fan, taxon symmetry, canonical cell weights, and normalized cell trees.")
    flush(stdout)

    degeneration_path = gf2_file("06_canonical_interior_degenerations.csv")
    records = csv_records(degeneration_path)
    require(sort([parse(Int, record["orbit"]) for record in records]) ==
            collect(eachindex(action.cell_orbits)), "degeneration orbit coverage is incomplete")
    object_hashes = Dict{String, String}()
    for (record_index, record) in enumerate(records)
        orbit_index = parse(Int, record["orbit"])
        orbit = action.cell_orbits[orbit_index]
        representative = first(orbit)
        require(parse(Int, record["representative_cell"]) == representative &&
                parse(Int, record["orbit_size"]) == length(orbit) &&
                integer_list(record["cells"]) == orbit, "degeneration orbit record differs")
        require(integer_list(record["rays"]) == normalized_data.cells[representative],
                "degeneration rays differ")
        require(rational.(split(record["canonical_weight"], ';')) == weights[representative],
                "degeneration canonical weight differs")
        require(record["monomial_free"] == "true", "a recorded degeneration is not tropical")
        initial_path = normpath(joinpath(dirname(degeneration_path), record["initial_ideal"]))
        saturation_path = normpath(joinpath(dirname(degeneration_path), record["torus_saturation"]))
        require(!isempty(record["torus_saturation"]) && isfile(initial_path) && isfile(saturation_path),
                "a saved degeneration ideal or saturation is missing")
        initial_ideal, saturated = load(initial_path), load(saturation_path)
        require_gf2(initial_ideal, "initial ideal for orbit $orbit_index")
        require_gf2(saturated, "torus saturation for orbit $orbit_index")
        require(!isone(initial_ideal) && !isone(saturated),
                "a saved degeneration has empty torus")
        initial_ring, saturation_ring = base_ring(initial_ideal), base_ring(saturated)
        require(symbols(initial_ring) == symbols(ring) == symbols(saturation_ring),
                "a saved degeneration has different coordinate names")
        inclusion_generators = [_parallel_trop_parse_polynomial(parallel_trop_gfan_polynomial(f),
                                  saturation_ring) for f in gens(initial_ideal)]
        require(all(iszero(normal_form(f, saturated)) for f in inclusion_generators),
                "saved saturation does not contain its initial ideal")
        saved_saturation = ideal(initial_ring,
            [_parallel_trop_parse_polynomial(parallel_trop_gfan_polynomial(f), initial_ring)
             for f in gens(saturated)])
        actual_saturation = saturation(initial_ideal, ideal(initial_ring, prod(gens(initial_ring))))
        require(actual_saturation == saved_saturation,
                "saved torus saturation differs from exact product saturation")
        require((saved_saturation == initial_ideal) == (record["torus_saturated"] == "true"),
                "saved torus-saturation flag differs")
        require(dim(initial_ideal) == parse(Int, record["dimension"]) == dim(defining_ideal),
                "saved initial-ideal dimension differs")
        object_hashes[relpath(initial_path, gf2_file())] = file_digest(initial_path)
        object_hashes[relpath(saturation_path, gf2_file())] = file_digest(saturation_path)
        if record_index % 8 == 0 || record_index == length(records)
            println("  Verified saved degeneration algebra $record_index/$(length(records)).")
            flush(stdout)
        end
    end
    println("Verified $(length(records)) saved GF(2) initial ideals and proper torus saturations.")
    flush(stdout)

    result = Dict{String, Any}(
        "schema_version" => 1, "verified" => true,
        "coefficient_field" => "GF(2)", "valuation" => "trivial", "convention" => "min",
        "verified_at_utc" => string(now(UTC)),
        "method" => "exact saved-artifact consistency, ideal membership, quotient fan equality, and complete whole-cone certificate coverage",
        "whole_cone_groebner_checks_recomputed" => false,
        "weighted_initial_ideals_recomputed" => false,
        "defining_generators" => length(gens(defining_ideal)),
        "ideal_dimension" => dim(defining_ideal),
        "tropical_basis_polynomials" => length(basis),
        "whole_cone_orbits" => length(whole),
        "certified_fine_maximal_cones" => length(certified_max.maximal_cones),
        "symmetry_group_order" => length(action.permutations),
        "canonical_degeneration_orbits" => length(records),
        "checks" => Dict("policy_hashes" => true, "generator_supports" => true,
            "basis_ideal_membership" => true, "basis_generates_defining_ideal" => true,
            "whole_cone_certificate_coverage" => true,
            "normalized_fan_equality" => true, "symmetry_actions" => true,
            "canonical_weight_actions" => true, "degeneration_orbit_coverage" => true,
            "saved_initial_ideals_gf2" => true, "saved_saturations_proper" => true,
            "saved_saturations_exact" => true,
            "initial_ideal_inclusions_in_saturations" => true),
        "fans" => Dict("prevariety" => fan_counts(prevariety), "raw" => fan_counts(raw),
                       "coarsest" => fan_counts(coarse), "normalized" => fan_counts(normalized)),
        "artifact_sha256" => Dict(name => file_digest(gf2_file(name)) for name in names),
        "degeneration_object_sha256" => object_hashes,
    )
    path = gf2_file("07_verification.toml")
    open(path, "w") do io
        TOML.print(io, result; sorted=true)
    end
    println("Verified GF(2) census: $path")
    return result
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    G2GF2CensusVerification.verify_gf2_census()
end
