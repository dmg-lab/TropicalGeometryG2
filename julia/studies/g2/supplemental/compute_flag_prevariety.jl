#!/usr/bin/env julia

"""
The G2 flag Dressian: the generator prevariety of the flag ideal, computed
with the census's Gfan and lifting routines. This is the `tropicalize` stage
of the `flag_dressian` target of run_census.jl; the certified tropical variety
is a separate artifact.
"""
module FlagDressianPrevariety

using Oscar
using TropicalCensus

include(joinpath(@__DIR__, "..", "ideals", "FlagG2Ideal.jl"))

function verify_flag_generator_reduction(original, data)
    names = FlagG2Ideal.original_variable_names()
    string.(symbols(base_ring(original))) == names || error("unexpected flag coordinates")
    ngens(original) == 78 || error("expected the 78 defining flag generators")
    signs = [name in ("p12", "p23", "p46") ? -1 : 1 for name in names]
    variables = gens(data.ring)
    substitution = hom(base_ring(original), data.ring,
        [signs[i] * variables[data.lift_indices[i]] for i in eachindex(names)])
    reduced = typeof(first(data.generators))[]
    eliminated = 0
    for polynomial in gens(original)
        image = substitution(polynomial)
        if iszero(image)
            length(polynomial) == 2 || error("a nonbinomial generator vanished")
            eliminated += 1
        else
            # No terms combine or cancel: their weights, including ties, are
            # preserved on the six binomial tropical hyperplanes.
            length(image) == length(polynomial) || error("reduction changed a tropical hypersurface")
            push!(reduced, image)
        end
    end
    eliminated == 6 || error("expected precisely six binomial equations")
    normalized(polynomials) = Set(FlagG2Ideal.gfan_polynomial(
        FlagG2Ideal._normalize_polynomial(f)) for f in polynomials)
    normalized(reduced) == normalized(data.generators) || error("provider generators disagree with the saved ideal")
    println("Verified: all 78 saved flag generators reduce without changing the prevariety.")
end

function verify_flag_prevariety_cones(original, data, native)
    # Keep distinct original terms even when binomial elimination gives equal
    # exponents. They still provide the required tropical tie.
    for (generator_index, polynomial) in enumerate(gens(original))
        supports = Vector{BigInt}[]
        for exponent in exponents(polynomial)
            reduced = zeros(BigInt, native.ambient_dimension)
            for i in eachindex(exponent)
                reduced[data.lift_indices[i]] += exponent[i]
            end
            push!(supports, reduced)
        end
        for lineality in native.lineality
            values = [sum(exponent .* lineality) for exponent in supports]
            all(==(first(values)), values) || error("generator $generator_index changes along lineality")
        end
        masks = UInt8[]
        for ray in native.rays
            values = [sum(exponent .* ray) for exponent in supports]
            maximum_value = maximum(values)
            mask = UInt8(0)
            for i in eachindex(values)
                values[i] == maximum_value && (mask |= UInt8(1) << (i - 1))
            end
            push!(masks, mask)
        end
        for (cone_index, cone) in enumerate(native.maximal_cones)
            mask = UInt8((1 << length(supports)) - 1)
            for ray_index in cone
                mask &= masks[ray_index + 1]
            end
            count_ones(mask) >= 2 || error("cone $cone_index fails flag generator $generator_index")
        end
    end
    println("Verified every maximal cone against all 78 original tropical equations.")
end

function compute_flag_prevariety(options=CensusOptions())
    data = FlagG2Ideal.build_flag_g2_ideal()
    ideal_path = data_path("g2", "flag", "00_ideal.mrdi")
    original = load(ideal_path)
    verify_flag_generator_reduction(original, data)

    directory = work_path("g2", "flag", "prevariety")
    mkpath(directory)
    input_path = joinpath(directory, "generators.gfan")
    checkpoint = joinpath(directory, "prevariety.gfan")
    output_path = data_path("g2", "flag_dressian", "01_prevariety_fan.mrdi")
    input_text = FlagG2Ideal.gfan_input(string.(symbols(data.ring)),
        FlagG2Ideal.gfan_polynomial.(data.generators))
    reuse = !options.force && isfile(input_path) && read(input_path, String) == input_text &&
        isfile(checkpoint) && TropicalCensus.prevariety_checkpoint_is_parseable(checkpoint)
    TropicalCensus.atomic_write(input_path, input_text)
    config = TropicalizationConfig(threads=options.threads, bits=options.bits, output=directory)
    isempty(config.gfan) && error("gfan_tropicalprevariety was not found")
    println("G2/B generator prevariety: $(length(data.generators)) equations in 22 reduced coordinates.")
    flush(stdout)
    # Gfan's threads order the rays and cones differently on every run, so the
    # fan is saved in canonical order, also a checkpoint saved without it.
    native = reuse ? TropicalCensus.canonicalize_gfan_fan(checkpoint) :
        TropicalCensus.run_prevariety(config, input_path, checkpoint, joinpath(directory, "gfan.log");
                                      canonical=true)
    verify_flag_prevariety_cones(original, data, native)
    TropicalCensus.save_lifted_fan(checkpoint, output_path, data.lift_indices)
    fan = load(output_path)
    ambient_dim(fan) == 28 || error("the saved flag fan must have 28 coordinates")
    n_maximal_cones(fan) == length(native.maximal_cones) || error("serialization changed the cones")
    lineality_dim(fan) == native.lineality_dimension || error("serialization changed lineality")
    f_vector_header = match(r"(?:^|\n)F_VECTOR\n([^\n]*)", native.text)
    isnothing(f_vector_header) && error("Gfan omitted the f-vector")
    quotient_f_vector = parse.(Int, split(only(f_vector_header.captures)))

    metadata = Dict{String, Any}(
        "component" => "initial_generator_prevariety",
        "convention" => "min",
        "original_generators" => ngens(original),
        "reduced_generators" => length(data.generators),
        "generator_reduction_verified" => true,
        "all_cones_verified_against_original_generators" => true,
        "ambient_dimension" => ambient_dim(fan),
        "dimension" => dim(fan),
        "lineality_dimension" => lineality_dim(fan),
        "f_vector_modulo_lineality" => quotient_f_vector,
        "rays" => length(native.rays),
        "maximal_cells" => length(native.maximal_cones),
        "pure" => is_pure(fan),
    )
    provenance = stage_provenance_path(data_path("g2", "flag_dressian"), tropicalize)
    write_provenance(provenance; study="G2", target="flag_dressian", stage=tropicalize,
        inputs=[ideal_path, joinpath(@__DIR__, "..", "ideals", "FlagG2Ideal.jl"), @__FILE__],
        outputs=[output_path, input_path, checkpoint],
        configuration=census_configuration(options), metadata=metadata)
    println("Saved $output_path")
    println("$(metadata["rays"]) rays, $(metadata["maximal_cells"]) maximal cones; dimension $(metadata["dimension"]), lineality $(metadata["lineality_dimension"]).")
    println("F-vector modulo lineality (including its zero cone): $quotient_f_vector")
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    FlagDressianPrevariety.compute_flag_prevariety(FlagDressianPrevariety.parse_census_options(
        vcat(["--target=flag_dressian", "--stage=tropicalize"], ARGS)))
end
