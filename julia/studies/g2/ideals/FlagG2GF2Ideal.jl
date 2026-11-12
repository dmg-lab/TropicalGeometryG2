module FlagG2GF2Ideal

using Oscar
using TropicalCensus
include("FlagG2Ideal.jl")
using .FlagG2Ideal

export build_full_flag_g2_gf2_ideal,
       G2_FLAG_GF2_PERMUTATIONS,
       flag_gf2_coordinate_images,
       verify_flag_tropical_purity

const G2_FLAG_GF2_PERMUTATIONS = G2_FLAG_PERMUTATIONS
const GF2_GRASSMANNIAN_TROPICAL_BASIS =
    TropicalCensus.data_path("g2", "grassmannian_gf2", "01_tropical_basis.txt")

"Signs in the point and exterior-square actions coincide in characteristic two."
flag_gf2_coordinate_images(permutation, variables) = vcat(
    [variables[permutation[index]] for index in 1:7],
    TropicalCensus.pluecker_coordinate_images(permutation, variables[8:28]),
)

"Parse a sum of monomials with coefficient one, as stored in the tropical basis."
function _gf2_polynomial(text, variables)
    return sum(split(text, " + ")) do monomial
        prod(split(strip(monomial), '*')) do factor
            name, power = occursin('^', factor) ? split(factor, '^') : (factor, "1")
            variables[String(name)]^parse(Int, power)
        end
    end
end

"""
    build_full_flag_g2_gf2_ideal()

The 78 defining equations of the flag ideal with coefficients reduced in
GF(2). The exported generators also contain the certified tropical basis of
the G2 Grassmannian over GF(2), whose polynomials lie in the flag ideal. Their
tropical hypersurfaces make the Pluecker projection of the generator
prevariety equal to the tropical Grassmannian over GF(2).
"""
function build_full_flag_g2_gf2_ideal()
    data = build_full_flag_g2_ideal(coefficient_field=GF(2))
    variables = Dict(string(symbol) => variable
                     for (symbol, variable) in zip(symbols(data.ring), gens(data.ring)))
    support(polynomial) = Tuple(sort!(Tuple.(collect(exponents(polynomial)))))
    known = Set(support.(data.generators))
    seeds = eltype(data.generators)[]
    for line in readlines(GF2_GRASSMANNIAN_TROPICAL_BASIS)
        isempty(strip(line)) && continue
        polynomial = _gf2_polynomial(line, variables)
        support(polynomial) in known && continue
        polynomial in data.ideal ||
            error("a Grassmannian tropical-basis polynomial is outside the flag ideal: $line")
        push!(known, support(polynomial))
        push!(seeds, polynomial)
    end
    return (
        ring=data.ring,
        ideal=data.ideal,
        generators=vcat(data.generators, seeds),
        defining_generators=data.generators,
    )
end

end # module
