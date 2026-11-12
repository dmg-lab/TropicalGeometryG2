module G2GrassmannianGF2Ideal

using Oscar
using TropicalCensus
include("G2GrassmannianIdeal.jl")
using .G2GrassmannianIdeal

export build_g2_grassmannian_gf2_ideal, G2_GF2_PERMUTATIONS, g2_gf2_coordinate_images

const FORM_SUPPORT = Set([(1, 2, 3), (4, 5, 6), (1, 4, 7), (2, 5, 7), (3, 6, 7)])
const G2_GF2_PERMUTATIONS = [
    permutation for permutation in
    [vcat(p, [7]) for p in all_permutations(collect(1:6))]
    if Set(Tuple(sort([permutation[i] for i in triple])) for triple in FORM_SUPPORT) == FORM_SUPPORT
]

"The 42 defining equations, with coefficients reduced in GF(2)."
build_g2_grassmannian_gf2_ideal() = build_g2_grassmannian_ideal(coefficient_field=GF(2))

"Signs in the exterior-square action coincide in characteristic two."
g2_gf2_coordinate_images(permutation, variables) =
    TropicalCensus.pluecker_coordinate_images(permutation, variables)

end
