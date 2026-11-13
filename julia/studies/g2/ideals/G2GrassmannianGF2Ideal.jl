module G2GrassmannianGF2Ideal

using Oscar
using TropicalGeometryG2
include("G2GrassmannianIdeal.jl")
using .G2GrassmannianIdeal
include("G2Symmetry.jl")
using .G2Symmetry: G2_PERMUTATIONS

export build_g2_grassmannian_gf2_ideal, G2_GF2_PERMUTATIONS, g2_gf2_coordinate_images

# The order-12 group of G2Symmetry.jl, under the name the census problems use.
const G2_GF2_PERMUTATIONS = G2_PERMUTATIONS

"The 42 defining equations, with coefficients reduced in GF(2)."
build_g2_grassmannian_gf2_ideal() = build_g2_grassmannian_ideal(coefficient_field=GF(2))

"Signs in the exterior-square action coincide in characteristic two."
g2_gf2_coordinate_images(permutation, variables) =
    TropicalGeometryG2.pluecker_coordinate_images(permutation, variables)

end
