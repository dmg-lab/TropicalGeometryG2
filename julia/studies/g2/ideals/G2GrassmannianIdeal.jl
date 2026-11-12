module G2GrassmannianIdeal

using Oscar

export build_g2_grassmannian_ideal

const PAIRS = [(i, j) for i in 1:6 for j in (i + 1):7]

"Build the Pluecker ideal of the G2/P2 linear section of Gr(2,7)."
function build_g2_grassmannian_ideal(; coefficient_field=QQ)
    variable_names = ["p$i$j" for (i, j) in PAIRS]
    ring, variables = polynomial_ring(coefficient_field, variable_names)
    p = Dict(pair => variable for (pair, variable) in zip(PAIRS, variables))
    linears = [
        p[(2, 3)] + p[(4, 7)],
        -p[(1, 3)] + p[(5, 7)],
        p[(1, 2)] + p[(6, 7)],
        p[(5, 6)] - p[(1, 7)],
        -p[(4, 6)] - p[(2, 7)],
        p[(4, 5)] - p[(3, 7)],
        p[(1, 4)] + p[(2, 5)] + p[(3, 6)],
    ]
    generators = vcat(linears, forget_grading.(gens(grassmann_pluecker_ideal(ring, 2, 7))))
    return (
        ring=ring,
        ideal=ideal(ring, generators),
        generators=generators,
    )
end

end
