using Oscar
using LinearAlgebra
using Combinatorics


function uniform_tropical_linear_space(r,n)
    ray_mat = matrix(QQ, I(n))
    lin_mat = matrix(QQ, ones(Int, 1, n))
    inc    = IncidenceMatrix(collect(Combinatorics.combinations(1:n, r)))
    fan = polyhedral_fan(inc, ray_mat, lin_mat; non_redundant=true)
    return fan
end


# Example of use
const OUTPUT = normpath(joinpath(@__DIR__, "..", "data", "fan_uniform_linear_space_16_21.mrdi"))
save(OUTPUT, uniform_tropical_linear_space(16, 21))
