# The Bergman fan B_4 of the uniform matroid of rank 17 on 21 elements (Method 1 of Section 10), min convention, in its
# coarse structure: lineality spanned by (1, ..., 1), one maximal cone per 16 of the 21 unit vectors. It is the input
# data/fan_uniform_linear_space_16_21.mrdi. Without --force an existing file is only compared, since a new save
# has new Oscar identifiers (and the stored file is the one of the paper and the Zenodo tarball).
#   julia --project=. studies/g2/curves/uniform_bergman_fan.jl [OUTPUT] [--force]
using Oscar
using LinearAlgebra
using Combinatorics

"The fan on the unit vectors of R^n with lineality (1, ..., 1) and one maximal cone per `r` of the unit vectors."
function uniform_bergman_fan(r, n)
    ray_mat = matrix(QQ, I(n))
    lin_mat = matrix(QQ, ones(Int, 1, n))
    inc    = IncidenceMatrix(collect(Combinatorics.combinations(1:n, r)))
    fan = polyhedral_fan(inc, ray_mat, lin_mat; non_redundant=true)
    return fan
end

"Whether two fans have the same rays, lineality and maximal cones, in the same order."
same_fan(F, G) = rays(F) == rays(G) && lineality_space(F) == lineality_space(G) &&
                 maximal_cones(IncidenceMatrix, F) == maximal_cones(IncidenceMatrix, G)

if abspath(PROGRAM_FILE) == @__FILE__
    positional = filter(a -> !startswith(a, "--"), ARGS)
    output = isempty(positional) ? normpath(joinpath(@__DIR__, "..", "..", "..", "data", "fan_uniform_linear_space_16_21.mrdi")) :
             abspath(only(positional))
    fan = uniform_bergman_fan(16, 21)
    if isfile(output) && !("--force" in ARGS)
        println(output, same_fan(fan, load(output)) ? ": the same fan, left as it is" : ": a DIFFERENT fan (pass --force to overwrite)")
    else
        mkpath(dirname(output))
        save(output, fan)
        println("wrote ", output)
    end
end
