using Oscar
using Combinatorics
import Oscar: pluecker_indices

function pluecker_indices(r::Int,n::Int)
    return [collect(s) for s in Combinatorics.combinations(1:n, r)]
end

function generate_uniform_tropical_linear_space(r, n, output_path)
    T = tropical_semiring(max)
    p = pluecker_indices(r,n)
    pVector = T.(zeros(Int64,length(p)))
    mkpath(dirname(output_path))
    save(output_path, polyhedral_complex(tropical_linear_space(p,pVector)))
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 3 || error("usage: generate_uniform_linear_space.jl R N OUTPUT")
    rank = parse(Int, ARGS[1])
    ground_set_size = parse(Int, ARGS[2])
    generate_uniform_tropical_linear_space(rank, ground_set_size, abspath(ARGS[3]))
end
