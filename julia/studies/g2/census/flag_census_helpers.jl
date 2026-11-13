# Flag-target helpers shared by the QQ census (run_census.jl) and the GF(2)
# census (run_census_gf2.jl). They read only the fan, never the coefficients.

"Weight coordinates moved by a taxon permutation: point block, then Pluecker block."
g2_flag_weight_coordinate_indices(permutation) =
    vcat(invperm(permutation), 7 .+ coordinate_permutation(permutation))

"The canonical cell weights of a flag fan, normalized as its rays are in stage 3."
flag_normalized_cell_trees(fan, weights) = normalized_cell_trees(
    weights, G2_GAUGE, Matrix{QQFieldElem}(matrix(QQ, rays_modulo_lineality(fan).lineality_basis)), 8:28;
    anchor_coordinate=7)

"Rays and normalized cell trees; vertical rays have no tree."
function write_flag_tree_census(path, rays, weights)
    open(path, "w") do io
        println(io, "kind,index,vertical,pendant_lengths,internal_splits")
        for (kind, vectors) in (("ray", eachrow(rays)), ("cell", weights))
            for (index, vector) in enumerate(vectors)
                pluecker = collect(vector)[8:28]
                vertical = all(iszero, pluecker)
                tree = vertical ? nothing : tree_from_min_pluecker(pluecker)
                pendant = vertical ? "" : join(TropicalGeometryG2.rational_key.(pendant_lengths(tree)), ';')
                splits = vertical ? "" : join(
                    (join(split.side, "") * "|" * join(split.complement, "") * ":" * TropicalGeometryG2.rational_key(split.length)
                     for split in internal_splits(tree)), ';')
                println(io, "$kind,$index,$vertical,\"$pendant\",\"$splits\"")
            end
        end
    end
end

"""
The subgroup of `candidates` acting on the flag fan, with its ray and cell
actions. A taxon permutation acts on the point block and on the Pluecker block.
"""
function flag_taxon_symmetry(fan, normalized_rays, candidates)
    cells = TropicalGeometryG2.fan_cells(fan)
    lineality = Matrix{QQFieldElem}(matrix(QQ, rays_modulo_lineality(fan).lineality_basis))
    quotient = Matrix{QQFieldElem}(kernel(matrix(QQ, lineality); side=:right))
    signatures = [quotient_ray_signature(row, quotient) for row in eachrow(normalized_rays)]
    signature_to_ray = Dict(signature => index for (index, signature) in enumerate(signatures))
    cell_index = Dict(Tuple(cell) => index for (index, cell) in enumerate(cells))
    permutations, ray_actions, cell_actions = Vector{Vector{Int}}(), Vector{Vector{Int}}(), Vector{Vector{Int}}()
    for permutation in candidates
        indices = g2_flag_weight_coordinate_indices(permutation)
        preserves_lineality(lineality, indices) || continue
        ray_action = Int[]
        valid = true
        for ray in eachrow(normalized_rays)
            signature = quotient_ray_signature(collect(ray)[indices], quotient)
            if !haskey(signature_to_ray, signature)
                valid = false
                break
            end
            push!(ray_action, signature_to_ray[signature])
        end
        valid || continue
        action = TropicalGeometryG2.induced_cell_action(ray_action, cells, cell_index)
        action === nothing && continue
        push!(permutations, permutation)
        push!(ray_actions, ray_action)
        push!(cell_actions, action)
    end
    return (
        permutations=permutations,
        ray_actions=ray_actions,
        cell_actions=cell_actions,
        ray_orbits=action_orbits(ray_actions, size(normalized_rays, 1)),
        cell_orbits=action_orbits(cell_actions, length(cells)),
        tested=length(candidates),
    )
end
