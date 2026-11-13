include(joinpath(@__DIR__, "atlas_common.jl"))
using TropicalGeometryG2: all_permutations, action_orbits, permutation_order, cycle_notation,
    generated_subgroup

const ORBIT_REPORT = joinpath(@__DIR__, "REPORT.md")

"""Topology-only rooted signature, used to cross-check triangle images."""
function rooted_topology_signature(tree::TreeMetric, permutation::Vector{Int}=collect(1:7))
    clusters = [sort([permutation[taxon] for taxon in edge.cluster])
                for edge in rooted_internal_edges(tree)]
    sort!(clusters; by=cluster -> (length(cluster), join(cluster, ",")))
    return join(cluster_key.(clusters), ";")
end

function unique_signature_index(trees, signature_function)
    result = Dict{String, Int}()
    for (index, tree) in enumerate(trees)
        signature = signature_function(tree)
        haskey(result, signature) &&
            error("signature collision between indices $(result[signature]) and $index: $signature")
        result[signature] = index
    end
    return result
end

function induced_ray_permutation(permutation, ray_trees, signature_to_ray)
    result = Int[]
    for tree in ray_trees
        signature = rooted_metric_signature(tree, permutation)
        haskey(signature_to_ray, signature) || return nothing
        push!(result, signature_to_ray[signature])
    end
    sort(result) == collect(eachindex(ray_trees)) || return nothing
    return result
end

function induced_triangle_permutation(ray_permutation, triangles, triangle_to_index)
    result = Int[]
    for triangle in triangles
        image = Tuple(sort(ray_permutation[triangle]))
        haskey(triangle_to_index, image) || return nothing
        push!(result, triangle_to_index[image])
    end
    sort(result) == collect(eachindex(triangles)) || return nothing
    return result
end

function write_group_csv(path, label_permutations, ray_actions, triangle_actions)
    open(path, "w") do io
        println(io, "element,label_images,cycle_notation,order,ray_permutation,triangle_permutation")
        ordering = sortperm(eachindex(label_permutations); by=index -> Tuple(label_permutations[index]))
        for (element, index) in enumerate(ordering)
            label_permutation = label_permutations[index]
            println(
                io,
                join([
                    string(element),
                    join(label_permutation, ";"),
                    cycle_notation(label_permutation),
                    string(permutation_order(label_permutation)),
                    join(ray_actions[index], ";"),
                    join(triangle_actions[index], ";"),
                ], ','),
            )
        end
    end
end

function write_ray_orbits_csv(path, orbits, ray_trees)
    open(path, "w") do io
        println(io, "orbit,size,rays,representative,rooted_metric_signature")
        for (orbit_index, orbit) in enumerate(orbits)
            representative = first(orbit)
            println(
                io,
                join([
                    string(orbit_index),
                    string(length(orbit)),
                    join(orbit, ";"),
                    string(representative),
                    rooted_metric_signature(ray_trees[representative]),
                ], ','),
            )
        end
    end
end

function fiber_type(tree::TreeMetric)
    internal_vertices = [vertex for vertex in keys(tree.adjacency) if vertex > tree.ntaxa]
    valencies = sort([length(tree.adjacency[vertex]) for vertex in internal_vertices])
    valencies == [3, 3, 3, 3, 3] && return "trivalent"
    valencies == [3, 3, 3, 4] && return "one_four_valent"
    return "valencies_" * join(valencies, "_")
end

function write_triangle_orbits_csv(path, orbits, triangles, triangle_trees)
    open(path, "w") do io
        println(io, "orbit,size,triangles,representative,representative_rays,fiber_type,rooted_topology")
        for (orbit_index, orbit) in enumerate(orbits)
            representative = first(orbit)
            types = unique(fiber_type(triangle_trees[index]) for index in orbit)
            length(types) == 1 || error("fiber type is not constant on label orbit $orbit_index")
            println(
                io,
                join([
                    string(orbit_index),
                    string(length(orbit)),
                    join(orbit, ";"),
                    string(representative),
                    join(triangles[representative], ";"),
                    only(types),
                    rooted_topology_signature(triangle_trees[representative]),
                ], ','),
            )
        end
    end
end

function refinement_description(coarse_orbits, fine_orbits)
    fine_sets = Set.(fine_orbits)
    result = Vector{Vector{Int}}()
    for coarse in coarse_orbits
        coarse_set = Set(coarse)
        pieces = [index for (index, fine_set) in enumerate(fine_sets) if fine_set <= coarse_set]
        union_of_pieces = reduce(union, (fine_sets[index] for index in pieces); init=Set{Int}())
        union_of_pieces == coarse_set || error("label orbits do not refine a combinatorial orbit")
        push!(result, pieces)
    end
    return result
end

function markdown_orbit_table(io, orbits)
    println(io, "| orbit | size | indices |")
    println(io, "|---:|---:|:---|")
    for (index, orbit) in enumerate(orbits)
        println(io, "| $index | $(length(orbit)) | `$(join(orbit, ", "))` |")
    end
end

function write_report(
    path;
    tested_permutations,
    ray_preserving_permutations,
    label_permutations,
    ray_orbits,
    triangle_orbits,
    abstract_group_order,
    abstract_ray_orbits,
    abstract_triangle_orbits,
    ray_refinement,
    triangle_refinement,
    expected_dihedral,
)
    open(path, "w") do io
        println(io, "# Label-induced orbits of the G2 tree fan\n")
        println(io, "## Result\n")
        println(io, "We rooted every exact metric tree at taxon `7`, tested all `$tested_permutations` permutations of taxa `1,...,6`, and fixed taxon `7`. A permutation was retained only when it\n")
        println(io, "1. permuted the 51 generator-ray metric signatures; and")
        println(io, "2. its induced ray permutation preserved all 251 maximal triangles.\n")
        println(io, "The ray test retained **$ray_preserving_permutations** permutations and the full fan-incidence test retained **$(length(label_permutations))**. The final group is $(expected_dihedral ? "exactly" : "not equal to") the expected hexagonal Weyl group\n")
        println(io, "```text")
        println(io, "rotation   (1 5 3 4 2 6)")
        println(io, "reflection (2 3)(5 6)")
        println(io, "taxon 7    fixed")
        println(io, "```\n")
        println(io, "Thus the label-induced group has order **$(length(label_permutations))**$(expected_dihedral ? " and is the dihedral Weyl group W(G2)." : ".")\n")
        println(io, "This is the permutation part visible after valuation. Algebraic Weyl representatives may also carry coordinate signs, but those signs disappear from the tropical tree data.\n")
        println(io, "The resulting action has **$(length(ray_orbits)) ray orbits** and **$(length(triangle_orbits)) triangle orbits**, with size distributions\n")
        println(io, "```text")
        println(io, "rays:      $(sort(length.(ray_orbits)))")
        println(io, "triangles: $(sort(length.(triangle_orbits)))")
        println(io, "```\n")

        println(io, "## Why this differs from the old atlas orbits\n")
        println(io, "The abstract simplicial automorphism group has order **$abstract_group_order**, with $(length(abstract_ray_orbits)) ray orbits and $(length(abstract_triangle_orbits)) triangle orbits. The label-induced Weyl group has index **$(abstract_group_order ÷ length(label_permutations))** in that group. It is the subgroup realized by permuting the six nonzero-weight taxa (with algebraic signs suppressed tropically). In particular, abstract triangle orbits of size 24 cannot be Weyl orbits because `|W(G2)|=12`.\n")
        println(io, "The finer label-orbit partition splits the abstract ray orbits as follows:\n")
        println(io, "| old abstract ray orbit | size | new label orbits |")
        println(io, "|---:|---:|:---|")
        for (index, pieces) in enumerate(ray_refinement)
            println(io, "| $index | $(length(abstract_ray_orbits[index])) | `$(join(pieces, ", "))` |")
        end
        println(io, "\nIt splits the abstract triangle orbits as follows:\n")
        println(io, "| old abstract orbit | size | new label orbits |")
        println(io, "|---:|---:|:---|")
        for (index, pieces) in enumerate(triangle_refinement)
            println(io, "| $index | $(length(abstract_triangle_orbits[index])) | `$(join(pieces, ", "))` |")
        end

        println(io, "\n## Ray orbits\n")
        markdown_orbit_table(io, ray_orbits)
        println(io, "\n## Triangle orbits\n")
        markdown_orbit_table(io, triangle_orbits)

        println(io, "\n## Files\n")
        println(io, "- `label_permutation_group.csv`: every label permutation and its induced actions;")
        println(io, "- `ray_label_orbits.csv`: exact rooted-metric ray orbits;")
        println(io, "- `triangle_label_orbits.csv`: triangle orbits, ray triples, fiber types, and rooted topologies.\n")
        println(io, "The paper ray catalog and the flag atlas read these files. Reproduce them from julia/ with\n")
        println(io, "```sh")
        println(io, "julia --project=. studies/g2/atlas/grassmannian/compute_label_orbits.jl")
        println(io, "```\n")
        println(io, "The maintained coarsest-fan symmetry artifacts are recomputed instead with\n")
        println(io, "```sh")
        println(io, "julia --project=. studies/g2/run_census.jl --target=base --stage=symmetry")
        println(io, "```")
    end
end

function compute_label_orbits()
    fan = load(DEFAULT_BALANCED_FAN_PATH)
    balanced_ray_matrix = load(DEFAULT_BALANCED_RAYS_PATH)
    @assert nrows(balanced_ray_matrix) == 51

    ray_trees = [tree_from_min_pluecker(row_as_vector(balanced_ray_matrix, index))
                 for index in 1:nrows(balanced_ray_matrix)]
    signature_to_ray = unique_signature_index(ray_trees, rooted_metric_signature)

    triangle_records = balanced_triangle_tree_data(fan, balanced_ray_matrix)
    triangles = [record.rays for record in triangle_records]
    triangle_trees = [record.tree for record in triangle_records]
    @assert length(triangles) == 251
    triangle_to_index = Dict(Tuple(triangle) => index for (index, triangle) in enumerate(triangles))
    @assert length(triangle_to_index) == 251

    tested = 0
    ray_preserving = 0
    label_permutations = Vector{Vector{Int}}()
    ray_actions = Vector{Vector{Int}}()
    triangle_actions = Vector{Vector{Int}}()

    for first_six in all_permutations(collect(1:6))
        tested += 1
        permutation = vcat(first_six, [7])
        ray_action = induced_ray_permutation(permutation, ray_trees, signature_to_ray)
        ray_action === nothing && continue
        ray_preserving += 1
        triangle_action = induced_triangle_permutation(ray_action, triangles, triangle_to_index)
        triangle_action === nothing && continue

        # The ray-incidence action must agree with directly relabelling the
        # rooted topology of a relative-interior tree in every triangle.
        for triangle_index in eachindex(triangle_trees)
            image_index = triangle_action[triangle_index]
            @assert rooted_topology_signature(triangle_trees[triangle_index], permutation) ==
                    rooted_topology_signature(triangle_trees[image_index])
        end

        push!(label_permutations, permutation)
        push!(ray_actions, ray_action)
        push!(triangle_actions, triangle_action)
    end

    @assert tested == factorial(6) == 720
    label_keys = Set(Tuple.(label_permutations))
    for left in label_permutations, right in label_permutations
        @assert Tuple(left[right]) in label_keys
    end

    rotation = [5, 6, 4, 2, 3, 1, 7]       # (1 5 3 4 2 6)
    reflection = [1, 3, 2, 4, 6, 5, 7]     # (2 3)(5 6)
    expected_group = generated_subgroup([rotation, reflection])
    expected_dihedral = Set(Tuple.(expected_group)) == label_keys

    ray_orbit_list = action_orbits(ray_actions, length(ray_trees))
    triangle_orbit_list = action_orbits(triangle_actions, length(triangles))

    complex = simplicial_complex(maximal_cones(IncidenceMatrix, fan))
    abstract_group = automorphism_group(complex)
    abstract_group_order = order(abstract_group)
    abstract_ray_orbit_list = ray_orbits(fan)
    abstract_triangle_orbit_list = triangle_orbits(fan)
    ray_refinement = refinement_description(abstract_ray_orbit_list, ray_orbit_list)
    triangle_refinement = refinement_description(abstract_triangle_orbit_list, triangle_orbit_list)

    write_group_csv(LABEL_GROUP_CSV, label_permutations, ray_actions, triangle_actions)
    write_ray_orbits_csv(RAY_ORBITS_CSV, ray_orbit_list, ray_trees)
    write_triangle_orbits_csv(
        TRIANGLE_ORBITS_CSV,
        triangle_orbit_list,
        triangles,
        triangle_trees,
    )
    write_report(
        ORBIT_REPORT;
        tested_permutations=tested,
        ray_preserving_permutations=ray_preserving,
        label_permutations=label_permutations,
        ray_orbits=ray_orbit_list,
        triangle_orbits=triangle_orbit_list,
        abstract_group_order=abstract_group_order,
        abstract_ray_orbits=abstract_ray_orbit_list,
        abstract_triangle_orbits=abstract_triangle_orbit_list,
        ray_refinement=ray_refinement,
        triangle_refinement=triangle_refinement,
        expected_dihedral=expected_dihedral,
    )

    println("tested label permutations: ", tested)
    println("ray-preserving permutations: ", ray_preserving)
    println("fan-preserving label group order: ", length(label_permutations))
    println("expected W(G2) dihedral group: ", expected_dihedral)
    println("label ray orbit sizes: ", sort(length.(ray_orbit_list)))
    println("label triangle orbit sizes: ", sort(length.(triangle_orbit_list)))
    println("abstract combinatorial group order: ", abstract_group_order)
    println("abstract ray orbit sizes: ", sort(length.(abstract_ray_orbit_list)))
    println("abstract triangle orbit sizes: ", sort(length.(abstract_triangle_orbit_list)))
    println("report: ", ORBIT_REPORT)
    return (
        label_permutations=label_permutations,
        ray_actions=ray_actions,
        triangle_actions=triangle_actions,
        ray_orbits=ray_orbit_list,
        triangle_orbits=triangle_orbit_list,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    compute_label_orbits()
end
