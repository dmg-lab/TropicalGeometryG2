include("utils.jl")
using TropicalGeometryG2: data_path

directory = data_path("g2", "grassmannian", "06_canonical_interior_degenerations_objects")

function characterize_from_directory(directory)
    for path in sort(readdir(directory; join=true))
        isfile(path) && endswith(path, ".poly") || continue
        # get the index
        m = match(r"^orbit_0(\d{2})_", basename(path))
        m === nothing && continue

        ij = parse(Int, m.captures[1])
        P= Oscar.load(path)
        result = RecursivePolytopeTypes.determine_polytope(P)
        println(ij, "| ", RecursivePolytopeTypes.describe_type(result))
        save_unrecognized(P, result, "orbit_$(lpad(ij, 3, '0'))_newton_okounkov_body")


        # Additional computations with P here
    end
end

characterize_from_directory(directory)