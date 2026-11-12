include("utils.jl")
import .RecursivePolytopeTypes
directory = joinpath("julia/data/g2/grassmannian/06_canonical_interior_degenerations_objects/")

function characterize_from_directory(directory)
    for path in sort(readdir(directory; join=true))
        isfile(path) && endswith(path, ".poly") || continue
        # get the index
        m = match(r"^orbit_0(\d{2})_", basename(path))
        m === nothing && continue

        ij = parse(Int, m.captures[1])
        P= Oscar.load(path)
        F= fiber_polytope(P,[4,5])
        result = RecursivePolytopeTypes.determine_polytope(F)
        println( ij, "| ", RecursivePolytopeTypes.describe_type(result))
        save_unrecognized(F, result, "orbit_$(lpad(ij, 3, '0'))_fiber")

        # Additional computations with P here
    end
end

characterize_from_directory(directory)