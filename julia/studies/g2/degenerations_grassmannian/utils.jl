# Function for computing the fiber of a polytope
using Oscar
import Polymake
using TropicalGeometryG2: RecursivePolytopeTypes

"""
Save every unrecognized part of the description `r` of `P` with `Oscar.save`,
as `<name>.mrdi` for the root and `<name>_<path>.mrdi` for a part below it.
"""
function save_unrecognized(P, r, name; directory=joinpath(@__DIR__, "unrecognized"))
    for (path, X) in RecursivePolytopeTypes.unrecognized_parts(P, r)
        mkpath(directory)
        file = isempty(path) ? name : "$(name)_$(path)"
        Oscar.save(joinpath(directory, file * ".mrdi"), convex_hull(QQ, X))
    end
end

function fiber_polytope(P::Polyhedron{QQFieldElem}, indices)
    n = ambient_dim(P)
    I = collect(indices)

    all(i -> i isa Integer && 1 <= i <= n, I) ||
        throw(ArgumentError("indices must be integers in 1:$n"))
    length(unique(I)) == length(I) ||
        throw(ArgumentError("indices must be distinct"))
    is_feasible(P) && is_bounded(P) ||
        throw(ArgumentError("P must be a nonempty polytope"))

    # Projection on homogeneous row coordinates.
    pi = zeros(Int, n + 1, length(I) + 1)
    pi[1, 1] = 1
    for (j, i) in enumerate(I)
        pi[i + 1, j + 1] = 1
    end

    C = Polymake.polytope.fiber_polytope(P.pm_polytope, Polymake.Matrix{Polymake.Rational}(pi))
    W = Matrix{Rational{BigInt}}(C.POINTS)
    F = convex_hull(QQ, W[:, 3:end])
    J = setdiff(1:ambient_dim(F), I)
    
    # Drop the output homogenizing coordinate and the constant
    # coordinate inherited from the source's homogeneous vertices.
    return convex_hull(QQ, [[v[j] for j in J] for v in vertices(F)])
end


function determine_3_polytope(
    P::Polyhedron{QQFieldElem}; check_johnson::Bool=false
)
    is_feasible(P) && is_bounded(P) ||
        throw(ArgumentError("P must be a nonempty bounded polytope."))
    dim(P) == 3 ||
        throw(ArgumentError("Polytope must be 3-dimensional."))

    fv = f_vector(P)
    fs = sort(collect(facet_sizes(P)))

    function matches(Q)
        f_vector(Q) == fv &&
        sort(collect(facet_sizes(Q))) == fs &&
        Polymake.polytope.isomorphic(P.pm_polytope, Q.pm_polytope)
    end

    # Platonic types, checked first to prefer familiar names.
    for (name, constructor) in (
        ("tetrahedron", () -> simplex(3)),
        ("cube",        () -> cube(3)),
        ("octahedron",  () -> cross_polytope(3)),
        ("dodecahedron", dodecahedron),
        ("icosahedron", icosahedron),
    )
        matches(constructor()) && return name
    end

    archimedean_names = (
        "truncated_tetrahedron", "cuboctahedron",
        "truncated_cube", "truncated_octahedron",
        "rhombicuboctahedron", "truncated_cuboctahedron",
        "snub_cube", "icosidodecahedron",
        "truncated_dodecahedron", "truncated_icosahedron",
        "rhombicosidodecahedron", "truncated_icosidodecahedron",
        "snub_dodecahedron",
    )

    for name in archimedean_names
        matches(archimedean_solid(name)) &&
            return replace(name, "_" => " ")
    end

    catalan_names = (
        "triakis_tetrahedron", "triakis_octahedron",
        "rhombic_dodecahedron", "tetrakis_hexahedron",
        "disdyakis_dodecahedron", "pentagonal_icositetrahedron",
        "pentagonal_hexecontahedron", "rhombic_triacontahedron",
        "triakis_icosahedron", "deltoidal_icositetrahedron",
        "pentakis_dodecahedron", "deltoidal_hexecontahedron",
        "disdyakis_triacontahedron",
    )

    for name in catalan_names
        matches(catalan_solid(name)) &&
            return replace(name, "_" => " ")
    end

    # The vertex count determines the possible polygon size.
    nv = n_vertices(P)

    n = nv - 1
    if n >= 3 && matches(pyramid(cyclic_polytope(2, n)))
        return "$n-gonal pyramid"
    end

    if iseven(nv)
        n = nv ÷ 2
        if n >= 3 && matches(prism(cyclic_polytope(2, n)))
            return "$n-gonal prism"
        end
    end

    n = nv - 2
    if n >= 3 && matches(bipyramid(cyclic_polytope(2, n)))
        return "$n-gonal bipyramid"
    end

    if check_johnson
        for j in 1:92
            matches(johnson_solid(j)) &&
                return "Johnson type J$j"
        end
    end

    return "Unrecognized in this catalog; f-vector = $(collect(fv))"
end

