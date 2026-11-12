using Oscar
# UniTriSat is not a project dependency, since it pulls in CUDA. Add it to the
# default environment: ] add https://github.com/krhuang/UniTriSat/
using UniTriSat

const G2_POLYTOPE = normpath(joinpath(
    @__DIR__, "..", "..", "..", "..",
    "data", "chevalley",
    "chevalley-G2-varpi_0-1-word_2-1-2-1-2.poly",
))

# Load the native polymake object and wrap it as an Oscar polyhedron.
P = polyhedron(Polymake.load(G2_POLYTOPE))
lecture_hall = convex_hull(matrix(QQ, [
    0 0 0 0 0
    1 0 0 0 0
    1 3 0 0 0
    1 3 2 0 0
    1 3 2 3 0
    1 3 2 3 1
]))

@assert P == lecture_hall
@assert normalized_volume(P) == 18
@assert length(lattice_points(P)) == 14

# Enumerate every unimodular triangulation and classify regular/flag solutions.
full_result = triangulate(
    P;
    find_all=true,
    unimodular=true,
    regular=true,
    flag_triangulation=true,
    return_triangulations="none",
    terminal_output="final",
    parallel_split_solving=false,
)
@assert full_result.total_number_of_triangulations_found == 418
@assert full_result.total_number_of_regular_triangulations_found == 418
@assert full_result.total_number_of_flag_triangulations_found == 180
@assert full_result.total_number_of_quadratic_triangulations_found == 180

# Project the x_1 = 1 facet to its four-dimensional lattice coordinates.
facet_vertices = [
    [QQ(v[i]) for i in 2:5]
    for v in vertices(P)
    if v[1] == 1
]
P_d4 = convex_hull(matrix(QQ, reduce(vcat, permutedims.(facet_vertices))))
facet_result = triangulate(
    P_d4;
    find_all=true,
    unimodular=true,
    regular=true,
    flag_triangulation=true,
    return_triangulations="none",
    terminal_output="final",
    parallel_split_solving=false,
)
@assert facet_result.total_number_of_triangulations_found == 418
@assert facet_result.total_number_of_regular_triangulations_found == 418
@assert facet_result.total_number_of_flag_triangulations_found == 180
@assert facet_result.total_number_of_quadratic_triangulations_found == 180

println("G2 Chevalley triangulation checks passed.")
