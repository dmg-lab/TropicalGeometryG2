###############################################################################
# A basis of QQ^21 for the random translations of the tropical linear space of
# the genus-10 curves. The coarsest tropical G2-Grassmannian spans the space W
# of dimension 15 where the coordinates of the six linear equations with two
# terms, such as p23 + p47, are equal; modulo its lineality space, W has
# dimension 12. The basis is: the weights of 12 maximal cells, independent
# modulo lineality, each the sum of the normalized generators of the cell, as
# its canonical weight is; an orthogonal basis of the orthogonal complement of
# W, the six normals e_b - e_a, on disjoint pairs of coordinates; and an integral
# basis of the lineality space. The translation of an experiment is a random
# integral combination of it, see run_experiments.jl.
#
# The cells are chosen by how often the part of the curve in them is a
# polytope, see shifted_bergman_cells.jl. First come the orbits without
# transverse bad subsets, whose part is a polytope for every translation
# outside finitely many proper subspaces, then the others by the share of
# CELL_BASIS_SAMPLES translations uniformly random in a cube giving a polytope.
# A cell is kept when its weight is independent of those kept before.
#
#   julia --project=. studies/g2/curves/shift_basis.jl [--output=FILE.mrdi]
#
# By default the basis, its rows in this order, goes to
# data/mukai_curves/genus_10_shift_basis.mrdi, with the cells kept, the
# numbers of normals and of lineality vectors, and for every orbit, in the
# order of the choice, its transverse bad subsets and its share of polytopes.
###############################################################################

include(joinpath(@__DIR__, "shifted_bergman_cells.jl"))
import TropicalGeometryG2

const CELL_BASIS_SAMPLES = 2000

"An orthogonal basis of the row space of `M`, by Gram-Schmidt over QQ, its rows primitive integral vectors."
function orthogonal_rows(M)
    rows = Vector{QQFieldElem}[]
    for i in 1:size(M, 1)
        v = M[i, :]
        foreach(u -> v -= (sum(v .* u) // sum(u .* u)) .* u, rows)
        push!(rows, v)
    end
    return primitive_rows(reduce(vcat, permutedims.(rows)))
end

"The rows of `M` scaled to primitive integral vectors."
primitive_rows(M) = reduce(vcat, permutedims.([TropicalGeometryG2.primitive_integer(M[i, :]) for i in 1:size(M, 1)]))

if abspath(PROGRAM_FILE) == @__FILE__
    outputs = [split(arg, '='; limit=2)[2] for arg in ARGS if startswith(arg, "--output=")]
    length(outputs) == length(ARGS) <= 1 || error("usage: shift_basis.jl [--output=FILE.mrdi]")
    output = isempty(outputs) ? data_path("mukai_curves", "genus_10_shift_basis.mrdi") : abspath(only(outputs))
    rng = Random.Xoshiro(2026)

    fan = load(data_path("g2", "grassmannian", "02_coarsest_fan.mrdi"))
    cones = collect(maximal_cones(fan))
    generators = Matrix{QQFieldElem}(load(data_path("g2", "grassmannian", "03_normalized_rays.mrdi")))
    incidence = maximal_cones(IncidenceMatrix, load(data_path("g2", "grassmannian", "03_normalized_coarsest_fan.mrdi")))
    weights = reduce(vcat, [sum(generators[collect(row(incidence, c)), :]; dims=1) for c in 1:nrows(incidence)])
    lineality = primitive_rows(Matrix{QQFieldElem}(matrix(QQ, lineality_space(fan))))
    # Two vectors spanning the lineality modulo the diagonal.
    diagonal = permutedims(fill(QQ(1), 21))
    torus = let (i, j) = first((i, j) for i in 1:3 for j in i+1:3 if rank(matrix(QQ, vcat(diagonal, lineality[[i, j], :]))) == 3)
        [lineality[i, :], lineality[j, :]]
    end
    orbits = map(readlines(data_path("g2", "grassmannian", "05_symmetry_cell_orbits.csv"))[2:end]) do line
        fields = split(line, ',')
        (orbit=parse(Int, fields[1]), representative=parse(Int, fields[3]),
         cells=parse.(Int, split(strip(fields[4], '"'), ';')))
    end

    cells = cached_shift_cells(orbits, cones, torus)
    transverse, share = Dict{Int,Int}(), Dict{Int,Float64}()
    for o in orbits
        transverse[o.orbit] = transverse_bad_subsets(cells[o.orbit], cell_rows(cones[o.representative], torus))
        share[o.orbit] = transverse[o.orbit] == 0 ? 1.0 : polytope_share(cells[o.orbit], CELL_BASIS_SAMPLES, rng)
        println("orbit $(o.orbit): $(transverse[o.orbit]) transverse bad subsets, polytope for a share $(share[o.orbit])")
        flush(stdout)
    end
    order = sort(orbits; by=o -> (transverse[o.orbit] > 0, -share[o.orbit], transverse[o.orbit], o.orbit))
    span = rank(matrix(QQ, vcat(lineality, weights))) - size(lineality, 1)
    kept = Int[]
    for o in order, c in o.cells
        length(kept) == span && break
        rank(matrix(QQ, vcat(lineality, weights[[kept; c], :]))) > size(lineality, 1) + length(kept) && push!(kept, c)
    end
    normals = orthogonal_rows(Matrix{QQFieldElem}(transpose(kernel(matrix(QQ, vcat(lineality, weights)); side=:right))))
    @req all(iszero, normals * transpose(vcat(lineality, weights))) "the normals are not orthogonal to the span"
    basis = matrix(QQ, vcat(weights[kept, :], normals, lineality))
    @req rank(basis) == 21 "the weights, the normals and the lineality are not a basis"

    println("\norbit | transverse bad subsets | share of polytopes | cells kept")
    for o in order
        println(o.orbit, " | ", transverse[o.orbit], " | ", share[o.orbit], " | ", join(filter(in(kept), o.cells), ", "))
    end
    pairs = [(i, j) for i in 1:7 for j in i+1:7]
    println("\nNormals of the span: ", join((join(("$(n[i] > 0 ? "+" : "-")p$(pairs[i][1])$(pairs[i][2])" for i in 1:21 if !iszero(n[i])), " ")
                                          for n in eachrow(normals)), ", "))
    println("Index of the lattice of the basis in ZZ^21: ", abs(det(basis)))
    mkpath(dirname(output))
    save(output, (basis=basis, cells=kept, normals=size(normals, 1), lineality=size(lineality, 1),
                  orbits=[o.orbit for o in order], transverse_bad_subsets=[transverse[o.orbit] for o in order],
                  polytope_shares=[share[o.orbit] for o in order]))
    println("Wrote $output")
end
