###############################################################################
# A curve as an Oscar polyhedral complex, from its one-dimensional pieces with
# their multiplicities, for the stage zenodo of the census
# (studies/g2/census/zenodo_census.jl), which writes the two genus-10
# curves of Section 10, data/mukai_curves/genus_10_method_{1,2}.mrdi.
###############################################################################

include(joinpath(@__DIR__, "intersect_complexes.jl"))
include(joinpath(@__DIR__, "decorated_curve_graph.jl"))

"""
The cells `multiplicity` (generator sets of `homogenized_cells`, each with its
multiplicity) after merging the two cells at every vertex that lies in exactly
two cells, when they are collinear through it and have the same multiplicity.
Such a vertex is not a vertex of the curve: the computation cut an edge there,
on the diagonal of a quadrilateral of the coarsest fan, a wall of the raw fan.
Returns the merged cells and the number of vertices removed.
"""
function merge_collinear_cells(multiplicity::Dict{Vector{Int},ZZRingElem}, generators)
    is_point(g) = !iszero(generators[g, 1])
    point(g) = generators[g, 2:end] ./ generators[g, 1]
    cells = copy(multiplicity)
    removed = 0
    while true
        around = Dict{Int,Vector{Vector{Int}}}()
        for cell in keys(cells), g in cell
            is_point(g) && push!(get!(around, g, Vector{Int}[]), cell)
        end
        merged = false
        for (g, (a, b)) in filter(entry -> length(last(entry)) == 2, around)
            cells[a] == cells[b] || continue
            ends = [only(setdiff(cell, [g])) for cell in (a, b)]
            any(is_point, ends) || continue
            d1, d2 = [is_point(h) ? point(h) .- point(g) : generators[h, 2:end] for h in ends]
            rank(matrix(QQ, [permutedims(d1); permutedims(d2)])) == 1 && sum(d1 .* d2) < 0 || continue
            m = pop!(cells, a)
            pop!(cells, b)
            cells[sort(ends)] = m
            removed += 1
            merged = true
            break
        end
        merged || return cells, removed
    end
end

"""
The polyhedral complex with the one-dimensional `pieces` as maximal cells,
built by `homogenized_cells` of intersect_complexes.jl, and the multiplicity
of each of its maximal cells. Equal pieces are merged, and so are two pieces
that continue each other straight through a vertex of valence 2 (see
`merge_collinear_cells`); they must have equal multiplicities.
"""
function curve_complex(pieces, multiplicities)
    keep = findall(piece -> dim(piece) == 1, pieces)
    pieces = Vector{Polyhedron{QQFieldElem}}(pieces[keep])
    multiplicities = isnothing(multiplicities) ? nothing : multiplicities[keep]
    generators, lineality, generator_sets = homogenized_cells(pieces)
    multiplicity = Dict{Vector{Int},ZZRingElem}()
    for (index, cell) in enumerate(generator_sets)
        value = isnothing(multiplicities) ? ZZ(0) : ZZ(multiplicities[index])
        get!(multiplicity, cell, value) == value || error("two equal pieces have different multiplicities")
    end
    merged, removed = merge_collinear_cells(multiplicity, generators)
    cells = sort!(collect(keys(merged)))
    used = sort!(unique!(reduce(vcat, cells)))
    number = Dict(g => k for (k, g) in enumerate(used))
    renumbered = [sort!([number[g] for g in cell]) for cell in cells]
    complex = dehomogenized_complex(renumbered, generators[used, :], lineality)
    dim(complex) == 1 || error("the curve is not one-dimensional")
    n_maximal_polyhedra(complex) == length(cells) || error("the complex merged or dropped cells")
    # Match the maximal cells of the complex to the cells by their generators.
    generator_vectors(rows, set) = Set(rows[g, :] for g in set)
    multiplicity_of = Dict(generator_vectors(generators, cell) => merged[cell] for cell in cells)
    maximal_generators, _, maximal_sets = homogenized_cells(collect(maximal_polyhedra(complex)))
    cell_multiplicities = map(maximal_sets) do set
        get(multiplicity_of, generator_vectors(maximal_generators, set)) do
            error("a maximal cell of the complex is not a cell of the curve")
        end
    end
    return (complex=complex, multiplicities=cell_multiplicities, removed=removed)
end

"""The first Betti number of the graph of bounded edges of a one-dimensional complex: the `genus` of its graph."""
bounded_betti_number(complex) = genus(decorated_graph(collect(maximal_polyhedra(complex)), nothing))

"""
Save to `output` the complex of the one-dimensional `pieces` of a genus-10
curve by `curve_complex`, with its multiplicities unless `multiplicities` is
`nothing`.
"""
function write_curve(output, pieces, multiplicities)
    result = curve_complex(collect(pieces), multiplicities)
    betti = bounded_betti_number(result.complex)
    betti == 10 || error("the curve of $output has genus $betti, not 10")
    println("  $(length(pieces)) pieces, $(n_maximal_polyhedra(result.complex)) maximal cells after removing ",
            "$(result.removed) vertices of valence 2, first Betti number $betti of the bounded part")
    mkpath(dirname(output))
    save(output, isnothing(multiplicities) ? (complex=result.complex,) :
                                             (complex=result.complex, multiplicities=result.multiplicities))
    println("  wrote $output")
end
