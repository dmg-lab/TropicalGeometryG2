###############################################################################
# When is the part of a genus-10 curve inside a maximal cell C of the coarsest
# tropical G2-Grassmannian a polytope? The curve is the intersection of the
# G2-Grassmannian with a translate B + t of the Bergman fan B of the uniform
# matroid of rank 17 on 21 elements, whose maximal cones are the cones B_T of
# the points whose minimum is attained on T, |T| = 5. Its part in C is
# C ∩ (B + t), the union of the pieces C ∩ (B_T + t). A nonempty piece is
# unbounded modulo the diagonal exactly when C ∩ B_T is larger than the
# diagonal: T is then a bad subset of C. So C ∩ (B + t) is a polytope exactly
# when t lies in none of the cones C - B_T of the bad subsets.
#
# For one cell of each orbit, this script finds the bad subsets, exactly, and
# tests translations in floating point: a point of every sector of the
# arrangement of the facets of C in its linear span, other points of the same
# sector, the same points plus random vectors orthogonal to the span, of
# relative lengths CELL_SHIFT_SCALES, and translations uniformly random in a
# cube, as in run_experiments.jl. It prints the table per orbit and
# writes it, with the sectors, to CSV files.
#
#   julia --project=. studies/g2/curves/shifted_bergman_cells.jl [--samples=N]
#
# The bad subsets are cached in work/curves/bad_subsets.mrdi and the
# tables go to results/curves/genus_10/.
###############################################################################

using Oscar
using Printf
import Combinatorics
import LinearAlgebra
import Random
using TropicalGeometryG2: data_path, project_path, work_path

const CELL_SHIFT_SCALES = (0.1, 1.0, 10.0)
const CELL_SHIFT_TRIALS = 20
const CELL_SHIFT_TOLERANCE = 1e-9

"""
A cell for the tests: its rays modulo lineality and two lineality vectors
besides the diagonal, as the rows of `G`, the first `k` the rays, and its bad
subsets, each with the matrices of the test of `piece_meets`.
"""
struct ShiftCell
    G::Matrix{Float64}
    k::Int
    bad::Vector{NamedTuple{(:T, :solve, :A, :E, :Q),
                           Tuple{Vector{Int}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}}
end

"The four rows `e_{T_1} - e_{T_i}` and the 16 rows `e_s - e_{T_1}`, `s` not in `T`."
function subset_rows(T)
    D = zeros(Int, 4, 21)
    for (i, t) in enumerate(T[2:end])
        D[i, T[1]] = 1
        D[i, t] = -1
    end
    S = setdiff(1:21, T)
    E = zeros(Int, length(S), 21)
    for (i, s) in enumerate(S)
        E[i, s] = 1
        E[i, T[1]] = -1
    end
    return D, E
end

"""
Whether the cone `{λ : Q λ >= 0}` of the plane or the line is larger than
the origin, `Q` exact with 1 or 2 columns, or with more by a linear program.
"""
function nonzero_cone(Q::Matrix{QQFieldElem})
    holds(λ) = !all(iszero, λ) && all(>=(0), Q * λ)
    d = size(Q, 2)
    d == 1 && return holds([QQ(1)]) || holds([QQ(-1)])
    if d == 2
        rows = [Q[i, :] for i in 1:size(Q, 1) if !all(iszero, Q[i, :])]
        isempty(rows) && return true
        return any(holds, [[-a[2], a[1]] for a in rows]) || any(holds, [[a[2], -a[1]] for a in rows])
    end
    # Q has no kernel, so the points with Q λ >= 0, sum(Q λ) = 1 are the nonzero ones up to scaling.
    return is_feasible(polyhedron((-Q, zeros(QQFieldElem, size(Q, 1))), (permutedims(vec(sum(Q; dims=1))), [QQ(1)])))
end

"""
The cell with rays `rays` modulo lineality and lineality spanned by the
diagonal and `torus`, with its bad subsets `T`: those for which some point of
the cell outside the diagonal has its minimum attained on `T`.
"""
function shift_cell(rays, torus)
    k = length(rays)
    G = matrix(QQ, reduce(vcat, permutedims.([rays; torus])))
    bad = fieldtype(ShiftCell, :bad)()
    for T in Combinatorics.combinations(1:21, 5)
        D, E = subset_rows(T)
        A = matrix(QQ, D) * transpose(G)
        K = Matrix{QQFieldElem}(kernel(A; side=:right))
        size(K, 2) == 0 && continue
        C = Matrix{QQFieldElem}(matrix(QQ, E) * transpose(G))
        Q = vcat(C * K, K[1:k, :])
        nonzero_cone(Q) || continue
        Af = Float64.(Matrix{QQFieldElem}(A))
        push!(bad, (T=T, solve=LinearAlgebra.pinv(Af), A=Af, E=Float64.(C), Q=Float64.(Q)))
    end
    return ShiftCell(Float64.(Matrix{QQFieldElem}(G)), k, bad)
end

"Whether some `λ` has `Q λ >= h` up to the tolerance, by Fourier-Motzkin elimination."
function fm_feasible(Q::AbstractMatrix{Float64}, h::AbstractVector{Float64})
    d = size(Q, 2)
    d == 0 && return all(<=(CELL_SHIFT_TOLERANCE), h)
    q = Q[:, d]
    upper = findall(>(CELL_SHIFT_TOLERANCE), q)
    lower = findall(<(-CELL_SHIFT_TOLERANCE), q)
    rows = [Q[i, 1:d-1] for i in eachindex(q) if abs(q[i]) <= CELL_SHIFT_TOLERANCE]
    rhs = [h[i] for i in eachindex(q) if abs(q[i]) <= CELL_SHIFT_TOLERANCE]
    for i in upper, j in lower
        push!(rows, Q[i, 1:d-1] ./ q[i] .- Q[j, 1:d-1] ./ q[j])
        push!(rhs, h[i] / q[i] - h[j] / q[j])
    end
    return fm_feasible(isempty(rows) ? zeros(0, d - 1) : reduce(vcat, permutedims.(rows)), rhs)
end

"""
Whether the piece of the bad subset `b` of `cell`, the points `x = u G` of
the cell with `x - t` in the cone of `b.T`, is nonempty: `A u = D t`, then
`u = u0 + K λ` with `E (x - t) >= 0` and the ray coefficients of `u` at least 0.
"""
function piece_meets(cell::ShiftCell, b, t)
    T = b.T
    Dt = [t[T[1]] - t[T[i]] for i in 2:5]
    u0 = b.solve * Dt
    LinearAlgebra.norm(b.A * u0 .- Dt) <= CELL_SHIFT_TOLERANCE * (1 + LinearAlgebra.norm(Dt)) || return false
    Et = [t[s] - t[T[1]] for s in 1:21 if !(s in T)]
    h = vcat(Et .- b.E * u0, -u0[1:cell.k])
    return fm_feasible(b.Q, h)
end

"Whether the part of the cell in the translate of the Bergman fan by `t` is a polytope."
bounded(cell::ShiftCell, t) = !any(b -> piece_meets(cell, b, t), cell.bad)

"""
The sectors of the arrangement of the facets of the cone `C` in its span,
each a representative unit vector of the span, the most central among
`samples` random combinations of the rays, and other points of the sector.
"""
function facet_sectors(C, rng; samples=4000)
    rays = [Float64.(collect(r)) for r in rays_modulo_lineality(C).rays_modulo_lineality]
    normals = Float64.(Matrix{QQFieldElem}(linear_inequality_matrix(facets(C))))
    sectors = Dict{Vector{Int},Vector{Tuple{Float64,Vector{Float64}}}}()
    for _ in 1:samples
        t = sum(randn(rng) .* r for r in rays)
        t ./= LinearAlgebra.norm(t)
        values = normals * t
        margin = minimum(abs.(values) ./ map(LinearAlgebra.norm, eachrow(normals)))
        push!(get!(sectors, Int.(sign.(values)), []), (margin, t))
    end
    return Dict(s => (centre=last(argmax(first, points)), points=last.(points)) for (s, points) in sectors)
end

"""
The number of bad subsets `T` of `cell` for which the cone `C - B_T` of the
translations meeting the piece is full-dimensional: those for which the span
of the cell, the rows of the exact `G` and the diagonal, and that of `B_T`
together span QQ^21. With none, the part of the curve in the cell is a
polytope for every translation outside finitely many proper subspaces.
"""
transverse_bad_subsets(cell::ShiftCell, G::Matrix{QQFieldElem}) =
    count(b -> rank(matrix(QQ, vcat(G, ones(QQFieldElem, 1, 21),
                                     [QQ(Int(i == s)) for s in setdiff(1:21, b.T), i in 1:21]))) == 21, cell.bad)

"The rays modulo lineality of the cone `C`, then the vectors `torus`, as the rows of an exact matrix."
cell_rows(C, torus) = reduce(vcat, permutedims.([[Vector{QQFieldElem}(collect(r)) for r in
                                                  rays_modulo_lineality(C).rays_modulo_lineality]; torus]))

"The share of `samples` translations uniformly random in the cube `[-1, 1]^21` for which the part of `cell` is a polytope."
polytope_share(cell::ShiftCell, samples, rng) = count(_ -> bounded(cell, 2 .* rand(rng, 21) .- 1), 1:samples) / samples

"""
The cells of `ShiftCell` of the representatives of `orbits` among `cones`,
read from the cache work/curves/bad_subsets.mrdi or computed and cached.
The cache holds vectors: Oscar saves no Dict of these cells once one of them
has no bad subsets.
"""
function cached_shift_cells(orbits, cones, torus)
    cache = work_path("curves", "bad_subsets.mrdi")
    if isfile(cache)
        stored = load(cache)
        return Dict(o => ShiftCell(c.G, c.k, c.bad) for (o, c) in zip(stored.orbits, stored.cells))
    end
    cells = Dict{Int,ShiftCell}()
    for o in orbits
        println("bad subsets of orbit $(o.orbit)")
        flush(stdout)
        rays = [Vector{QQFieldElem}(collect(r)) for r in rays_modulo_lineality(cones[o.representative]).rays_modulo_lineality]
        cells[o.orbit] = shift_cell(rays, torus)
    end
    mkpath(dirname(cache))
    save(cache, (orbits=collect(keys(cells)), cells=[(G=c.G, k=c.k, bad=c.bad) for c in values(cells)]))
    return cells
end

"An orthonormal basis of the complement of the span of the rays, the lineality and the diagonal of `C`."
function orthogonal_complement(cell::ShiftCell)
    return LinearAlgebra.nullspace(vcat(cell.G, ones(1, 21)))
end

if abspath(PROGRAM_FILE) == @__FILE__
    sample_options = [parse(Int, split(arg, '=')[2]) for arg in ARGS if startswith(arg, "--samples=")]
    length(sample_options) == length(ARGS) <= 1 || error("usage: shifted_bergman_cells.jl [--samples=N]")
    samples = isempty(sample_options) ? 500 : only(sample_options)
    rng = Random.Xoshiro(2026)

    fan = load(data_path("g2", "grassmannian", "02_coarsest_fan.mrdi"))
    cones = collect(maximal_cones(fan))
    lineality = Matrix{QQFieldElem}(matrix(QQ, lineality_space(fan)))
    diagonal = permutedims(fill(QQ(1), 21))
    i, j = first((i, j) for i in 1:3 for j in i+1:3 if rank(matrix(QQ, vcat(diagonal, lineality[[i, j], :]))) == 3)
    torus = [lineality[i, :], lineality[j, :]]
    orbits = map(readlines(data_path("g2", "grassmannian", "05_symmetry_cell_orbits.csv"))[2:end]) do line
        fields = split(line, ',')
        (orbit=parse(Int, fields[1]), size=parse(Int, fields[2]), representative=parse(Int, fields[3]))
    end
    types = Dict(parse(Int, split(l, ',')[3]) => parse(Int, split(l, ',')[1])
                 for l in readlines(data_path("g2", "maximal_cell_tree_metrics.csv"))[2:end] if !isempty(split(l, ',')[3]))

    cells = cached_shift_cells(orbits, cones, torus)

    output = project_path("results", "curves", "genus_10")
    mkpath(output)
    summary = open(joinpath(output, "shifted_bergman_cells.csv"), "w")
    println(summary, "orbit,type,cells,bad_subsets,transverse_bad_subsets,sectors,sectors_bounded,sector_points_agree,",
            join(("orthogonal_$(s)" for s in CELL_SHIFT_SCALES), ","), ",uniform")
    detail = open(joinpath(output, "shifted_bergman_sectors.csv"), "w")
    println(detail, "orbit,sector,bounded_at_centre,bounded_points_of_sector,",
            join(("bounded_orthogonal_$(s)" for s in CELL_SHIFT_SCALES), ","))
    println("| orbit | Type | cells | bad subsets | transverse | sectors | sectors giving a polytope | ",
            join(("orthogonal part $(s)" for s in CELL_SHIFT_SCALES), " | "), " | uniform shifts |")
    println("|---:|---:|---:|---:|---:|---:|---:|", repeat("---:|", length(CELL_SHIFT_SCALES)), "---:|")
    for o in orbits
        cell, C = cells[o.orbit], cones[o.representative]
        sectors = facet_sectors(C, rng)
        complement = orthogonal_complement(cell)
        bounded_sectors, agree = 0, true
        orthogonal = zeros(length(CELL_SHIFT_SCALES))
        for (sign_vector, sector) in sort!(collect(sectors); by=first)
            centre = bounded(cell, sector.centre)
            bounded_sectors += centre
            points = sector.points[1:min(end, CELL_SHIFT_TRIALS)]
            inside = count(t -> bounded(cell, t), points)
            agree &= inside in (0, length(points))
            fractions = map(CELL_SHIFT_SCALES) do scale
                count(1:CELL_SHIFT_TRIALS) do _
                    v = complement * randn(rng, size(complement, 2))
                    bounded(cell, sector.centre .+ scale .* v ./ LinearAlgebra.norm(v))
                end / CELL_SHIFT_TRIALS
            end
            orthogonal .+= fractions ./ length(sectors)
            println(detail, "$(o.orbit),", join(s > 0 ? "+" : "-" for s in sign_vector), ",$centre,",
                    "$inside/$(length(points)),", join(fractions, ","))
        end
        uniform = polytope_share(cell, samples, rng)
        transverse = transverse_bad_subsets(cell, cell_rows(C, torus))
        @printf("| %d | %d | %d | %d | %d | %d | %d | %s | %.2f |\n", o.orbit, types[o.orbit], o.size, length(cell.bad), transverse,
                length(sectors), bounded_sectors, join((@sprintf("%.2f", f) for f in orthogonal), " | "), uniform)
        println(summary, "$(o.orbit),$(types[o.orbit]),$(o.size),$(length(cell.bad)),$transverse,$(length(sectors)),$bounded_sectors,$agree,",
                join(orthogonal, ","), ",$uniform")
        flush(stdout)
    end
    close(summary)
    close(detail)
    println("\nWrote $(joinpath(output, "shifted_bergman_cells.csv")) and $(joinpath(output, "shifted_bergman_sectors.csv"))")
end
