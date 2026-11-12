using Dates
using Random
using LinearAlgebra

# Computes the certificate hexagon_certificate.mrdi of Theorem thm:completing_hexagons (README.md): for each of the 56
# coordinates f of Lemma lem:complete_hexagon, a subdivision of every maximal cone C of Trop = trop(linespace) ×
# trop(pointspace) into closed cells σ, each with a tropically faithful representative f_σ = x^(-u) g ∈ f + I. Uses only
# the definitions (polynomial_expressions.jl, tropical_complete_hexagon.jl) and the algorithm for tropically faithful
# representatives (tropically_faithful_representatives.jl). Min convention throughout.
#
# Stages, each resumable (a rerun continues), writing below the work folder (--work=DIR, default certificate/work/):
#
#   cones     The dimensions of the 6017 cones of Trop, in cones.csv (the maximal ones are subdivided by cover).
#   setups   Steps 1-3 of the algorithm on every cone of Trop (6017: the 547 faces of the raw fan of trop(linespace) times
#             the 11 cones of the tropical quadric), at the point of the cone: the reduced Gröbner basis G of I, in_w(G),
#             J = in_w(<G>) : (x_1 ⋯ x_28)^∞, its minimal primes over QQ and the exponent m of the monomial of the loop.
#             One file setups/<face>_<quadric>.mrdi per cone.
#   greedy    Starting lists R_f, first part: three sample points w on every cone (its point and two random points of its
#             relative interior). At each w the representatives found so far are tested, and if none is faithful the loop
#             of the algorithm gives a new one. Writes greedy/samples_<part>.csv (the value val(f(x)) of a generic lift x
#             of w for each coordinate) and greedy/candidates_<part>.mrdi. The parts are contiguous blocks of cones.
#   merge     Starting lists R_f, second part: a greedy set cover of all samples by the candidates. Writes greedy/R.mrdi.
#   cover     The lazy subdivision of each maximal cone C for each coordinate f: C is cut by the normal fan of the first
#             representative h_1 in R_f; cells on which h_1 is faithful are kept, the others are cut by h_2, and so on. A
#             cell on which no representative of R_f is faithful gets a new one by the algorithm at a point of it, appended
#             to R_f. Faithfulness on a cell is tested at the sum of its rays, a point of its relative interior: in_w(g) is
#             constant there. Writes cover/trees_<coordinate>.csv and cover/representatives_<coordinate>.mrdi.
#   assemble  Writes hexagon_certificate.mrdi from cover/, or from --from=DIR (a folder with these files, or with trees.csv
#             and representatives.mrdi, the format of an earlier certificate).
#
# From the `julia` folder (run_compute_overnight.sh runs all stages):
#     julia --project=. studies/hexagon/certificate/compute_certificate.jl --stage=cones
#     julia --project=. studies/hexagon/certificate/compute_certificate.jl --stage=setups [--part=i --parts=n]
#     julia --project=. studies/hexagon/certificate/compute_certificate.jl --stage=greedy --part=i --parts=4    (i = 1..4)
#     julia --project=. studies/hexagon/certificate/compute_certificate.jl --stage=merge
#     julia --project=. studies/hexagon/certificate/compute_certificate.jl --stage=cover [--coordinates=1,2] [--part=i --parts=n]
#     julia --project=. studies/hexagon/certificate/compute_certificate.jl --stage=assemble [--from=DIR]
# `--cones=N` restricts greedy and cover to the first N cones of their part (a test). The parts of greedy change the
# starting lists (the random points of part i use the seed i); --parts=4 is what the certificate of 2026-11-07 used. The
# parts of setups and cover only distribute the work. Coordinates are numbered 1..56: x2_1..7, x3_1..7, p23_1..21, p34_1..21.

include(joinpath(@__DIR__, "..", "tropical_complete_hexagon.jl"))

const CERTIFICATE_PATH = joinpath(@__DIR__, "hexagon_certificate.mrdi")
const GREEDY_MULTIPLICATIONS = (0, 1, 2, 4, 8, 16, 32, 64)
const SAMPLES_PER_CONE = 3

# --- The coordinates and the cones of Trop ------------------------------------------------------------------

"The 56 coordinates f of Lemma lem:complete_hexagon, as (name, f): x2_1, ..., x2_7, x3_1, ..., p34_21."
function certificate_coordinates()
    H = hexagon_expressions()
    return [("$(name)_$c", H[name][c]) for name in (:x2, :x3, :p23, :p34) for c in eachindex(H[name])]
end

"A polynomial with the variables of `EXPRESSION_RING` as an element of it (each `load` creates a new ring)."
function to_expression_ring(p)
    B = MPolyBuildCtx(EXPRESSION_RING)
    foreach((c, e) -> push_term!(B, QQ(c), e), AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
    return finish(B)
end

"The cones of the tropical quadric: the sets S of at least two terms among x1x4, x2x5, x3x6, x7^2 of minimal weight."
quadric_cones() = [S for d in 2:4 for S in subsets(collect(1:4), d)]

"All cones of Trop, as (face of the raw fan, cone of the quadric)."
trop_cones() = [(k, S) for k in eachindex(raw_fan_faces()) for S in quadric_cones()]

"A closed cone {E w = 0, A w ≤ 0} with its integer rays modulo lineality, a lineality basis, and the sum of the rays."
struct Region
    E
    A
    rays::Vector{Vector{Int}}
    lineality::Vector{Vector{Int}}
    point::Vector{Int}
end

function Region(E, A)
    R = rays_modulo_lineality(cone_from_inequalities(A, E))
    rays = [integer_weight(r) for r in R[1]]
    lineality = [integer_weight(l) for l in R[2]]
    return Region(E, A, rays, lineality, isempty(rays) ? zeros(Int, ncols(E)) : sum(rays))
end

"The weights of the terms x1x4, x2x5, x3x6, x7^2 of q as linear forms on R^28."
function quadric_term_form(t)
    v = zeros(Int, 28)
    t <= 3 ? (v[21 + t] = 1; v[24 + t] = 1) : (v[28] = 2)
    return v
end

"The cone of Trop of the face `k` of the raw fan and the cone `S` of the quadric."
function product_region(k, S)
    face = raw_fan_faces()[k]
    padded(M) = hcat(M, zero_matrix(QQ, nrows(M), 7))
    equations = [permutedims(quadric_term_form(t) .- quadric_term_form(S[1])) for t in S[2:end]]
    inequalities = [permutedims(quadric_term_form(S[1]) .- quadric_term_form(t)) for t in setdiff(1:4, S)]
    E = vcat(padded(face.E), matrix(QQ, reduce(vcat, equations)))
    A = isempty(inequalities) ? padded(face.A) : vcat(padded(face.A), matrix(QQ, reduce(vcat, inequalities)))
    return Region(E, A)
end

"""
The maximal cones of Trop, as (face, quadric string), sorted. Their dimensions are cached in cones.csv of the work folder.
"""
function maximal_trop_cones(work)
    path = joinpath(work, "cones.csv")
    if !isfile(path)
        span_dim(R) = length(independent_rays(vcat(R.rays, R.lineality)))
        rows = ["$k,$(join(S)),$(span_dim(product_region(k, S)))" for (k, S) in trop_cones()]
        write(path * ".tmp", "face,quadric,dim\n" * join(rows, "\n") * "\n")
        mv(path * ".tmp", path; force=true)
    end
    rows = [split(l, ",") for l in Iterators.drop(eachline(path), 1)]
    top = maximum(parse(Int, r[3]) for r in rows)
    return sort([(parse(Int, r[1]), String(r[2])) for r in rows if parse(Int, r[3]) == top])
end

parse_quadric(Sstring) = [parse(Int, ch) for ch in Sstring]

function append_csv(path, row)
    isfile(path) || write(path, join(string.(keys(row)), ",") * "\n")
    open(io -> println(io, join(string.(values(row)), ",")), path, "a")
end

# --- Stage setups -------------------------------------------------------------------------------------------

setup_path(work, k, Sstring) = joinpath(work, "setups", "$(k)_$(Sstring).mrdi")

"Steps 1-3 of the algorithm at the point of the cone (k, S), saved to the work folder."
function save_setup(work, k, S)
    T = faithful_setup(line_point_ideal(), vcat(raw_fan_faces()[k].point, quadric_cone_point(S)); saturate=false)
    wi = integer_weight(T.weight)
    data = Dict{String,Any}(
        "face" => k, "quadric" => S, "weight" => wi, "order_weight" => (maximum(wi) + 1) .- wi,
        "groebner_basis" => collect(T.groebner_basis), "initial_ideal" => gens(T.initial_ideal),
        "saturation" => gens(T.saturation), "k" => T.k, "primes" => [gens(P) for P in T.primes], "multiplier" => T.multiplier)
    path = setup_path(work, k, join(S))
    save(path * ".tmp", data)
    mv(path * ".tmp", path; force=true)
end

"The setup of the cone (k, Sstring) from the work folder, as a `FaithfulSetup` in `EXPRESSION_RING`."
function load_setup(work, k, Sstring)
    loaded = load(setup_path(work, k, Sstring))
    G = to_expression_ring.(loaded["groebner_basis"])
    S = EXPRESSION_RING
    ordering = weight_ordering(Int.(loaded["order_weight"]), degrevlex(S))
    return FaithfulSetup(ideal(S, G), QQ.(loaded["weight"]), ordering, G, ideal(S, to_expression_ring.(loaded["initial_ideal"])),
                         ideal(S, to_expression_ring.(loaded["saturation"])),
                         [ideal(S, to_expression_ring.(P)) for P in loaded["primes"]], loaded["k"], Int.(loaded["multiplier"]),
                         tropical_semiring_map(QQ), Dict())
end

function setups_stage(work, part, parts; log)
    mkpath(joinpath(work, "setups"))
    mine = [c for (n, c) in enumerate(trop_cones()) if (n - 1) % parts == part - 1 && !isfile(setup_path(work, c[1], join(c[2])))]
    log("setups, part $part: $(length(mine)) cones to do")
    for (n, (k, S)) in enumerate(mine)
        save_setup(work, k, S)
        n % 100 == 0 && log("setups, part $part: $n/$(length(mine))")
        n % 200 == 0 && GC.gc(true)
    end
    log("setups, part $part finished")
end

# --- Representatives -----------------------------------------------------------------------------------------

"Whether x^(-u) g is faithful at w for the setup T: in_w(g) lies in no minimal prime of J."
faithful_at(T, g, w) = (h = initial(g, T.tropical_map, QQ.(w)); !any(P -> h in P, T.primes))

"trop(x^(-u) g)(w) for integer weights."
tropical_at(g, u, w) = minimum(sum(e .* w) for e in AbstractAlgebra.exponent_vectors(g)) - sum(u .* w)

"""
A representative x^(-u) g of f faithful at w by the loop of the algorithm with the setup T, as (g, u) with
g = NF(x^u f) and u = K m, trying K = 0, 1, 2, 4, ..., or `nothing`.
"""
function faithful_representative_at(T, f, w)
    S = parent(f)
    for K in GREEDY_MULTIPLICATIONS
        (K > 0 && iszero(T.multiplier)) && break
        g = reduce_modulo(T, prod(gens(S) .^ (K .* T.multiplier)) * f)
        iszero(g) && error("f lies in I")
        faithful_at(T, g, w) && return (g, K .* T.multiplier)
    end
    return nothing
end

save_representatives(path, reps) = save(path, Dict{String,Any}("g" => first.(reps), "u" => last.(reps)))
load_representatives(path) = (d = load(path); [(to_expression_ring(g), Int.(u)) for (g, u) in zip(d["g"], d["u"])])

# --- Stages greedy and merge ---------------------------------------------------------------------------------

"Random points of the relative interior of a region: positive combinations of its rays plus combinations of its lineality."
random_point(rng, R) = sum(rand(rng, 1:20) * v for v in R.rays; init=zeros(Int, 28)) .+
                       sum(rand(rng, -20:20) * v for v in R.lineality; init=zeros(Int, 28))

function greedy_stage(work, part, parts; cones=nothing, log)
    dir = joinpath(work, "greedy")
    mkpath(dir)
    coordinates = certificate_coordinates()
    all_cones = trop_cones()
    # contiguous parts, so that consecutive cones are neighbours and share representatives
    chunk = cld(length(all_cones), parts)
    mine = all_cones[(part - 1) * chunk + 1:min(part * chunk, length(all_cones))]
    cones === nothing || (mine = first(mine, cones))
    samples_path = joinpath(dir, "samples_$part.csv")
    candidates_path = joinpath(dir, "candidates_$part.mrdi")
    done = isfile(samples_path) ?
        Set((parse(Int, split(l, ",")[1]), split(l, ",")[2]) for l in Iterators.drop(eachline(samples_path), 1)) : Set()
    # candidates[i]: the representatives (g, u) found for coordinate i, most recently successful first
    candidates = [Tuple{Any,Vector{Int}}[] for _ in coordinates]
    if isfile(candidates_path)
        data = load(candidates_path)
        for (i, g, u) in zip(data["coordinate"], data["g"], data["u"])
            push!(candidates[i], (to_expression_ring(g), Int.(u)))
        end
    end
    function save_candidates()
        is = [i for i in eachindex(candidates) for _ in candidates[i]]
        save(candidates_path * ".tmp", Dict{String,Any}("coordinate" => is, "g" => [c[1] for cs in candidates for c in cs],
                                                         "u" => [c[2] for cs in candidates for c in cs]))
        mv(candidates_path * ".tmp", candidates_path; force=true)
    end
    rng = MersenneTwister(part)
    log("greedy, part $part: $(length(mine)) cones, $(length(done)) done before, $(sum(length, candidates)) candidates")
    new_since_save = 0
    for (n, (k, S)) in enumerate(mine)
        # the random points are drawn before the check, so that a rerun continues with the same points
        R = product_region(k, S)
        T = (k, join(S)) in done ? nothing : load_setup(work, k, join(S))
        points = vcat([T === nothing ? Int[] : Int.(T.weight)], [random_point(rng, R) for _ in 2:SAMPLES_PER_CONE])
        T === nothing && continue
        values = zeros(Int, length(points), length(coordinates))
        for (i, (_, f)) in enumerate(coordinates), (j, w) in enumerate(points)
            hit = findfirst(c -> faithful_at(T, c[1], w), candidates[i])
            if hit === nothing
                r = faithful_representative_at(T, f, w)
                r === nothing && error("cone ($k, $(join(S))): no decision at $w for coordinate $i")
                pushfirst!(candidates[i], r)
                new_since_save += 1
            elseif hit > 1
                pushfirst!(candidates[i], popat!(candidates[i], hit))   # move to front
            end
            values[j, i] = tropical_at(candidates[i][1]..., w)
        end
        open(samples_path, "a") do io
            filesize(samples_path) == 0 && println(io, "face,quadric,sample,point,", join(first.(coordinates), ","))
            for (j, w) in enumerate(points)
                println(io, "$k,$(join(S)),$j,$(join(w, " ")),", join(values[j, :], ","))
            end
        end
        if new_since_save >= 20 || n == length(mine)
            save_candidates()
            new_since_save = 0
        end
        n % 100 == 0 && log("greedy, part $part: $n/$(length(mine)) cones, $(sum(length, candidates)) candidates")
        n % 50 == 0 && GC.gc(true)
    end
    save_candidates()
    log("greedy, part $part finished: candidates per coordinate $(length.(candidates))")
end

"Exponent matrix (terms × 28) of g as Float64, for fast evaluation of trop(g) at many points (entries are small integers)."
exponent_matrix(g) = Float64.(reduce(vcat, [permutedims(e) for e in AbstractAlgebra.exponent_vectors(g)]))

function merge_stage(work; log)
    dir = joinpath(work, "greedy")
    coordinates = certificate_coordinates()
    sample_points = Vector{Int}[]
    sample_cones = Tuple{Int,String}[]
    values = Vector{Int}[]
    for file in filter(startswith("samples_"), readdir(dir))
        for l in Iterators.drop(eachline(joinpath(dir, file)), 1)
            r = split(l, ",")
            push!(sample_cones, (parse(Int, r[1]), r[2]))
            push!(sample_points, parse.(Int, split(r[4])))
            push!(values, parse.(Int, r[5:end]))
        end
    end
    W = Float64.(reduce(hcat, sample_points))          # 28 × samples
    V = reduce(hcat, values)                              # coordinates × samples
    log("merge: $(size(W, 2)) samples")
    candidates = [Tuple{Any,Vector{Int}}[] for _ in coordinates]
    for file in filter(f -> startswith(f, "candidates_") && endswith(f, ".mrdi"), readdir(dir))
        data = load(joinpath(dir, file))
        for (i, g, u) in zip(data["coordinate"], data["g"], data["u"])
            push!(candidates[i], (to_expression_ring(g), Int.(u)))
        end
    end
    chosen = Vector{Tuple{Any,Vector{Int}}}[]
    for (i, (name, f)) in enumerate(coordinates)
        covers = map(candidates[i]) do (g, u)
            values_at = vec(minimum(exponent_matrix(g) * W; dims=1)) .- vec(Float64.(u)' * W)
            BitVector(values_at .== V[i, :])
        end
        uncovered = trues(size(W, 2))
        R = Tuple{Any,Vector{Int}}[]
        while any(uncovered) && !isempty(covers)
            gains = [count(cv .& uncovered) for cv in covers]
            best = argmax(gains)
            gains[best] == 0 && break
            push!(R, candidates[i][best])
            uncovered .&= .!covers[best]
        end
        # samples covered by no candidate (e.g. after an interrupted part): the algorithm there
        for j in findall(uncovered)
            r = faithful_representative_at(load_setup(work, sample_cones[j]...), f, sample_points[j])
            r === nothing && error("no decision at sample $j for coordinate $name")
            push!(R, r)
            gv = vec(minimum(exponent_matrix(r[1]) * W; dims=1)) .- vec(Float64.(r[2])' * W)
            uncovered .&= .!(gv .== V[i, :])
        end
        push!(chosen, R)
        log("merge: $name: $(length(candidates[i])) candidates, $(length(R)) chosen")
    end
    save(joinpath(dir, "R.mrdi"), Dict{String,Any}("coordinate" => first.(coordinates),
                                                    "g" => [first.(R) for R in chosen], "u" => [last.(R) for R in chosen]))
    log("merge finished")
end

# --- Stage cover ---------------------------------------------------------------------------------------------

"Independent rows among `rays`."
function independent_rays(rays)
    B = Vector{Int}[]
    for r in rays
        candidate = vcat(B, [r])
        rank(matrix(QQ, reduce(vcat, [permutedims(b) for b in candidate]))) == length(candidate) && (B = candidate)
    end
    return B
end

"""
The maximal cells of the normal fan of the Newton polytope of a polynomial with exponent vectors `exps` on the closed
region `R`, as inequality matrices A with {E w = 0, A w <= 0} (E = R.E): one for each vertex of the Newton polyhedron
conv(exps) + R^∨ in coordinates on the span of R. A single cell `R.A` if the polynomial is a monomial on R.
"""
function maximal_cells(R::Region, exps)
    B = independent_rays(vcat(R.rays, R.lineality))
    groups = Dict{Vector{QQFieldElem},Vector{Vector{Int}}}()
    for e in exps
        push!(get!(groups, QQ.([dot(e, b) for b in B]), Vector{Int}[]), e)
    end
    forms = collect(keys(groups))
    length(forms) == 1 && return [R.A]
    as_matrix(rows) = isempty(rows) ? zero_matrix(QQ, 0, length(B)) : matrix(QQ, reduce(vcat, rows))
    M = as_matrix(permutedims.(B))
    coordinates(v) = permutedims(collect(solve(M, QQ.(v); side=:left)))
    dual = polarize(positive_hull(as_matrix(coordinates.(R.rays)), as_matrix(coordinates.(R.lineality))))
    N = convex_hull(as_matrix(permutedims.(forms)), as_matrix([permutedims(QQFieldElem.(collect(r))) for r in rays(dual)]))
    vertex_forms = [QQFieldElem.(collect(v)) for v in vertices(N)]
    length(vertex_forms) == 1 && return [R.A]
    return map(vertex_forms) do v
        base = groups[v][1]
        vcat(R.A, matrix(QQ, reduce(vcat, [permutedims(base .- groups[u][1]) for u in vertex_forms if u != v])))
    end
end

"""
The tree of the subdivision of the cell `region`, which is cut by the normal fan of the j-th representative (appended to
`reps` by the algorithm at the point of the cell if there is none). The children that are not faithful are subdivided
further in reverse order: the order in which a stack of cells was worked off for the certificate of 2026-11-07, which
decides where new representatives are computed.
"""
function subdivision_tree(T, f, reps, region, j; log, where)
    if j > length(reps)
        r = faithful_representative_at(T, f, region.point)
        r === nothing && error("$where: no decision at $(region.point)")
        push!(reps, r)
        log("  $where: new representative $(length(reps)) at $(region.point)")
    end
    g = reps[j][1]
    exps = collect(AbstractAlgebra.exponent_vectors(g))
    sorted = sort(exps)
    children = Tuple{Int,Any}[]
    open_cells = Tuple{Int,Region}[]
    for A in maximal_cells(region, exps)
        cell = Region(region.E, A)
        # the vertex of the cell: the smallest exponent of minimal weight at a point of its relative interior
        at = [dot(e, cell.point) for e in sorted]
        t = findfirst(==(minimum(at)), at)
        push!(children, (t, "."))
        faithful_at(T, g, cell.point) || push!(open_cells, (length(children), cell))
    end
    for (n, cell) in reverse(open_cells)
        children[n] = (children[n][1], subdivision_tree(T, f, reps, cell, j + 1; log, where))
    end
    allunique(first.(children)) || error("$where: two cells with the same vertex")
    return "$j(" * join(["$t:$c" for (t, c) in sort(children; by=first)], ",") * ")"
end

cover_trees_path(work, name) = joinpath(work, "cover", "trees_$name.csv")
cover_reps_path(work, name) = joinpath(work, "cover", "representatives_$name.mrdi")

function cover_coordinate(work, i; cones=nothing, log)
    name, f = certificate_coordinates()[i]
    path, rpath = cover_trees_path(work, name), cover_reps_path(work, name)
    reps = isfile(rpath) ? load_representatives(rpath) : (d = load(joinpath(work, "greedy", "R.mrdi"));
        [(to_expression_ring(g), Int.(u)) for (g, u) in zip(d["g"][i], d["u"][i])])
    isfile(rpath) || save_representatives(rpath, reps)
    done = isfile(path) ? Set((parse(Int, split(l, ",")[1]), split(l, ",")[2]) for l in Iterators.drop(eachline(path), 1)) : Set()
    maximal = maximal_trop_cones(work)
    cones === nothing || (maximal = first(maximal, cones))
    log("cover $name: $(length(reps)) representatives, $(length(maximal)) maximal cones, $(length(done)) done before")
    started = time()
    for (n, (k, Sstring)) in enumerate(maximal)
        (k, Sstring) in done && continue
        before = length(reps)
        tree = subdivision_tree(load_setup(work, k, Sstring), f, reps, product_region(k, parse_quadric(Sstring)), 1;
                                log, where="$name, cone ($k, $Sstring)")
        length(reps) > before && save_representatives(rpath, reps)   # before the tree, so that its indices are saved
        isfile(path) || write(path, "face,quadric,tree\n")
        open(io -> println(io, "$k,$Sstring,$tree"), path, "a")
        n % 100 == 0 && (log("cover $name: $n/$(length(maximal)) cones, $(length(reps)) representatives"); GC.gc(true))
    end
    log("cover $name: done, $(length(reps)) representatives, $(round(time() - started; digits=1)) s")
end

function cover_stage(work, coordinates, part, parts; cones=nothing, log)
    mkpath(joinpath(work, "cover"))
    foreach(i -> cover_coordinate(work, i; cones, log), coordinates[part:parts:end])
    log("cover, part $part finished")
end

# --- Stage assemble ------------------------------------------------------------------------------------------

const CERTIFICATE_DESCRIPTION = """
Certificate for Theorem thm:completing_hexagons, see README.md. "coordinate": the names of the 56 coordinates f of Lemma
lem:complete_hexagon, "f": the polynomials f, "g" and "u": for each coordinate its representatives h_j = x^(-u_j) g_j,
"face" and "quadric": the maximal cones C of Trop, "tree": for each maximal cone and each coordinate the tree of the
subdivision of C, node := "." | j "(" t ":" node ("," t ":" node)* ")"."""

"""
Writes hexagon_certificate.mrdi from the folder `from`: either trees_<coordinate>.csv and representatives_<coordinate>.mrdi
of the stage cover, or trees.csv and representatives.mrdi of the certificate of 2026-11-07.
"""
function assemble_stage(from; log)
    coordinates = certificate_coordinates()
    names = first.(coordinates)
    trees = Dict{Tuple{Int,String},Dict{String,String}}()
    add_tree(k, Sstring, name, tree) = (get!(trees, (parse(Int, k), String(Sstring)), Dict{String,String}())[name] = tree)
    if isfile(joinpath(from, "trees.csv"))
        for l in Iterators.drop(eachline(joinpath(from, "trees.csv")), 1)
            add_tree(split(l, ","; limit=4)...)
        end
        d = load(joinpath(from, "representatives.mrdi"))
        collect(d["coordinate"]) == names || error("representatives.mrdi: the coordinates are not $(join(names, " "))")
        reps = [[(to_expression_ring(g), Int.(u)) for (g, u) in zip(gs, us)] for (gs, us) in zip(d["g"], d["u"])]
    else
        for name in names, l in Iterators.drop(eachline(joinpath(from, "trees_$name.csv")), 1)
            k, Sstring, tree = split(l, ","; limit=3)
            add_tree(k, Sstring, name, tree)
        end
        reps = [load_representatives(joinpath(from, "representatives_$name.mrdi")) for name in names]
    end
    cones = sort(collect(keys(trees)))
    for c in cones
        Set(keys(trees[c])) == Set(names) || error("cone $c: trees for $(length(trees[c])) coordinates")
    end
    log("assemble: $(length(cones)) maximal cones, $(sum(length, reps)) representatives")
    save(CERTIFICATE_PATH * ".tmp", Dict{String,Any}(
        "description" => CERTIFICATE_DESCRIPTION,
        "coordinate" => names, "f" => last.(coordinates),
        "g" => [first.(r) for r in reps], "u" => [last.(r) for r in reps],
        "face" => first.(cones), "quadric" => last.(cones),
        "tree" => [[trees[c][name] for name in names] for c in cones]))
    mv(CERTIFICATE_PATH * ".tmp", CERTIFICATE_PATH; force=true)
    log("wrote $CERTIFICATE_PATH")
end

# --- Main ----------------------------------------------------------------------------------------------------

function compute_main(args)
    option(name) = (a = findfirst(startswith("--$name="), args); a === nothing ? nothing : split(args[a], "="; limit=2)[2])
    log(s) = (println("[$(now())] $s"); flush(stdout))
    work = option("work") === nothing ? joinpath(@__DIR__, "work") : abspath(option("work"))
    part = option("part") === nothing ? 1 : parse(Int, option("part"))
    parts = option("parts") === nothing ? 1 : parse(Int, option("parts"))
    cones = option("cones") === nothing ? nothing : parse(Int, option("cones"))
    stage = option("stage")
    mkpath(work)
    if stage == "cones"
        log("$(length(maximal_trop_cones(work))) maximal cones of Trop, in $(joinpath(work, "cones.csv"))")
    elseif stage == "setups"
        setups_stage(work, part, parts; log)
    elseif stage == "greedy"
        greedy_stage(work, part, parts; cones, log)
    elseif stage == "merge"
        merge_stage(work; log)
    elseif stage == "cover"
        coordinates = option("coordinates") === nothing ? collect(1:56) : parse.(Int, split(option("coordinates"), ","))
        cover_stage(work, coordinates, part, parts; cones, log)
    elseif stage == "assemble"
        assemble_stage(option("from") === nothing ? joinpath(work, "cover") : abspath(option("from")); log)
    else
        error("--stage=cones|setups|greedy|merge|cover|assemble expected")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    compute_main(ARGS)
end
