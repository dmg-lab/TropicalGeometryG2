using Dates

# Independent verification of the certificate hexagon_certificate.mrdi for Theorem thm:completing_hexagons. Only the
# definitions are taken from the other files: the ideal I of linespace × pointspace and the 56 coordinates f of Lemma
# lem:complete_hexagon (polynomial_expressions.jl), and Trop = trop(linespace) × trop(pointspace) as the products of the
# faces of the raw fan of trop(linespace) with the cones of the tropical quadric (tropical_complete_hexagon.jl). Nothing
# else computed before is used. Min convention throughout.
#
# The certificate (README.md, written by compute_certificate.jl): the coordinates f, for each a list of representatives
# h_j = x^(-u_j) g_j, and for each maximal cone C of Trop and each f a tree. The root is C; a node σ split by h_j has as
# children the cells σ ∩ N_v, for v the vertices of the normal fan of Newt(h_j) on σ; the leaves are the cells of the
# cover, with f_σ = the h_j of their parent. The cells are computed here from the tree.
#
# Checked:
#   (a) ∪ σ = Trop. Every cone of Trop lies in a maximal one (Trop is pure), there is a tree for every maximal cone, and at
#       each split the children are exactly the vertices v of conv(exponents of h_j) + σ^∨ (in coordinates on the span of
#       σ). These are the v whose N_v meets σ in a full-dimensional cone; as the N_v cover R^28, the children cover σ, so
#       the leaves cover C.
#   (b) f_σ ∈ f + I: the normal form of x^u f - g modulo I is 0. The f of the certificate are those of the definitions.
#   (c) f_σ is tropically faithful on relint σ: σ = σ' ∩ N_v for its parent σ', so in_w(g) is the same for all
#       w ∈ relint σ, and at w = the sum of the rays of σ (a point of relint σ) in_w(g) ∉ J_C. Here
#       J_C = in_w(I) : (x_1 ⋯ x_28)^∞ for w ∈ relint C, computed from a reduced Gröbner basis G of I at a point of
#       relint C and certified to be the same on all of relint C: in_w(g) is constant on relint C for every g ∈ G.
#   (d) J is prime on every cone of Trop (all 6017, also those that are not maximal, for the faces of the cells): J has
#       a single minimal prime P over QQ, and P ⊆ J. The Gröbner cone of (c) is certified for every cone as well.
# With J_C prime over QQ, in_w(g) ∉ J_C means that in_w(g) vanishes on no geometric component of V(in_w I) in the torus:
# these components are Galois conjugate and in_w(g) has rational coefficients.
#
# The cones are split into parts (--part=i --parts=n); each part writes to verification/:
#   part_<i>_cones.csv:           per cone of Trop: face, quadric, dim, maximal, point_in_relint, in_maximal_cone (a),
#                                 groebner_cone, J_prime (c, d), trees (coordinates checked), cells, failures (a, c)
#   part_<i>_representatives.csv: per representative: coordinate, j, in_f_plus_I (b)
#   part_<i>_finished.csv:        written at the end of the part
#   part_<i>_failures.log:        what failed, if anything
# A rerun continues where a part stopped. `--summary` reads all parts and states whether everything is verified.
#
# From the `julia` folder (run_verification_overnight.sh runs all of it):
#     julia --project=. studies/hexagon/certificate/verify_certificate.jl --part=1 --parts=3     (and parts 2, 3)
#     julia --project=. studies/hexagon/certificate/verify_certificate.jl --summary
# `--limit=k` checks only k cones of each kind and k representatives (a timing test), `--out=dir` writes elsewhere.

include(joinpath(@__DIR__, "..", "tropical_complete_hexagon.jl"))

const VERIFY_CERTIFICATE_PATH = joinpath(@__DIR__, "hexagon_certificate.mrdi")
const VERIFICATION_DIR = joinpath(@__DIR__, "verification")

# --- The certificate ----------------------------------------------------------------------------------------

"The 56 coordinates of Lemma lem:complete_hexagon, as (name, polynomial)."
function verification_coordinates()
    H = hexagon_expressions()
    return [("$(name)_$c", H[name][c]) for name in (:x2, :x3, :p23, :p34) for c in eachindex(H[name])]
end

"A polynomial with the variables of `EXPRESSION_RING` as an element of it (each `load` creates a new ring)."
function to_expression_ring_verify(p)
    B = MPolyBuildCtx(EXPRESSION_RING)
    foreach((c, e) -> push_term!(B, QQ(c), e), AbstractAlgebra.coefficients(p), AbstractAlgebra.exponent_vectors(p))
    return finish(B)
end

"""
The certificate: for each coordinate (in the order of `verification_coordinates`) its representatives (g, u),
h_j = x^(-u) g, and the trees as a Dict (face, quadric) => the 56 trees. Fails unless its coordinates and their
polynomials f are those of the definitions.
"""
function load_certificate()
    d = load(VERIFY_CERTIFICATE_PATH)
    coordinates = verification_coordinates()
    collect(d["coordinate"]) == first.(coordinates) || error("the coordinates are not $(join(first.(coordinates), " "))")
    all(to_expression_ring_verify(f) == f0 for (f, (_, f0)) in zip(d["f"], coordinates)) ||
        error("the polynomials f are not those of polynomial_expressions.jl")
    reps = [[(to_expression_ring_verify(g), Int.(u)) for (g, u) in zip(gs, us)] for (gs, us) in zip(d["g"], d["u"])]
    trees = Dict((Int(k), String(S)) => collect(String, t) for (k, S, t) in zip(d["face"], d["quadric"], d["tree"]))
    return reps, trees
end

"""
A tree of the certificate, `node := "." | j "(" child ("," child)* ")"` with `child := t ":" node`: `nothing` for a leaf,
else (j, children) with children a vector of (t, node).
"""
function parse_tree(s)
    pos = Ref(1)
    peek() = pos[] <= ncodeunits(s) ? s[pos[]] : '\0'
    expect(c) = (peek() == c || error("tree: '$c' expected at $(pos[])"); pos[] += 1)
    function number()
        start = pos[]
        while isdigit(peek())
            pos[] += 1
        end
        pos[] > start || error("tree: number expected at $start")
        return parse(Int, s[start:pos[] - 1])
    end
    function node()
        peek() == '.' && (pos[] += 1; return nothing)
        j = number()
        expect('(')
        children = Tuple{Int,Any}[]
        while true
            t = number()
            expect(':')
            push!(children, (t, node()))
            peek() == ',' || break
            pos[] += 1
        end
        expect(')')
        return (j=j, children=children)
    end
    tree = node()
    pos[] == ncodeunits(s) + 1 || error("tree: unexpected text at $(pos[])")
    return tree
end

rows_matrix(vs) = isempty(vs) ? zero_matrix(QQ, 0, 28) : matrix(QQ, reduce(vcat, [permutedims(QQ.(v)) for v in vs]))

# --- The cones of Trop --------------------------------------------------------------------------------------

"The weights of the terms x1x4, x2x5, x3x6, x7^2 of q as linear forms on R^28."
function quadric_term_weight(t)
    v = zeros(Int, 28)
    t <= 3 ? (v[21 + t] = 1; v[24 + t] = 1) : (v[28] = 2)
    return v
end

"""
The closed cone of Trop of the face `k` of the raw fan and the set `S` of terms of q of minimal weight, as
{E w = 0, A w ≤ 0}, with a point of its relative interior (where A w < 0).
"""
function trop_cone(k, S)
    face = raw_fan_faces()[k]
    pad(M) = hcat(M, zero_matrix(QQ, nrows(M), 7))
    q = quadric_term_weight.(1:4)
    E = vcat(pad(face.E), rows_matrix([q[t] .- q[S[1]] for t in S[2:end]]))
    A = vcat(pad(face.A), rows_matrix([q[S[1]] .- q[t] for t in setdiff(1:4, S)]))
    C = cone_from_inequalities(A, E)
    return (cone=C, E=E, A=A, point=vcat(face.point, quadric_cone_point(S)), dim=dim(C))
end

in_closed_cone(c, w) = (v = matrix(QQ, 28, 1, QQ.(w)); all(iszero, c.E * v) && all(<=(0), c.A * v))
in_open_cone(c, w) = (v = matrix(QQ, 28, 1, QQ.(w)); all(iszero, c.E * v) && all(<(0), c.A * v))

"The rays modulo lineality and a basis of the lineality space of a cone, as integer vectors."
integer_rays(C) = ((R, L) = rays_modulo_lineality(C); ([integer_weight(r) for r in R], [integer_weight(l) for l in L]))

"in_w(g): the terms of g of minimal weight <e, w>."
function initial_form(g, w)
    exps = collect(AbstractAlgebra.exponent_vectors(g))
    coeffs = collect(AbstractAlgebra.coefficients(g))
    at = [sum(e .* w) for e in exps]
    B = MPolyBuildCtx(parent(g))
    foreach(i -> push_term!(B, coeffs[i], exps[i]), findall(==(minimum(at)), at))
    return finish(B)
end

"""
Whether in_w(g) = in_{w0}(g) for every g in G and every w in relint C, for w0 in relint C: the terms of minimal weight at
w0 have equal weight on the span of C, all others weight at least as large on the rays (hence larger on relint C).
Then G, a reduced Gröbner basis for the order refining w0, is one for the order refining w with the same tie-break, and
in_w(I) = <in_w(G)> is constant on relint C.
"""
function groebner_cone_certified(G, w0, C)
    rays_, lineality_ = integer_rays(C)
    for g in G
        exps = collect(AbstractAlgebra.exponent_vectors(g))
        at = [sum(e .* w0) for e in exps]
        e0 = exps[argmin(at)]
        for (e, a) in zip(exps, at)
            all(l -> sum((e .- e0) .* l) == 0, lineality_) || return false
            for r in rays_
                d = sum((e .- e0) .* r)
                (a == minimum(at) ? d == 0 : d >= 0) || return false
            end
        end
    end
    return true
end

"""
For a cone c of Trop: whether its point is in its relative interior, whether the Gröbner cone of its point is certified
on all of relint c, and J_c and whether it is prime.
"""
function check_trop_cone(I, c)
    S = base_ring(I)
    w = c.point
    G = without_gc(() -> collect(groebner_basis(I; ordering=min_weight_ordering(S, w), complete_reduction=true)))
    inG = ideal(S, [initial_form(g, w) for g in G])
    J = saturation(inG, ideal(S, [prod(gens(S))]))
    P = minimal_primes(J)
    prime = length(P) == 1 && all(p -> p in J, gens(only(P)))
    return (point_in_relint=in_open_cone(c, w), groebner_cone=groebner_cone_certified(G, w, c.cone), J=J, J_prime=prime)
end

# --- The trees ----------------------------------------------------------------------------------------------

"Independent vectors among `vs`, a basis of their span."
function span_basis(vs)
    B = Vector{Int}[]
    for v in vs
        rank(rows_matrix(vcat(B, [v]))) == length(B) + 1 && push!(B, v)
    end
    return B
end

"""
The checks (a) and (c) below the node `nd` = (j, children) on the cell σ (with its integer rays and lineality): the
children must be the vertices v of conv(forms) + σ^∨, where in coordinates c on a basis B of the span of σ the exponent
e of g_j is the linear form (<e, b>)_b; a vertex is named by the position of the smallest exponent with its form in
`sorted[j]`, the sorted exponents of g_j. The cell of v is σ ∩ N_v = σ ∩ {<v - u, w> ≤ 0 for all vertices u}: on σ the
minimum over all exponents is attained at a vertex. A leaf is tested for faithfulness at the sum of its rays.
Returns the number of leaves; failures are pushed to `failures`.
"""
function check_node(σ, rays_, lineality_, nd, J, reps, sorted, failures, where)
    1 <= nd.j <= length(reps) || (push!(failures, "$where: no representative $(nd.j)"); return 0)
    g = reps[nd.j][1]
    exps = sorted[nd.j]
    B = span_basis(vcat(rays_, lineality_))
    groups = Dict{Vector{QQFieldElem},Int}()   # form => position of its smallest exponent
    for (i, e) in enumerate(exps)
        get!(groups, QQ.([sum(e .* b) for b in B]), i)
    end
    M = rows_matrix(B)
    coordinates(v) = collect(solve(M, QQ.(v); side=:left))
    coordinate_rows(vs) = isempty(vs) ? zero_matrix(QQ, 0, length(B)) : matrix(QQ, reduce(vcat, [permutedims(coordinates(v)) for v in vs]))
    dual_rays = [QQFieldElem.(collect(r)) for r in rays(polarize(positive_hull(coordinate_rows(rays_), coordinate_rows(lineality_))))]
    forms = collect(keys(groups))
    newton = convex_hull(matrix(QQ, reduce(vcat, permutedims.(forms))),
                         isempty(dual_rays) ? zero_matrix(QQ, 0, length(B)) : matrix(QQ, reduce(vcat, permutedims.(dual_rays))))
    V = [QQFieldElem.(collect(v)) for v in vertices(newton)]
    vertex = Dict(groups[v] => v for v in V)
    sort(first.(nd.children)) == sort(collect(keys(vertex))) ||
        (push!(failures, "$where: the children $(sort(first.(nd.children))) are not the vertices $(sort(collect(keys(vertex))))"); return 0)
    Aσ = linear_inequality_matrix(facets(σ))
    Eσ = linear_equation_matrix(linear_span(σ))
    cells = 0
    for (t, child) in nd.children
        v = vertex[t]
        cell = cone_from_inequalities(vcat(Aσ, rows_matrix([exps[t] .- exps[groups[u]] for u in V if u != v])), Eσ)
        dim(cell) == dim(σ) || (push!(failures, "$where, child $t: not full-dimensional"); continue)
        cell_rays, cell_lineality = integer_rays(cell)
        if child === nothing
            initial_form(g, sum(cell_rays)) in J && push!(failures, "$where, child $t: not faithful")
            cells += 1
        else
            cells += check_node(cell, cell_rays, cell_lineality, child, J, reps, sorted, failures, "$where, $t")
        end
    end
    return cells
end

# --- Main ---------------------------------------------------------------------------------------------------

function append_csv(path, row)
    isfile(path) || write(path, join(string.(keys(row)), ",") * "\n")
    open(io -> println(io, join(string.(values(row)), ",")), path, "a")
end

done_keys(path, columns) = isfile(path) ? Set(join(split(l, ",")[columns], ",") for l in Iterators.drop(eachline(path), 1)) : Set{String}()

function verify_part(part, parts; dir, limit, log)
    mkpath(dir)
    I = line_point_ideal()
    S = base_ring(I)
    coordinates = verification_coordinates()
    names = first.(coordinates)
    reps, trees = load_certificate()
    sorted = [[sort(collect(AbstractAlgebra.exponent_vectors(g))) for (g, _) in r] for r in reps]
    failures_path = joinpath(dir, "part_$(part)_failures.log")
    report(s) = (open(io -> println(io, s), failures_path, "a"); log("FAILED: $s"))

    # (b)
    reps_path = joinpath(dir, "part_$(part)_representatives.csv")
    done = done_keys(reps_path, 1:2)
    pairs = [(i, j) for i in eachindex(coordinates) for j in eachindex(reps[i])]
    mine = [p for (n, p) in enumerate(pairs) if (n - 1) % parts == part - 1 && !("$(names[p[1]]),$(p[2])" in done)]
    limit === nothing || (mine = first(mine, limit))
    log("part $part: $(length(mine)) representatives to check")
    started = time()
    for (i, j) in mine
        g, u = reps[i][j]
        ok = iszero(normal_form(prod(gens(S) .^ u) * coordinates[i][2] - g, I))
        ok || report("representative $j of $(names[i]) is not in f + I")
        append_csv(reps_path, (coordinate=names[i], j=j, in_f_plus_I=ok))
    end
    log("part $part: representatives done in $(round(time() - started; digits=1)) s")

    # (a), (c), (d)
    all_cones = [(k, T) for k in eachindex(raw_fan_faces()) for T in [T for d in 2:4 for T in subsets(collect(1:4), d)]]
    data = Dict((k, join(T)) => trop_cone(k, T) for (k, T) in all_cones)
    top = maximum(c.dim for c in values(data))
    maximal = sort([key for (key, c) in data if c.dim == top])
    index = Dict((k, join(T)) => n for (n, (k, T)) in enumerate(all_cones))
    is_mine(key) = (index[key] - 1) % parts == part - 1
    cones_path = joinpath(dir, "part_$(part)_cones.csv")
    done = done_keys(cones_path, 1:2)
    todo_count = count(key -> is_mine(key) && !("$(key[1]),$(key[2])" in done), keys(data))
    log("part $part: $(length(all_cones)) cones, $(length(maximal)) maximal (dimension $top), $todo_count to check")
    function cone_check(key)
        c = data[key]
        t = check_trop_cone(I, c)
        t.point_in_relint || report("cone $key: its point is not in its relative interior")
        t.groebner_cone || report("cone $key: the Gröbner cone is not certified")
        t.J_prime || report("cone $key: J is not prime")
        in_maximal = c.dim == top || any(m -> in_closed_cone(data[m], c.point), maximal)
        in_maximal || report("cone $key: in no maximal cone")
        return (; t..., in_maximal)
    end
    cone_row(key, t, trees, cells, nfailures) =
        (face=key[1], quadric=key[2], dim=data[key].dim, maximal=data[key].dim == top, point_in_relint=t.point_in_relint,
         in_maximal_cone=t.in_maximal, groebner_cone=t.groebner_cone, J_prime=t.J_prime, trees=trees, cells=cells,
         failures=nfailures)
    checked = 0
    # the cones that are not maximal: (a) and (d) only
    for (k, T) in all_cones
        key = (k, join(T))
        (data[key].dim == top || !is_mine(key) || "$k,$(join(T))" in done) && continue
        limit !== nothing && checked >= limit && break
        append_csv(cones_path, cone_row(key, cone_check(key), 0, 0, 0))
        checked += 1
        checked % 100 == 0 && (log("part $part: $checked cones"); GC.gc(true))
    end
    # the maximal cones, with the trees of all coordinates
    checked_maximal = 0
    for key in maximal
        (!is_mine(key) || "$(key[1]),$(key[2])" in done) && continue
        limit !== nothing && checked_maximal >= limit && break
        t0 = time()
        t = cone_check(key)
        failures = String[]
        cone_trees = get(trees, key, String[])
        length(cone_trees) == length(names) || push!(failures, "$(length(cone_trees)) trees instead of $(length(names))")
        C = data[key].cone
        C_rays, C_lineality = integer_rays(C)
        cells = 0
        for (i, s) in enumerate(cone_trees)
            i > length(names) && break
            tree = try
                parse_tree(s)
            catch err
                push!(failures, "$(names[i]): $(sprint(showerror, err))")
                continue
            end
            tree === nothing && (push!(failures, "$(names[i]): the root is a leaf"); continue)
            cells += check_node(C, C_rays, C_lineality, tree, t.J, reps[i], sorted[i], failures, names[i])
        end
        foreach(s -> report("cone $key: $s"), failures)
        append_csv(cones_path, cone_row(key, t, length(cone_trees), cells, length(failures)))
        checked_maximal += 1
        log("part $part: maximal cone $key: $cells cells, $(length(failures)) failures, $(round(time() - t0; digits=1)) s")
        checked_maximal % 20 == 0 && GC.gc(true)
    end
    if limit === nothing
        leftover = !issubset(keys(trees), Set(maximal))
        leftover && report("trees of cones that are not maximal cones of Trop")
        append_csv(joinpath(dir, "part_$(part)_finished.csv"), (part=part, parts=parts, leftover_rows=leftover))
    end
    log("part $part finished")
end

function verify_summary(; dir, log)
    files(suffix) = filter(f -> startswith(f, "part_") && endswith(f, suffix), readdir(dir))
    rows(file) = [split(l, ",") for l in Iterators.drop(eachline(joinpath(dir, file)), 1)]
    cone_rows = reduce(vcat, rows.(files("_cones.csv")); init=Vector{SubString{String}}[])
    rep_rows = reduce(vcat, rows.(files("_representatives.csv")); init=Vector{SubString{String}}[])
    finished = reduce(vcat, rows.(files("_finished.csv")); init=Vector{SubString{String}}[])
    n_coordinates = length(verification_coordinates())
    n_reps = sum(length, first(load_certificate()))
    problems = String[]
    keys_ = [(r[1], r[2]) for r in cone_rows]
    n_cones = 11 * length(raw_fan_faces())
    length(keys_) == length(Set(keys_)) == n_cones || push!(problems, "$(length(Set(keys_))) cones checked (of $n_cones), $(length(keys_)) rows")
    for r in cone_rows
        all(==("true"), r[5:8]) || push!(problems, "cone $(r[1]) $(r[2]): $(join(r[5:8], " "))")
        r[4] == "true" && !(r[9] == string(n_coordinates) && r[11] == "0") && push!(problems, "cone $(r[1]) $(r[2]): trees $(r[9]), failures $(r[11])")
    end
    rep_keys = [(r[1], r[2]) for r in rep_rows]
    length(rep_keys) == length(Set(rep_keys)) == n_reps || push!(problems, "$(length(Set(rep_keys))) representatives checked (of $n_reps)")
    foreach(r -> r[3] == "true" || push!(problems, "representative $(r[2]) of $(r[1]) not in f + I"), rep_rows)
    parts = isempty(finished) ? 0 : parse(Int, finished[1][2])
    (parts > 0 && sort(unique(parse.(Int, first.(finished)))) == 1:parts && all(r -> r[3] == "false", finished)) ||
        push!(problems, "not all parts finished, or trees of cones that are not maximal")
    maximal = count(r -> r[4] == "true", cone_rows)
    log("$(length(keys_)) cones ($maximal maximal), $(sum(parse(Int, r[10]) for r in cone_rows; init=0)) cells, $(length(rep_keys)) representatives")
    isempty(problems) ? log("VERIFIED: (a) the cells cover Trop, (b) f_σ ∈ f + I, (c) f_σ is faithful on relint σ, (d) J is prime on all cones") :
                        foreach(p -> log("PROBLEM: $p"), problems)
end

function verify_main(args)
    option(name) = (a = findfirst(startswith("--$name="), args); a === nothing ? nothing : split(args[a], "=")[2])
    log(s) = (println("[$(now())] $s"); flush(stdout))
    dir = option("out") === nothing ? VERIFICATION_DIR : option("out")
    if "--summary" in args
        verify_summary(; dir, log)
    else
        part = option("part") === nothing ? 1 : parse(Int, option("part"))
        parts = option("parts") === nothing ? 1 : parse(Int, option("parts"))
        limit = option("limit") === nothing ? nothing : parse(Int, option("limit"))
        verify_part(part, parts; dir, limit, log)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    verify_main(ARGS)
end
