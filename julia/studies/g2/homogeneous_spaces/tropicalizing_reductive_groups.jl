#=
Theorem thm:modification of main.tex (Section 8).

Compares the Pluecker tropicalization trop(G) in R^21/R1 with the tropicalization in the weight basis,
obtained by forgetting the coordinate p_36, using the census data of the G2-Grassmannian.  Sections 1-5
compare the two, with exact linear algebra over QQ.
Section 6 computes, for every non-vertical maximal cell, the initial ideal of the weight-basis ideal
at an interior weight, and tests it for being toric with the criterion of canonical_initial_ideal in
julia/src/degenerations.jl.

Usage, from the repository root:
    julia --project=julia julia/studies/g2/homogeneous_spaces/tropicalizing_reductive_groups.jl \
        julia/data/g2/grassmannian > julia/studies/g2/homogeneous_spaces/tropicalizing_reductive_groups.txt
Reads:  00_ideal.mrdi, 04_trees.csv, 02_coarsest_cells.csv, 06_canonical_interior_degenerations.csv,
        and, for the paper's labels, julia/studies/g2/atlas/flag/local_flag_cell_atlas.csv
        (cell types) and julia/studies/g2/atlas/grassmannian/ray_label_orbits.csv (ray orbits).
Census orbits are printed as paper cell types, e.g. type 1 is census orbit 29.
This is a standalone script; it is not part of the TropicalCensus module.
=#

using Oscar

# ---------------------------------------------------------------- CSV and printing

"Split one CSV line into fields, honoring double quotes."
function csv_fields(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    inquote = false
    chars = collect(line)
    i = 1
    while i <= length(chars)
        c = chars[i]
        if inquote
            if c == '"'
                if i < length(chars) && chars[i + 1] == '"'
                    write(buf, '"')
                    i += 1
                else
                    inquote = false
                end
            else
                write(buf, c)
            end
        elseif c == '"'
            inquote = true
        elseif c == ','
            push!(fields, String(take!(buf)))
        else
            write(buf, c)
        end
        i += 1
    end
    push!(fields, String(take!(buf)))
    return fields
end

"Read a CSV file with a header line into a vector of rows."
function read_csv(path::AbstractString)
    lines = filter(!isempty, readlines(path))
    header = csv_fields(lines[1])
    return [Dict(zip(header, csv_fields(l))) for l in lines[2:end]]
end

pyrepr(s::AbstractString) = "'" * s * "'"
pyrepr(b::Bool) = b ? "True" : "False"
pyrepr(x::Integer) = string(x)
pyrepr(v::AbstractVector) = "[" * join(pyrepr.(v), ", ") * "]"
pyrepr(ps::AbstractVector{<:Pair}) = "{" * join(["$(pyrepr(k)): $(pyrepr(v))" for (k, v) in ps], ", ") * "}"

"Count the elements of a collection, keeping the order of first occurrence (like Python's Counter)."
function ordered_count(xs)
    keys_ = Any[]
    counts = Dict{Any, Int}()
    for x in xs
        haskey(counts, x) || push!(keys_, x)
        counts[x] = get(counts, x, 0) + 1
    end
    return [k => counts[k] for k in keys_]
end

# ---------------------------------------------------------------- exact linear algebra with Oscar

rowmat(vs) = matrix(QQ, reduce(vcat, [permutedims(v) for v in vs]))
colmat(vs) = transpose(rowmat(vs))

"The row span of the given vectors, as the nonzero rows of their reduced row echelon form."
function span_key(vectors, coords=nothing)
    A = rowmat(coords === nothing ? vectors : [v[coords] for v in vectors])
    r, R = rref(A)
    return [[R[i, j] for j in 1:ncols(R)] for i in 1:r]
end

"A solution of M x = b over QQ, or nothing if there is none."
function solve_q(M, b)
    ok, x = can_solve_with_solution(M, matrix(QQ, length(b), 1, b); side=:right)
    return ok ? [x[i, 1] for i in 1:nrows(x)] : nothing
end

# ---------------------------------------------------------------- data

const D = ARGS[1]
const PAIRS = [(i, j) for i in 1:6 for j in (i + 1):7]
const K36 = findfirst(==((3, 6)), PAIRS)
const KEEP = [i for i in 1:21 if i != K36]
pidx(ij) = findfirst(==(ij), PAIRS)

function tree_metric(row)
    s = parse.(Int, split(row["pendant_lengths"], ';'))
    splits = Tuple{Set{Int}, Int}[]
    for t in filter(!isempty, split(row["internal_splits"], ';'))
        a, l = split(t, ':')
        A = Set(parse(Int, string(c)) for c in split(a, '|')[1])
        push!(splits, (A, parse(Int, l)))
    end
    return QQFieldElem[s[i] + s[j] + sum((l for (A, l) in splits if (i in A) != (j in A)); init=0)
                       for (i, j) in PAIRS]
end

trees = read_csv(joinpath(D, "04_trees.csv"))
ray_order = [r["index"] for r in trees if r["kind"] == "ray"]
rays = Dict(r["index"] => tree_metric(r) for r in trees if r["kind"] == "ray")
cell_point = Dict(r["index"] => tree_metric(r) for r in trees if r["kind"] == "cell")
coarsest = read_csv(joinpath(D, "02_coarsest_cells.csv"))
cell_order = [c["coarsest_cell"] for c in coarsest]
cell_rays = Dict(c["coarsest_cell"] => String.(split(c["extremal_original_rays"], ';')) for c in coarsest)
ATLAS = joinpath(D, "..", "..", "..", "studies", "g2", "atlas")
cell_type = Dict(r["qq_orbit"] => r["label"]
                 for r in read_csv(joinpath(ATLAS, "flag", "local_flag_cell_atlas.csv")) if !isempty(r["qq_orbit"]))
# article labels of the ray orbits, as in julia/to_tex/build_styled_ray_orbit_tables.jl
RAY_LABELS = ["beta_1", "beta_2", "alpha_1", "alpha_3", "alpha_2", "alpha_4", "gamma_1", "gamma_2",
              "alpha_5", "gamma_3", "delta"]
ray_type = Dict{String, String}()
for r in read_csv(joinpath(ATLAS, "grassmannian", "ray_label_orbits.csv"))
    for i in split(r["rays"], ';')
        ray_type[String(i)] = RAY_LABELS[parse(Int, r["orbit"])]
    end
end
toric = Dict{String, Bool}()
orbit = Dict{String, String}()
for r in read_csv(joinpath(D, "06_canonical_interior_degenerations.csv"))
    for c in split(r["cells"], ';')
        toric[String(c)] = r["toric"] == "true"
        orbit[String(c)] = cell_type[r["orbit"]]
    end
end

# lineality: 1 and the two star trees of eq:lineality
STARS = [[1, -1, 0, -1, 1, 0, 0], [-1, 0, 1, 1, 0, -1, 0]]
LIN = [fill(QQ(1), 21), [QQFieldElem[l[i] + l[j] for (i, j) in PAIRS] for l in STARS]...]
E36 = zeros(QQ, 21)
E36[K36] = 1

println("== 1. the vertical ray -e_36")
M = colmat([LIN; [E36]])
for idx in ray_order
    x = solve_q(M, rays[idx])
    x === nothing && continue
    c = Int(x[4])
    println("ray $idx (orbit $(ray_type[idx])) = $(c >= 0 ? "+" : "")$c e_36 mod lineality")
end
VR = "30"

println("== 2. cells whose span contains e_ij, for ij in 14, 25, 36")
for ij in [(1, 4), (2, 5), (3, 6)]
    e = zeros(QQ, 21)
    e[pidx(ij)] = 1
    n = count(cell_order) do c
        B = [[rays[r] for r in cell_rays[c]]; LIN]
        rank(rowmat([B; [e]])) == rank(rowmat(B))
    end
    println("($(ij[1]), $(ij[2])) $n of $(length(cell_order))")
end
vertical = [c for c in cell_order if rank(rowmat([[rays[r] for r in cell_rays[c]]; LIN; [E36]])) == 6]
println("cells containing e_36 in their span: $(length(vertical))  all contain ray $VR : ",
        pyrepr(all(VR in cell_rays[c] for c in vertical)),
        "  cells containing ray $VR : $(count(c -> VR in cell_rays[c], cell_order))")

println("== 3. classification of maximal cells by (d14, d25, d36), with toric status")
cat = Dict{String, String}()
for c in cell_order
    d = cell_point[c]
    a, b, e = d[pidx((1, 4))], d[pidx((2, 5))], d[K36]
    cat[c] = a == b == e ? "d14=d25=d36" :
             e == max(a, b) ? "d36=max(d14,d25)" :
             c in vertical ? "d36<d14=d25, vertical" : "d36<d14=d25, not vertical"
end
cnt = ordered_count((cat[c], toric[c]) for c in cell_order)
for (k, n) in sort(cnt; by=first)
    println(rpad(k[1], 28), " toric=", rpad(pyrepr(k[2]), 5), " ", n)
end
println("paper cell types (with multiplicity) by class:")
for k in sort(unique(values(cat)))
    types = sort(ordered_count(orbit[c] for c in cell_order if cat[c] == k); by=t -> parse(Int, t.first))
    println("  ", k, " ", pyrepr(types))
end

println("== 4. images of the non-vertical cells after forgetting p_36")
group_keys = Any[]
groups = Dict{Any, Vector{String}}()
for c in cell_order
    c in vertical && continue
    k = span_key([[rays[r] for r in cell_rays[c]]; LIN], KEEP)
    haskey(groups, k) || push!(group_keys, k)
    push!(get!(groups, k, String[]), c)
end
pairs = [groups[k] for k in group_keys if length(groups[k]) > 1]
println("non-vertical cells $(sum(length(g) for g in values(groups)))  distinct spans $(length(groups))",
        "  pairs with equal span $(length(pairs))  groups larger than 2: ",
        pyrepr([g for g in pairs if length(g) > 2]))
keep_lin = [l[KEEP] for l in LIN]
n_convex = 0
for g in pairs
    A, B = cell_rays[g[1]], cell_rays[g[2]]
    sh = [r for r in A if r in B]
    c1 = [r for r in A if !(r in B)]
    c2 = [r for r in B if !(r in A)]
    # c2 = a*sh0 + b*sh1 + z*c1 modulo lineality.  The cells lie on opposite sides iff z < 0, and then
    # their union is convex iff the segment from c1 to c2 meets cone(sh0, sh1), that is, a, b >= 0.
    a, b, z = solve_q(colmat([[rays[sh[1]][KEEP], rays[sh[2]][KEEP], rays[c1[1]][KEEP]]; keep_lin]),
                      rays[c2[1]][KEEP])[1:3]
    convex = z < 0 && a >= 0 && b >= 0
    global n_convex += convex
    println("  cells ", lpad(g[1], 3), ",", lpad(g[2], 3), " (types $(orbit[g[1]]),$(orbit[g[2]]), toric ",
            pyrepr(toric[g[1]]), ",", pyrepr(toric[g[2]]), ", $(cat[g[1]])): share $(length(sh)) rays, ",
            z < 0 ? "opposite sides" : "SAME SIDE", ", ", convex ? "convex union" : "NOT CONVEX")
end
println("pairs with convex union: $n_convex of $(length(pairs))")
both_toric = [g for g in pairs if toric[g[1]] && toric[g[2]]]
println("pairs with both halves toric: $(length(both_toric))  classes: ",
        pyrepr(sort(unique([cat[c] for g in both_toric for c in g]))))

println("== 5. the link of the vertical ray")
adj_keys = String[]
adj = Dict{String, Set{String}}()
for c in vertical
    a, b = [r for r in cell_rays[c] if r != VR]
    for (x, y) in ((a, b), (b, a))
        haskey(adj, x) || push!(adj_keys, x)
        push!(get!(adj, x, Set{String}()), y)
    end
end
seen = Set{String}()
comps = 0
for v in adj_keys
    v in seen && continue
    global comps += 1
    stack = [v]
    while !isempty(stack)
        x = pop!(stack)
        x in seen && continue
        push!(seen, x)
        append!(stack, adj[x])
    end
end
nv, ne = length(adj), length(vertical)
println("vertices $nv, edges $ne, components $comps, first Betti number $(ne - nv + comps)")
println("vertex degrees ", pyrepr(ordered_count(length(adj[v]) for v in adj_keys)))
println("vertices (census ray indices) ", pyrepr(sort(adj_keys; by=x -> parse(Int, x))))
println("their ray orbits ", pyrepr(ordered_count(ray_type[v] for v in adj_keys)))
println("distinct images of these rays after forgetting p_36: ",
        length(unique(span_key([[rays[r]]; LIN], KEEP) for r in adj_keys)))

# ---------------------------------------------------------------- 6. initial ideals in the weight basis

"""
Toric test of canonical_initial_ideal in julia/src/degenerations.jl: the initial ideal is toric if
it is prime, has a binomial reduced Groebner basis, contains no monomial, and the exponent-difference
lattice of that basis is saturated.  The binomial test runs first, since it is the cheapest.
"""
function is_toric_initial(J)
    reduced = gens(groebner_basis(J; complete_reduction=true))
    all(f -> length(f) <= 2, reduced) || return false
    S = base_ring(J)
    isone(saturation(J, ideal(S, prod(gens(S))))) && return false
    is_prime(J) || return false
    differences = [e[1] - e[2] for e in (collect(exponents(f)) for f in reduced) if length(e) == 2]
    isempty(differences) && return true
    L = snf(matrix(ZZ, reduce(vcat, permutedims.(differences))))
    return prod((abs(L[i, i]) for i in 1:min(size(L)...) if !iszero(L[i, i])); init=ZZ(1)) == 1
end

println("== 6. initial ideals in the weight basis, at an interior point of each non-vertical cell")
I = load(joinpath(D, "00_ideal.mrdi"))
R = base_ring(I)
p = gens(R)
names_w = [string(x) for (k, x) in enumerate(p) if k != K36]
Rw, pw = polynomial_ring(QQ, names_w)
images = Vector{elem_type(Rw)}(undef, 21)
for (k, idx) in enumerate(KEEP)
    images[idx] = pw[k]
end
images[K36] = -pw[findfirst(==(pidx((1, 4))), KEEP)] - pw[findfirst(==(pidx((2, 5))), KEEP)]
phi = hom(R, Rw, images)
Iw = ideal(Rw, filter(!iszero, phi.(gens(I))))
# Apart from the 6 linear forms cutting out L in the weight basis, the ideal has no linear forms.
println("linear forms in the weight-basis ideal (6 expected, those cutting out L): ",
        count(f -> total_degree(f) == 1, gens(groebner_basis(Iw; complete_reduction=true))))
nu = tropical_semiring_map(QQ, min)
initial_w = Dict{String, Any}()
toric_w = Dict{String, Bool}()
t0 = time()
for c in cell_order
    c in vertical && continue
    # max-convention tree metric d; Oscar's min convention uses the weight -d
    w = [-cell_point[c][k] for k in KEEP]
    J = initial(Iw, nu, w)
    initial_w[c] = J
    toric_w[c] = is_toric_initial(J)
end
println("computed $(length(toric_w)) initial ideals in $(round(Int, time() - t0)) s")
for (k, n) in sort(ordered_count((cat[c], toric[c], toric_w[c]) for c in cell_order if !(c in vertical)); by=first)
    println(rpad(k[1], 28), " toric=", rpad(pyrepr(k[2]), 5), " weight-basis toric=", rpad(pyrepr(k[3]), 5), " ", n)
end
mismatch = [c for c in keys(toric_w) if toric_w[c] != (toric[c] && cat[c] != "d36<d14=d25, not vertical")]
println("cells where the weight-basis toric status differs from the prediction of thm:modification: ",
        pyrepr(sort(mismatch; by=x -> parse(Int, x))))
println("toric cells in the weight basis: $(count(values(toric_w)))")
println("merged pairs whose halves have equal weight-basis initial ideals: ",
        count(g -> initial_w[g[1]] == initial_w[g[2]], pairs), " of $(length(pairs))")
distinct_toric = Any[]
for c in cell_order
    get(toric_w, c, false) || continue
    any(J -> J == initial_w[c], distinct_toric) || push!(distinct_toric, initial_w[c])
end
println("distinct toric initial ideals in the weight basis: $(length(distinct_toric))")
