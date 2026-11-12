###############################################################################
# Independent check of the Christophersen--Ilten paragraph of Section 6
# (CI, Math. Z. 278 (2014), Theorem 3.1, case g = 10, and Figures 1 and 3).
#
# 1. CI's 8 term orders, one for each resolution of their three comparisons,
#    give 8 squarefree quadratic initial ideals, each the Stanley--Reisner
#    ideal of K * Delta_2 with Delta_2 = {u, v, r} and K a 2-sphere on 11
#    vertices; K is T_11, built by CI's starring rule, for exactly two of them.
# 2. A cell of trop is below an order when its canonical weight lies in the
#    closed Groebner cone: in_<(in_w I) = in_<(I), tested in our 21 variables
#    with the lift of CI's order read off their matrix.
# 3. For toric cells this is tested again without Groebner bases: the regular
#    subdivision of the point configuration of in_w(I) under CI's weight is
#    a unimodular triangulation with the facets of K * Delta_2 (Sturmfels,
#    Theorem 8.3). Both tests must agree on every toric cell.
#
# Run with julia --project=julia from any directory; writes verify.txt next
# to this file. Characteristic 0 only.
###############################################################################

using TropicalCensus
using Oscar

const DATA = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "data", "g2", "grassmannian"))
const TORIC_TABLE = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "..", "tables",
                                      "toric_degeneration_table.tex"))
const OUTPUT = open(joinpath(@__DIR__, "verify.txt"), "w")
say(parts...) = (println(parts...); println(OUTPUT, parts...))

pfaffian(M, i, j, k, l) = M[i, j] * M[k, l] - M[i, k] * M[j, l] + M[i, l] * M[j, k]
const QUADRUPLES = [(i, j, k, l) for i in 1:7 for j in i+1:7 for k in j+1:7 for l in k+1:7]

# ------------------------------------------------------------------ 1. CI's side
# The matrix of the case g = 10, copied from the published paper.
const NAMES = ["u", "v", "w", "r", "x00", "x01", "x10", "x11", "y00", "y01", "y10", "y11", "z0", "z1"]
S, V = polynomial_ring(QQ, NAMES)
var(name) = V[findfirst(==(name), NAMES)]
(u, v, w, r, x00, x01, x10, x11, y00, y01, y10, y11, z0, z1) = V
A = S[ 0    -x10  x11  w    y11  y10   u;
       x10   0   -v    y00  r    z0    x00;
      -x11   v    0    y01  z1  -w-r   x01;
      -w    -y00 -y01  0    x01 -x00   v;
      -y11  -r   -z1  -x01  0    u     x11;
      -y10  -z0   w+r  x00 -u    0     x10;
      -u    -x00 -x01 -v   -x11 -x10   0 ]
@assert all(A[i, j] == -A[j, i] for i in 1:7, j in 1:7)
const PFAFFIANS = [pfaffian(A, q...) for q in QUADRUPLES]
J = ideal(S, PFAFFIANS)
@assert dim(J) == 6 && degree(J) == 18

"""
CI's conditions u, v < r, w, x_ij < y_ij, z_i and middle * middle > right * left,
with the three open comparisons resolved by the signs s: s_1 for x00 x11 vs
x01 x10, s_2 for x00 y01 vs x01 y00, s_3 for x10 y11 vs x11 y10 (+1: left wins).
"""
function ci_weight(s)
    base = Dict("u" => 0, "v" => 0, "w" => 10, "r" => 10, "x00" => 10, "x01" => 10, "x10" => 10,
                "x11" => 10, "y00" => 15, "y01" => 15, "y10" => 15, "y11" => 15, "z0" => 15, "z1" => 15)
    base["x00"] += s[1]; base["y01"] += 2s[2]; base["y11"] += s[3]
    return [base[name] for name in NAMES]
end
const SIGNS = [(s1, s2, s3) for s1 in (1, -1) for s2 in (1, -1) for s3 in (1, -1)]
for s in SIGNS
    ω = ci_weight(s)
    group(name) = name in ("u", "v") ? 1 : (name[1] in "yz" ? 3 : 2)
    @assert all(ω[a] < ω[b] for a in 1:14, b in 1:14 if group(NAMES[a]) < group(NAMES[b]))
    @assert all(ω[a] + ω[b] > ω[c] + ω[d] for a in 1:14, b in 1:14, c in 1:14, d in 1:14
                if group(NAMES[a]) == group(NAMES[b]) == 2 && group(NAMES[c]) == 3 && group(NAMES[d]) == 1)
    lhs = ((ω[5] + ω[8], ω[6] + ω[7]), (ω[5] + ω[10], ω[6] + ω[9]), (ω[7] + ω[12], ω[8] + ω[11]))
    @assert all(sign(a - b) == si for ((a, b), si) in zip(lhs, s))
end

"The squarefree supports of the minimal generators of a monomial ideal."
supports(M) = Set(Tuple(findall(!iszero, first(exponents(g)))) for g in gens(M))

# ---- triangulated 2-spheres as sets of vertex triples; T_n by CI's starring rule
edges_of(T) = unique(sort([minmax(a, b) for t in T for (a, b) in ((t[1], t[2]), (t[1], t[3]), (t[2], t[3]))]))
vertices_of(T) = sort(unique(vcat(collect.(T)...)))
valency(T, x) = count(e -> x in e, edges_of(T))
edge_link(T, (a, b)) = [only(setdiff(collect(t), [a, b])) for t in T if a in t && b in t]
function star_edge(T, (a, b))
    new = maximum(vertices_of(T)) + 1
    starred = [t for t in T if !(a in t && b in t)]
    for c in edge_link(T, (a, b)), x in (a, b)
        push!(starred, Tuple(sort([x, c, new])))
    end
    return starred
end
function edge_graph(T)
    labels = vertices_of(T)
    graph = Graph{Undirected}(length(labels))
    for (a, b) in edges_of(T)
        add_edge!(graph, findfirst(==(a), labels), findfirst(==(b), labels))
    end
    return graph
end
"Star a vertex into an edge chosen by `rule`; all choices must give one sphere up to isomorphism."
function star_by_rule(T, rule)
    choices = [star_edge(T, e) for e in edges_of(T) if rule(map(x -> valency(T, x), edge_link(T, e)))]
    @assert !isempty(choices) && all(is_isomorphic(edge_graph(first(choices)), edge_graph(c)) for c in choices)
    return first(choices)
end
# Read literally, CI's rule fails at T_5 -> T_6: starring an edge of the bipyramid
# with link of valency 4 gives valencies 3 3 4 4 5 5, not the octahedron they
# call T_6. Start from the octahedron; from there on the rule is as stated.
T = [(a, c, e) for a in (1, 2) for c in (3, 4) for e in (5, 6)]                  # T_6
for n in 7:10
    global T = star_by_rule(T, link -> count(==(4), link) == 2)
    @assert length(vertices_of(T)) == n && all(valency(T, x) in (4, 5) for x in vertices_of(T))
end
const T11 = star_by_rule(T, link -> count(==(4), link) == 1)
"Cliques of size 3 and 4 in the edge graph of a sphere."
cliques(T, k) = [c for c in Iterators.product(ntuple(_ -> vertices_of(T), k)...)
                 if issorted(c) && allunique(c) && all(minmax(a, b) in edges_of(T) for a in c for b in c if a < b)]
@assert length(T11) == 18 && length(cliques(T11, 3)) == 18 && isempty(cliques(T11, 4))   # T_11 is flag
say("T_11 from CI's rule: 11 vertices, 27 edges, 18 triangles, flag, valencies ",
    join(sort([valency(T11, x) for x in vertices_of(T11)]), " "))

"The flag complex of a squarefree quadratic monomial ideal: cone points and the 2-sphere K."
function flag_complex(nonedges)
    used = sort(unique(vcat(collect.(collect(nonedges))...)))
    edges = [(a, b) for a in used for b in used if a < b && !((a, b) in nonedges)]
    triangles = [(a, b, c) for (a, b) in edges for c in used if c > b &&
                 (a, c) in edges && (b, c) in edges]
    return (cone=setdiff(1:14, used), sphere=triangles, edges=edges)
end
"Is a list of edges one cycle?"
function is_cycle(edges)
    ends = unique(vcat(collect.(edges)...))
    length(edges) == length(ends) && all(y -> count(e -> y in e, edges) == 2, ends) || return false
    seen, frontier = Set([first(ends)]), [first(ends)]
    while !isempty(frontier)
        y = pop!(frontier)
        for e in edges, z in e
            y in e && !(z in seen) && (push!(seen, z); push!(frontier, z))
        end
    end
    return length(seen) == length(ends)
end
"A closed 2-manifold with Euler characteristic 2."
is_two_sphere(K) = all(count(t -> a in t && b in t, K) == 2 for (a, b) in edges_of(K)) &&
    all(x -> is_cycle([Tuple(setdiff(collect(t), [x])) for t in K if x in t]), vertices_of(K)) &&
    length(vertices_of(K)) - length(edges_of(K)) + length(K) == 2

CI_IDEALS = []
for s in SIGNS
    orders = [weight_ordering(ci_weight(s), degrevlex(S)), weight_ordering(ci_weight(s), lex(S))]
    leading = [leading_ideal(J; ordering=o) for o in orders]
    @assert leading[1] == leading[2]                                    # tie-breaks do not matter
    pfaffian_terms = ideal(S, [leading_term(f; ordering=orders[1]) for f in PFAFFIANS])
    @assert pfaffian_terms == leading[1]                               # CI's Sturmfels--Zelevinsky step
    nonedges = supports(leading[1])
    @assert all(length(e) == 2 for e in nonedges) && length(nonedges) == 28
    complex = flag_complex(nonedges)
    @assert NAMES[complex.cone] == ["u", "v", "r"] && is_two_sphere(complex.sphere) &&
            length(vertices_of(complex.sphere)) == 11
    t11 = is_isomorphic(edge_graph(complex.sphere), edge_graph(T11))
    push!(CI_IDEALS, (signs=s, weight=ci_weight(s), nonedges=nonedges, complex=complex, t11=t11))
end
@assert length(unique(c.nonedges for c in CI_IDEALS)) == 8 && count(c -> c.t11, CI_IDEALS) == 2
for (n, c) in enumerate(CI_IDEALS)
    say("order ", n, " signs ", c.signs, ": 28 squarefree quadrics, cone over {u,v,r} of a 2-sphere",
        c.t11 ? " isomorphic to T_11" : "")
end

# ------------------------------------------------------------------ 2. our side, Groebner test
const PAIRS = [(i, j) for i in 1:7 for j in i+1:7]
R, P = polynomial_ring(QQ, ["p$(i)$(j)" for (i, j) in PAIRS])
stored = load(joinpath(DATA, "00_ideal.mrdi"))
I = ideal(R, hom(base_ring(stored), R, P).(gens(stored)))
# p_ij -> A[i, j] maps our ideal onto CI's: their coordinates are our Pluecker coordinates.
substitution = hom(R, S, [A[i, j] for (i, j) in PAIRS])
@assert ideal(S, substitution.(gens(I))) == J

"Each CI variable is A[i, j] (up to sign) for one or two pairs; p36 maps to -w-r."
const IMAGES = Dict(name => [k for (k, (i, j)) in enumerate(PAIRS) if A[i, j] in (var(name), -var(name))]
                    for name in NAMES)
const KEPT = Dict(name => first(IMAGES[name]) for name in NAMES)
const P36 = findfirst(==((3, 6)), PAIRS)
@assert all(1 <= length(IMAGES[n]) <= 2 for n in NAMES) && A[3, 6] == -w - r &&
        sort(vcat(values(IMAGES)...)) == setdiff(1:21, [P36])

"CI's order on our 21 variables: each eliminated variable leads its linear form."
function lifted_weight(ω)
    lifted = zeros(Int, 21)
    for (n, name) in enumerate(NAMES), k in IMAGES[name]
        lifted[k] = 4ω[n] + (k == KEPT[name] ? 0 : 1)
    end
    lifted[P36] = max(lifted[KEPT["w"]], lifted[KEPT["r"]]) + 1
    return lifted
end
for c in CI_IDEALS
    expected = ideal(R, vcat([P[KEPT[NAMES[a]]] * P[KEPT[NAMES[b]]] for (a, b) in c.nonedges],
                             [P[k] for name in NAMES for k in IMAGES[name] if k != KEPT[name]], [P[P36]]))
    @assert leading_ideal(I; ordering=weight_ordering(lifted_weight(c.weight), degrevlex(R))) == expected
end

# The census stores weights in the min convention; here the largest weight leads.
weights = load(joinpath(DATA, "04_canonical_cell_weights.mrdi"))
cell_weights = [[-Int(ZZ(weights[cell, k])) for k in 1:21] for cell in 1:nrows(weights)]
orbit_of = zeros(Int, length(cell_weights))
for line in readlines(joinpath(DATA, "05_symmetry_cell_orbits.csv"))[2:end]
    fields = split(line, ','; limit=4)
    orbit_of[parse.(Int, split(strip(fields[4], '"'), ';'))] .= parse(Int, fields[1])
end
toric_orbit = Dict{Int, Bool}()
for line in readlines(joinpath(DATA, "06_canonical_interior_degenerations.csv"))[2:end]
    fields = split(line, ',')
    toric_orbit[parse(Int, fields[1])] = fields[15] == "true"
end
@assert count(values(toric_orbit)) == 17
type_of_orbit = Dict{Int, Int}()
for line in readlines(TORIC_TABLE)
    m = match(r"\\oldstylenums\{(\d+)\}.*% QQ orbit (\d+)", line)
    isnothing(m) || (type_of_orbit[parse(Int, m.captures[2])] = parse(Int, m.captures[1]))
end
@assert sort(collect(keys(type_of_orbit))) == sort([o for (o, t) in toric_orbit if t])
label(cell) = haskey(type_of_orbit, orbit_of[cell]) ? "type $(type_of_orbit[orbit_of[cell]])" :
              "non-toric orbit $(orbit_of[cell])"

function cells_below(c)
    order = weight_ordering(lifted_weight(c.weight), degrevlex(R))
    theirs = leading_ideal(I; ordering=order)
    return Set(cell for cell in eachindex(cell_weights) if leading_ideal(I; ordering=weight_ordering(
        cell_weights[cell] .- minimum(cell_weights[cell]) .+ 1, order)) == theirs)
end
below_groebner = [cells_below(c) for c in CI_IDEALS]

# ------------------------------------------------------------------ 3. polyhedral test
# Lower hull of (points, heights) is the subdivision for the largest-weight-leads
# convention: x1 x3 - x2^2 with heights (0, 1, 0) has initial term x2^2.
@assert collect.(maximal_cells(subdivision_of_points(QQ[0; 1; 2], [0, 1, 0]))) == [[1, 3]]

"Points of the toric ideal in_w(I), one row per variable, in a lattice basis of Z^5."
function point_configuration(cell)
    initial_ideal = initial(I, tropical_semiring_map(QQ, max), QQ.(cell_weights[cell]))
    differences = [e[1] - e[2] for e in (collect(exponents(g)) for g in gens(groebner_basis(initial_ideal)))]
    @assert all(length(collect(exponents(g))) == 2 for g in gens(groebner_basis(initial_ideal)))
    K = kernel(matrix(ZZ, reduce(vcat, permutedims.(differences))); side=:right)   # 21 x 6, saturated
    @assert ncols(K) == 6
    shifted = matrix(ZZ, [K[k, j] - K[1, j] for k in 1:21, j in 1:6])
    H = hnf(shifted)
    basis = H[1:5, :]
    @assert iszero(H[6:end, :])
    coordinates(k) = (x = solve(basis, shifted[k:k, :]; side=:left); [Int(x[1, j]) for j in 1:5])
    return reduce(vcat, [permutedims(coordinates(k)) for k in 1:21])
end

function below_polyhedral(points, c)
    rows = [Vector(points[KEPT[name], :]) for name in NAMES]
    allunique(rows) && all(Vector(points[k, :]) in rows for k in 1:21) || return false
    cells = collect.(maximal_cells(subdivision_of_points(QQ.(reduce(vcat, permutedims.(rows))), c.weight)))
    all(length(cell) == 6 for cell in cells) || return false
    unimodular = all(abs(det(matrix(ZZ, [rows[cell[i]][j] - rows[cell[1]][j] for i in 2:6, j in 1:5]))) == 1
                     for cell in cells)
    facets = Set(sort(vcat(c.complex.cone, collect(t))) for t in c.complex.sphere)
    return unimodular && Set(sort.(cells)) == facets
end

toric_cells = [cell for cell in eachindex(cell_weights) if toric_orbit[orbit_of[cell]]]
configurations = Dict(cell => point_configuration(cell) for cell in toric_cells)
below_polytope = [Set(cell for cell in toric_cells if below_polyhedral(configurations[cell], c))
                  for c in CI_IDEALS]
for n in 1:8
    @assert below_polytope[n] == intersect(below_groebner[n], Set(toric_cells))
end
say("Groebner and polyhedral tests agree on all ", length(toric_cells), " toric cells")

# The point configuration of a cell is the Newton--Okounkov body of its type.
for cell in union(below_polytope...)
    body = load(joinpath(DATA, "06_canonical_interior_degenerations_objects",
                         "orbit_$(lpad(orbit_of[cell], 3, '0'))_newton_okounkov_body.poly"))
    @assert affine_unimodular_equivalent(convex_hull(QQ, configurations[cell]), body)
end

for (n, c) in enumerate(CI_IDEALS)
    say("order ", n, c.t11 ? " (T_11)" : "", ": cells ",
        join(("$cell ($(label(cell)))" for cell in sort(collect(below_groebner[n]))), ", "))
end
types(cells) = sort(unique(type_of_orbit[orbit_of[c]] for c in cells if haskey(type_of_orbit, orbit_of[c])))
say("toric types below some order: ", join(types(union(below_groebner...)), ", "))
say("toric types below a T_11 order: ",
    join(types(union((below_groebner[n] for n in 1:8 if CI_IDEALS[n].t11)...)), ", "))
say("non-toric orbits below some order: ",
    join(sort(unique(orbit_of[c] for c in union(below_groebner...) if !toric_orbit[orbit_of[c]])), ", "))
close(OUTPUT)
