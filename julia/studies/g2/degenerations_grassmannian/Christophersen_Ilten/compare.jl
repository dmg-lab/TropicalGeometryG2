###############################################################################
# The Stanley--Reisner degenerations of the G2-Grassmannian of Christophersen
# and Ilten (arXiv:1102.4521, Theorem 3.1, case g = 10) inside our Groebner fan.
#
# Their 7x7 skew matrix is our Pluecker matrix on the linear space (eq:7linear),
# so their 14 variables are Pluecker coordinates. Their conditions on the term
# order leave three comparisons open, which gives eight squarefree initial
# ideals; two of them are the ideal of T11 * Delta_2. For each of the eight we
# list the coarsest cells of the tropical G2-Grassmannian whose canonical
# weight lies in the closure of the Groebner cone, that is, the initial
# degenerations of ours that degenerate further to theirs.
#
# Run with julia --project=julia from any directory; writes compare.txt next
# to this file. Characteristic 0 only.
###############################################################################

using Oscar, Random, LinearAlgebra

const DATA = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "data", "g2", "grassmannian"))
const TORIC_TABLE = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "..", "tables",
                                      "toric_degeneration_table.tex"))
const OUTPUT = open(joinpath(@__DIR__, "compare.txt"), "w")
say(parts...) = (println(parts...); println(OUTPUT, parts...))

pfaffian(M, i, j, k, l) = M[i, j] * M[k, l] - M[i, k] * M[j, l] + M[i, l] * M[j, k]
const QUADRUPLES = [(i, j, k, l) for i in 1:7 for j in i+1:7 for k in j+1:7 for l in k+1:7]

# ------------------------------------------------------------------ their 14 variables
const NAMES = ["u", "v", "w", "r", "x00", "x01", "x10", "x11", "y00", "y01", "y10", "y11", "z0", "z1"]
const POSITION = Dict(name => index for (index, name) in enumerate(NAMES))
S, V = polynomial_ring(QQ, NAMES)
(u, v, w, r, x00, x01, x10, x11, y00, y01, y10, y11, z0, z1) = V
z = zero(S)
A = S[ z    -x10  x11  w    y11  y10   u;
       x10   z   -v    y00  r    z0    x00;
      -x11   v    z    y01  z1  -w-r   x01;
      -w    -y00 -y01  z    x01 -x00   v;
      -y11  -r   -z1  -x01  z    u     x11;
      -y10  -z0   w+r  x00 -u    z     x10;
      -u    -x00 -x01 -v   -x11 -x10   z ]
@assert all(A[i, j] == -A[j, i] for i in 1:7, j in 1:7)
J = ideal(S, [pfaffian(A, q...) for q in QUADRUPLES])
say("ideal of the 35 Pfaffians in their coordinates: dimension ", dim(J), ", degree ", degree(J))

# u, v < r, w, x_ij < y_ij, z_i, and middle * middle > right * left.
const LEFT, MIDDLE, RIGHT = 1:2, 3:8, 9:14
function admissible_weight(rng)
    weight = zeros(Int, 14)
    weight[LEFT] .= 1; weight[MIDDLE] .= 1001; weight[RIGHT] .= 1501
    return weight .+ rand(rng, 0:60, 14)
end

"""
The leading exponents of the reduced Groebner basis for the weight order
(largest weight leads, ties broken by degrevlex), and the inequalities
`weight * row >= 0` of the closure of its Groebner cone.
"""
function initial_data(ring, ideal_value, weight)
    ordering = weight_ordering(weight, degrevlex(ring))
    basis = gens(groebner_basis(ideal_value; ordering=ordering, complete_reduction=true))
    leading = [first(exponents(leading_monomial(g; ordering=ordering))) for g in basis]
    inequalities = Vector{Int}[]
    for (g, e0) in zip(basis, leading), e in exponents(g)
        e == e0 || push!(inequalities, e0 .- e)
    end
    return (key=sort(leading), inequalities=unique(inequalities), weight=weight)
end

rng = MersenneTwister(20261010)
found = Dict{Any, Any}()
for _ in 1:160
    data = initial_data(S, J, admissible_weight(rng))
    haskey(found, data.key) || (found[data.key] = data)
end
say("distinct initial ideals for 160 random admissible term orders: ", length(found))

"The diagonals of the three quadrangles of their Figure 3 that are edges of the complex."
function diagonals(data)
    nonedge(a, b) = (e = zeros(Int, 14); e[POSITION[a]] = 1; e[POSITION[b]] = 1; e in data.key)
    return join((nonedge("x00", "x11") ? "x01-x10" : "x00-x11",
                 nonedge("x00", "y01") ? "x01-y00" : "x00-y01",
                 nonedge("x10", "y11") ? "x11-y10" : "x10-y11"), ", ")
end

function graph_of(edges, vertices)
    graph = Graph{Undirected}(vertices)
    for (a, b) in edges
        add_edge!(graph, a, b)
    end
    return graph
end

"T11: the gyroelongated square bipyramid T10 with a vertex starred into an edge of a square."
function t11_graph()
    apex, other_apex, new = 1, 10, 11
    s(i) = 1 + mod1(i, 4); t(i) = 5 + mod1(i, 4)
    edges = Set{Tuple{Int, Int}}()
    for i in 1:4, edge in ((apex, s(i)), (other_apex, t(i)), (s(i), s(i + 1)), (t(i), t(i + 1)),
                           (s(i), t(i)), (s(i + 1), t(i)))
        push!(edges, minmax(edge...))
    end
    delete!(edges, minmax(s(1), s(2)))
    for a in (s(1), s(2), apex, t(1))
        push!(edges, minmax(new, a))
    end
    return graph_of(edges, 11)
end
const T11 = t11_graph()

"The flag complex of a squarefree quadratic monomial ideal: cone points and edge graph."
function describe(data)
    key = data.key
    @assert all(sum(e) == 2 && maximum(e) == 1 for e in key)
    used = findall(i -> any(e[i] == 1 for e in key), 1:14)
    nonedges = Set(Tuple(findall(==(1), e)) for e in key)
    edges = [(i, j) for i in used for j in used if i < j && !((i, j) in nonedges)]
    relabel = Dict(a => index for (index, a) in enumerate(used))
    graph = graph_of([(relabel[i], relabel[j]) for (i, j) in edges], length(used))
    return (simplex=NAMES[setdiff(1:14, used)], vertices=length(used), edges=length(edges),
            valencies=sort([count(e -> i in e, edges) for i in used]),
            is_t11=length(used) == 11 && length(edges) == 27 && is_isomorphic(graph, T11))
end

ideals = sort(collect(values(found)); by=diagonals)
for (n, data) in enumerate(ideals)
    info = describe(data)
    say("M$n  diagonals {", diagonals(data), "}  generators ", length(data.key),
        "  simplex ", join(info.simplex, " "), "  sphere: ", info.vertices, " vertices, ",
        info.edges, " edges, valencies ", join(info.valencies, " "), info.is_t11 ? "   <== T11" : "")
end

# ------------------------------------------------------------------ our 21 variables
const PAIRS = [(i, j) for i in 1:7 for j in i+1:7]
const INDEX = Dict(pair => index for (index, pair) in enumerate(PAIRS))
R, P = polynomial_ring(QQ, ["p$(i)$(j)" for (i, j) in PAIRS])
p(i, j) = P[INDEX[(i, j)]]
linear_forms = [p(2,3) + p(4,7), p(1,3) - p(5,7), p(1,2) + p(6,7), p(1,7) - p(5,6),
                p(2,7) + p(4,6), p(3,7) - p(4,5), p(1,4) + p(2,5) + p(3,6)]
pluecker = [i == j ? zero(R) : (i < j ? p(i, j) : -p(j, i)) for i in 1:7, j in 1:7]
I = ideal(R, vcat(linear_forms, [pfaffian(pluecker, q...) for q in QUADRUPLES]))
stored = load(joinpath(DATA, "00_ideal.mrdi"))
@assert ideal(R, hom(base_ring(stored), R, P).(gens(stored))) == I

# Their matrix is our Pluecker matrix on the linear space (eq:7linear).
@assert all(iszero, hom(R, S, [A[i, j] for (i, j) in PAIRS]).(linear_forms))
const KEPT = Dict("u" => (1,7), "v" => (4,7), "w" => (1,4), "r" => (2,5), "x00" => (2,7),
                  "x01" => (3,7), "x10" => (6,7), "x11" => (1,3), "y00" => (2,4), "y01" => (3,4),
                  "y10" => (1,6), "y11" => (1,5), "z0" => (2,6), "z1" => (3,5))
const ELIMINATED = [((5,6), (1,7)), ((2,3), (4,7)), ((4,6), (2,7)), ((4,5), (3,7)),
                    ((1,2), (6,7)), ((5,7), (1,3))]
say("dictionary: ", join(("$name = p$(KEPT[name][1])$(KEPT[name][2])" for name in NAMES), ", "),
    ", -w-r = p36")

# The census stores weights in the min convention; here the largest weight leads.
weights = load(joinpath(DATA, "04_canonical_cell_weights.mrdi"))
cell_weights = [[-Int(ZZ(weights[cell, k])) for k in 1:21] for cell in 1:nrows(weights)]

orbit_of = zeros(Int, length(cell_weights))
for line in readlines(joinpath(DATA, "05_symmetry_cell_orbits.csv"))[2:end]
    fields = split(line, ','; limit=4)
    orbit_of[parse.(Int, split(strip(fields[4], '"'), ';'))] .= parse(Int, fields[1])
end
# Index of the toric orbits in the table of the paper, read off its generated comments.
table_index = Dict{Int, String}()
for line in readlines(TORIC_TABLE)
    m = match(r"\\oldstylenums\{(\d+)\}.*% QQ orbit (\d+)", line)
    m === nothing || (table_index[parse(Int, m.captures[2])] = m.captures[1])
end
label(cell) = string(cell, " (", get(table_index, orbit_of[cell], "non-toric orbit $(orbit_of[cell])"), ")")

"Their chart eliminates p36: the weight must attain the maximum of p14, p25, p36 at p36."
function in_their_chart(weight)
    @assert all(weight[INDEX[a]] == weight[INDEX[b]] for (a, b) in ELIMINATED)
    return weight[INDEX[(3,6)]] == max(weight[INDEX[(1,4)]], weight[INDEX[(2,5)]])
end
restricted(weight) = [weight[INDEX[KEPT[name]]] for name in NAMES]
say("coarsest cells: ", length(cell_weights), ", of which ", count(in_their_chart, cell_weights),
    " have the maximum of p14, p25, p36 attained at p36")

below = [[cell for cell in eachindex(cell_weights) if in_their_chart(cell_weights[cell]) &&
          all(dot(restricted(cell_weights[cell]), row) >= 0 for row in data.inequalities)]
         for data in ideals]

# ------------------------------------------------------------------ check in 21 variables
"Their term order on our 21 variables: each eliminated variable leads its linear form."
function lifted_weight(weight)
    lifted = zeros(Int, 21)
    for name in NAMES
        lifted[INDEX[KEPT[name]]] = 4 * weight[POSITION[name]]
    end
    for (a, b) in ELIMINATED
        lifted[INDEX[a]] = lifted[INDEX[b]] + 1
    end
    lifted[INDEX[(3,6)]] = max(lifted[INDEX[(1,4)]], lifted[INDEX[(2,5)]]) + 1
    return lifted
end
leading_key(ordering) = sort([first(exponents(leading_monomial(g; ordering=ordering)))
    for g in gens(groebner_basis(I; ordering=ordering, complete_reduction=true))])
function lifted_key(data)
    key = Vector{Int}[]
    for e in data.key
        lifted = zeros(Int, 21)
        for name in NAMES
            lifted[INDEX[KEPT[name]]] = e[POSITION[name]]
        end
        push!(key, lifted)
    end
    for pair in vcat(first.(ELIMINATED), [(3,6)])
        lifted = zeros(Int, 21); lifted[INDEX[pair]] = 1
        push!(key, lifted)
    end
    return sort(key)
end

# in_<(in_w(I)) = in_<(I) for the weight w of the cell refined by their order <.
for (n, data) in enumerate(ideals)
    theirs = weight_ordering(lifted_weight(data.weight), degrevlex(R))
    @assert leading_key(theirs) == lifted_key(data)
    direct = [cell for cell in eachindex(cell_weights) if
              leading_key(weight_ordering(cell_weights[cell] .- minimum(cell_weights[cell]) .+ 1, theirs)) ==
              lifted_key(data)]
    @assert direct == below[n]
    say("M$n: cells ", join(label.(below[n]), ", "))
end
say("cells below all eight: ", join(label.(intersect(below...)), ", "))
say("the test in 14 variables and the direct test in 21 variables agree on all cells")
close(OUTPUT)
