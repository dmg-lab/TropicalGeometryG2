###############################################################################
# The G2 polytope of Feigin, Fourier and Littelmann (Transform. Groups 22
# (2017), Definition 11.5) at lambda = omega_2, against the Newton--Okounkov
# polytopes of the toric degenerations of the G2-Grassmannian in the census.
#
# Their roots are beta_1, ..., beta_6 = 3a1+2a2, 3a1+a2, 2a1+a2, a1+a2, a2, a1
# with a1 short, f_beta_1 > ... > f_beta_6, and the homogeneous reverse
# lexicographic order on PBW monomials (Section 11). V(omega_2) is the adjoint
# module with cyclic vector the highest root vector, so (k, l) = (0, 1).
#
# 1. The essential multi-exponents of V(omega_2) and of V(2 omega_2), the
#    Cartan component in V(omega_2)^{(x)2}, computed from the Lie bracket, are
#    the lattice points of Definition 11.5 for (k, l) = (0, 1) and (0, 2).
# 2. The polytope for (0, 1), on p_6 = 0, is lattice equivalent to the body
#    of exactly one type of Table tab:g2-grassmannian-newton-okounkov.
# 3. For the flag, (k, l) = (1, 1): the essential multi-exponents of
#    V(omega_1) and V(omega_1 + omega_2), the Cartan component in
#    V(omega_1) (x) V(omega_2), are the lattice points for (1, 0) and (1, 1),
#    and P(omega_1 + omega_2) is compared with the body of every toric orbit
#    of the tropical flag variety, both alone and as the pair of summands
#    (P(omega_1), P(omega_2)) under one lattice map, via Cayley polytopes
#    (../flag_census.jl).
#
# Run with julia --project=julia from any directory; writes verify.txt next
# to this file. Characteristic 0 only.
###############################################################################

using TropicalCensus
using Oscar
include(joinpath(@__DIR__, "..", "flag_census.jl"))

const DATA = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "data", "g2"))
const TORIC_TABLE = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "..", "tables",
                                      "toric_degeneration_table.tex"))
const OUTPUT = open(joinpath(@__DIR__, "verify.txt"), "w")
say(parts...) = (println(parts...); println(OUTPUT, parts...))

# ------------------------------------------------------------------ the adjoint module
L = lie_algebra(QQ, :G, 2)
R = root_system(L)
x, y, _ = chevalley_basis(L)
# Oscar lists the positive roots as a1, a2, a1+a2, 2a1+a2, 3a1+a2, 3a1+2a2.
@assert [Int.(vec(Matrix(coefficients(r)))) for r in positive_roots(R)] ==
        [[1, 0], [0, 1], [1, 1], [2, 1], [3, 1], [3, 2]]
const F = [y[7 - i] for i in 1:6]        # F[i] = f_{beta_i}
const AD = [matrix(QQ, reduce(hcat, [Vector{QQFieldElem}(coefficients(f * b)) for b in basis(L)]))
            for f in F]                  # ad f_{beta_i} on the basis of L, acting on columns
const HIGHEST = matrix(QQ, 14, 1, Vector{QQFieldElem}(coefficients(x[6])))

"f_1^p_1 ... f_6^p_6 on a vector, f_6 first. `act(i, v)` applies f_i."
function monomial_action(act, p, vector)
    for i in 6:-1:1, _ in 1:p[i]
        vector = act(i, vector)
    end
    return vector
end

"Homogeneous reverse lexicographic order with f_1 > ... > f_6: is f^a < f^b?"
revlex_less(a, b) = sum(a) == sum(b) ?
    (i = findlast(!iszero, a .- b); !isnothing(i) && a[i] > b[i]) : sum(a) < sum(b)

compositions(degree, parts) = parts == 1 ? [[degree]] :
    [vcat(first, rest) for first in 0:degree for rest in compositions(degree - first, parts - 1)]

as_row(X) = matrix(QQ, 1, nrows(X) * ncols(X), [X[i, j] for j in 1:ncols(X) for i in 1:nrows(X)])

"Essential multi-exponents: f^p v not in the span of the f^q v with q < p."
function essential_exponents(act, cyclic, dimension; less=revlex_less)
    essential, span = Vector{Int}[], zero_matrix(QQ, 0, nrows(cyclic) * ncols(cyclic))
    for degree in 0:typemax(Int)
        for p in sort(compositions(degree, 6); lt=less)
            vector = monomial_action(act, p, cyclic)
            extended = vcat(span, as_row(vector))
            rank(extended) > nrows(span) || continue
            push!(essential, p)
            span = extended
            nrows(span) == dimension && return essential
        end
    end
end

# On V^{(x)2}, a 14 x 14 matrix X stands for sum X_ab e_a (x) e_b.
adjoint_action(i, v) = AD[i] * v
square_action(i, X) = AD[i] * X + X * transpose(AD[i])
# The cyclic modules are V(omega_2) and V(2 omega_2), of dimensions 14 and 77.
@assert dim_of_simple_module(L, [0, 1]) == 14 && dim_of_simple_module(L, [0, 2]) == 77
es1 = essential_exponents(adjoint_action, HIGHEST, 14)
es2 = essential_exponents(square_action, HIGHEST * transpose(HIGHEST), 77)

# ------------------------------------------------------------------ Definition 11.5
"P(k omega_1 + l omega_2) in R^6, as printed in Definition 11.5."
function ffl_polytope(k, l)
    A = [-1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 -1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0; 0 0 0 0 0 -1;
         0 0 0 0 1 0; 0 0 0 0 0 1;
         0 1 1 0 0 1; 0 0 1 1 0 1; 0 0 0 1 1 1;
         1 1 1 1 1 0; 0 1 1 1 1 1]
    b = [0, 0, 0, 0, 0, 0, l, k, k + l, k + l, k + l, k + 2l, k + 2l]
    return polyhedron(A, b)
end
lattice_set(P) = Set(Vector{Int}(Int.(collect(q))) for q in lattice_points(P))

S1, S2 = lattice_set(ffl_polytope(0, 1)), lattice_set(ffl_polytope(0, 2))
say("V(omega_2): ", length(es1), " essential exponents; Definition 11.5 at (0,1): ",
    length(S1), " lattice points; equal: ", Set(es1) == S1)
say("V(2 omega_2): ", length(es2), " essential exponents; Definition 11.5 at (0,2): ",
    length(S2), " lattice points; equal: ", Set(es2) == S2,
    "; S(2 omega_2) = S(omega_2) + S(omega_2): ", S2 == Set(a + b for a in S1 for b in S1))
deglex_less(a, b) = sum(a) == sum(b) ? (i = findfirst(!iszero, a .- b); !isnothing(i) && a[i] < b[i]) :
                                       sum(a) < sum(b)
say("homogeneous lexicographic order instead: equal to Definition 11.5: ",
    Set(essential_exponents(adjoint_action, HIGHEST, 14; less=deglex_less)) == S1)
say("essential exponents (p_1..p_6) of V(omega_2): ", join((join(p) for p in sort(es1)), " "))
@assert Set(es1) == S1 && Set(es2) == S2

# p_6 = 0 on P(omega_2); dropping it is a lattice isomorphism onto Z^5.
P6 = ffl_polytope(0, 1)
@assert all(v[6] == 0 for v in vertices(P6))
FFL = convex_hull(QQ, matrix(QQ, reduce(vcat, [permutedims(collect(v)[1:5]) for v in vertices(P6)])))
say("FFL polytope: dimension ", dim(FFL), ", f-vector ", Int.(f_vector(FFL)), ", ",
    length(lattice_points(FFL)), " lattice points, normalized volume ",
    normalized_volume(FFL), ", lattice polytope: ", is_lattice_polytope(FFL))

# ------------------------------------------------------------------ the census bodies
body_path(field, orbit) = joinpath(DATA, field == "QQ" ? "grassmannian" : "grassmannian_gf2",
    "06_canonical_interior_degenerations_objects",
    "orbit_$(lpad(orbit, 3, '0'))_newton_okounkov_body.poly")
types = Tuple{Int, String, Int}[]
for line in readlines(TORIC_TABLE)
    m = match(r"\\oldstylenums\{(\d+)\}.*% (?:QQ orbit (\d+), )?GF\(2\) orbit (\d+)", line)
    isnothing(m) && continue
    field, orbit = isnothing(m.captures[2]) ? ("GF(2)", m.captures[3]) : ("QQ", m.captures[2])
    push!(types, (parse(Int, m.captures[1]), field, parse(Int, orbit)))
end
@assert first.(types) == 1:20

matches = Int[]
for (type, field, orbit) in types
    body = load(body_path(field, orbit))
    equivalent = affine_unimodular_equivalent(FFL, body)
    equivalent && push!(matches, type)
    say("type ", lpad(type, 2), " (", field, " orbit ", lpad(orbit, 2), "): f-vector ",
        Int.(f_vector(body)), ", lattice equivalent to FFL: ", equivalent)
end
say("types lattice equivalent to the FFL polytope P(omega_2): ", join(matches, ", "))

# ------------------------------------------------------------------ the flag, (k, l) = (1, 1)
V1 = simple_module(L, [1, 0])
module_matrix(e) = matrix(QQ, reduce(hcat, [Vector{QQFieldElem}(coefficients(e * v)) for v in basis(V1)]))
const AD1 = [module_matrix(f) for f in F]  # f_{beta_i} on V(omega_1), acting on columns
const HIGHEST1 = kernel(reduce(vcat, [module_matrix(e) for e in x]); side=:right)
@assert size(HIGHEST1) == (7, 1)
# On V(omega_1) (x) V(omega_2), a 7 x 14 matrix X stands for sum X_ab v_a (x) e_b.
standard_action(i, v) = AD1[i] * v
tensor_action(i, X) = AD1[i] * X + X * transpose(AD[i])
@assert dim_of_simple_module(L, [1, 0]) == 7 && dim_of_simple_module(L, [1, 1]) == 64
es10 = essential_exponents(standard_action, HIGHEST1, 7)
es11 = essential_exponents(tensor_action, HIGHEST1 * transpose(HIGHEST), 64)

S10, S11 = lattice_set(ffl_polytope(1, 0)), lattice_set(ffl_polytope(1, 1))
say("V(omega_1): ", length(es10), " essential exponents; Definition 11.5 at (1,0): ",
    length(S10), " lattice points; equal: ", Set(es10) == S10)
say("V(omega_1+omega_2): ", length(es11), " essential exponents; Definition 11.5 at (1,1): ",
    length(S11), " lattice points; equal: ", Set(es11) == S11,
    "; S(omega_1+omega_2) = S(omega_1) + S(omega_2): ", S11 == Set(a + b for a in S10 for b in S1))
@assert Set(es10) == S10 && Set(es11) == S11
FFL11 = ffl_polytope(1, 1)
say("FFL polytope P(omega_1+omega_2): dimension ", dim(FFL11), ", f-vector ", Int.(f_vector(FFL11)),
    ", ", length(S11), " lattice points, normalized volume ", normalized_volume(FFL11),
    ", lattice polytope: ", is_lattice_polytope(FFL11))

# The FFL points: degree (1,0) for S(omega_1), (0,1) for S(omega_2).
FFL_CAYLEY = pair_cayley_polytope(sort(collect(S10)), sort(collect(S1)))
for field in ("QQ", "GF(2)")
    report_flag_census(say, field, flag_census_matches(field, FFL11, FFL_CAYLEY), "P(omega_1+omega_2)")
end
close(OUTPUT)
