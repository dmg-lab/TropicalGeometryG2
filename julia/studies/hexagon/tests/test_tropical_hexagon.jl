using Oscar
using Random

# Tests for the hexagon code of Section 9 of main.tex. Min convention.
#
# From the `julia` folder:  julia --project=. studies/hexagon/tests/test_tropical_hexagon.jl
# runs `run_all_tests()`. The first half of the file has the membership tests for tropical hexagons:
#
#   is_tropical_hexagon(H.flags)              # Oscar membership in the fans, slow
#   check_tropical_incidences(H.flags)        # the same with details, fast after the first call

include(joinpath(@__DIR__, "..", "tropical_complete_hexagon.jl"))

const LINE_FAN_PATH = joinpath(@__DIR__, "..", "..", "..", "data", "g2", "grassmannian", "03_normalized_coarsest_fan.mrdi")

# --- Membership in the tropical point space, line space and flag space -------------------------------------

"Is `x` a point of the tropical quadric? The minimum of x1+x4, x2+x5, x3+x6, 2*x7 must be attained at least twice."
is_tropical_point(x) = length(quadric_cone(x)) >= 2

"Does the vector `v` lie in some maximal cone of the fan `F`?"
in_fan(v, F) = any(C -> QQ.(v) in C, maximal_cones(F))

"Is `p` (21 Plücker coordinates) a point of the tropical G2-Grassmannian?"
is_tropical_line(p) = in_fan(p, load(LINE_FAN_PATH))

"Is `(x, p)` a point of the tropical G2 flag variety, i.e. a tropical point on a tropical line?"
is_tropical_flag(x, p) = in_fan(vcat(x, p), tropical_flag_fan())

"""
Is `flags = [(x_1, p_1), ..., (x_6, p_6)]` a tropical hexagon? Every (x_k, p_k) and every (x_k, p_{k-1}), with p_0 = p_6,
must be a point of the tropical flag variety: each point lies on its own line and on the previous one.
"""
function is_tropical_hexagon(flags)
    length(flags) == 6 || return false
    return all(is_tropical_flag(flags[k][1], flags[k][2]) && is_tropical_flag(flags[k][1], flags[mod1(k - 1, 6)][2]) for k in 1:6)
end

const _FLAG_CONES = Ref{Any}(nothing)

"""
Every maximal cone of the fan of trop(flagspace) as matrices `(E, A)`, the cone being `E*v = 0`, `A*v <= 0`. Computed on
first use (about a minute) and cached, as the membership of Oscar recomputes the facets of the cones on every call.
"""
function flag_cones()
    if _FLAG_CONES[] === nothing
        _FLAG_CONES[] = [(E=linear_equation_matrix(linear_span(C)), A=linear_inequality_matrix(facets(C)))
                         for C in maximal_cones(tropical_flag_fan())]
    end
    return _FLAG_CONES[]
end

"The index of a maximal cone of the fan of trop(flagspace) containing `(x, p)`, or `nothing`."
function flag_cone_index(x, p)
    v = matrix(QQ, 28, 1, QQ.(vcat(x, p)))
    return findfirst(c -> all(iszero, c.E * v) && all(<=(0), c.A * v), flag_cones())
end

"""
The fast version of `is_tropical_hexagon`, with details: `ok` (all 12 incidences hold), `own` and `previous` (for
k = 1..6 the cone containing (x_k, p_k), resp. (x_k, p_{k-1}), or `nothing`), and `failed`, the failing incidences.
"""
function check_tropical_incidences(flags)
    own = [flag_cone_index(flags[k][1], flags[k][2]) for k in 1:6]
    previous = [flag_cone_index(flags[k][1], flags[mod1(k - 1, 6)][2]) for k in 1:6]
    failed = vcat([(k, :own) for k in 1:6 if own[k] === nothing], [(k, :previous) for k in 1:6 if previous[k] === nothing])
    return (ok=isempty(failed), own=own, previous=previous, failed=failed)
end

# --- Lemma 9.1 and Remark 9.2 ------------------------------------------------------------------------------

"The first coordinate of x^(3) as printed in Remark 9.2 of main.tex."
function remark_x3_1()
    p(i, j) = pl(P_VARS, i, j)
    x1, x2, x3, x4, x5, x6, x7 = X_VARS
    return p(1, 2) * x5 * x7 + p(1, 3) * x6 * x7 + p(1, 4) * x1 * x7 + p(1, 5) * x2 * x7 + p(1, 5) * x4 * x6 +
           p(1, 6) * x3 * x7 - p(1, 6) * x4 * x5 + p(1, 7) * x1 * x4 - 2 * p(1, 7) * x7^2 + p(2, 5) * x5 * x6 -
           p(2, 6) * x5^2 + p(2, 7) * x1 * x5 + p(3, 5) * x6^2 - p(3, 6) * x5 * x6 + p(3, 7) * x1 * x6 +
           p(4, 5) * x1 * x6 - p(4, 6) * x1 * x5 + p(4, 7) * x1^2 - p(5, 6) * x2 * x5 - p(5, 6) * x3 * x6 +
           p(5, 7) * x1 * x2 + 2 * p(5, 7) * x6 * x7 + p(6, 7) * x1 * x3 - 2 * p(6, 7) * x5 * x7
end

"The expressions have the bidegrees of Remark 9.2 and integer coefficients, x^(3)_1 is as printed there."
function test_expressions()
    H = hexagon_expressions()
    for (name, degree) in ((:x2, (1, 1)), (:x3, (1, 2)), (:p23, (2, 3)), (:p34, (1, 3)))
        @assert unique(filter(!isnothing, bidegree.(H[name]))) == [degree] "bidegree of $name"
        @assert all(isone(denominator(c)) for f in H[name] for c in coefficients(f)) "denominators in $name"
    end
    @assert H.x3[1] == remark_x3_1() "x^(3)_1 differs from Remark 9.2"
    (x, p), y = random_flag(QQ), random_point(QQ)
    A = evaluate_expressions(p, y)
    v = vcat(p, y)
    @assert A.p23 == [evaluate(f, v) for f in H.p23] && A.p34 == [evaluate(f, v) for f in H.p34]
    println("expressions: bidegrees and x^(3)_1 as in Remark 9.2, integer coefficients")
end

"""
`complete_hexagon` over several fields, in particular characteristic 2: its checks pass, x^(2) is orthogonal to x^(4), and
completing the opposite flags (x^(k), p^(k,k+1)), (x^(k+3), p^(k+3,k+4)) gives the same hexagon.
"""
function test_complete_hexagon(; trials=20)
    Kt, t = rational_function_field(QQ, "t")
    for K in (QQ, GF(2), GF(2, 2), GF(2, 4), GF(3), GF(5), GF(7, 2), Kt)
        count = 0
        for _ in 1:trials
            f1, f4 = random_flag(K), random_flag(K)
            H = try
                complete_hexagon(f1, f4)
            catch err
                err isa ArgumentError && occursin("not generic", err.msg) ? continue : rethrow()
            end
            count += 1
            @assert all(iszero(polar(H.points[k], H.points[mod1(k + 2, 6)])) for k in 1:6) "orthogonality over $K"
            H2 = try
                complete_hexagon(H.flags[2], H.flags[5])
            catch err
                err isa ArgumentError && occursin("not generic", err.msg) ? continue : rethrow()
            end
            for k in 1:6
                @assert proportional(H2.points[k], H.points[mod1(k + 1, 6)]) "rotated hexagon over $K"
                @assert proportional(H2.lines[k], H.lines[mod1(k + 1, 6)]) "rotated hexagon over $K"
            end
        end
        @assert count > 0 "no generic flags over $K"
        println("complete_hexagon: $count of $trials random pairs over $K")
    end
    @assert (try complete_hexagon(random_flag(QQ)[1], random_flag(QQ)[2], random_flag(QQ)...); false catch err; err isa ArgumentError end)
end

# --- Algorithm 2 ------------------------------------------------------------------------------------------

"""
Algorithm 2 on small examples: Example ex:tropically_faithful has no faithful representative; a monomial is faithful; on
x1 x2 + x1 x3 + x4^2 the initial ideal is not saturated, and one iteration of the loop gives x2 + x3 = -x4^2 / x1.
"""
function test_algorithm()
    S, (x1, x2, x3) = polynomial_ring(QQ, [:x1, :x2, :x3])
    I = ideal(S, [x1 * x2 + x2^2 - x3^2])
    w = [1, 0, 0]
    T = faithful_setup(I, w)
    @assert length(T.primes) == 2 && T.k == 0
    r = tropically_faithful_representative(x2 + x3, T)
    @assert r.status == :no_representative && r.value == 0 "Example ex:tropically_faithful"
    @assert tropically_faithful_representative(I, w, x1 * x2^2).status == :faithful
    r = tropically_faithful_representative(x2^2 - x3^2, T)
    @assert r.status == :faithful && r.value == 1 "x2^2 - x3^2 = -x1 x2 on V(I)"
    r = tropically_faithful_representative(x1 * x2 + x2^2 - x3^2, T)
    @assert r.status == :in_ideal

    S, (a, b, c, d) = polynomial_ring(QQ, [:a, :b, :c, :d])
    I = ideal(S, [a * b + a * c + d^2])
    T = faithful_setup(I, [0, 0, 0, 1])
    @assert T.k == 1 && T.saturation == ideal(S, [b + c])
    r = tropically_faithful_representative(b + c, T)
    @assert r.status == :faithful && r.steps == 1 && r.value == 2 && r.shift == [1, 0, 0, 0]
    @assert (try faithful_setup(I, [1, 0, 0, 0]); false catch err; err isa ArgumentError end) "w is not in the tropical variety"
    println("Algorithm 2: Example ex:tropically_faithful, and a non-saturated initial ideal")
end

# --- Tropical hexagons --------------------------------------------------------------------------------------

"""
A random element of QQ(t) for lifts: a random large integer times t^a, plus a term of higher t-order, with `a` in `orders`.
"""
function random_lift_element(K; orders=0:40, coeff=10^4)
    t = gen(K)
    a = rand(orders)
    c = rand(setdiff(-coeff:coeff, 0))
    return K(c) * t^a + K(rand(-coeff:coeff)) * t^(a + 1)
end

"The t-adic valuations of a classical hexagon over QQ(t), every vector shifted to minimum 0."
tropicalized_hexagon(H) = (points=shift_to_zero.(tropicalize.(H.points)), lines=shift_to_zero.(tropicalize.(H.lines)))

"""
`tropical_complete_hexagon` against the classical hexagon of random flags over QQ(t), which is generic: all valuations
agree. The output of random points of the fan passes all 12 incidences.
"""
function test_tropical_complete_hexagon(; lifts=10, fan_pairs=10)
    K, t = rational_function_field(QQ, "t")
    for n in 1:lifts
        sample = () -> random_lift_element(K)
        H = complete_hexagon(random_flag(K; sample), random_flag(K; sample))
        Tr = tropicalized_hexagon(H)
        T = tropical_complete_hexagon(tropicalize(H.points[1]), tropicalize(H.lines[1]),
                                      tropicalize(H.points[4]), tropicalize(H.lines[4]))
        @assert all(isempty, T.no_representative) "lift $n: no faithful representative"
        @assert T.points == Tr.points "lift $n: points differ"
        @assert T.lines == Tr.lines "lift $n: lines differ"
    end
    println("tropical_complete_hexagon: agrees with the valuations of $lifts lifts over QQ(t), $(length(_CONE_SETUPS)) Gröbner cones so far")
    for n in 1:fan_pairs
        H = tropical_complete_hexagon(random_tropical_flag(), random_tropical_flag())
        result = check_tropical_incidences(H.flags)
        @assert result.ok "random flags $n: failing incidences $(result.failed)"
    end
    println("tropical_complete_hexagon: all incidences hold for $fan_pairs random pairs of fan points, $(length(_CONE_SETUPS)) Gröbner cones so far")
end

"""
A cone of type 33, where the initial ideal of the ideal of Remark 9.2 is not saturated: the normal form of the coordinate
x^(3)_6 has an initial form in the prime, so its tropical value is too small, and the loop of Algorithm 2 raises it.
"""
function test_type33()
    @assert is_prime(line_point_ideal()) "the ideal of Remark 9.2 is not prime"
    p = [-10, -10, -10, 0, -10, -10, -12, 0, -10, -12, -12, -12, -10, -10, -10, -10, -12, -12, -10, -10, -10]
    x = [-6, -8, 22, 18, 8, 4, 0]
    w = vcat(p, x)
    T = cone_setup(w)
    f = hexagon_expressions().x3[6]
    r = tropically_faithful_representative(f, T; w)
    @assert T.k == 1 && r.status == :faithful && r.steps == 1
    @assert tropical_value(reduce_modulo(T, f), QQ.(w)) < r.value
    println("Algorithm 2 on a cone of type 33: x^(3)_6 is raised by the loop")
end

function run_all_tests()
    Random.seed!(1)
    test_expressions()
    test_complete_hexagon()
    test_algorithm()
    test_type33()
    test_tropical_complete_hexagon()
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_all_tests()
end
