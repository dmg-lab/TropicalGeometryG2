using Test
using Oscar
using Random

# The incidences of Lemma 9.1 (lem:complete_hexagon) of main.tex for `hexagon_expressions` of
# studies/hexagon/polynomial_expressions.jl: from p^(12) and x^(4) it constructs x^(2), x^(3), p^(23), p^(34) with
#
#     p^(12) ∋ x^(2) ∈ p^(23) ∋ x^(3) ∈ p^(34) ∋ x^(4).
#
# The study code is loaded into its own module, as it defines top-level constants such as `EXPRESSION_RING`.
#
# From the `julia` folder:  julia --project=. studies/hexagon/tests/hexagon_expressions.jl

module HexagonExpressionsTest
include(joinpath(@__DIR__, "..", "complete_hexagon.jl"))
end

const HX = HexagonExpressionsTest

@testset "hexagon_expressions: incidences of Lemma 9.1 on linespace × pointspace" begin
    H = HX.hexagon_expressions()
    p12, x4 = HX.P_VARS, HX.X_VARS
    I = HX.line_point_ideal()
    in_I(f) = iszero(normal_form(f, I))

    # x^(2) lies on p^(12) and is orthogonal to x^(4); x^(2), x^(3) lie on the quadric.
    @test all(in_I, HX.incidence_relations(H.x2, p12))
    @test in_I(HX.polar(H.x2, x4))
    @test in_I(HX.quadric(H.x2))
    @test in_I(HX.quadric(H.x3))

    # p^(23) and p^(34) lie in linespace: the 7 linear forms modulo I, the Plücker relations identically.
    @test all(in_I, HX.linear_forms(H.p23))
    @test all(in_I, HX.linear_forms(H.p34))
    @test all(iszero, HX.plucker_relations(H.p34))

    # x^(2), x^(3) lie on p^(23) and x^(3), x^(4) on p^(34), identically.
    @test all(iszero, HX.incidence_relations(H.x2, H.p23))
    @test all(iszero, HX.incidence_relations(H.x3, H.p23))
    @test all(iszero, HX.incidence_relations(H.x3, H.p34))
    @test all(iszero, HX.incidence_relations(x4, H.p34))
end

@testset "hexagon_expressions: incidences of Lemma 9.1 at random points over several fields" begin
    Random.seed!(1)
    for K in (QQ, GF(2), GF(2, 4), GF(3), GF(7, 2))
        count = 0
        for _ in 1:20
            _, p12 = HX.random_flag(K)
            x4 = HX.random_point(K)
            A = HX.evaluate_expressions(p12, x4)
            # p^(12) ⊆ (x^(4))^⊥ is excluded by the lemma, and then x^(2) = 0.
            all(iszero, A.x2) && continue
            count += 1
            @test HX.is_flag(A.x2, p12)
            @test HX.is_flag(A.x2, A.p23)
            @test HX.is_flag(A.x3, A.p23)
            @test HX.is_flag(A.x3, A.p34)
            @test HX.is_flag(x4, A.p34)
            @test iszero(HX.polar(A.x2, x4))
            @test !HX.proportional(A.x2, A.x3) && !HX.proportional(A.x3, x4) && !HX.proportional(A.x2, x4)
            @test !HX.proportional(p12, A.p23) && !HX.proportional(A.p23, A.p34)
        end
        @test count > 0
    end
end

@testset "hexagon_expressions: x^(2) = 0 if p^(12) ⊆ (x^(4))^⊥" begin
    Random.seed!(2)
    for K in (QQ, GF(2), GF(3))
        # A point on p^(12) is orthogonal to all of p^(12).
        x, p = HX.random_flag(K)
        @test all(iszero, HX.evaluate_expressions(p, x).x2)
    end
end
