using Oscar

# Completing a classical hexagon over any field, by the polynomial functions of Lemma 9.1 (polynomial_expressions.jl).
#
# A hexagon has points x^(1), ..., x^(6) and lines p^(12), p^(23), ..., p^(61). Two flags (x^(1), p^(12)) and
# (x^(4), p^(45)) with x^(1), x^(4) opposite determine it: Lemma 9.1 applied to (p^(12), x^(4)) gives
# x^(2), x^(3), p^(23), p^(34), and applied to (p^(45), x^(1)) it gives x^(5), x^(6), p^(56), p^(61).

include(joinpath(@__DIR__, "polynomial_expressions.jl"))

"Whether x lies on the line with Plücker vector p, i.e. x ∧ p = 0."
on_line(x, p) = all(iszero, incidence_relations(x, p))

"Whether p is a point of the G2-Grassmannian: the Plücker relations and the linear space (eq:7linear)."
in_linespace(p) = any(!iszero, p) && all(iszero, linear_forms(p)) && all(iszero, plucker_relations(p))

"Whether (x, p) is a point of the flag variety: x on the quadric, p in the G2-Grassmannian, x on p."
is_flag(x, p) = any(!iszero, x) && iszero(quadric(x)) && in_linespace(p) && on_line(x, p)

"Whether two vectors are equal as points of projective space."
function proportional(x, y)
    n = length(x)
    return all(iszero(x[i] * y[j] - x[j] * y[i]) for i in 1:n for j in (i + 1):n)
end

"""
A vector over a field up to scaling, in a canonical form: over QQ and QQ(t) entries are coprime and without denominators
(over QQ(t) polynomials in t with coprime integer coefficients); over other fields the first nonzero entry is 1.
"""
function projective_normal_form(v::Vector{QQFieldElem})
    d = reduce(lcm, denominator.(v); init=ZZ(1))
    g = reduce(gcd, numerator.(v .* d); init=ZZ(0))
    return v .* (d // g)
end

function projective_normal_form(v::Vector{<:AbstractAlgebra.Generic.RationalFunctionFieldElem})
    K = parent(v[1])
    d = reduce(lcm, denominator.(v))
    n = numerator.(v .* d)
    g = reduce(gcd, n)
    content_lcm = reduce(lcm, [denominator(c) for f in n for c in coefficients(divexact(f, g))]; init=ZZ(1))
    content_gcd = reduce(gcd, [numerator(c * content_lcm) for f in n for c in coefficients(divexact(f, g))]; init=ZZ(0))
    return v .* (d // g) .* K(content_lcm // content_gcd)
end

function projective_normal_form(v::Vector{<:FieldElem})
    return v ./ v[findfirst(!iszero, v)]
end

"""
    complete_hexagon(x1, p12, x4, p45; check=true, normalize=true)
    complete_hexagon((x1, p12), (x4, p45); ...)

The hexagon through the flags (x^(1), p^(12)) and (x^(4), p^(45)), with points as 7-vectors and lines as 21-vectors of
Plücker coordinates (lexicographic in i < j), all with entries in one field. Returns

- `points`: `[x1, ..., x6]`,
- `lines`: `[p12, p23, p34, p45, p56, p61]`, where `lines[k]` joins `points[k]` and `points[k + 1]`,
- `flags`: `[(points[k], lines[k]) for k in 1:6]`, so `flags[1]` and `flags[4]` are the input flags.

With `normalize` all vectors are returned in `projective_normal_form`. With `check`, an `ArgumentError` is thrown if an
input is not a flag or the result is not a non-degenerate hexagon.
"""
function complete_hexagon(x1, p12, x4, p45; check::Bool=true, normalize::Bool=true)
    if check
        is_flag(x1, p12) || throw(ArgumentError("(x1, p12) is not a flag"))
        is_flag(x4, p45) || throw(ArgumentError("(x4, p45) is not a flag"))
    end
    A = evaluate_expressions(p12, x4)
    B = evaluate_expressions(p45, x1)
    all(iszero, A.x2) && throw(ArgumentError("not generic: p12 is contained in the orthogonal space of x4"))
    all(iszero, B.x2) && throw(ArgumentError("not generic: p45 is contained in the orthogonal space of x1"))
    points = [x1, A.x2, A.x3, x4, B.x2, B.x3]
    lines = [p12, A.p23, A.p34, p45, B.p23, B.p34]
    if check
        for k in 1:6
            is_flag(points[k], lines[k]) || throw(ArgumentError("not generic: (x$k, line $k) is not a flag"))
            on_line(points[k], lines[mod1(k - 1, 6)]) || throw(ArgumentError("not generic: x$k is not on line $(mod1(k - 1, 6))"))
            for l in (k + 1):6
                proportional(points[k], points[l]) && throw(ArgumentError("not generic: x$k and x$l coincide"))
                proportional(lines[k], lines[l]) && throw(ArgumentError("not generic: lines $k and $l coincide"))
            end
        end
    end
    if normalize
        points, lines = projective_normal_form.(points), projective_normal_form.(lines)
    end
    return (points=points, lines=lines, flags=[(points[k], lines[k]) for k in 1:6])
end

complete_hexagon(f1::Tuple, f4::Tuple; kwargs...) = complete_hexagon(f1..., f4...; kwargs...)

# --- Random flags over any field -------------------------------------------------------------------------

"A random element for test data: small integers over QQ, uniform over a finite field."
random_element(K::QQField) = K(rand(-9:9))
random_element(K::FinField) = rand(K)

"""
A random element of K = QQ(t): `random_element(QQ) t^a` plus a term of higher t-order, where `a` is drawn from `orders`.
"""
function random_element(K::AbstractAlgebra.Generic.RationalFunctionField; orders=0:40)
    t = gen(K)
    a = rand(orders)
    return K(rand(-9:9)) * t^a + K(rand(-9:9)) * t^(a + 1)
end

"A random point of the quadric: all coordinates random, except for x4 which is solved for."
function random_point(K; sample=() -> random_element(K))
    while true
        x = [sample() for _ in 1:7]
        iszero(x[1]) && continue
        x[4] = (x[7]^2 - x[2] * x[5] - x[3] * x[6]) / x[1]
        return x
    end
end

"""
    random_flag(K; sample=() -> random_element(K))

A random flag (x, p) over the field K: a random point x and the line through x and a random point y with x ∧ y in the
linear space (eq:7linear). These y form the plane of all special lines through x, the kernel of a 7×7 matrix.
"""
function random_flag(K; sample=() -> random_element(K))
    while true
        x = random_point(K; sample)
        M = matrix(K, [linear_forms(wedge(x, basis_vector(j, K)))[i] for i in 1:7, j in 1:7])
        N = kernel(M; side=:right)
        c = N * matrix(K, ncols(N), 1, [sample() for _ in 1:ncols(N)])
        y = [c[i, 1] for i in 1:7]
        p = wedge(x, y)
        is_flag(x, p) && return (x, p)
    end
end
