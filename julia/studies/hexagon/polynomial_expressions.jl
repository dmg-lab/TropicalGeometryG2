using Oscar

# The polynomial functions of Lemma 9.1 (lem:complete_hexagon) of main.tex. Given a line p = p^(12) in Plücker
# coordinates p_ij and a point x = x^(4) = (x1 : ... : x7), we construct, in this order,
#
#     x^(2) = the point of p^(12) orthogonal to x^(4),
#     x^(3) = x^(2) × x^(4),
#     p^(23) = x^(2) ∧ x^(3),   p^(34) = x^(3) ∧ x^(4),
#
# as polynomials in K[p_ij, x_k]. Conventions are those of the paper:
#
#     q(x)    = x1 x4 + x2 x5 + x3 x6 - x7^2                                    (eq:quadric)
#     <x,y>   = x1 y4 + x2 y5 + x3 y6 + x4 y1 + x5 y2 + x6 y3 - 2 x7 y7          (polar form of q)
#     Φ       = e123 + e456 + e147 + e257 + e367                                 (eq:trilinear)
#     x × y   = (x ∧ y) ⌟ Φ, identified with a vector via <.,.>, i.e. <x × y, z> = Φ(x, y, z)
#
# and p_ii = 0, p_ij = -p_ji. The cross product involves 1/2, but it cancels in x^(3): all expressions have integer
# coefficients, so they define the same functions over every field, also in characteristic 2.
#
# From the `julia` folder:  julia --project=. studies/hexagon/polynomial_expressions.jl [--lines]
# prints the coordinates of x^(2), x^(3) as in Remark 9.2 (with --lines also p^(23), p^(34), several thousand terms).

const PAIRS = [(i, j) for i in 1:7 for j in (i + 1):7]
const P_INDEX = Dict(pair => k for (k, pair) in enumerate(PAIRS))

const EXPRESSION_RING, _vars = polynomial_ring(QQ, vcat(["p$i$j" for (i, j) in PAIRS], ["x$k" for k in 1:7]))
const P_VARS = _vars[1:21]
const X_VARS = _vars[22:28]

"Entry (i, j) of the antisymmetric matrix of a Plücker vector `L` (21 entries, lexicographic in i < j)."
pl(L, i, j) = i == j ? zero(L[1]) : i < j ? L[P_INDEX[(i, j)]] : -L[P_INDEX[(j, i)]]

"The quadratic form q of (eq:quadric)."
quadric(x) = x[1] * x[4] + x[2] * x[5] + x[3] * x[6] - x[7]^2

"The polar form <x, y> of q."
polar(x, y) = x[1] * y[4] + x[2] * y[5] + x[3] * y[6] + x[4] * y[1] + x[5] * y[2] + x[6] * y[3] - 2 * x[7] * y[7]

det3(x, y, z, (i, j, k)) = x[i] * (y[j] * z[k] - y[k] * z[j]) - x[j] * (y[i] * z[k] - y[k] * z[i]) +
                           x[k] * (y[i] * z[j] - y[j] * z[i])

"The trilinear form Φ = e123 + e456 + e147 + e257 + e367."
Phi(x, y, z) = sum(det3(x, y, z, t) for t in ((1, 2, 3), (4, 5, 6), (1, 4, 7), (2, 5, 7), (3, 6, 7)))

basis_vector(k, K=QQ) = [i == k ? one(K) : zero(K) for i in 1:7]

"""
The cross product x × y: the vector z with <z, e_k> = Φ(x, y, e_k) for all k. Since <z, e_k> is z_{k+3},
z_{k-3}, -2 z_7 for k ≤ 3, 4 ≤ k ≤ 6, k = 7, this is a permutation of the covector, with z_7 = -Φ(x, y, e_7)/2.
"""
function cross(x, y)
    phi = [Phi(x, y, basis_vector(k)) for k in 1:7]
    return vcat(phi[4:6], phi[1:3], [-QQ(1, 2) * phi[7]])
end

"Plücker vector of the line through x and y."
wedge(x, y) = [x[i] * y[j] - x[j] * y[i] for (i, j) in PAIRS]

# --- The construction of Lemma 9.1, step by step ---------------------------------------------------------

"Step 1. x^(2)_i = sum_j <x^(4), e_j> p_ji = x4 p_1i + x5 p_2i + x6 p_3i + x1 p_4i + x2 p_5i + x3 p_6i - 2 x7 p_7i."
step_x2(p, x) = [sum(polar(x, basis_vector(j)) * pl(p, j, i) for j in 1:7) for i in 1:7]

"Step 2. x^(3) = x^(2) × x^(4)."
step_x3(x2, x) = cross(x2, x)

"Step 3. The lines p^(23) through x^(2), x^(3) and p^(34) through x^(3), x^(4)."
step_lines(x2, x3, x) = (wedge(x2, x3), wedge(x3, x))

const _HEXAGON_EXPRESSIONS = Ref{Any}(nothing)

"""
    hexagon_expressions()

The coordinates of `x2`, `x3`, `p23`, `p34` of Lemma 9.1 as polynomials in K[p_ij, x_k] (`EXPRESSION_RING`, over QQ with
integer coefficients). Built on first use and cached.
"""
function hexagon_expressions()
    if _HEXAGON_EXPRESSIONS[] === nothing
        x2 = step_x2(P_VARS, X_VARS)
        x3 = step_x3(x2, X_VARS)
        p23, p34 = step_lines(x2, x3, X_VARS)
        _HEXAGON_EXPRESSIONS[] = (x2=x2, x3=x3, p23=p23, p34=p34)
    end
    return _HEXAGON_EXPRESSIONS[]
end

"The value of a polynomial with integer coefficients at a point `v` with entries in any field."
function evaluate_integral(f, v)
    K = parent(v[1])
    S, _ = polynomial_ring(K, symbols(parent(f)); cached=true)
    return evaluate(map_coefficients(c -> K(ZZ(c)), f; parent=S), v)
end

"""
    evaluate_expressions(p, x)

The polynomial functions of Lemma 9.1 at a line `p` (21 entries) and a point `x` (7 entries) over any field: returns
`(x2, x3, p23, p34)`. The polynomials `p23`, `p34` are by definition wedge products, so they are evaluated as such.
"""
function evaluate_expressions(p, x)
    H = hexagon_expressions()
    v = vcat(p, x)
    x2 = [evaluate_integral(f, v) for f in H.x2]
    x3 = [evaluate_integral(f, v) for f in H.x3]
    p23, p34 = step_lines(x2, x3, x)
    return (x2=x2, x3=x3, p23=p23, p34=p34)
end

# --- The ideal of linespace × pointspace (Remark 9.2) ---------------------------------------------------

"The 7 linear forms of (eq:7linear)."
linear_forms(p) = [pl(p, 2, 3) + pl(p, 4, 7), pl(p, 1, 3) - pl(p, 5, 7), pl(p, 1, 2) + pl(p, 6, 7),
                   pl(p, 1, 7) - pl(p, 5, 6), pl(p, 2, 7) + pl(p, 4, 6), pl(p, 3, 7) - pl(p, 4, 5),
                   pl(p, 1, 4) + pl(p, 2, 5) + pl(p, 3, 6)]

"The 35 Plücker relations."
plucker_relations(p) = [pl(p, i, j) * pl(p, k, l) - pl(p, i, k) * pl(p, j, l) + pl(p, i, l) * pl(p, j, k)
                        for i in 1:7 for j in (i + 1):7 for k in (j + 1):7 for l in (k + 1):7]

"The 35 relations of x lying on the line p, i.e. x ∧ p = 0."
incidence_relations(x, p) = [x[i] * pl(p, j, k) - x[j] * pl(p, i, k) + x[k] * pl(p, i, j)
                             for i in 1:7 for j in (i + 1):7 for k in (j + 1):7]

"The ideal I of Remark 9.2, cutting out linespace × pointspace in `EXPRESSION_RING`."
line_point_ideal() = ideal(EXPRESSION_RING, vcat(linear_forms(P_VARS), plucker_relations(P_VARS), [quadric(X_VARS)]))

# --- Readable output -------------------------------------------------------------------------------------

"Bidegree (degree in the p_ij, degree in the x_k) of a nonzero bihomogeneous polynomial, or `nothing`."
function bidegree(f)
    iszero(f) && return nothing
    degs = unique((sum(e[1:21]), sum(e[22:28])) for e in AbstractAlgebra.exponent_vectors(f))
    return length(degs) == 1 ? only(degs) : error("not bihomogeneous: $degs")
end

latex_variable(k) = k <= 21 ? "p_{$(PAIRS[k][1])$(PAIRS[k][2])}" : "x_$(k - 21)"

function latex_monomial(e)
    parts = [latex_variable(k) * (e[k] > 1 ? "^$(e[k])" : "") for k in eachindex(e) if e[k] > 0]
    return isempty(parts) ? "1" : join(parts)
end

function latex_coefficient(c, first::Bool)
    sign = c < 0 ? (first ? "-" : " - ") : (first ? "" : " + ")
    a = abs(c)
    return sign * (isone(a) ? "" : string(numerator(a)))
end

"""
The polynomial `f` as LaTeX, as in Remark 9.2: terms in lex order (p_12 > ... > p_67 > x_1 > ... > x_7),
`terms_per_line` terms per line separated by `\\\\`.
"""
function latex_expression(f; terms_per_line=8)
    iszero(f) && return "0"
    out = IOBuffer()
    for (n, (c, e)) in enumerate(zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f)))
        n > 1 && (n - 1) % terms_per_line == 0 && print(out, " \\\\\n")
        print(out, latex_coefficient(c, n == 1), latex_monomial(e))
    end
    return String(take!(out))
end

"Print the coordinates of `v` (a point, 7 entries, or a line, 21 Plücker entries) with term counts and bidegrees."
function print_coordinates(name, v)
    labels = length(v) == 7 ? string.(1:7) : ["$i$j" for (i, j) in PAIRS]
    for (label, f) in zip(labels, v)
        println("($name)_$label   [$(length(f)) terms, bidegree $(bidegree(f))]")
        println(latex_expression(f), "\n")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    H = hexagon_expressions()
    print_coordinates("x^(2)", H.x2)
    print_coordinates("x^(3)", H.x3)
    if "--lines" in ARGS
        print_coordinates("p^(23)", H.p23)
        print_coordinates("p^(34)", H.p34)
    end
end
