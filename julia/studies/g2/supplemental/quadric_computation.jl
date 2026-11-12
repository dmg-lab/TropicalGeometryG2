using Oscar

# Coefficient field QQ(x1, ..., x7).
S, polynomial_variables = polynomial_ring(QQ, [Symbol("x$i") for i in 1:7])
K = fraction_field(S)
x1, x2, x3, x4, x5, x6, x7 = K.(polynomial_variables)

# Exterior algebra K[e1, ..., e7].
E, (e1, e2, e3, e4, e5, e6, e7) =
    exterior_algebra(K, [Symbol("e$i") for i in 1:7])

phi =
    e1 * e2 * e3 +
    e4 * e5 * e6 +
    e1 * e4 * e7 +
    e2 * e5 * e7 +
    e3 * e6 * e7

C =
    x1 * (e2 * e3 + e4 * e7) -
    x2 * (e1 * e3 - e5 * e7) +
    x3 * (e1 * e2 + e6 * e7) -
    x4 * (e1 * e7 - e5 * e6) -
    x5 * (e2 * e7 + e4 * e6) -
    x6 * (e3 * e7 - e4 * e5) +
    x7 * (e1 * e4 + e2 * e5 + e3 * e6)

result = C * C * phi

quadric = 6 * (x1 * x4 + x2 * x5 + x3 * x6 - x7^2)
volume_form = e1 * e2 * e3 * e4 * e5 * e6 * e7
@assert result == quadric * volume_form

println("C*C*phi = ", result)
println("Top exterior coefficient = ", quadric)
