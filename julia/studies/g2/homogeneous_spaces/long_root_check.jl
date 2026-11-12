#=
Section 8 of main.tex: is there a basis of K^7 with a Pluecker coordinate that restricts to a long root line on the
G2-Grassmannian? Such a coordinate would be a decomposable 2-form in s(p14 - p25) + span(L1, ..., L7), where
L1, ..., L7 are the seven linear forms eq:7linear of main.tex. Decomposable means that all 4x4 Pfaffians vanish.
  long root eps1 - eps2, the form p14 - p25: the Pfaffians generate the unit ideal, so there is no such basis;
  control, the short root eps1, the form p14: the Pfaffians generate (a1, ..., a7), the single solution a = 0.

Usage, from the repository root (a few seconds after loading Oscar):
    julia --project=julia julia/studies/g2/homogeneous_spaces/long_root_check.jl
=#

using Oscar

R, a = polynomial_ring(QQ, ["a$i" for i in 1:7])

"The 4x4 Pfaffians of the 2-form with the entries w14 and w25 at (1,4) and (2,5), in the coefficients a of the L_i."
function pfaffians(w14, w25)
    w = zero_matrix(R, 7, 7)
    for ((i, j), v) in [(1, 2) => a[3], (1, 3) => a[2], (1, 4) => w14, (1, 7) => a[4],
                        (2, 3) => a[1], (2, 5) => w25, (2, 7) => a[5],
                        (3, 6) => a[7], (3, 7) => a[6],
                        (4, 5) => -a[6], (4, 6) => a[5], (4, 7) => a[1],
                        (5, 6) => -a[4], (5, 7) => -a[2],
                        (6, 7) => a[3]]
        w[i, j] = v
    end
    return ideal(R, [w[i, j] * w[k, l] - w[i, k] * w[j, l] + w[i, l] * w[j, k]
                     for i in 1:7 for j in i+1:7 for k in j+1:7 for l in k+1:7])
end

long_root = pfaffians(1 + a[7], -1 + a[7])
short_root = pfaffians(1 + a[7], a[7])
println("long root p14 - p25: ", gens(groebner_basis(long_root)))
println("short root p14: ", gens(groebner_basis(short_root)))
is_one(long_root) || error("the Pfaffians of the long root do not generate the unit ideal")
short_root == ideal(R, a) || error("the Pfaffians of the short root do not generate (a1, ..., a7)")
println("no Pluecker coordinate restricts to a long root line; for the short root only a = 0")
