###############################################################################
# The Chevalley polytopes of G2/B at varpi = omega_1 + omega_2 (Spacek--Wang,
# arXiv:2411.10276), one for each of the two reduced words of w_0, against the
# Newton--Okounkov polytopes of the toric degenerations of the G2 flag variety
# in the census.
#
# The valuation of a word is defined on all of C[G2/B]. Its values on the
# sections of omega_1 and of omega_2 are the two summands; their sums must be
# its values on omega_1 + omega_2. The comparison with the census is that of
# ../Feigin_Fourier_Littelmann/verify.jl, through ../flag_census.jl.
#
# Run from julia/ with
#   julia --project=. studies/g2/degenerations_grassmannian/Spacek_Wang/compare_flag.jl
# It writes compare_flag.txt next to this file. Characteristic 0 only.
###############################################################################

using TropicalGeometryG2
using Oscar
include(joinpath(@__DIR__, "..", "flag_census.jl"))

say, OUTPUT = output_file(@__DIR__, "compare_flag.txt")

"""
The values of the valuation of the paper word `word` of G2/B on the sections
of `varpi`, computed as in `chevalley_data`, which only accepts the word of
the parabolic of `varpi`.
"""
function valuation_values(varpi, word)
    L = lie_algebra(QQ, :G, 2)
    V = simple_module(L, varpi)
    expansion = torus_expansion(L, V, only(findall(==(varpi), basis_weights(L, V))), word)
    values = adapted_valuations(expansion.vector, expansion.ordering)
    @assert length(values) == dim(V)
    return values
end

invariants(P) = string("dimension ", dim(P), ", f-vector ", Int.(f_vector(P)), ", ",
                       length(lattice_points(P)), " lattice points, normalized volume ", normalized_volume(P))
"The convex hull of integer points in the affine lattice of their differences."
intrinsic_hull(points, dimension) = TropicalGeometryG2.affine_lattice_body(
    matrix(ZZ, reduce(vcat, permutedims.(points))); expected_dimension=dimension).polytope
# The Chevalley polytopes of G2/P1 and G2/P2, each for its only word.
const PARABOLIC = Dict(varpi => chevalley_polytope(chevalley_data(:G, 2, varpi)) for varpi in ([1, 0], [0, 1]))

R = root_system(:G, 2)
words = all_paper_words(minimal_coset_representative(weyl_group(R), parabolic_indices([1, 1]));
                        up_to_commutation=true)
@assert length(words) == 2
for word in words
    D = chevalley_data(:G, 2, [1, 1]; word=word)
    first_values, second_values, values = valuation_values.(([1, 0], [0, 1], [1, 1]), Ref(word))
    @assert values == D.valuation_image
    # The polytopes in the lattice of their own points, as the census bodies are.
    body = intrinsic_hull(values, 6)
    summands = intrinsic_hull.((first_values, second_values), 5)
    say("word ", join(word), ": Chevalley polytope P(omega_1+omega_2): ", invariants(body),
        "; values on omega_1 + omega_2 = values on omega_1 + values on omega_2: ",
        Set(values) == Set(a + b for a in first_values for b in second_values))
    for (summand, varpi, space) in zip(summands, ([1, 0], [0, 1]), ("G2/P1", "G2/P2"))
        say("  summand ", varpi, ": ", invariants(summand), "; lattice equivalent to the Chevalley polytope of ",
            space, ": ", affine_unimodular_equivalent(summand, PARABOLIC[varpi]))
    end
    cayley = pair_cayley_polytope(first_values, second_values)
    for field in ("QQ", "GF(2)")
        report_flag_census(say, field, flag_census_matches(field, body, cayley),
                           "the Chevalley polytope of word " * join(word))
    end
end
close(OUTPUT)
