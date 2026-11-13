const G2_CHEVALLEY = Ref{Any}(nothing)

"""
Whether the Newton-Okounkov body of a toric degeneration of G2/P2 is unimodularly
equivalent to the Chevalley polytope of Spacek and Wang (cell type 1 of main.tex),
P_{G2/P2, varpi_2, s} for s = s_2 s_1 s_2 s_1 s_2, the only reduced word of w^P.
Shared by the G2 census drivers.
"""
function g2_chevalley_simplex(body)
    if isnothing(G2_CHEVALLEY[])
        G2_CHEVALLEY[] = chevalley_polytope(chevalley_data(:G, 2, [0, 1]; word=[2, 1, 2, 1, 2]))
    end
    return affine_unimodular_equivalent(body.polytope, G2_CHEVALLEY[])
end
