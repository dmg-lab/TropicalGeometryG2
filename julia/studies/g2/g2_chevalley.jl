"""
Whether the Newton-Okounkov body of a toric degeneration of G2/P2 is unimodularly
equivalent to the Chevalley polytope of Spacek and Wang (cell type 1 of main.tex).
Shared by the G2 census drivers.
"""
function g2_chevalley_simplex(body)
    source = matrix(QQ, matrix(ZZ, vertices(body.polytope)))
    size(source) == (6, 5) || return false
    target = matrix(QQ, [
        0 0 0 0 0; 1 0 0 0 0; 1 3 0 0 0;
        1 3 2 0 0; 1 3 2 3 0; 1 3 2 3 1
    ])
    source_edges = matrix(QQ, 5, 5, [source[i, j] - source[1, j] for i in 2:6 for j in 1:5])
    iszero(det(source_edges)) && return false
    for target_base in 1:6
        remaining = [index for index in 1:6 if index != target_base]
        for order in all_permutations(remaining)
            target_edges = matrix(QQ, 5, 5, [target[i, j] - target[target_base, j] for i in order for j in 1:5])
            linear = inv(source_edges) * target_edges
            all(denominator(linear[i, j]) == 1 for i in 1:5 for j in 1:5) || continue
            abs(det(linear)) == 1 || continue
            translation = [target[target_base, j] - sum(source[1, k] * linear[k, j] for k in 1:5) for j in 1:5]
            mapped = Set(Tuple(sum(source[i, k] * linear[k, j] for k in 1:5) + translation[j] for j in 1:5) for i in 1:6)
            mapped == Set(Tuple(target[i, j] for j in 1:5) for i in 1:6) && return true
        end
    end
    return false
end
