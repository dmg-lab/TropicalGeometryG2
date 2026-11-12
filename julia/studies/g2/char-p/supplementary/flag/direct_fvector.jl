# f-vector and Weyl-orbit counts of the directly constructed tropical flag variety
# (direct-*.csv from direct.jl), independently of the census fan.
# Rays are stored projected orthogonally to the flag lineality space (W-invariant), so the
# taxon permutations act on them by permuting coordinates: x_i -> x_s(i), p_ij -> p_s(i)s(j).
using Oscar
cells = Vector{Vector{BigInt}}[]
for f in filter(s -> occursin(r"^direct-\d+\.csv$", s), readdir(@__DIR__))
  for line in eachline(joinpath(@__DIR__, f))
    (startswith(line, "tau") || startswith(line, "#")) && continue
    rs = split(split(line, "\"")[4], ";")
    push!(cells, [parse.(BigInt, split(r)) for r in rs])
  end
end
rays = unique(vcat(cells...))
rindex = Dict(r => k for (k, r) in enumerate(rays))
inc = IncidenceMatrix([[rindex[r] for r in C] for C in cells])
D = joinpath(@__DIR__, "../../../../../data/g2")
FF = load(joinpath(D, "flag/02_coarsest_fan.mrdi"))
L = Matrix{QQFieldElem}(generator_matrix(lineality_space(FF)))
PF = polyhedral_fan(inc, matrix(QQ, reduce(vcat, [permutedims(QQ.(r)) for r in rays])), matrix(QQ, L))
println("constructed cells: ", length(cells), ", rays: ", length(rays))
println("f-vector (mod lineality): ", f_vector(PF), "   census: (141, 1352, 3584, 2714)")
println("maximal cells by number of rays: ", sort(collect(Dict(k => count(C -> length(C) == k, cells) for k in unique(length.(cells))))),
        "   census: 2564 tetrahedra, 108 square pyramids (5 rays), 42 triangular prisms (6 rays)")
# Weyl group: the taxon permutations of the census (Grassmannian, fixing taxon 7)
perms = load(joinpath(D, "grassmannian/05_symmetry_taxon_permutations.mrdi"))
PAIRS = [(i, j) for i in 1:6 for j in (i + 1):7]; IDX = Dict(p => k for (k, p) in enumerate(PAIRS))
function act(s, r)
  v = similar(r)
  for i in 1:7; v[s[i]] = r[i]; end
  for (k, (i, j)) in enumerate(PAIRS)
    a, b = minmax(s[i], s[j]); v[7 + IDX[(a, b)]] = r[7 + k]
  end
  v
end
sl = [[Int(perms[i, j]) for j in 1:7] for i in 1:nrows(perms)]
# combinatorial types: (number of rays, number of facets) of each maximal cell
types = Dict{Tuple{Int,Int},Int}()
for C in maximal_cones(PF)
  t = (length(rays_modulo_lineality(C).rays_modulo_lineality), n_facets(C)); types[t] = get(types, t, 0) + 1
end
println("maximal cells by (rays, facets): ", sort(collect(types)),
        "   [tetrahedron (4,4), square pyramid (5,5), triangular bipyramid (5,6), triangular prism (6,5), octahedron (6,8)]")
println("number of taxon permutations: ", length(sl))
cellset = Set(Set(C) for C in cells)
orbit_rep(C) = minimum(sort([sort([rindex[act(s, r)] for r in C]) for s in sl]))
@assert all(Set(Set(act(s, r) for r in C)) in cellset for C in cells for s in sl)
println("Weyl orbits on maximal cells: ", length(unique(orbit_rep.(cells))), "   census: 270")
println("Weyl orbits on rays: ", length(unique([minimum(rindex[act(s, r)] for s in sl) for r in rays])), "   census: 27")
