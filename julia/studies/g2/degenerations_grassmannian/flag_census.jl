###############################################################################
# A candidate Newton--Okounkov polytope of G2/B at omega_1 + omega_2, given
# with its two summands, against the toric orbits of the flag census.
#
# The summands are compared through Cayley polytopes: an affine unimodular map
# of Cayley polytopes keeping each summand in its layer is an isomorphism of
# the bigraded semigroups, i.e. of the toric varieties in P^6 x P^13.
#
# Included by Feigin_Fourier_Littelmann/verify.jl and Spacek_Wang/compare_flag.jl.
###############################################################################

using TropicalCensus
using Oscar

vertex_set(P) = Set(Vector{QQFieldElem}(collect(v)) for v in vertices(P))

"""
The Cayley polytope of the degree-one points, rows of `values`, in the affine
lattice of their differences, and the summand (1 or 2) of each of its vertices.
"""
function cayley_polytope(values, summand)
    body = TropicalCensus.affine_lattice_body(values; expected_dimension=7)
    point_summand = Dict([Int(body.points[i, j]) for j in 1:7] => summand[i] for i in 1:nrows(values))
    return body.polytope, [point_summand[Int.(collect(v))] for v in vertices(body.polytope)]
end

"The Cayley polytope of integer points of degree (1,0) (omega_1) and (0,1) (omega_2)."
function pair_cayley_polytope(first_points, second_points)
    rows = vcat([[[1, 0]; p] for p in first_points], [[[0, 1]; q] for q in second_points])
    return cayley_polytope(matrix(ZZ, reduce(vcat, permutedims.(rows))),
                           [fill(1, length(first_points)); fill(2, length(second_points))])
end

"""
Is there an affine unimodular map between Cayley polytopes keeping the vertices
of each summand in that summand? Same search as `affine_unimodular_equivalent`:
one face-lattice isomorphism followed by every combinatorial automorphism of
the target.
"""
function summand_preserving_equivalent((source_polytope, source_summand), (target_polytope, target_summand))
    f_vector(source_polytope) == f_vector(target_polytope) || return false
    source, target = matrix(QQ, vertices(source_polytope)), matrix(QQ, vertices(target_polytope))
    found = Oscar.Polymake.polytope.find_facet_vertex_permutations(
        Oscar.pm_object(source_polytope), Oscar.pm_object(target_polytope))
    isnothing(found) && return false
    permutation = Int.(collect(last(found))) .+ 1
    target_facets = TropicalCensus.facet_vertex_sets(target_polytope)
    initial = first(filter(seed -> all(facet -> sort!(seed[facet]) in target_facets,
                                       TropicalCensus.facet_vertex_sets(source_polytope)),
                           [permutation, invperm(permutation)]))
    for automorphism in automorphism_group(target_polytope; action=:on_vertices)
        mapping = [Int(automorphism(index)) for index in 1:nrows(target)][initial]
        source_summand == target_summand[mapping] &&
            TropicalCensus.affine_map_is_unimodular(source, target, mapping) && return true
    end
    return false
end

"""
The toric orbits of the flag census over `field` ("QQ" or "GF(2)") against
`body` and its Cayley polytope `cayley`: their number, the orbits with the
f-vector of `body`, and for each orbit whose body is lattice equivalent to
`body` its census row and whether the Cayley polytopes are equivalent.
Coordinates 1:7 of the census are the vector coordinates (omega_1), 8:28 the
Pluecker coordinates (omega_2).
"""
function flag_census_matches(field, body, cayley)
    directory = data_path("g2", field == "QQ" ? "flag" : "flag_gf2")
    fan = load(joinpath(directory, "03_normalized_coarsest_fan.mrdi"))
    normalized_rays = load(joinpath(directory, "03_normalized_rays.mrdi"))
    lineality = Matrix{QQFieldElem}(matrix(QQ, rays_modulo_lineality(fan).lineality_basis))
    toric = filter(row -> row["toric"] == "true",
                   TropicalCensus.csv_rows(joinpath(directory, "06_canonical_interior_degenerations.csv")))
    same_f_vector, matches = Int[], NamedTuple[]
    for row in toric
        census_body = load(joinpath(directory, row["newton_okounkov_body"]))
        f_vector(census_body) == f_vector(body) && push!(same_f_vector, parse(Int, row["orbit"]))
        affine_unimodular_equivalent(body, census_body) || continue
        rays = parse.(Int, split(row["rays"], ';'))
        values = TropicalCensus.cell_values(lineality, normalized_rays, rays)
        # The 28 values give back the stored body: Minkowski sum of the two slices.
        @assert vertex_set(TropicalCensus.intrinsic_cell_body(
            lineality, normalized_rays, rays; gradings=[1:7, 8:28]).polytope) == vertex_set(census_body)
        push!(matches, (orbit=parse(Int, row["orbit"]), row=row,
                        pair=summand_preserving_equivalent(cayley, cayley_polytope(values, [fill(1, 7); fill(2, 21)]))))
    end
    return (toric=length(toric), same_f_vector=same_f_vector, matches=matches)
end

"Print `flag_census_matches` for the polytope called `name`."
function report_flag_census(say, field, census, name)
    for match in census.matches
        say(field, " flag orbit ", match.orbit, " (", match.row["orbit_size"], " cells, body ",
            match.row["body_type"], ", slices ", match.row["body_slice_types"],
            "): lattice equivalent to ", name, "; as the pair of summands under one map: ", match.pair)
    end
    say(field, " flag: ", census.toric, " toric orbits, ", length(census.same_f_vector),
        " with the f-vector of ", name, " (", join(census.same_f_vector, ", "),
        "), lattice equivalent to it: ", join((match.orbit for match in census.matches), ", "))
end
