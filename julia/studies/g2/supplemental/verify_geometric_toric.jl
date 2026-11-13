#!/usr/bin/env julia

# Verify the stage-6 toric labels after extending scalars from QQ to CC.
# A binomial prime avoiding the coordinate hyperplanes is geometrically prime
# exactly when its exponent-difference lattice is saturated.

using Oscar
using TropicalGeometryG2: csv_rows, binomial_lattice_index

const DATA = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "g2"))

"The stage-6 census of `target`, one dictionary per orbit."
function census_rows(target)
    path = joinpath(DATA, target, "06_canonical_interior_degenerations.csv")
    rows = csv_rows(path)
    needed = ("orbit", "initial_ideal", "prime", "binomial", "monomial_free", "toric")
    all(haskey(row, name) for row in rows for name in needed) || error("unexpected census header: $path")
    return rows
end

"Index of the exponent-difference lattice of the reduced Groebner basis in its saturation."
lattice_index(ideal_value) =
    binomial_lattice_index(gens(groebner_basis(ideal_value; complete_reduction=true)))

for target in ("grassmannian", "flag")
    rows = census_rows(target)
    toric = filter(row -> row["toric"] == "true", rows)
    all(row["monomial_free"] == "true" for row in rows) ||
        error("a canonical weight is not tropical for $target")
    count(row -> row["prime"] == "true" && row["binomial"] == "true", rows) == length(toric) ||
        error("toric labels disagree with prime binomial labels for $target")
    for row in toric
        row["prime"] == "true" && row["binomial"] == "true" ||
            error("inconsistent toric label for $target orbit $(row["orbit"])")
        ideal_path = joinpath(DATA, target, row["initial_ideal"])
        index = lattice_index(load(ideal_path))
        index == 1 ||
            error("$target orbit $(row["orbit"]) splits over CC: lattice index $index")
    end
    println("$target: $(length(toric)) geometrically prime binomial ideals / $(length(rows)) orbits")
end

exception = only(filter(row -> row["prime"] == "true" && row["binomial"] == "false",
                        census_rows("grassmannian")))
exception["orbit"] == "31" || error("unexpected prime nonbinomial orbit")
ideal_value = load(joinpath(DATA, "grassmannian", exception["initial_ideal"]))
is_prime(ideal_value) && !is_binomial(ideal_value) || error("exceptional QQ classification changed")
decomposition = absolute_primary_decomposition(ideal_value)
length(decomposition) == 1 || error("unexpected absolute primary decomposition")
_, _, component, extension_degree = only(decomposition)
extension_degree == 2 || error("unexpected number of geometric components")
!is_binomial(component) || error("geometric component unexpectedly binomial")
Oscar.degree(component) == 9 || error("unexpected geometric component degree")
field_polynomial = defining_polynomial(coefficient_ring(base_ring(component)))
println("grassmannian orbit 31: prime nonbinomial over QQ; ",
        "two nonbinomial conjugate components over CC, each of degree 9; ",
        "field polynomial = $field_polynomial")
