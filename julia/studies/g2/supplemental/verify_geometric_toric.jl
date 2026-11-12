#!/usr/bin/env julia

# Verify the stage-6 toric labels after extending scalars from QQ to CC.
# A binomial prime avoiding the coordinate hyperplanes is geometrically prime
# exactly when its exponent-difference lattice is saturated.

using Oscar

const DATA = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "g2"))

# Split one CSV line, honouring double-quoted fields such as `body_type`.
function csv_fields(line)
    characters = collect(line)
    fields = String[]
    buffer = IOBuffer()
    quoted = false
    index = 1
    while index <= length(characters)
        character = characters[index]
        if character == '"'
            if quoted && index < length(characters) && characters[index + 1] == '"'
                print(buffer, '"')
                index += 1
            else
                quoted = !quoted
            end
        elseif character == ',' && !quoted
            push!(fields, String(take!(buffer)))
        else
            print(buffer, character)
        end
        index += 1
    end
    quoted && error("unterminated CSV quote")
    push!(fields, String(take!(buffer)))
    return fields
end

function census_rows(target)
    path = joinpath(DATA, target, "06_canonical_interior_degenerations.csv")
    lines = readlines(path)
    isempty(lines) && error("empty census: $path")
    header = split(first(lines), ',')
    columns = Dict(name => index for (index, name) in enumerate(header))
    needed = ("orbit", "initial_ideal", "prime", "binomial", "monomial_free", "toric")
    all(haskey(columns, name) for name in needed) || error("unexpected census header: $path")
    return [
        begin
            fields = csv_fields(line)
            length(fields) == length(header) || error("unexpected CSV column count in $path")
            Dict(name => fields[index] for (name, index) in columns)
        end
        for line in lines[2:end]
    ]
end

function exponent_lattice_torsion(ideal_value)
    basis = gens(groebner_basis(ideal_value; complete_reduction=true))
    all(length(polynomial) == 2 for polynomial in basis) ||
        error("a purported binomial prime has a nonbinomial reduced Groebner basis")
    differences = [
        let exponents_of_terms = collect(exponents(polynomial))
            exponents_of_terms[1] - exponents_of_terms[2]
        end
        for polynomial in basis
    ]
    nvariables = ngens(base_ring(ideal_value))
    lattice = matrix(ZZ, length(differences), nvariables, reduce(vcat, differences))
    smith = snf(lattice)
    nontrivial = [
        smith[index, index]
        for index in 1:min(nrows(smith), ncols(smith))
        if !iszero(smith[index, index]) && abs(smith[index, index]) != 1
    ]
    return nontrivial
end

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
        torsion = exponent_lattice_torsion(load(ideal_path))
        isempty(torsion) ||
            error("$target orbit $(row["orbit"]) splits over CC: Smith factors $torsion")
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
