# Pairwise affine unimodular equivalence of the Newton--Okounkov bodies of the toric orbits of the
# tropical G2-Grassmannian over QQ, for the table of Newton--Okounkov polytopes in main.tex.
# Row k is the k-th toric orbit of 06_canonical_interior_degenerations.csv.
#
#   julia --project=julia julia/studies/g2/supplemental/toric_body_equivalence.jl \
#     > julia/studies/g2/reports/toric_body_equivalence.txt
using Oscar
using TropicalCensus: affine_unimodular_equivalent

const DIR = joinpath(@__DIR__, "..", "..", "..", "data", "g2", "grassmannian")

lines = readlines(joinpath(DIR, "06_canonical_interior_degenerations.csv"))
header = split(lines[1], ',')
column(name) = findfirst(==(name), header)
rows = [split(line, ',') for line in lines[2:end]]
toric = [row for row in rows if row[column("toric")] == "true"]
orbits = [parse(Int, row[column("orbit")]) for row in toric]
bodies = [load(joinpath(DIR, row[column("newton_okounkov_body")])) for row in toric]

pairs = [(i, j) for i in eachindex(bodies) for j in (i + 1):length(bodies)]
equivalent = [(i, j) for (i, j) in pairs if affine_unimodular_equivalent(bodies[i], bodies[j])]
println("toric orbits (rows 1-", length(orbits), "): ", join(orbits, ", "))
println("pairs tested: ", length(pairs))
println("equivalent pairs of rows: ", isempty(equivalent) ? "none" :
        join(("$i and $j (orbits $(orbits[i]) and $(orbits[j]))" for (i, j) in equivalent), "; "))
