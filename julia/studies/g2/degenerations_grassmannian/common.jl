###############################################################################
# Shared by the comparisons with the literature in this folder: the output
# file, the census orbit of each cell and the paper's type of each toric
# orbit of the tropical G2-Grassmannian.
###############################################################################

using TropicalGeometryG2
using Oscar

"""
Open the file `name` in `directory`. Returns `say`, which prints its arguments
to stdout and to the file, and the file, to close at the end.
"""
function output_file(directory, name)
    io = open(joinpath(directory, name), "w")
    say(parts...) = (println(parts...); println(io, parts...))
    return say, io
end

"The census directory of the Grassmannian over `field`, \"QQ\" or \"GF(2)\"."
grassmannian_census(field) = data_path("g2", field == "QQ" ? "grassmannian" : "grassmannian_gf2")

"The orbit of each cell of the Grassmannian census over QQ (stage 5)."
function cell_orbits()
    orbit_of = Dict(parse(Int, cell) => parse(Int, row["orbit"])
                    for row in csv_rows(joinpath(grassmannian_census("QQ"), "05_symmetry_cell_orbits.csv"))
                    for cell in split(row["cells"], ';'))
    return [orbit_of[cell] for cell in 1:length(orbit_of)]
end

"""
The paper's type of each toric orbit of the Grassmannian census over `field`,
by orbit: the types 1-20 of ../tables/toric_degeneration_table.tex, read from
data/g2/grassmannian_orbits.csv.
"""
function toric_types(field)
    column = field == "QQ" ? "qq_orbit" : "gf2_orbit"
    toric = Set(row["orbit"] for row in csv_rows(joinpath(grassmannian_census(field),
                                                          "06_canonical_interior_degenerations.csv"))
                if row["toric"] == "true")
    return Dict(parse(Int, row[column]) => parse(Int, row["orbit_type"])
                for row in csv_rows(data_path("g2", "grassmannian_orbits.csv"))
                if row["kind"] == "cell" && row[column] in toric)
end
