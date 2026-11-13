###############################################################################
# The numbers of the orbits of maximal cells of the tropical G2-Grassmannian
# in characteristic 0 and in characteristic 2, by which the tables of main.tex
# label them: the rows of the stage-6 censuses of both fields, and the labels,
# their order, the matching of the orbits of both fields and the support cell
# of each orbit from the G2-tree atlas local_flag_cell_atlas.csv. Only CSV
# files are read.
#
# Included by the tables of to_tex, through toric_degeneration_table.jl,
# and by cell_orbit_trees.jl, which adds the tree of each support cell.
###############################################################################

using TropicalGeometryG2: csv_fields, csv_rows

const TORIC_TABLE_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "..", ".."))
const TORIC_TABLE_CENSUS = (
    qq="julia/data/g2/grassmannian/06_canonical_interior_degenerations.csv",
    gf2="julia/data/g2/grassmannian_gf2/06_canonical_interior_degenerations.csv",
)
const TORIC_TABLE_ATLAS = "julia/studies/g2/atlas/flag/local_flag_cell_atlas.csv"
const TORIC_TABLE_FIELDS = (:qq, :gf2)
const TORIC_TABLE_FIELD_NAMES = (qq="QQ", gf2="GF(2)")

# The CSV reader of the package; toric_table_csv_fields serves the scripts that include this file.
toric_table_csv_fields(line) = csv_fields(line)

"""The rows of a CSV file of the repository, each a dictionary by column name."""
toric_table_csv_rows(path) = csv_rows(joinpath(TORIC_TABLE_ROOT, path))

"""The rows of the census of each field, by the number of their orbit."""
function toric_table_census()
    return map(TORIC_TABLE_CENSUS) do path
        Dict(parse(Int, row["orbit"]) => row for row in toric_table_csv_rows(path))
    end
end

"""The number of the atlas `label` of an orbit, by which the paper labels it."""
toric_table_label_number(label::AbstractString) = parse(Int, label)

"""
The support orbits of both fields in atlas order that are toric in some
field or, with `toric=false`, that are not toric in some field. A row has the
atlas `label`, the `census` `orbits` of these fields, their `orbit_size`, and
the atlas rays of the `support` cell that the atlas draws.
"""
function toric_table_support_rows(census=toric_table_census(); toric::Bool=true)
    pages = Dict{String,String}[]
    for row in toric_table_csv_rows(TORIC_TABLE_ATLAS)
        all(field -> field in ("qq", "gf2"), split(row["fields"], '+')) || continue
        any(page -> page["label"] == row["label"], pages) || push!(pages, row)
    end
    rows = []
    for page in pages
        orbits = Pair{Symbol,Int}[]
        for field in TORIC_TABLE_FIELDS
            orbit = page["$(field)_orbit"]
            (string(field) in split(page["fields"], '+')) == !isempty(orbit) ||
                error("atlas label $(page["label"]) has inconsistent fields")
            isempty(orbit) && continue
            (census[field][parse(Int, orbit)]["toric"] == "true") == toric &&
                push!(orbits, field => parse(Int, orbit))
        end
        isempty(orbits) && continue
        sizes = unique!([parse(Int, census[field][orbit]["orbit_size"])
                         for (field, orbit) in orbits])
        length(sizes) == 1 ||
            error("the orbit of atlas label $(page["label"]) differs between the fields")
        push!(rows, (
            label=page["label"], orbits=orbits, orbit_size=only(sizes),
            support=parse.(Int, split(page["support"], ';')),
        ))
    end
    for field in TORIC_TABLE_FIELDS
        expected = sort!([orbit for (orbit, row) in census[field]
                          if (row["toric"] == "true") == toric])
        listed = sort!([orbit for row in rows for (other, orbit) in row.orbits if other == field])
        expected == listed || error(
            "$TORIC_TABLE_ATLAS does not list the $(toric ? "toric" : "other") orbits of the " *
            "$(TORIC_TABLE_FIELD_NAMES[field]) census once each; rebuild the atlas",
        )
    end
    return rows
end
